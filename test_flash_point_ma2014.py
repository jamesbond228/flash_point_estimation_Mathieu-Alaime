"""Unit tests for flash_point_ma2014."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from flash_point_ma2014 import (
    DEFAULT_PARAMETER_PATH,
    NegativeTermPolicy,
    ParameterError,
    analyze_fragments,
    convert_temperature,
    load_parameters,
    predict_flash_point,
)
from rdkit import Chem


@pytest.fixture(scope="module")
def model():
    return load_parameters()


def _manual_tf(counts: dict[str, int], model) -> float:
    """Independent re-computation of Eq.(2) for positive t_i."""
    f = sum(
        n * model.parameters[k].t_k ** (1 / model.alpha)
        for k, n in counts.items()
    )
    return f**model.alpha


# ---- parameters ----------------------------------------------------------
def test_parameters_load_21_entries(model):
    assert len(model.parameters) == 21
    assert model.alpha == pytest.approx(0.4)


def _write_variant(tmp_path: Path, mutate) -> Path:
    raw = json.loads(DEFAULT_PARAMETER_PATH.read_text(encoding="utf-8"))
    mutate(raw)
    out = tmp_path / "p.json"
    out.write_text(json.dumps(raw), encoding="utf-8")
    return out


def test_missing_parameter_fails(tmp_path):
    path = _write_variant(tmp_path, lambda r: r["parameters"].pop("C4-3"))
    with pytest.raises(ParameterError):
        load_parameters(path)


def test_missing_source_fails(tmp_path):
    def mutate(raw):
        raw["parameters"]["C4-3"]["source"] = ""

    with pytest.raises(ParameterError):
        load_parameters(_write_variant(tmp_path, mutate))


def test_duplicate_key_fails(tmp_path):
    text = DEFAULT_PARAMETER_PATH.read_text(encoding="utf-8")
    text = text.replace('"C4-3":', '"C4-2":', 1)
    out = tmp_path / "dup.json"
    out.write_text(text, encoding="utf-8")
    with pytest.raises(ParameterError):
        load_parameters(out)


# ---- fragments -----------------------------------------------------------
@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CCO", {"C4-3": 1, "C4-2": 1, "O2-1": 1}),
        ("CC#N", {"C4-3": 1, "CN": 1}),
        ("CC(C)=O", {"C4-3": 2, "C3-0": 1, "O1-0": 1}),
        ("COC", {"C4-3": 2, "O2-0": 1}),
        ("c1ccccc1", {"C3-1": 6, "R6a": 1}),
        ("N#Cc1ccccn1", {"C3-1": 4, "C3-0": 1, "N2-0": 1, "CN": 1, "R6a": 1}),
        ("CCS", {"C4-3": 1, "C4-2": 1, "S2-1": 1}),
        ("CCCl", {"C4-3": 1, "C4-2": 1, "Cl1-0": 1}),
    ],
)
def test_fragment_counts(smiles, expected):
    result = analyze_fragments(Chem.MolFromSmiles(smiles))
    assert result.unsupported == []
    assert result.counts == expected


def test_unsupported_atoms():
    assert analyze_fragments(Chem.MolFromSmiles("CC#C")).unsupported
    assert analyze_fragments(Chem.MolFromSmiles("c1ccccc1[N+](=O)[O-]")).unsupported
    assert analyze_fragments(Chem.MolFromSmiles("CBr")).unsupported


def test_fused_or_mixed_ring_not_guessed():
    assert analyze_fragments(Chem.MolFromSmiles("c1ccc2CCCc2c1")).unsupported


# ---- model ---------------------------------------------------------------
def test_paper_worked_example_2_cyanopyridine():
    res = predict_flash_point("N#Cc1ccccn1", output_unit="K")
    assert res.quality["status"] == "COMPUTED"
    assert res.prediction["value"] == pytest.approx(356.5, abs=0.1)


def test_ethanol_matches_independent_formula(model):
    res = predict_flash_point("CCO", output_unit="K")
    expected = _manual_tf({"C4-3": 1, "C4-2": 1, "O2-1": 1}, model)
    assert res.prediction["value"] == pytest.approx(expected, rel=1e-12)


def test_contribution_fractions_sum_to_one():
    res = predict_flash_point("CCO")
    total = sum(c["fraction_of_f"] for c in res.contributions.values())
    assert total == pytest.approx(1.0)


def test_ring_r6_rejected_by_default():
    res = predict_flash_point("C1CCCCC1")
    assert res.quality["status"] == "UNDEFINED_PARAMETER_TERM"
    assert res.prediction["value"] is None


def test_signed_power_is_opt_in_and_warned():
    res = predict_flash_point(
        "C1CCCCC1", negative_term_policy=NegativeTermPolicy.SIGNED_POWER
    )
    assert res.quality["status"] == "COMPUTED"
    assert any("UNVERIFIED" in w for w in res.quality["warnings"])


def test_low_occurrence_warning():
    res = predict_flash_point("CCCS")  # S2-1 n_occ=8 -> no warning
    assert not any("S2-1" in w for w in res.quality["warnings"])
    res = predict_flash_point("CSC")  # S2-0 n_occ=4 -> warning
    assert any("S2-0" in w for w in res.quality["warnings"])


# ---- invalid input / units ----------------------------------------------
@pytest.mark.parametrize("bad", ["", "   ", "not_a_smiles", "C(C"])
def test_invalid_input(bad):
    res = predict_flash_point(bad)
    assert res.quality["status"] == "INVALID_INPUT"
    assert res.prediction["value"] is None


def test_mixture_rejected():
    res = predict_flash_point("CCO.CC")
    assert res.quality["status"] == "UNSUPPORTED_STRUCTURE"


def test_unit_conversion():
    assert convert_temperature(300.0, "K") == 300.0
    assert convert_temperature(273.15, "degC") == pytest.approx(0.0)
    assert convert_temperature(273.15, "degF") == pytest.approx(32.0)
    with pytest.raises(ValueError):
        convert_temperature(300.0, "R")


# ---- applicability domain -------------------------------------------------
def test_ad_in_domain_partial_for_ethanol():
    ad = predict_flash_point("CCO").applicability_domain
    assert ad["status"] == "IN_DOMAIN_PARTIAL"
    assert ad["details"]["not_evaluated_criteria"]


def test_ad_out_for_low_mw_and_amide():
    assert (
        predict_flash_point("CO").applicability_domain["status"]
        == "OUT_OF_DOMAIN"
    )
    res = predict_flash_point("CC(=O)N")
    assert res.applicability_domain["status"] == "OUT_OF_DOMAIN"
    assert res.quality["status"] == "COMPUTED"
    assert any("applicability" in w for w in res.quality["warnings"])
