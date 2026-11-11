"""Reference implementation of Mathieu & Alaime (2014) flash point model.

Mathieu, D.; Alaime, T. J. Hazard. Mater. 2014, 267, 169-174.
DOI: 10.1016/j.jhazmat.2013.12.047

Model (Eq.1-2, p.170)::

    T_F = f ** alpha,   f = sum_i n_i * t_i ** (1 / alpha),   alpha = 2/5

Anything not specified in the paper is NOT guessed: such structures return
a non-computed status with the reason attached. Pure-substance screening
aid only; not a substitute for standard flash point testing.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors

MODEL_NAME = "Mathieu_Alaime_2014"
DEFAULT_PARAMETER_PATH = Path(__file__).with_name("source_parameters.json")

RING_KEYS = frozenset({"R5", "R6", "R6a"})
CYANO_KEY = "CN"
EXPECTED_KEYS = frozenset(
    {
        "C3-0", "C3-1", "C3-2", "C4-0", "C4-1", "C4-2", "C4-3",
        "N2-0", "N3-0", "N3-1", "N3-2", "O1-0", "O2-0", "O2-1",
        "S2-0", "S2-1", "Cl1-0", CYANO_KEY, "R5", "R6", "R6a",
    }
)
ATOMIC_KEYS = EXPECTED_KEYS - RING_KEYS - {CYANO_KEY}

# Applicability domain, Section 2 (p.171): molecular weight bounds in g/mol.
AD_MW_MIN = 45.0
AD_MW_MAX = 370.0
LOW_OCCURRENCE_THRESHOLD = 5  # Section 3: Nocc < 5 -> possibly ill-defined

_CYANO = Chem.MolFromSmarts("[CX2]#[NX1]")
_AMIDE = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
_PEROXIDE = Chem.MolFromSmarts("[OX2][OX2]")
_NITRO = Chem.MolFromSmarts("[$([NX3](=O)=O),$([N+](=O)[O-])]")


class ParameterError(ValueError):
    """Raised when the parameter file is incomplete or inconsistent."""


class UndefinedTermError(ValueError):
    """Raised when t_i ** (1/alpha) is undefined for a real number."""


class NegativeTermPolicy(str, Enum):
    """Handling of negative t_i (R5, R6) in t_i ** (1/alpha).

    The paper does not state how a negative t_i is raised to 5/2.
    REJECT (default): molecules with R5/R6 are not computed.
    SIGNED_POWER: -(|t|) ** (1/alpha). UNVERIFIED assumption.
    """

    REJECT = "reject"
    SIGNED_POWER = "signed_power"


class AdStatus(str, Enum):
    OUT_OF_DOMAIN = "OUT_OF_DOMAIN"
    IN_DOMAIN_PARTIAL = "IN_DOMAIN_PARTIAL"
    NOT_EVALUATED = "NOT_EVALUATED"


# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Parameter:
    key: str
    t_k: float
    sd_k: float
    n_occ: int
    source: str


@dataclass(frozen=True)
class ModelParameters:
    model_version: str
    alpha: float
    alpha_source: str
    equation_source: str
    parameters: dict[str, Parameter]


def _reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ParameterError(f"duplicate key in parameter file: {key}")
        result[key] = value
    return result


def load_parameters(path: Optional[Path] = None) -> ModelParameters:
    """Load and validate the parameter file (fail-fast on any defect)."""
    path = Path(path) if path else DEFAULT_PARAMETER_PATH
    with path.open(encoding="utf-8") as fh:
        raw = json.load(fh, object_pairs_hook=_reject_duplicates)

    if raw.get("model") != MODEL_NAME:
        raise ParameterError("model name mismatch")
    if raw.get("temperature_unit") != "K":
        raise ParameterError("temperature_unit must be 'K'")

    alpha = raw["alpha"]["value"]
    if not (isinstance(alpha, (int, float)) and 0.0 < alpha < 1.0):
        raise ParameterError("alpha must satisfy 0 < alpha < 1")

    entries = raw.get("parameters", {})
    keys = set(entries)
    if keys != EXPECTED_KEYS:
        raise ParameterError(
            f"missing={sorted(EXPECTED_KEYS - keys)} "
            f"unexpected={sorted(keys - EXPECTED_KEYS)}"
        )

    params: dict[str, Parameter] = {}
    for key, entry in entries.items():
        t_k, sd_k = entry.get("t_k"), entry.get("sd_k")
        n_occ, source = entry.get("n_occ"), entry.get("source")
        if not all(
            isinstance(v, (int, float)) and math.isfinite(v)
            for v in (t_k, sd_k)
        ):
            raise ParameterError(f"{key}: t_k/sd_k must be finite numbers")
        if not isinstance(n_occ, int) or n_occ < 1:
            raise ParameterError(f"{key}: invalid n_occ")
        if not isinstance(source, str) or not source.strip():
            raise ParameterError(f"{key}: source is required")
        params[key] = Parameter(key, float(t_k), float(sd_k), n_occ, source)

    return ModelParameters(
        model_version=str(raw["model_version"]),
        alpha=float(alpha),
        alpha_source=str(raw["alpha"]["source"]),
        equation_source=str(raw["equation_source"]),
        parameters=params,
    )


# --------------------------------------------------------------------------
# Fragment recognition
# --------------------------------------------------------------------------
@dataclass
class FragmentAnalysis:
    counts: dict[str, int]
    unsupported: list[str] = field(default_factory=list)


def analyze_fragments(mol: Chem.Mol) -> FragmentAnalysis:
    """Count geometrical fragments X<coordination>-<nH> and ring terms.

    Coordination number is the total number of neighbours including
    hydrogens (methyl C = C4-3, nitrile C = C2-0). Bond orders are ignored.
    The cyano group (C and N atoms) is counted once as ``CN``.
    Atoms/rings outside Table 1 are reported in ``unsupported``.
    """
    mol = Chem.RemoveHs(mol)
    counts: Counter[str] = Counter()
    unsupported: list[str] = []

    cyano_atoms: set[int] = set()
    for c_idx, n_idx in mol.GetSubstructMatches(_CYANO):
        cyano_atoms.update((c_idx, n_idx))
        counts[CYANO_KEY] += 1

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in cyano_atoms:
            continue
        label = (
            f"{atom.GetSymbol()}{atom.GetTotalDegree()}"
            f"-{atom.GetTotalNumHs()}"
        )
        if atom.GetFormalCharge() != 0 or atom.GetNumRadicalElectrons():
            unsupported.append(f"atom {idx} ({label}): charged or radical")
        elif label in ATOMIC_KEYS:
            counts[label] += 1
        else:
            unsupported.append(f"atom {idx} ({label}): not in Table 1")

    for ring in mol.GetRingInfo().AtomRings():
        size = len(ring)
        flags = [mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring]
        if all(flags) and size == 6:
            counts["R6a"] += 1
        elif not any(flags) and size in (5, 6):
            counts[f"R{size}"] += 1
        else:
            unsupported.append(
                f"ring size={size} aromatic_atoms={sum(flags)}/{size}: "
                "ring correction not defined in paper"
            )

    return FragmentAnalysis(dict(counts), unsupported)


# --------------------------------------------------------------------------
# Model
# --------------------------------------------------------------------------
def _term(t_k: float, alpha: float, policy: NegativeTermPolicy) -> float:
    """Return t ** (1/alpha) following the selected negative-term policy."""
    exponent = 1.0 / alpha
    if t_k > 0.0:
        return t_k**exponent
    if policy is NegativeTermPolicy.SIGNED_POWER:
        return -((-t_k) ** exponent)
    raise UndefinedTermError(
        f"t={t_k} K is negative: t**(1/alpha) is undefined in the paper"
    )


def compute_flash_point_k(
    counts: dict[str, int],
    model: ModelParameters,
    policy: NegativeTermPolicy = NegativeTermPolicy.REJECT,
) -> tuple[float, dict[str, dict[str, float]]]:
    """Evaluate Eq.(1)-(2). Returns (T_F in K, per-fragment contributions)."""
    contributions: dict[str, dict[str, float]] = {}
    f_total = 0.0
    for key, n in counts.items():
        t_k = model.parameters[key].t_k
        term = n * _term(t_k, model.alpha, policy)
        contributions[key] = {"count": n, "t_K": t_k, "n_times_t_pow": term}
        f_total += term

    if not math.isfinite(f_total) or f_total <= 0.0:
        raise UndefinedTermError(f"f={f_total} is not a positive finite value")
    t_f = f_total**model.alpha
    if not math.isfinite(t_f):
        raise UndefinedTermError("T_F is not finite")

    for item in contributions.values():
        item["fraction_of_f"] = item["n_times_t_pow"] / f_total
    return t_f, contributions


def convert_temperature(value_k: float, unit: str) -> float:
    """Convert from Kelvin. Supported units: K, degC, degF."""
    if unit == "K":
        return value_k
    if unit == "degC":
        return value_k - 273.15
    if unit == "degF":
        return (value_k - 273.15) * 9.0 / 5.0 + 32.0
    raise ValueError(f"unsupported unit: {unit!r} (use K, degC, degF)")


# --------------------------------------------------------------------------
# Applicability domain (Section 2, p.171)
# --------------------------------------------------------------------------
def evaluate_applicability(mol: Chem.Mol) -> dict[str, Any]:
    """Evaluate the implementable parts of the paper's AD definition.

    Paper criteria: (a) MW 45-370 g/mol; (b) no nitro, peroxide, amide;
    (c) no rings/groups absent from S0; (d) no strong intramolecular
    H-bond or dipole-dipole interaction between geminal groups.
    (d) is not operationally defined -> reported as not evaluated.
    Amide SMARTS (includes urea/carbamate) and the ring rule in (c) are
    this implementation's interpretation, flagged in ``details``.
    """
    mw = Descriptors.MolWt(mol)
    violations: list[str] = []
    if not AD_MW_MIN <= mw <= AD_MW_MAX:
        violations.append(f"molecular_weight={mw:.2f} outside 45-370")
    for name, patt in (
        ("nitro", _NITRO), ("peroxide", _PEROXIDE), ("amide", _AMIDE)
    ):
        if mol.HasSubstructMatch(patt):
            violations.append(f"contains {name}")

    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        aromatic = all(a.GetIsAromatic() for a in atoms)
        hetero = any(a.GetSymbol() != "C" for a in atoms)
        if len(ring) not in (5, 6):
            violations.append(f"ring size {len(ring)} (interpretation)")
        elif not aromatic and hetero:
            violations.append(
                "saturated heterocycle (interpretation of 'not in S0')"
            )

    status = AdStatus.OUT_OF_DOMAIN if violations else AdStatus.IN_DOMAIN_PARTIAL
    return {
        "status": status.value,
        "method": "Mathieu_Alaime_2014 Section 2 (partial implementation)",
        "details": {
            "molecular_weight": mw,
            "violations": violations,
            "not_evaluated_criteria": [
                "strong intramolecular hydrogen bonding",
                "dipole-dipole interaction between geminal groups",
            ],
        },
    }


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------
@dataclass
class FlashPointResult:
    model: str
    model_version: str
    input: dict[str, Any]
    prediction: dict[str, Any]
    fragments: dict[str, int]
    contributions: dict[str, Any]
    applicability_domain: dict[str, Any]
    quality: dict[str, Any]
    provenance: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def predict_flash_point(
    smiles: str,
    output_unit: str = "degC",
    model: Optional[ModelParameters] = None,
    negative_term_policy: NegativeTermPolicy = NegativeTermPolicy.REJECT,
) -> FlashPointResult:
    """Predict the flash point of a pure organic compound from SMILES."""
    convert_temperature(0.0, output_unit)  # validate the unit early
    model = model or load_parameters()
    warnings: list[str] = []

    def build(
        status: str,
        value: Optional[float] = None,
        value_k: Optional[float] = None,
        counts: Optional[dict[str, int]] = None,
        contribs: Optional[dict[str, Any]] = None,
        ad: Optional[dict[str, Any]] = None,
        reasons: Optional[list[str]] = None,
    ) -> FlashPointResult:
        return FlashPointResult(
            model=MODEL_NAME,
            model_version=model.model_version,
            input={"smiles": smiles},
            prediction={
                "value": value,
                "unit": output_unit,
                "value_K": value_k,
            },
            fragments=counts or {},
            contributions=contribs or {},
            applicability_domain=ad
            or {
                "status": AdStatus.NOT_EVALUATED.value,
                "method": None,
                "details": None,
            },
            quality={
                "status": status,
                "reasons": reasons or [],
                "warnings": warnings,
            },
            provenance={
                "equation_source": model.equation_source,
                "parameter_source": "Table 1, p.170",
                "alpha_source": model.alpha_source,
                "negative_term_policy": negative_term_policy.value,
            },
        )

    if not isinstance(smiles, str) or not smiles.strip():
        return build("INVALID_INPUT", reasons=["empty or non-string SMILES"])
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return build("INVALID_INPUT", reasons=["RDKit failed to parse SMILES"])

    ad = evaluate_applicability(mol)
    if len(Chem.GetMolFrags(mol)) > 1:
        return build(
            "UNSUPPORTED_STRUCTURE",
            ad=ad,
            reasons=["multiple fragments (mixture/salt) not supported"],
        )

    analysis = analyze_fragments(mol)
    if analysis.unsupported:
        return build(
            "UNSUPPORTED_STRUCTURE",
            counts=analysis.counts,
            ad=ad,
            reasons=analysis.unsupported,
        )

    try:
        t_f_k, contribs = compute_flash_point_k(
            analysis.counts, model, negative_term_policy
        )
    except UndefinedTermError as exc:
        return build(
            "UNDEFINED_PARAMETER_TERM",
            counts=analysis.counts,
            ad=ad,
            reasons=[str(exc)],
        )

    if any(model.parameters[k].t_k < 0 for k in analysis.counts):
        warnings.append(
            "negative ring term evaluated with UNVERIFIED signed-power rule"
        )
    for key in analysis.counts:
        if model.parameters[key].n_occ < LOW_OCCURRENCE_THRESHOLD:
            warnings.append(
                f"{key}: n_occ={model.parameters[key].n_occ} (<5) in "
                "training set; parameter may be numerically ill-defined"
            )
    if ad["status"] != AdStatus.IN_DOMAIN_PARTIAL.value:
        warnings.append("outside applicability domain: low reliability")

    return build(
        "COMPUTED",
        value=convert_temperature(t_f_k, output_unit),
        value_k=t_f_k,
        counts=analysis.counts,
        contribs=contribs,
        ad=ad,
        reasons=[],
    )


if __name__ == "__main__":
    print(json.dumps(predict_flash_point("N#Cc1ccccn1").to_dict(), indent=2))
