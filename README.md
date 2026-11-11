# Mathieu & Alaime (2014) 引火点推算モデル

Mathieu, D.; Alaime, T. *Insight into the contribution of individual functional groups to the flash point of organic compounds*. J. Hazard. Mater. 2014, 267, 169-174. DOI: [10.1016/j.jhazmat.2013.12.047](https://doi.org/10.1016/j.jhazmat.2013.12.047)

純物質の引火点を SMILES から推算する Reference Implementation。論文に記載のない仕様は推測で実装せず、未確認事項は未実装として明示する。

## モデル

```
T_F = f^alpha
f   = sum_i n_i * t_i^(1/alpha)      alpha = 2/5
```

- `n_i`: フラグメント i の個数、`t_i`: Table 1 のパラメータ（単位 K）。
- 21パラメータ: 原子17種、シアノ基 CN、環補正 R5 / R6 / R6a。
- 出典: Eq.(1)-(2) および Table 1（p.170）、alpha は p.171。

## 構成

| ファイル | 内容 |
| --- | --- |
| `flash_point_ma2014.py` | パラメータ読込、フラグメント認識、推算式、適用領域判定、公開API |
| `source_parameters.json` | 21パラメータと alpha（出典付き） |
| `test_flash_point_ma2014.py` | pytest（28件） |

## 要件

- Python 3.10 以上
- RDKit
- pytest（テスト時）

```bash
pip install rdkit pytest
```

## 使い方

```python
from flash_point_ma2014 import predict_flash_point

result = predict_flash_point("N#Cc1ccccn1", output_unit="K")
print(result.prediction["value"])   # 356.5 付近（論文 Eq.3 の例）
print(result.to_dict())             # JSON 化可能な辞書
```

- `output_unit`: `"K"`、`"degC"`（既定）、`"degF"`。内部計算は常に K で、変換は出力時のみ。
- 丸めは行わない。
- `model=load_parameters(path)` で別のパラメータファイルを指定できる。

### 出力

| キー | 内容 |
| --- | --- |
| `prediction` | `value`（出力単位）、`unit`、`value_K` |
| `fragments` | フラグメント個数 |
| `contributions` | 各項の `count`、`t_K`、`n_times_t_pow`、`fraction_of_f` |
| `applicability_domain` | `status`、`method`、`details` |
| `quality` | `status`、`reasons`、`warnings` |
| `provenance` | 式・パラメータ・alpha の出典、負項の扱い |

### `quality.status`

| 値 | 意味 |
| --- | --- |
| `COMPUTED` | 計算済み |
| `INVALID_INPUT` | SMILES が空、または RDKit で解釈不能 |
| `UNSUPPORTED_STRUCTURE` | Table 1 外の原子、電荷・ラジカル、アルキン、未定義の環、混合物・塩 |
| `UNDEFINED_PARAMETER_TERM` | 負の t_i（R5/R6）を含み、既定設定では計算不能 |

非 `COMPUTED` の場合、`prediction.value` は `null`。

### `applicability_domain.status`

| 値 | 意味 |
| --- | --- |
| `OUT_OF_DOMAIN` | 判定可能な基準のいずれかに違反 |
| `IN_DOMAIN_PARTIAL` | 実装済み基準は全て満足。ただし未評価基準が残る |
| `NOT_EVALUATED` | 判定していない（計算不能時の初期値） |

`IN_DOMAIN_PARTIAL` は完全な域内を意味しない。

## フラグメント定義

- 原子は `X<配位数>-<H数>`。配位数は H を含む総結合数（メチル C = `C4-3`、ケトン C = `C3-0`）。結合次数は無視する。
- シアノ基の C と N は `CN` として1回だけ加算する（`C2-0`、`N1-0` は単独では加算しない）。
- 環補正: 完全芳香族の6員環は `R6a`、非芳香族の5員環・6員環は `R5` / `R6`。

## 適用領域（論文 Section 2, p.171）

実装済み:
- 分子量 45-370 g/mol
- ニトロ基、過酸化物、アミドの除外

実装者の解釈（`details` に明示）:
- アミド SMARTS は尿素・カルバメートを含む。
- 5員環・6員環以外の環、飽和ヘテロ環は域外とする。

未評価（判定基準が論文に明記されていない）:
- 強い分子内水素結合
- ジェミナル基間の双極子相互作用

## 既知の制限と未確認事項

1. **負の t_i（R5 = -113.0 K、R6 = -81.9 K）**: t_i^(5/2) の実数値は未定義で、論文に扱いの記載がない。
   - 既定 `NegativeTermPolicy.REJECT` では環系分子を計算しない。
   - `NegativeTermPolicy.SIGNED_POWER`（-|t|^(5/2)）は未検証の仮定であり、使用時は警告を出す。
2. **縮合環・部分芳香族環**（インダン、テトラリン等）と芳香族5員環: 環補正の割当規則が論文にないため未サポート。
3. **論文データでの検証（MA14-009）**: 補足資料（1457化合物の実測値・計算値）が未入手のため未実施。検証済みなのは論文本文の計算例（2-cyanopyridine）のみ。
4. 係数の再学習・混合物・溶液は対象外。
5. 論文の注記: Nocc < 5 のパラメータ（N3-0、S2-0、CN、R5）は数値的に不安定な可能性があり、使用時に警告を出す。
6. 論文の記載差: S1+S2 の AAE が本文 13.4 K、Table 2 で 13.3 K。実装には影響しない。

## 論文報告の性能（参考値）

| テストセット | 化合物数 | AAE (K) |
| --- | --- | --- |
| 訓練 S0 | 287 | 6.7 |
| 外部 S1（適用領域内） | 488 | 8.6 |
| 外部 S2 | 682 | 16.6 |

これらは本実装での再現結果ではない。

## テスト

```bash
python -m pytest -q
```

内容: パラメータ検証、フラグメント個数、論文の計算例（356.5 K）、独立再計算との一致、環の負項、AD、入力不備、単位変換。

## 運用上の注意

- 推算値は実測値ではない。危険物該当性の判定、設備設計、運転条件の決定に推算値のみを使用しない。
- 安全性評価には、標準試験法による実測を優先する。
- 適用領域外・未評価の結果を、適用領域内と同等に扱わない。
- 本モデルは純物質用で、混合物・溶液・不純物を含む試料には適用できない。
- 密閉式・開放式など測定法の異なる実測値は、区別して比較する。
