# AGENTS.md — Downstream Analysis Rules for In-Silico Genetic Variance

This file defines **strict downstream column usage** for gene-level parquet outputs.
Follow these rules exactly to avoid mixing sanity/permutation artifacts with biological null models.

---

## 1) Dataset Types and Intended Use

You will typically analyze four sources:

| source | Meaning | AF basis |
|---|---|---|
| `background` | Real gnomAD background genes | observed AF |
| `clingen` | Real gnomAD ClinGen HI genes | observed AF |
| `background_null` | Synthetic/null background genes | permuted AF (null design) |
| `clingen_null` | Synthetic/null ClinGen genes | permuted AF (null design) |

---

## 2) CRITICAL Rule: `vg_predicted_perm_sanity` is **NOT** a downstream null

`vg_predicted_perm_sanity` (and any sanity-derived real-file columns) are generated only as an internal QC/sanity check of within-gene AF shuffling.

- Use it to verify pipeline behavior/QC.
- **Do NOT** use it as biological null in depletion/enrichment inference.
- For biological obs-vs-null comparisons, use `*_NULL/*_Synth_*` datasets.

---

## 3) Which Vg Columns to Use

## Real datasets (`background`, `clingen`)
Use:
- `vg_predicted`
- `vg_common`, `vg_rare`
- `vg_distal_upstream`, `vg_proximal_upstream`, `vg_promoter_core`, `vg_down_proximal`, `vg_down_distal`

Do **not** use as null:
- `vg_predicted_perm_sanity`
- any `_perm_sanity`-derived depletion/weighted columns (if present)

## Synthetic/null datasets (`background_null`, `clingen_null`)
Use:
- `vg_predicted_perm` (null expectation for total Vg)
- `vg_common_perm`, `vg_rare_perm`
- `vg_distal_upstream_perm`, `vg_proximal_upstream_perm`, `vg_promoter_core_perm`, `vg_down_proximal_perm`, `vg_down_distal_perm`

Do **not** use:
- non-perm `vg_*` columns in null files for null inference (these are not the target null signal).

---

## 4) Correct Obs/Null Comparison Patterns

### Total Vg depletion
- Observed: real sources, column `vg_predicted`
- Null: null sources, column `vg_predicted_perm`

### AF-bin depletion
- Observed: `vg_common`, `vg_rare` from real
- Null: `vg_common_perm`, `vg_rare_perm` from null

### Spatial-window depletion
- Observed: `vg_{window}` from real
- Null: `vg_{window}_perm` from null

Join real/null by gene and mapped source:
- `background_null -> background`
- `clingen_null -> clingen`

---

## 5) Confidence Intervals: What They Mean and How to Use Them

CI columns are Monte Carlo AF-resampling uncertainty summaries, not the primary biological null comparison.

Available CI columns (when `--calc-ci` was used):
- `vg_predicted_CI_mean`, `vg_predicted_CI_p05`, `vg_predicted_CI_p95`
- `vg_common_CI_mean`, `vg_common_CI_p05`, `vg_common_CI_p95`
- `vg_rare_CI_mean`, `vg_rare_CI_p05`, `vg_rare_CI_p95`
- `vg_distal_upstream_CI_mean`, `vg_distal_upstream_CI_p05`, `vg_distal_upstream_CI_p95`
- `vg_proximal_upstream_CI_mean`, `vg_proximal_upstream_CI_p05`, `vg_proximal_upstream_CI_p95`
- `vg_promoter_core_CI_mean`, `vg_promoter_core_CI_p05`, `vg_promoter_core_CI_p95`
- `vg_down_proximal_CI_mean`, `vg_down_proximal_CI_p05`, `vg_down_proximal_CI_p95`
- `vg_down_distal_CI_mean`, `vg_down_distal_CI_p05`, `vg_down_distal_CI_p95`

Downstream usage guidance:
- Use CI columns for uncertainty/context around per-gene Vg.
- Do not replace real-vs-null dataset comparisons with CI-only inference.
- Keep CI interpretation separate from synthetic null depletion tests.

---

## 6) Frequently Misused Columns — Do NOT Use This Way

```python
# WRONG: using real-file sanity column as null
ratio = df_real["vg_predicted"] / df_real["vg_predicted_perm_sanity"]

# WRONG: using real-file sanity/depletion columns as biological null evidence
signal = df_real["depletion_promoter_core"]

# WRONG: comparing real vg_window to real perm_sanity window (if any appear)
ratio = df_real["vg_promoter_core"] / df_real["vg_promoter_core_perm_sanity"]