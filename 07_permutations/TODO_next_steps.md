# Next Steps — AF Cutoff & Vg Distribution Analysis
> PI feedback debrief · 2026-02-25
> Notebook: `06_AF_cutoffs/01_genetic_variance_distribution copy.ipynb`

---

## Task 1 — Empirical p-value via permutation distribution

### Hypothesis
Observed Vg differences between ClinGen HI and background genes reflect genuine purifying selection, not sampling noise. This can be formalised through an empirical p-value derived from a permutation null.

### Execution

**Formula**

```python
p_empirical = (np.sum(perm_distribution >= obs_value) + 1) / len(perm_distribution)
```

- `obs_value` — observed statistic (e.g. median Vg, MW test statistic) for the real gene set
- `perm_distribution` — vector of the same statistic computed over permuted assignments;  
  `n_perm = len(perm_distribution)` (no fixed size; use however many permutation replicates are generated)
- `+1` continuity correction to avoid p = 0

**Implementation steps**
- [ ] Implement `empirical_pval(obs, perm_vec)` utility in `utils/stats_utils.py`
- [ ] Apply across all four source comparisons: `background`, `background_null`, `clingen`, `clingen_null`
- [ ] Apply per spatial region separately (promoter core, proximal/distal upstream, proximal/distal downstream)
- [ ] Report p-values alongside current Mann-Whitney U results for cross-validation

---

## Task 2 — Standardise Vg log-scale display via config

### Hypothesis / Motivation
`vg_predicted` spans multiple orders of magnitude; linear scale misleads visual comparisons. Centralising this as a config flag enforces consistency across all notebooks.

### Execution

**In `config.py`**
- [ ] Add a `VG_LOG_SCALE: bool = True` flag (and optionally `VG_MIN_CLIP: float = 1e-12` for zero-value handling)

**In each plot building Vg axes**
- [ ] Import `VG_LOG_SCALE` and conditionally call `ax.set_yscale('log')` / `ax.set_xscale('log')`
- [ ] Applies to: `vg_predicted`, `vg_common`, `vg_rare`, all `vg_region_*` columns

---

## Task 3 — Vg-value permutation null (1000 iterations, gnomAD)

This is the main new analysis block. Currently, permutation is done by shuffling AF vectors and re-assigning those AFs to synthetic variants. The new approach goes further: **permute the Vg contribution values themselves** to create a richer null distribution.

---

### 3a — Build the joint effect-size vector and permutation distribution

#### Hypothesis
If ClinGen HI genes are truly depleted of regulatory variation, their observed Vg should sit in the tail of a permutation distribution built from the combined gnomAD + synthetic pool.

#### Execution
- [ ] Construct a joint vector of per-variant Vg contributions from `background` gnomAD + `clingen` gnomAD (observed AFs)
- [ ] Run **1000 permutation iterations**: randomly shuffle gene labels (or resample variant-to-gene assignment) and recompute per-gene Vg
- [ ] Store the resulting null distribution for each gene set and spatial window
- [ ] Compute `p_empirical` (Task 1 formula) for each gene

---

### 3b — MW test distribution and scatter validation

#### Hypothesis
Permuted-label Mann-Whitney U statistics should cluster around the identity line vs. the real MW statistic. Real ClinGen HI–vs–background MW statistics should appear as outliers, confirming genuine depletion.

#### Execution
- [ ] For each permutation replicate, compute MW U-statistic between permuted HI and permuted background
- [ ] Build scatter plot: **x = permuted MW stat, y = observed MW stat** for each spatial region  
  _Expected: dots should cluster below the identity diagonal (observed stat > permuted stats)_
- [ ] Add `within_observed_variance` reference band (±1 SD of permuted distribution)

---

### 3c — QQ plot of empirical p-values

#### Hypothesis
Under the null (random label assignment), p-values should follow Uniform(0,1). Deviation toward small p-values in ClinGen HI promoter windows indicates non-random depletion.

#### Execution
- [ ] Plot QQ: observed –log₁₀(p) vs. expected –log₁₀(p) for each gene set × spatial region
- [ ] Highlight ClinGen HI promoter core as primary signal track
- [ ] Add genomic inflation factor λ annotation

---

### 3d — Stacked cumulative Vg plot by AF bin

#### Hypothesis
Rare variants (ultra-low AF) disproportionately drive Vg in constrained genes. A cumulative plot along the AF axis should show a steeper initial rise for ClinGen HI genes.

#### Execution
- [ ] Per gene, sort variants by AF (ascending)
- [ ] Compute cumulative sum of `vg_contribution` as AF increases (x-axis = AF)
- [ ] Stack / overlay per-gene cumulative curves, grouped by source (`background`, `clingen`)
- [ ] Optionally normalise by gene-total Vg to compare shape (not magnitude)

---

### 3e — Spatial Vg partition: single-gene deep-dive + log-mean summary

#### Hypothesis
HI genes are more constrained specifically in the **promoter core** (TSS ± 200 bp). Partitioning Vg by spatial region within individual genes should reveal significantly lower promoter-core Vg in HI genes relative to their own distal-upstream window.

#### Execution
- [ ] For a representative gene (e.g. an established HI gene), produce a **bar / waterfall chart** partitioning total Vg into the five spatial regions
- [ ] Compute **log-mean Vg** per spatial region across all genes in each source  
  `log_mean = exp(mean(log(vg_region + ε)))`
- [ ] Compare log-means: `clingen` vs `background` per region — expected: largest effect in promoter core
- [ ] Statistical test: paired Wilcoxon signed-rank on log Vg per region (ClinGen HI vs background, matched by genomic-length quantile)

---

## Notes on analysis structure

> "New structure must always follow: **hypothesis → execution (biological and/or technical)**"

Each analysis block above is written in that form. When adding new cells to the notebook, start with a markdown cell stating:
1. The biological question being answered
2. The expected direction of effect
3. The technical approach used to test it

---

## Column reference (AGENTS.md)

| Use case | Real sources | Null sources |
|---|---|---|
| Total Vg | `vg_predicted` | `vg_predicted_perm` |
| Common AF Vg | `vg_common` | `vg_common_perm` |
| Rare AF Vg | `vg_rare` | `vg_rare_perm` |
| Spatial regions | `vg_{window}` | `vg_{window}_perm` |

Do **not** use `vg_predicted_perm_sanity` as a biological null (QC use only).
