# Repository Map

This repository contains both the thesis-facing analysis layer and earlier exploratory work. The cleanup policy is to preserve exploratory work for traceability while making the public entrypoints clear.

## Thesis-Facing

| Path | Role |
| --- | --- |
| `10_THESIS/` | Main thesis notebooks, thesis figure generation, and thesis-specific plotting utilities. This is the preferred entrypoint for the final analysis narrative. |
| `10_THESIS/5.manuscript_figure_polish.py` | Script used to polish selected manuscript figures after notebook exploration. |
| `10_THESIS/thesis_style.py` | Shared thesis plotting style and figure-export helpers. |
| `docs/research-summary.md` | Public summary of thesis methods, results, caveats, and figure provenance. |
| `docs/figures/` | Curated GitHub-renderable figures derived from the thesis `Figures` directory. |
| `utils/thesis_model_validation.py` | Thesis-era validation helpers, when present in the worktree. |
| `utils/eqtl_susie_alpha_genome.py` | eQTL/SuSiE helper code, when present in the worktree. |

## Active or Recent Exploratory Work

| Path | Role |
| --- | --- |
| `09_3103_biweekly/` | eQTL/SuSiE validation, singleton/haplotype exploration, and late-stage comparison analyses. Some files here are generated outputs and should not be committed without curation. |
| `07_permutations/` | Permutation, null-model, and allele-frequency sanity analyses. |
| `08_figure_outputs_to_start_somewhere/` | Earlier figure-development notebooks that informed the final thesis figures. |

## Earlier Analysis Layers

| Path | Role |
| --- | --- |
| `01_qc/` and `01_initial_qc_biweekly_0701/` | Early quality-control notebooks for variant- and gene-level outputs. |
| `02_vg_perm/` | Early genetic-variance and permutation development. |
| `03_variance_coupled/` | Earlier coupled/uncoupled variance analyses and ESHG-oriented statistics. |
| `04_ESHG_abstract/` | Analyses prepared for ESHG abstract development. |
| `04_cis_regulatory_constraint_pipeline/` | Earlier cis-regulatory constraint exploration and candidate tables. |
| `05_Hyper_Buffered_Distal_Regulation_Of_Haploinsufficient_Genes/` | Earlier distal/proximal architecture exploration; the directory name may contain legacy ligature characters in the filesystem. |
| `06_AF_cutoffs/` | Earlier allele-frequency cutoff analyses. |

## Data and Generated Outputs

| Path | Role |
| --- | --- |
| `experiments_data/` | Local analysis data, small reference inputs, and generated parquet outputs. Treat large parquet files as local or derived artifacts unless explicitly curated. |
| `figures/` | Generated local figure exports. The root figure directory is intentionally ignored by Git; curated public figures live under `docs/figures/`. |
| `results/` | Local derived result tables and diagnostics. Commit only curated, thesis-facing summaries. |

## Public Documentation Boundary

The public documentation should describe the final thesis analysis and its limitations. It should not promote every exploratory notebook as a stable workflow. When an exploratory folder contains useful history, describe it here rather than moving or deleting it during documentation cleanup.

## Commit Boundary

The `cleanup-audit` branch should stage only curated public-facing files:

- root `README.md`
- `docs/*.md`
- curated `docs/figures/*.png`
- `.gitignore`
- `CITATION.cff`

Do not stage local kernels, caches, temporary AlphaGenome outputs, bulk parquet files, duplicate draft figures, or notebook execution churn unless those files are intentionally selected for a later reproducibility release.
