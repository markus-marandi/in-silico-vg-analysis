# Research Summary

This repository contains the downstream analyses for **Deep-Learning Prediction of Human Gene Expression Variability from Genetic Variation**. The thesis evaluates whether AlphaGenome-predicted variant effects can be aggregated into a gene-level measure of cis-regulatory genetic variance (`V_G`) that reflects dosage constraint in haploinsufficient genes.

## Scientific Motivation

Haploinsufficient genes are sensitive to reduced dosage: one functional copy is not enough to maintain normal physiology. Existing constraint metrics such as LOEUF, pLI, pTriplo, ncGERP, ncCADD, and Episcore capture important dimensions of coding, structural, conservation, and epigenomic constraint. The thesis addresses a complementary question: whether predicted cis-regulatory expression variance is depleted in dosage-sensitive genes.

The analysis uses the genetic-variance formulation:

```text
V_G = sum_i 2 p_i (1 - p_i) beta_i^2
```

where `p_i` is allele frequency and `beta_i` is the predicted expression effect size for variant `i`.

## Analysis Design

The study compares:

| Gene set | Count | Definition |
| --- | ---: | --- |
| ClinGen haploinsufficient | 316 | Highest-confidence ClinGen Dosage Sensitivity Map tier, after transcript/expression/promoter-overlap filtering |
| Background coding | 349 | Protein-coding genes from the gnomAD v4.1 constraint table, excluding high-confidence ClinGen haploinsufficient genes |

Variants were analyzed in +/-10 kb windows around MANE Select v1.4 transcription start sites on GRCh38. AlphaGenome used its 1 Mb sequence context during inference, but the tested variant set was proximal to the TSS. The primary expression signal was skeletal-muscle polyA+ RNA-seq, matching the GTEx skeletal-muscle validation context.

## Notebook-to-Thesis Map

| Thesis section | Main notebook(s) | Purpose |
| --- | --- | --- |
| 4.1 Model validation | `10_THESIS/4.1.model_validation_and_baseline_characterization.ipynb`, `10_THESIS/4.1b_eqtl_v10_concordance.ipynb` | Compare predicted `V_G` with eQTL-derived `V_G`, LOEUF, and other constraint metrics |
| 4.2 Total `V_G` comparison | `10_THESIS/4.2.Total_Vg_global_comparison.ipynb`, `10_THESIS/4.2.Vg_depletion_in_haploinsufficient_genes.ipynb` | Test total predicted regulatory variance in ClinGen haploinsufficient versus background genes |
| 4.3 Spatial architecture | `10_THESIS/4.3.Spatial_architecture_of_regulatory_constraint.ipynb` | Decompose `V_G` by TSS-relative windows and cluster regional architecture |
| 4.4 Allele-frequency architecture | `10_THESIS/4.4.Allele_frequency_architecture_of_constraint.ipynb` | Split predicted `V_G` by common/rare variants and allele-count bins |
| Supplementary analyses | `10_THESIS/4.Supplementary_robustness_and_subclass_figures.ipynb` | Robustness checks and subclass figures |

Exploratory and historical notebooks remain in other numbered folders. They are preserved for traceability but are not the clean public entrypoint for thesis reproduction.

## Key Results

| Question | Thesis result | Interpretation |
| --- | --- | --- |
| Per-variant eQTL concordance | `rho = 0.019`, `P = 0.784`, `n = 213` matched variants | No detectable aggregate per-variant directional concordance in the matched SuSiE set |
| Gene-level eQTL concordance | `rho = 0.205`, `P = 4.65e-6`, `n = 434` | Aggregation into gene-level `V_G` recovers a significant but modest empirical signal |
| LOEUF association | ClinGen HI `rho = 0.40`, `P = 2.3e-13`; background `rho = 0.31`, `P = 1.7e-8` | Genes under stronger coding constraint tend to have lower predicted regulatory variance |
| Total predicted `V_G` | ClinGen HI/background median ratio `0.68x`, Mann-Whitney `P = 1.12e-3` | ClinGen haploinsufficient genes carry lower predicted regulatory variance |
| Large-effect tail | HI retention `0.384` at `|Delta| >= 1.0`, uncorrected `P = 4.79e-2` | Hypothesis-generating; not a corrected definitive result |
| Spatial clustering | Distal-upstream cluster OR `1.88`, `P = 4.9e-3`; promoter-core cluster OR `0.54`, `P = 3.6e-3` | Cluster-level architecture differs between gene sets; individual-window contrasts are more tentative |

## Curated Figures

The public figure set in `docs/figures/` is derived from the thesis `Figures` directory:

| File | Thesis source | What it shows |
| --- | --- | --- |
| `gene-level-eqtl-vg-concordance.png` | `Figures/4_1b_gene_level_vg_eqtl_v10.pdf` | Predicted `V_G` versus GTEx eQTL-derived `V_G` |
| `total-vg-clingen-background.png` | `Figures/4_2_total_vg_global_comparison_4_2_total_vg_observed_only_22052026_1245.pdf` | Lower total predicted `V_G` in ClinGen haploinsufficient genes |
| `allele-frequency-vg-depletion.png` | `Figures/4_4_vg_depletion_violin_publication.pdf` | Common/rare partition of predicted `V_G` |
| `spatial-kmeans-architecture.png` | `Figures/4_3_spatial_kmeans_clusters_publication.pdf` | Regional `V_G` architecture clusters |
| `constraint-metric-correlations.png` | `Figures/4.1.4_constraint_metric_correlations.pdf` | Relationship to established constraint metrics |

## Limitations

Predicted `V_G` is a gene-level research screening metric, not clinical evidence. Per-variant accuracy was not supported in the matched eQTL validation set. The synthetic null is useful but incomplete: the ClinGen haploinsufficient/background difference persists in the null, implying that sequence composition and CpG-rich promoter architecture are not fully separated from selective constraint. All biological interpretations should therefore be treated as computational and hypothesis-generating until experimentally validated.

## Acknowledgements and Funding

This work was supervised by Philipp Rentzsch and Tuuli Lappalainen at the Lappalainen Lab, KTH / Science for Life Laboratory.

Computational resources were provided by NAISS at PDC Center for High Performance Computing, KTH Royal Institute of Technology, allocation `NAISS 2024/6-322`, and at NSC, Linköping University, allocation `Berzelius-2025-176`.

Funding support came from a Wallenberg Scholar award 2024 to T. Lappalainen, Knut and Alice Wallenberg Foundation grant `KAW 2023.0337`, and a Göran Gustavsson award 2023 to T. Lappalainen.
