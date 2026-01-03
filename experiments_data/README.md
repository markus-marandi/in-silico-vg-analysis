# Data Catalog: Variant and Gene Level Summaries

This directory contains the results of variant scoring and gene-level aggregation for three primary datasets: ClinGen (pathogenic/likely pathogenic), Background (gnomAD-derived), and matched NULL sets.

## Directory Structure

```text
experiments_data/
├── dataset3_ClinGen
│   ├── dataset3_clingen_gene_level_summary.csv
│   └── dataset3_clingen_variant_level_summary.parquet
├── dataset4_background
│   └── dataset4_background_gene_level_summary.csv
└── dataset5_NULL
    ├── dataset5_Background_NULL_variant_level_summary.parquet
    ├── dataset5_ClinGen_NULL_variant_level_summary.parquet
    ├── dataset5_Clingen_NULL_gene_level_summary.csv
    └── dataset5_background_NULL_gene_level_summary.csv

```

## The "Matched NULL" Sets

Dataset 5 (NULL) represents a "matched null" control set. For each gene (TSS ±10kb), we constructed a synthetic variant set that mirrors the mutational opportunities of the observed variants.

**The Goal:** To determine if differences in AlphaGenome LFC scores between ClinGen and Background gene sets persist when mutation-spectrum artifacts (e.g., CpG methylation, context bias) are held constant.

**How it was built:**

1. **Enumeration:** We enumerated all possible candidate SNVs within the TSS ±10kb windows.
2. **Exclusion:** We removed any sites already present in gnomAD or previously scored variant tables.
3. **Weighted Sampling:** We sampled ~1.8M variants per group using weights to match the **6-class mutation spectrum** and **methylation-bin spectrum** of the observed gnomAD variants.
4. **Inference:** By holding the mutational spectrum fixed, any remaining differences in LFC distributions are interpretable as true differences in predicted regulatory impact rather than mutational bias.

---

## Scoring Logic: AlphaGenome & Exon Mask

The scores in these files were generated using the **AlphaGenome Variant Scorer**, which quantifies the impact of a variant on overall gene transcript abundance.

### Gene Exon Mask

The primary scorer calculates the **Log Fold Change (LFC)** of gene expression between the ALT and REF alleles. To isolate the signal relevant to specific transcripts, we apply a **Exon Mask**.

A spatial mask defines the genomic region of interest (e.g., exons, TSS, or gene body) based on GTF annotations. Values outside of this mask are discarded during the prediction process to ensure the LFC reflects the impact on the target gene's coding and regulatory sequences.

*Source: [AlphaGenome Documentation*](https://www.alphagenomedocs.com/variant_scoring.html)

---

## Column Definitions (Gene Level Summaries)

| Column Group | Description |
| --- | --- |
| **Identifiers** | `gene_id`, `min_variant_id`, `max_variant_id` |
| **Genetic Variance ()** | `vg_predicted_sum`: Calculated as , where  is AF and  is the raw LFC score. |
| **Exon Stats** | Counts (`n_variants_exon`) and distribution of absolute/signed LFC scores (`mean`, `median`, `std`, `min`, `max`, `q90`). |
| **Spatial Bins** | Aggregated stats for variants in specific windows: `up2kb`, `up10kb`, `up100kb`, `down2kb`, and `promo` (promoter). |
| **Enrichment & Decay** | `enrich_up_vs_down_2kb`, `enrich_up_vs_intron_2kb`, and `distance_decay_slope` (LFC magnitude vs. TSS distance). |
| **External Metrics** | Population/Biological constraints: `loeuf_score`, `pLI`, `ncRVIS`, `ncGERP`, `pHaplo`, `pTriplo`, and `Episcore`. |
| **Expression** | `tpm`, `median_tpm`, `log10_tpm1` (Muscle skeletal context). |
| **Architecture** | `genomic_length`, `exonic_length`, `intronic_length`, `utr_length`, `utr5_length`. |

---

## Calculation Notes

* **Variant Level :** For an individual variant, the contribution is .
* **Gene Level :** The `vg_predicted_sum` is the sum of all variant-level  contributions within the gene's defined window.
* **Missing AF:** For variants without an observed Allele Frequency (common in ClinGen or synthetic NULL sets), AF is handled via backfilling or excluded from the variance sum depending on the specific experiment config.