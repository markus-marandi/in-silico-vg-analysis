For a gene g with variants i \in g, let the per-variant contribution to genetic variance be
v_i = 2 p_i(1-p_i)\beta_i^2
and total gene variance V_g = \sum_i v_i. Sort variants in descending v_i.

Add new parameters per gene:
	1.	top1_frac — fraction of V_g explained by the single top variant (largest v_i).
	2.	top3_frac — fraction of V_g explained by the top 3 variants (or fewer if gene has <3 variants).
	3.	top10_frac — fraction of V_g explained by the top 10 variants (or fewer).
	4.	frac_variants_for_90 — minimal fraction of variants (k / n_variants) needed to reach ≥ 90% of V_g.
	5.	frac_variants_for_99 — minimal fraction of variants needed to reach ≥ 99% of V_g.
	6.	count_variants_for_90, count_variants_for_99 — the raw counts k (optional, but often useful).


	•	Variant A: AC = 1 → singleton
	•	Variant B: AC = 2 → doubleton
	•	Variant C: AC = 3 → tripleton
---
Gene Annotation Specifications: Variance Architecture and Population Rarity

Overview

Each gene is annotated along two complementary axes:
	1.	Variance architecture — how genetic variance (V_g) is distributed across variants within a gene.
	2.	Population rarity composition — how rare the gene’s variants are in an external reference population.

These annotations capture distinct biological properties:
	•	Variance architecture reflects regulatory effect concentration.
	•	Rarity composition reflects population-level selective pressure.

The two axes are orthogonal by construction and are analyzed jointly to assess whether genes with constrained variance architectures are enriched for rare variants.

⸻

Data Sources

External reference datasets
	•	gnomAD
Used to annotate population allele counts (AC) and allele frequencies (AF) for all variants.
	•	AC defines singleton/doubleton/tripleton status.
	•	AF is used to compute expected genetic variance.

Internal model outputs
	•	Exon-mask regulatory scores
Variant-level predicted effect sizes (β) derived from exon-masked regulatory models.
	•	Variants are defined as ±10 kb around genes, filtered by exon mask.
	•	Absolute effect size is denoted as |\beta|.

----

Variant-Level Quantities

For each variant i:

Population rarity (external)
	•	Allele count AC_i from gnomAD:
	•	AC_i = 1 → singleton
	•	AC_i = 2 → doubleton
	•	AC_i = 3 → tripleton

Genetic variance contribution (internal + external)

v_i = 2 \cdot p_i \cdot (1 - p_i) \cdot \beta_i^2
where:
	•	p_i is the allele frequency from gnomAD,
	•	\beta_i is the exon-mask regulatory effect size.

⸻

Gene-Level Annotations

All gene-level metrics are computed using variants assigned to that gene under the exon-mask definition.

⸻

Axis 1: Variance Architecture Metrics (Internal)

These metrics describe how concentrated or diffuse genetic variance is within a gene.

Let V_g = \sum_i v_i be total gene-level genetic variance, and variants be ranked in descending order of v_i.

Core metrics
	1.	top1_frac
	•	Fraction of V_g explained by the single largest-variance variant.
\text{top1\_frac} = \frac{v_{(1)}}{V_g}
	2.	top3_frac
	•	Fraction of V_g explained by the top 3 variants (or fewer if n<3).
	3.	top10_frac
	•	Fraction of V_g explained by the top 10 variants (or fewer if n<10).
	4.	frac_variants_for_90
	•	Minimal fraction of variants required to explain ≥90% of V_g.
\frac{k_{90}}{n_{\text{variants}}}
	5.	frac_variants_for_99
	•	Minimal fraction of variants required to explain ≥99% of V_g.
	6.	count_variants_for_90 / count_variants_for_99
	•	Raw number of variants k required to reach 90% / 99% of V_g.

Interpretation
	•	High top1_frac or low frac_variants_for_90 indicates variance concentration.
	•	Low top1_frac and high frac_variants_for_90 indicates polygenic architecture.

⸻

Axis 2: Population Rarity Composition (External)

These annotations summarize how rare the gene’s variants are in the population.

Variant-level definitions (from gnomAD)
	•	Singleton: AC = 1
	•	Doubleton: AC = 2
	•	Tripleton: AC = 3

Gene-level counts
For each gene g:
	•	n_singleton: number of exon-mask variants with AC = 1
	•	n_doubleton: number with AC = 2
	•	n_tripleton: number with AC = 3
	•	n_variants: total number of exon-mask variants

Optional gene categories
Genes may additionally be categorized as:
	•	has_singleton
	•	singleton_only
	•	singleton_dominant (e.g., >50% of variants are singletons)
	•	mixed_rare
	•	no_singletons

These categories are derived solely from AC, independent of effect size.

⸻

Complementarity of the Two Annotation Axes

The two annotation axes answer different but related biological questions:

Axis	Captures	Derived from
Variance architecture	How few or many variants drive gene variance	Model effects + AF
Population rarity	How rare the gene’s variants are	gnomAD AC

Key principle

Rarity does not imply dominance, and dominance does not imply rarity.

	•	A singleton variant may or may not dominate gene variance.
	•	A dominant variance-driving variant may or may not be a singleton.

⸻

Downstream Correlation Analysis

Joint analysis of these annotations enables questions such as:
	•	Are genes with high variance concentration enriched for singleton-dominated variant sets?
	•	Do ClinGen genes show a stronger coupling between rarity and variance architecture than background genes?
	•	Are genes under strong constraint characterized by both:
	•	low observed V_g (via permutation tests), and
	•	high variance concentration in rare variants?

Correlation analyses may include:
	•	Spearman or Kendall correlation between top1_frac and n_singleton
	•	Stratified comparisons of frac_variants_for_90 across rarity-defined gene categories
	•	Regression models predicting constraint rank from both axes jointly

⸻

Summary
	•	Variance architecture metrics quantify how genetic variance is distributed within a gene.
	•	Rarity annotations quantify how rare the gene’s variants are in the population.
	•	Using exon-mask scores strengthens biological interpretability without conflating the two axes.
	•	Together, these annotations enable a principled assessment of selective constraint and regulatory architecture.
