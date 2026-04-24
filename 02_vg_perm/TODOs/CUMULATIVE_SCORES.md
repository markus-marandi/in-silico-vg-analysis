Cumulative scoring — from per-variant absolute scores to gene-level cumulative scores
⸻

Which “cumulative” do you mean? (pick one or more)
	1.	Simple sum of effects
S_g = \sum_{i\in g} |\beta_i|
— easy measure of total absolute effect per gene.
	2.	Sum of squared effects (proportional to variance if AF ignored)
Q_g = \sum_{i\in g} \beta_i^2
— useful when effects combine additively in variance.
	3.	Allele-frequency weighted (genetic variance) — your existing V_g
V_g = \sum_{i\in g} 2 p_i (1-p_i) \beta_i^2
— best when you want expected per-individual variance.
	4.	Cumulative contribution by rank within a gene
Sort a gene’s variants by |β| (or by per-variant variance contribution) and compute cumulative sums to see how many variants explain X% of the gene’s total.

⸻


⸻

Polars + NumPy code snippets

Assume a Parquet with columns: gene_id, abs_effect (|β|), beta (signed), af (allele frequency), variant_id.

1) Compute basic gene-level aggregates (polars)

import polars as pl
import numpy as np
from typing import Any

# lazy load and compute aggregates
lf = pl.scan_parquet("path/to/variants.parquet").select(
    ["gene_id", "abs_effect", "beta", "af", "variant_id"]
)

agg = (
    lf
    .filter(pl.col("abs_effect").is_not_null())
    .groupby("gene_id")
    .agg([
        pl.sum("abs_effect").alias("sum_abs_effect"),
        pl.sum(pl.col("abs_effect")**2).alias("sum_abs_effect_sq"),
        pl.count().alias("n_variants"),
        pl.mean("abs_effect").alias("mean_abs_effect"),
        # Vg per gene using observed AF:
        pl.sum(2 * pl.col("af") * (1 - pl.col("af")) * (pl.col("abs_effect")**2)).alias("observed_vg"),
    ])
    .collect()
)

print(agg.head())

2) Compute per-variant variance and get per-gene cumulative (in-memory with NumPy)

This computes, for each gene, the ranked cumulative fraction explained by variants.

import polars as pl
import numpy as np

df = pl.read_parquet("path/to/variants.parquet")  # use collect if small enough
# ensure necessary columns exist
df = df.filter(pl.col("abs_effect").is_not_null())

# compute per-variant contribution to Vg
df = df.with_columns(
    (2 * pl.col("af") * (1 - pl.col("af")) * (pl.col("abs_effect") ** 2)).alias("vg_i")
)

# group into python dict of arrays (fast if gene-wise arrays are moderate)
groups = df.groupby("gene_id").agg(
    [
        pl.col("variant_id").list().alias("variants"),
        pl.col("abs_effect").list().alias("abs_effects"),
        pl.col("vg_i").list().alias("vg_i_list")
    ]
).to_dict(as_series=False)

# now compute within-gene cumulative curves
cumulative_results = []
for gene, variants, abs_effects, vg_i_list in zip(
    groups["gene_id"], groups["variants"], groups["abs_effects"], groups["vg_i_list"]
):
    vg_arr = np.array(vg_i_list, dtype=float)
    if vg_arr.size == 0:
        continue
    # rank variants by contribution (descending)
    order = np.argsort(vg_arr)[::-1]
    vg_sorted = vg_arr[order]
    cumsum = np.cumsum(vg_sorted)
    total = cumsum[-1]
    frac_explained = cumsum / total  # fraction of gene Vg explained by top k variants

    cumulative_results.append({
        "gene_id": gene,
        "n_variants": vg_arr.size,
        "total_vg": float(total),
        "top1_frac": float(frac_explained[0]),
        "top3_frac": float(frac_explained[min(2, vg_arr.size-1)]),
        "top10_frac": float(frac_explained[min(9, vg_arr.size-1)]),
        # optionally keep the full curve (but memory heavy)
        # "frac_curve": frac_explained.tolist()
    })

This gives top1_frac, top3_frac, etc. — a compact measure of whether a gene’s variance is dominated by a few big variants.

3) Genome-wide cumulative curve (Lorenz / Pareto)

Compute per-gene totals, sort genes, and plot cumulative fraction of the total Vg across genes.

import polars as pl
import numpy as np
import matplotlib.pyplot as plt

df = pl.read_parquet("path/to/gene_aggregates.parquet")  # contains gene_id, total_vg
arr = df.select(["gene_id", "total_vg"]).to_numpy()

# sort by total_vg descending
sorted_idx = np.argsort(arr[:,1].astype(float))[::-1]
v_sorted = arr[sorted_idx,1].astype(float)
cum = np.cumsum(v_sorted)
cum_frac = cum / cum[-1]
gene_frac = np.arange(1, len(v_sorted)+1) / len(v_sorted)

plt.plot(gene_frac, cum_frac)
plt.xlabel("Fraction of genes (ranked by Vg)")
plt.ylabel("Cumulative fraction of total Vg")
plt.title("Genome-level cumulative contribution of genes")
plt.grid(True)
plt.show()

This visual immediately shows whether variance is concentrated in a few genes (steep initial rise) or distributed.

⸻

Practical recommendations / choices
	•	If you care about biological impact: use the AF-weighted Vg (option 3). It’s interpretable as expected variance in individuals.
	•	If you care about effect magnitude only (independent of AF): use sum(|β|) or sum(β^2).
	•	To describe sparsity vs polygenicity per gene: compute within-gene cumulative curves (option 4) and report top1, top3, top10 fractions.
	•	For comparing genes of different sizes: normalize by n_variants or by gene length or use percentile ranks.
	•	For visualization: use Lorenz (genome-wide) and per-gene cumulative plots (for top genes).

⸻

Integration with your permutation engine

If you want cumulative metrics in permutations:
	•	For each permutation iteration you already compute vg_per_gene_shuffled. Keep it.
	•	You can also compute, per iteration, the top-k fraction for each gene by:
	•	computing per-variant vg_i with the sampled AFs,
	•	for each gene, sort its vg_i and compute cumsum/top_total.
	•	Memory constraint: instead of storing curves for all genes and iterations, store only scalar summaries per gene per iteration (e.g., top1_frac, top3_frac, top10_frac). Then compute percentile ranks of those summaries vs observed.

Example: store top3_frac per gene per iteration in your shuffled_vg_accumulator-style matrix, then compute (count(shuffled_top3_frac > observed_top3_frac)/n_iter).

⸻

Suggested small API for your codebase

def gene_totals_from_variants(df_variants: pl.DataFrame) -> pl.DataFrame:
    """Return per-gene sums: sum_abs, sum_sq, observed_vg, n_variants"""
    ...

def gene_within_cumulative(df_variants: pl.DataFrame, top_k: int = 3) -> pl.DataFrame:
    """Return per-gene fraction explained by top_k variants (using vg_i)"""
    ...

I can generate these two helper functions in Polars/NumPy style and add them to your existing permutation_vg.py so they are computed for both observed and shuffled data (and included in outputs). Want me to add those functions and wire them into the engine so permutation outputs include top1_frac, top3_frac, and top10_frac?