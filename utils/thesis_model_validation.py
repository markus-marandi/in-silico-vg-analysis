from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns
from scipy import stats

from config import DISPLAY_NAMES, GENE_PATHS, SOURCE_PALETTE, VARIANT_PATHS


PROJECT_ROOT = Path(__file__).resolve().parents[1]
THESIS_FIGURE_DIR = PROJECT_ROOT / "10_THESIS" / "FIGURES_FOR_THESIS"
AG_BORZOI_ROOT = Path(
    "/Users/markus/master-thesis-lappalainen/seq2expr-variance/data/output/dataset1/AG-Borzoi_correlations"
)
AG_BORZOI_MATCHED_TISSUES = AG_BORZOI_ROOT / "matched_tissues.txt"
AG_BORZOI_BY_TISSUE = AG_BORZOI_ROOT / "by_tissue_correlations.tsv"

REAL_DATASETS = ("clingen", "background")
POOLED_REAL_DATASET = "pooled_real"
REAL_LABELS = {
    key: DISPLAY_NAMES[key].replace("\n", " ") for key in (*REAL_DATASETS, POOLED_REAL_DATASET)
}
SMALL_EFFECT_THRESHOLD = 0.02
VARIANT_SIZE_EPS = 1e-12


def _ensure_output_dir() -> Path:
    THESIS_FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    return THESIS_FIGURE_DIR


def _dataset_order(value: str) -> int:
    return REAL_DATASETS.index(value)


def _sample_dataframe(frame: pl.DataFrame, sample_n: int, seed: int = 42) -> pl.DataFrame:
    if frame.height <= sample_n:
        return frame
    return frame.sample(n=sample_n, shuffle=True, seed=seed)


def _safe_log10(values: np.ndarray, minimum: float = 1e-12) -> np.ndarray:
    return np.log10(np.clip(values, minimum, None))


def _compute_log_pearson(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    lx = _safe_log10(x)
    ly = _safe_log10(y)
    result = stats.pearsonr(lx, ly)
    return float(result.statistic), float(result.pvalue)


def _replace_zero_for_log(values: np.ndarray, fallback: float = 1e-6) -> np.ndarray:
    positive = values[values > 0]
    replacement = max(float(np.min(positive)) / 2.0, fallback) if positive.size else fallback
    return np.where(values > 0, values, replacement)


def _scale_marker_sizes(values: np.ndarray) -> np.ndarray:
    safe = _safe_log10(values, minimum=VARIANT_SIZE_EPS)
    lo, hi = np.quantile(safe, [0.05, 0.95])
    if np.isclose(hi, lo):
        return np.full(values.shape, 36.0)
    scaled = np.clip((safe - lo) / (hi - lo), 0.0, 1.0)
    return 12.0 + 156.0 * scaled


def _correlation_record(
    dataset_key: str,
    x_metric: str,
    y_metric: str,
    x_values: np.ndarray,
    y_values: np.ndarray,
) -> dict[str, float | str | int]:
    spearman = stats.spearmanr(x_values, y_values)
    pearson_raw = stats.pearsonr(x_values, y_values)
    pearson_log_r, pearson_log_p = _compute_log_pearson(x_values, y_values)
    return {
        "dataset_key": dataset_key,
        "dataset": REAL_LABELS[dataset_key],
        "x_metric": x_metric,
        "y_metric": y_metric,
        "n_genes": int(len(x_values)),
        "spearman_rho": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
        "pearson_raw_r": float(pearson_raw.statistic),
        "pearson_raw_p": float(pearson_raw.pvalue),
        "pearson_log10_r": pearson_log_r,
        "pearson_log10_p": pearson_log_p,
    }


@lru_cache(maxsize=1)
def load_real_gene_tables() -> dict[str, pl.DataFrame]:
    return {dataset: pl.read_parquet(GENE_PATHS[dataset]) for dataset in REAL_DATASETS}


def _load_pooled_real_gene_table(columns: list[str]) -> pd.DataFrame:
    return pl.concat([table.select(columns) for table in load_real_gene_tables().values()]).to_pandas()


@lru_cache(maxsize=None)
def load_real_variant_table(dataset: str, columns: tuple[str, ...] | None = None) -> pl.DataFrame:
    if dataset not in REAL_DATASETS:
        raise ValueError(f"Unsupported dataset: {dataset}")

    lazy = pl.scan_parquet(VARIANT_PATHS[dataset])
    if columns:
        requested = [column for column in columns if column in lazy.collect_schema().names()]
        lazy = lazy.select(requested)
    return lazy.collect()


@lru_cache(maxsize=None)
def choose_representative_genes() -> pd.DataFrame:
    rows: list[dict[str, float | str | int]] = []
    for dataset in REAL_DATASETS:
        gene_df = (
            load_real_gene_tables()[dataset]
            .select(["gene_id", "gene_symbol", "N90", "n_variants", "vg_predicted"])
            .filter(pl.col("gene_symbol").is_not_null() & (pl.col("gene_symbol") != ""))
        )

        median_n90 = float(gene_df["N90"].median())
        median_variants = float(gene_df["n_variants"].median())
        median_vg = float(gene_df["vg_predicted"].median())

        representative = (
            gene_df.with_columns(
                [
                    ((pl.col("N90") - median_n90).abs() / max(median_n90, 1.0)).alias("n90_dev"),
                    ((pl.col("n_variants") - median_variants).abs() / max(median_variants, 1.0)).alias(
                        "variant_dev"
                    ),
                    ((pl.col("vg_predicted") - median_vg).abs() / max(median_vg, 1e-12)).alias("vg_dev"),
                ]
            )
            .with_columns((pl.col("n90_dev") + pl.col("variant_dev") + pl.col("vg_dev")).alias("rep_score"))
            .sort("rep_score")
            .row(0, named=True)
        )

        rows.append(
            {
                "dataset_key": dataset,
                "dataset": REAL_LABELS[dataset],
                "gene_id": representative["gene_id"],
                "gene_symbol": representative["gene_symbol"],
                "N90": int(representative["N90"]),
                "n_variants": int(representative["n_variants"]),
                "vg_predicted": float(representative["vg_predicted"]),
                "selection_score": float(representative["rep_score"]),
                "selection_rule": "closest to dataset medians for N90, variant count, and vg_predicted",
            }
        )

    return pd.DataFrame(rows).sort_values("dataset_key", key=lambda col: col.map(_dataset_order)).reset_index(
        drop=True
    )


def summarize_variant_effects() -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries: list[dict[str, float | str | int]] = []
    gene_tables = load_real_gene_tables()
    representatives = choose_representative_genes()

    for dataset in REAL_DATASETS:
        variant_df = load_real_variant_table(dataset, ("gene_id", "raw_score", "AF"))
        abs_scores = variant_df["raw_score"].abs()
        gene_df = gene_tables[dataset]
        representative = representatives.loc[representatives["dataset_key"] == dataset].iloc[0]

        summaries.append(
            {
                "dataset": REAL_LABELS[dataset],
                "n_variants": int(variant_df.height),
                "n_genes": int(variant_df["gene_id"].n_unique()),
                "prop_|raw_score|<0.02": float(
                    variant_df.select((pl.col("raw_score").abs() < SMALL_EFFECT_THRESHOLD).mean()).item()
                ),
                "median_|raw_score|": float(abs_scores.median()),
                "q90_|raw_score|": float(abs_scores.quantile(0.90)),
                "q99_|raw_score|": float(abs_scores.quantile(0.99)),
                "median_N90": float(gene_df["N90"].median()),
                "iqr_low_N90": float(gene_df["N90"].quantile(0.25)),
                "iqr_high_N90": float(gene_df["N90"].quantile(0.75)),
                "representative_gene": str(representative["gene_symbol"]),
            }
        )

    summary_df = pd.DataFrame(summaries)
    return summary_df, representatives.copy()


def _representative_gene_variants(dataset: str, gene_id: str) -> pd.DataFrame:
    frame = (
        pl.scan_parquet(VARIANT_PATHS[dataset])
        .filter(pl.col("gene_id") == gene_id)
        .select(["variant_id", "raw_score", "AF", "perm_AF"])
        .collect()
        .with_columns(
            [
                pl.col("raw_score").abs().alias("abs_score"),
                (2.0 * pl.col("AF") * (1.0 - pl.col("AF")) * pl.col("raw_score").pow(2)).alias("variant_vg"),
            ]
        )
    )

    pdf = frame.to_pandas()
    pdf["AF_plot"] = _replace_zero_for_log(pdf["AF"].to_numpy())
    return pdf


def plot_variant_effect_distributions(
    sample_n: int = 200_000,
    representatives: pd.DataFrame | None = None,
) -> plt.Figure:
    """Plot the variant effect distribution histogram for both real gene sets."""
    _ensure_output_dir()

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    for dataset in REAL_DATASETS:
        sample = _sample_dataframe(load_real_variant_table(dataset, ("raw_score",)), sample_n=sample_n)
        values = _safe_log10(np.abs(sample["raw_score"].to_numpy()), minimum=1e-6)
        ax.hist(
            values,
            bins=60,
            density=True,
            alpha=0.45,
            color=SOURCE_PALETTE[dataset],
            label=REAL_LABELS[dataset],
        )
    ax.set_xlabel("log10(|effect size|)")
    ax.set_ylabel("Density")
    ax.set_title("Variant effect distribution")
    ax.legend(frameon=True)
    return fig


def plot_representative_gene_scatter(
    dataset: str,
    representative: pd.Series | None = None,
) -> plt.Figure:
    """Plot effect size vs allele frequency for a single representative gene.

    args:
        dataset (str): dataset key, one of 'clingen' or 'background'
        representative (pd.Series | None): row from choose_representative_genes(); chosen automatically if None

    returns:
        plt.Figure: single-panel scatter figure
    """
    _ensure_output_dir()

    if representative is None:
        reps = choose_representative_genes()
        representative = reps.loc[reps["dataset_key"] == dataset].iloc[0]

    panel_labels = {
        "clingen": "C. Representative ClinGen HI gene",
        "background": "D. Representative background gene",
    }

    plot_df = _representative_gene_variants(dataset, str(representative["gene_id"]))
    sizes = _scale_marker_sizes(plot_df["variant_vg"].to_numpy())

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.scatter(
        plot_df["AF_plot"],
        plot_df["raw_score"],
        s=sizes,
        alpha=0.35,
        color=SOURCE_PALETTE[dataset],
        edgecolor="none",
    )

    top_points = plot_df.nlargest(min(12, len(plot_df)), "variant_vg")
    if not top_points.empty:
        ax.scatter(
            top_points["AF_plot"],
            top_points["raw_score"],
            s=_scale_marker_sizes(top_points["variant_vg"].to_numpy()) + 28.0,
            facecolors="none",
            edgecolors="black",
            linewidths=0.7,
        )

    ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.4)
    ax.set_xscale("log")
    ax.set_xlabel("Allele frequency")
    ax.set_ylabel("Effect size")
    ax.set_title(panel_labels.get(dataset, f"Representative {dataset} gene"))
    ax.text(
        0.03,
        0.97,
        (
            f"{representative['gene_symbol']}\n"
            f"n={int(representative['n_variants']):,} variants\n"
            f"N90={int(representative['N90'])}"
        ),
        transform=ax.transAxes,
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )

    plt.tight_layout()
    return fig


def summarize_ag_borzoi_correlations() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not AG_BORZOI_BY_TISSUE.exists():
        raise FileNotFoundError(f"Missing AG/Borzoi correlation table: {AG_BORZOI_BY_TISSUE}")

    by_tissue = pl.read_csv(AG_BORZOI_BY_TISSUE, separator="\t")
    gene_spearman = (
        by_tissue.filter((pl.col("level") == "gene") & (pl.col("method") == "spearman"))
        .sort("rho", descending=True)
        .to_pandas()
    )

    spearman_level_summary = (
        by_tissue.filter(pl.col("method") == "spearman")
        .group_by("level")
        .agg(
            [
                pl.len().alias("n_tissues"),
                pl.col("rho").mean().alias("mean_rho"),
                pl.col("rho").median().alias("median_rho"),
                pl.col("rho").min().alias("min_rho"),
                pl.col("rho").max().alias("max_rho"),
            ]
        )
        .sort("level")
        .to_pandas()
    )

    top = gene_spearman.iloc[0]
    bottom = gene_spearman.iloc[-1]
    overview = pd.DataFrame(
        [
            {
                "n_tissues": int(len(gene_spearman)),
                "positive_fraction": float((gene_spearman["rho"] > 0).mean()),
                "median_gene_spearman": float(gene_spearman["rho"].median()),
                "top_tissue": str(top["tissue"]),
                "top_rho": float(top["rho"]),
                "bottom_tissue": str(bottom["tissue"]),
                "bottom_rho": float(bottom["rho"]),
                "matched_genes_per_tissue": int(gene_spearman["n"].median()),
            }
        ]
    )

    return overview, spearman_level_summary, gene_spearman


def plot_ag_borzoi_correlations(
    gene_spearman: pd.DataFrame | None = None,
    spearman_level_summary: pd.DataFrame | None = None,
) -> plt.Figure:
    _ensure_output_dir()
    overview, level_summary, gene_level = summarize_ag_borzoi_correlations()
    gene_spearman = gene_level if gene_spearman is None else gene_spearman
    spearman_level_summary = level_summary if spearman_level_summary is None else spearman_level_summary

    by_tissue = pl.read_csv(AG_BORZOI_BY_TISSUE, separator="\t")
    spearman_all_levels = by_tissue.filter(pl.col("method") == "spearman").to_pandas()

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(16, 9),
        gridspec_kw={"width_ratios": [2.3, 1.2]},
        constrained_layout=True,
    )

    ax = axes[0]
    colors = np.where(gene_spearman["rho"] >= 0, "#0F766E", "#B45309")
    ax.barh(gene_spearman["tissue"], gene_spearman["rho"], color=colors)
    ax.axvline(0.0, color="black", linewidth=1.0, alpha=0.6)
    ax.invert_yaxis()
    ax.set_xlabel("Spearman rho")
    ax.set_ylabel("GTEx tissue code")
    ax.set_title("A. Gene-level AlphaGenome vs Borzoi correlation")
    ax.text(
        0.02,
        0.02,
        (
            f"{overview.loc[0, 'positive_fraction']:.0%} tissues > 0\n"
            f"Top: {overview.loc[0, 'top_tissue']} ({overview.loc[0, 'top_rho']:.2f})\n"
            f"Bottom: {overview.loc[0, 'bottom_tissue']} ({overview.loc[0, 'bottom_rho']:.2f})"
        ),
        transform=ax.transAxes,
        va="bottom",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )

    ax = axes[1]
    level_order = ["gene", "variant", "row"]
    sns.boxplot(
        data=spearman_all_levels,
        x="rho",
        y="level",
        order=level_order,
        ax=ax,
        color="#D1D5DB",
        showfliers=False,
        orient="h",
    )
    sns.stripplot(
        data=spearman_all_levels,
        x="rho",
        y="level",
        order=level_order,
        ax=ax,
        hue="level",
        palette={"gene": "#0F766E", "variant": "#2563EB", "row": "#7C3AED"},
        size=4,
        alpha=0.55,
        orient="h",
        dodge=False,
    )
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()
    ax.axvline(0.0, color="black", linewidth=1.0, alpha=0.6)
    ax.set_xlabel("Spearman rho")
    ax.set_ylabel("")
    ax.set_title("B. Correlation spread by aggregation level")
    return fig


def summarize_eqtl_concordance() -> pd.DataFrame:
    rows: list[dict[str, float | str | int]] = []
    pairs = (
        ("vg_predicted", "vg_eqtl"),
        ("sum_sq_raw_score", "vg_eqtl"),
        ("vg_predicted", "vgh"),
    )

    pdf = _load_pooled_real_gene_table(["vg_predicted", "vg_eqtl", "vgh", "sum_sq_raw_score"])

    for x_metric, y_metric in pairs:
        subset = pdf[[x_metric, y_metric]].dropna()
        subset = subset[(subset[x_metric] > 0) & (subset[y_metric] > 0)]
        rows.append(
            _correlation_record(
                POOLED_REAL_DATASET,
                x_metric,
                y_metric,
                subset[x_metric].to_numpy(),
                subset[y_metric].to_numpy(),
            )
        )

    summary = pd.DataFrame(rows)
    summary["comparison"] = summary["x_metric"] + " vs " + summary["y_metric"]
    return summary[
        [
            "dataset",
            "comparison",
            "n_genes",
            "spearman_rho",
            "pearson_log10_r",
            "pearson_raw_r",
            "spearman_p",
        ]
    ]


def _add_scatter_panel(ax: plt.Axes, x_metric: str, y_metric: str, title: str, x_label: str, y_label: str) -> None:
    plot_df = _load_pooled_real_gene_table([x_metric, y_metric])
    subset = plot_df[[x_metric, y_metric]].dropna()
    subset = subset[(subset[x_metric] > 0) & (subset[y_metric] > 0)]
    x_values = subset[x_metric].to_numpy()
    y_values = subset[y_metric].to_numpy()
    log_x = _safe_log10(x_values)
    log_y = _safe_log10(y_values)
    slope, intercept = np.polyfit(log_x, log_y, deg=1)
    line_x = np.geomspace(x_values.min(), x_values.max(), 200)
    line_y = 10 ** (intercept + slope * _safe_log10(line_x))

    ax.scatter(
        x_values,
        y_values,
        s=26,
        alpha=0.55,
        color=SOURCE_PALETTE[POOLED_REAL_DATASET],
        label=REAL_LABELS[POOLED_REAL_DATASET],
    )
    ax.plot(
        line_x,
        line_y,
        color="#111827",
        linewidth=1.6,
        label="log-log fit",
    )

    spearman = stats.spearmanr(x_values, y_values)
    pearson_log_r, _ = _compute_log_pearson(x_values, y_values)
    metric_line = (
        f"{REAL_LABELS[POOLED_REAL_DATASET]}: ρ={float(spearman.statistic):.2f}, "
        f"r_log={pearson_log_r:.2f}, n={len(x_values)}"
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(frameon=True, fontsize=9)
    ax.text(
        0.03,
        0.03,
        metric_line,
        transform=ax.transAxes,
        va="bottom",
        fontsize=8.5,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )


def plot_eqtl_concordance() -> plt.Figure:
    _ensure_output_dir()

    fig, ax = plt.subplots(figsize=(5.5, 5), constrained_layout=True)
    _add_scatter_panel(
        ax,
        x_metric="vg_predicted",
        y_metric="vg_eqtl",
        title="Predicted Vg vs GTEx eQTL variance",
        x_label="Predicted Vg",
        y_label="GTEx eQTL variance",
    )
    return fig


__all__ = [
    "AG_BORZOI_BY_TISSUE",
    "AG_BORZOI_MATCHED_TISSUES",
    "AG_BORZOI_ROOT",
    "POOLED_REAL_DATASET",
    "REAL_DATASETS",
    "REAL_LABELS",
    "THESIS_FIGURE_DIR",
    "choose_representative_genes",
    "load_real_gene_tables",
    "load_real_variant_table",
    "plot_ag_borzoi_correlations",
    "plot_eqtl_concordance",
    "plot_representative_gene_scatter",
    "plot_variant_effect_distributions",
    "summarize_ag_borzoi_correlations",
    "summarize_eqtl_concordance",
    "summarize_variant_effects",
]
