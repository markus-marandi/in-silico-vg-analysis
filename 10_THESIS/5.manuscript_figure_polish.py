from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

ROOT = Path(__file__).resolve().parents[1]
THESIS_ROOT = ROOT.parent / "2026SS-thesis-medical-engineering-msc"
MAIN_FIGURE_DIR = THESIS_ROOT / "Figures"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import GENE_PATHS, SOURCE_PALETTE, VG_MIN_CLIP  # noqa: E402
from thesis_style import apply_thesis_style, save_figure  # noqa: E402


apply_thesis_style()


def clip_for_log(values: np.ndarray, clip: float = VG_MIN_CLIP) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return np.where(values <= 0, clip, values)


def save_manuscript_figure(fig: plt.Figure, name: str) -> list[Path]:
    paths: list[Path] = []
    paths.extend(
        save_figure(
            fig,
            name=name,
            outdir=ROOT / "10_THESIS" / "FIGURES_FOR_THESIS",
            formats=("pdf",),
            tight=False,
        )
    )
    paths.extend(
        save_figure(
            fig,
            name=name,
            outdir=MAIN_FIGURE_DIR,
            formats=("pdf",),
            tight=False,
        )
    )
    return paths


def load_gene_tables() -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    frames = [
        pl.scan_parquet(path).with_columns(pl.lit(source).alias("source"))
        for source, path in GENE_PATHS.items()
    ]
    df_all = pl.concat(frames, how="diagonal_relaxed").collect()
    df_real = df_all.filter(~pl.col("source").str.contains("null"))
    df_synth = df_all.filter(pl.col("source").str.contains("null"))
    df_hi = df_real.filter(pl.col("source") == "clingen")
    df_bg = df_real.filter(pl.col("source") == "background")
    return df_all, df_real, df_synth, df_hi, df_bg


def build_regional_boxplot(df_all: pl.DataFrame) -> None:
    source_order = ["background", "background_null", "clingen", "clingen_null"]
    source_labels = {
        "background": "Background\n(gnomAD)",
        "background_null": "Background\n(Synthetic)",
        "clingen": "ClinGen HI\n(gnomAD)",
        "clingen_null": "ClinGen HI\n(Synthetic)",
    }

    required_cols = [
        "source",
        "vg_promoter_core",
        "vg_promoter_core_perm",
        "vg_proximal_upstream",
        "vg_proximal_upstream_perm",
        "vg_distal_upstream",
        "vg_distal_upstream_perm",
        "vg_down_proximal",
        "vg_down_proximal_perm",
        "vg_down_distal",
        "vg_down_distal_perm",
    ]
    plot_df = (
        df_all.select(required_cols)
        .with_columns(
            [
                pl.when(pl.col("source").str.contains("null"))
                .then(pl.col("vg_promoter_core_perm"))
                .otherwise(pl.col("vg_promoter_core"))
                .alias("vg_region_promoter_core"),
                pl.when(pl.col("source").str.contains("null"))
                .then(pl.col("vg_proximal_upstream_perm"))
                .otherwise(pl.col("vg_proximal_upstream"))
                .alias("vg_region_proximal_upstream"),
                pl.when(pl.col("source").str.contains("null"))
                .then(pl.col("vg_distal_upstream_perm"))
                .otherwise(pl.col("vg_distal_upstream"))
                .alias("vg_region_distal_upstream"),
                pl.when(pl.col("source").str.contains("null"))
                .then(pl.col("vg_down_proximal_perm"))
                .otherwise(pl.col("vg_down_proximal"))
                .alias("vg_region_down_proximal"),
                pl.when(pl.col("source").str.contains("null"))
                .then(pl.col("vg_down_distal_perm"))
                .otherwise(pl.col("vg_down_distal"))
                .alias("vg_region_down_distal"),
            ]
        )
        .to_pandas()
    )

    plot_df["source_label"] = plot_df["source"].map(source_labels)
    region_df = plot_df.copy()
    region_df["vg_upstream_total"] = region_df["vg_region_proximal_upstream"] + region_df["vg_region_distal_upstream"]
    region_df["vg_downstream_total"] = region_df["vg_region_down_proximal"] + region_df["vg_region_down_distal"]

    region_long = region_df.melt(
        id_vars=["source_label"],
        value_vars=["vg_region_promoter_core", "vg_upstream_total", "vg_downstream_total"],
        var_name="region_group",
        value_name="region_vg",
    )
    region_long = region_long[region_long["region_vg"] > 0].copy()
    region_labels = {
        "vg_region_promoter_core": "Promoter core\n(TSS ± 200 bp)",
        "vg_upstream_total": "Upstream\n(-200 bp to -10 kb)",
        "vg_downstream_total": "Downstream\n(+200 bp to +10 kb)",
    }
    region_order = list(region_labels.values())
    region_long["region_group"] = region_long["region_group"].map(region_labels)

    region_palette = {
        region_labels["vg_region_promoter_core"]: "#2563EB",
        region_labels["vg_upstream_total"]: "#F59E0B",
        region_labels["vg_downstream_total"]: "#10B981",
    }

    fig, ax = plt.subplots(figsize=(11.2, 6.7), constrained_layout=True)
    sns.boxplot(
        data=region_long,
        x="source_label",
        y="region_vg",
        hue="region_group",
        order=[source_labels[key] for key in source_order],
        hue_order=region_order,
        palette=region_palette,
        showfliers=True,
        flierprops={"marker": "o", "markersize": 2.4, "alpha": 0.55},
        medianprops={"color": "white", "linewidth": 2.0},
        width=0.76,
        ax=ax,
    )
    ax.set_yscale("log")
    ax.set_xlabel("Gene set")
    ax.set_ylabel(r"Predicted genetic variance, $V_G$")
    ax.set_title("Predicted regulatory variance by genomic region")
    ax.tick_params(axis="x", rotation=0)
    ax.grid(axis="y", alpha=0.28)
    ax.grid(axis="x", visible=False)
    ax.legend(
        title="Region",
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.13),
        frameon=False,
    )
    sns.despine(ax=ax)
    save_manuscript_figure(fig, "4_3_regional_vg_boxplot_publication")
    plt.close(fig)


def compute_cluster_summary(df_real: pl.DataFrame) -> tuple[pd.DataFrame, float]:
    cluster_cols = [
        "vg_distal_upstream",
        "vg_proximal_upstream",
        "vg_promoter_core",
        "vg_down_proximal",
        "vg_down_distal",
    ]
    share_cols = [f"{col}_share" for col in cluster_cols]
    cluster_region_names = {
        "vg_distal_upstream": "Distal upstream",
        "vg_proximal_upstream": "Proximal upstream",
        "vg_promoter_core": "Promoter core",
        "vg_down_proximal": "Proximal downstream",
        "vg_down_distal": "Distal downstream",
    }

    cluster_df = (
        df_real.select(["gene_id", "gene_symbol", "source", *cluster_cols])
        .drop_nulls(["gene_id", "source", *cluster_cols])
        .to_pandas()
    )

    cluster_values = np.clip(cluster_df[cluster_cols].to_numpy(dtype=float), a_min=0.0, a_max=None)
    totals = cluster_values.sum(axis=1)
    cluster_df = cluster_df.loc[totals > 0].copy()
    cluster_values = cluster_values[totals > 0]
    cluster_df["regional_total_vg"] = cluster_values.sum(axis=1)
    cluster_shares = cluster_values / cluster_df["regional_total_vg"].to_numpy()[:, None]
    cluster_df[share_cols] = cluster_shares

    candidate_ks = list(range(2, min(6, len(cluster_df) - 1) + 1))
    silhouette_rows = []
    for k in candidate_ks:
        model = KMeans(n_clusters=k, n_init=50, random_state=11)
        labels = model.fit_predict(cluster_shares)
        silhouette_rows.append(
            {
                "k": k,
                "silhouette": silhouette_score(
                    cluster_shares,
                    labels,
                    sample_size=min(4000, len(cluster_df)),
                    random_state=11,
                ),
            }
        )
    silhouette_df = pd.DataFrame(silhouette_rows).sort_values(
        ["silhouette", "k"],
        ascending=[False, True],
    )
    best_k = int(silhouette_df.iloc[0]["k"])

    kmeans_model = KMeans(n_clusters=best_k, n_init=100, random_state=11)
    raw_cluster_ids = kmeans_model.fit_predict(cluster_shares)

    centroid_df = pd.DataFrame(kmeans_model.cluster_centers_, columns=share_cols)
    centroid_df["center_of_mass"] = centroid_df[share_cols].to_numpy() @ np.arange(len(share_cols))
    centroid_df["dominant_idx"] = np.argmax(centroid_df[share_cols].to_numpy(), axis=1)

    cluster_rank = centroid_df.sort_values(["center_of_mass", "dominant_idx"]).index.to_list()
    cluster_map = {old_id: new_id for new_id, old_id in enumerate(cluster_rank, start=1)}

    cluster_df["cluster_id"] = pd.Series(raw_cluster_ids, index=cluster_df.index).map(cluster_map)
    centroid_df["cluster_id"] = centroid_df.index.map(cluster_map)
    centroid_df = centroid_df.sort_values("cluster_id").reset_index(drop=True)
    centroid_df["dominant_region"] = centroid_df[share_cols].idxmax(axis=1).str.replace("_share", "", regex=False)
    centroid_df["dominant_region_label"] = centroid_df["dominant_region"].map(cluster_region_names)

    overall_hi_fraction = float((cluster_df["source"] == "clingen").mean())

    rows = []
    for cluster_id in sorted(cluster_df["cluster_id"].unique()):
        sub = cluster_df[cluster_df["cluster_id"] == cluster_id].copy()
        hi_n = int((sub["source"] == "clingen").sum())
        bg_n = int((sub["source"] == "background").sum())
        other_hi = int(((cluster_df["cluster_id"] != cluster_id) & (cluster_df["source"] == "clingen")).sum())
        other_bg = int(((cluster_df["cluster_id"] != cluster_id) & (cluster_df["source"] == "background")).sum())
        odds_ratio, fisher_p = stats.fisher_exact([[hi_n, bg_n], [other_hi, other_bg]], alternative="two-sided")

        centroid_row = centroid_df.loc[centroid_df["cluster_id"] == cluster_id].iloc[0]
        rows.append(
            {
                "cluster_id": cluster_id,
                "cluster_label": f"Cluster {cluster_id}",
                "n_genes": len(sub),
                "dominant_region_key": centroid_row["dominant_region"],
                "dominant_region": centroid_row["dominant_region_label"],
                "background_n": bg_n,
                "clingen_n": hi_n,
                "background_fraction": bg_n / len(sub),
                "clingen_fraction": hi_n / len(sub),
                "fisher_odds_ratio": odds_ratio,
                "fisher_odds_ratio_label": "inf" if np.isinf(odds_ratio) else f"{odds_ratio:.2f}",
                "fisher_p": fisher_p,
                **{share_col: float(centroid_row[share_col]) for share_col in share_cols},
            }
        )

    cluster_summary = pd.DataFrame(rows).sort_values("cluster_id").reset_index(drop=True)
    return cluster_summary, overall_hi_fraction


def build_cluster_figure(df_real: pl.DataFrame) -> None:
    cluster_summary, overall_hi_fraction = compute_cluster_summary(df_real)
    share_cols = [
        "vg_distal_upstream_share",
        "vg_proximal_upstream_share",
        "vg_promoter_core_share",
        "vg_down_proximal_share",
        "vg_down_distal_share",
    ]
    axis_labels = [
        "Distal\nupstream",
        "Proximal\nupstream",
        "Promoter\ncore",
        "Proximal\ndownstream",
        "Distal\ndownstream",
    ]
    region_profile_colors = {
        "vg_distal_upstream": "#2563EB",
        "vg_proximal_upstream": "#60A5FA",
        "vg_promoter_core": "#D97706",
        "vg_down_proximal": "#F97316",
        "vg_down_distal": "#DC2626",
    }

    angles = np.linspace(0, 2 * np.pi, len(axis_labels), endpoint=False)
    angles = np.concatenate([angles, [angles[0]]])
    radar_max = float(cluster_summary[share_cols].to_numpy().max())
    radar_max = max(0.30, np.ceil(radar_max / 0.05) * 0.05)

    fig = plt.figure(figsize=(13.2, 10.8), constrained_layout=True)
    outer = fig.add_gridspec(2, 3, wspace=0.18, hspace=0.24)

    for idx, row in cluster_summary.iterrows():
        panel = outer[idx // 3, idx % 3].subgridspec(2, 1, height_ratios=[3.0, 1.45], hspace=0.02)
        profile_color = region_profile_colors[row["dominant_region_key"]]

        radar_ax = fig.add_subplot(panel[0], projection="polar")
        profile = row[share_cols].to_numpy(dtype=float)
        profile = np.concatenate([profile, [profile[0]]])
        radar_ax.plot(angles, profile, color=profile_color, lw=2.4)
        radar_ax.fill(angles, profile, color=profile_color, alpha=0.18)
        radar_ax.set_xticks(angles[:-1])
        radar_ax.set_xticklabels(axis_labels, fontsize=8)
        radar_ax.set_ylim(0, radar_max)
        yticks = np.linspace(radar_max / 3, radar_max, 3)
        radar_ax.set_yticks(yticks)
        radar_ax.set_yticklabels([f"{tick:.0%}" for tick in yticks], fontsize=8)
        radar_ax.grid(alpha=0.35)
        radar_ax.spines["polar"].set_alpha(0.22)
        radar_ax.set_title(
            f"{row['cluster_label']} | {row['dominant_region']} dominant\nn = {row['n_genes']}",
            fontsize=11,
            pad=16,
        )
        radar_ax.text(
            0.5,
            1.08,
            f"HI share {row['clingen_fraction']:.0%} | OR {row['fisher_odds_ratio_label']} | p = {row['fisher_p']:.1e}",
            transform=radar_ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="#334155",
        )

        mix_ax = fig.add_subplot(panel[1])
        bg_pct = row["background_fraction"] * 100
        hi_pct = row["clingen_fraction"] * 100
        overall_pct = overall_hi_fraction * 100

        mix_ax.bar([0], [bg_pct], color=SOURCE_PALETTE["background"], width=0.56)
        mix_ax.bar([0], [hi_pct], bottom=[bg_pct], color=SOURCE_PALETTE["clingen"], width=0.56)
        mix_ax.axhline(overall_pct, color="#0F172A", lw=1.1, ls="--", alpha=0.8)
        mix_ax.set_ylim(0, 100)
        mix_ax.set_xticks([])
        mix_ax.set_yticks([0, 50, 100])
        if idx % 3 == 0:
            mix_ax.set_ylabel("Cluster composition (%)")
        else:
            mix_ax.set_yticklabels([])
        mix_ax.grid(axis="y", alpha=0.25)
        mix_ax.spines[["top", "right"]].set_visible(False)
        mix_ax.text(
            0,
            bg_pct / 2,
            f"BG\n{row['background_n']}",
            ha="center",
            va="center",
            fontsize=8.5,
            fontweight="bold",
            color="white" if bg_pct > 18 else "#111827",
        )
        mix_ax.text(
            0,
            bg_pct + hi_pct / 2,
            f"HI\n{row['clingen_n']}",
            ha="center",
            va="center",
            fontsize=8.5,
            fontweight="bold",
            color="white" if hi_pct > 18 else "#111827",
        )

    legend_ax = fig.add_subplot(outer[1, 2])
    legend_ax.axis("off")
    radar_handles = [
        Line2D([0], [0], color=color, lw=2.4, label=label)
        for label, color in [
            ("Distal upstream", region_profile_colors["vg_distal_upstream"]),
            ("Proximal upstream", region_profile_colors["vg_proximal_upstream"]),
            ("Promoter core", region_profile_colors["vg_promoter_core"]),
            ("Proximal downstream", region_profile_colors["vg_down_proximal"]),
            ("Distal downstream", region_profile_colors["vg_down_distal"]),
        ]
    ]
    mix_handles = [
        Patch(facecolor=SOURCE_PALETTE["clingen"], label="ClinGen HI"),
        Patch(facecolor=SOURCE_PALETTE["background"], label="Background"),
        Line2D([0], [0], color="#0F172A", lw=1.1, ls="--", label=f"Overall HI fraction ({overall_hi_fraction:.1%})"),
    ]
    legend_ax.text(
        0.0,
        0.98,
        "Figure guide",
        fontsize=12,
        fontweight="bold",
        va="top",
    )
    legend_ax.legend(
        handles=radar_handles + mix_handles,
        loc="upper left",
        bbox_to_anchor=(0, 0.92),
        frameon=False,
        fontsize=9,
    )
    legend_ax.text(
        0.0,
        0.22,
        "Each card shows a cluster centroid as a radar profile\n"
        "and the corresponding Background / ClinGen HI mix below.\n"
        "The 3-column layout enlarges each cluster enough to stay\n"
        "readable at thesis page width.",
        fontsize=9,
        color="#334155",
        va="top",
    )

    save_manuscript_figure(fig, "4_3_spatial_kmeans_clusters_publication")
    plt.close(fig)


def build_af_overview(df_hi: pl.DataFrame, df_bg: pl.DataFrame) -> None:
    overview_cols = [
        ("vg_predicted", "Total Vg"),
        ("vg_common", "Common variants\n(AF ≥ 0.05)"),
        ("vg_rare", "Rare variants\n(AF < 0.05)"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 5.8), constrained_layout=True)

    rng = np.random.default_rng(42)
    for ax, (col, label) in zip(axes, overview_cols):
        hi_values = df_hi.select(col).to_pandas()[col].dropna().to_numpy()
        bg_values = df_bg.select(col).to_pandas()[col].dropna().to_numpy()
        hi_log = np.log10(clip_for_log(hi_values))
        bg_log = np.log10(clip_for_log(bg_values))
        combined_log = np.concatenate([bg_log, hi_log])

        violin = ax.violinplot([bg_log, hi_log], positions=[0, 1], showmedians=True, showextrema=False)
        for patch, color in zip(violin["bodies"], [SOURCE_PALETTE["background"], SOURCE_PALETTE["clingen"]]):
            patch.set_facecolor(color)
            patch.set_alpha(0.72)
        violin["cmedians"].set_colors(["white", "white"])
        violin["cmedians"].set_linewidth(2.5)

        for idx, (values, color) in enumerate(
            [(bg_log, SOURCE_PALETTE["background"]), (hi_log, SOURCE_PALETTE["clingen"])]
        ):
            sample_size = min(250, len(values))
            sample = rng.choice(values, sample_size, replace=False)
            jitter = rng.uniform(-0.08, 0.08, size=sample_size)
            ax.scatter(idx + jitter, sample, s=4, alpha=0.20, color=color, zorder=2)

        _, p_dep = stats.mannwhitneyu(bg_values, hi_values, alternative="greater")
        delta = float(np.median(bg_log) - np.median(hi_log))
        fold = 10 ** delta
        sig = "***" if p_dep < 0.001 else ("**" if p_dep < 0.01 else ("*" if p_dep < 0.05 else "ns"))

        if col == "vg_common":
            lower = float(np.percentile(combined_log, 1.0) - 0.35)
            upper = float(np.percentile(combined_log, 99.7) + 0.75)
            ax.text(
                0.03,
                0.03,
                "lower tail clipped for display",
                transform=ax.transAxes,
                fontsize=7.5,
                color="#475569",
                va="bottom",
            )
        else:
            lower = float(combined_log.min() - 0.25)
            upper = float(np.percentile(combined_log, 99.7) + 0.75)

        y_bracket = upper - 0.30
        ax.plot([0, 0, 1, 1], [y_bracket - 0.13, y_bracket, y_bracket, y_bracket - 0.13], color="black", lw=1.2)
        ax.text(0.5, y_bracket + 0.06, sig, ha="center", fontsize=13, fontweight="bold")
        ax.text(
            0.5,
            y_bracket + 0.32,
            f"p = {p_dep:.2e} | Δ median = {delta:+.2f} log10 ({fold:.1f}×)",
            ha="center",
            fontsize=8,
            color="#444",
        )

        ax.set_ylim(lower, upper)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Background\n(n={len(bg_values)})", f"ClinGen HI\n(n={len(hi_values)})"], fontsize=9)
        ax.set_ylabel(r"$\log_{10}(V_G)$", fontsize=10)
        ax.set_title(label, fontsize=11, fontweight="bold", pad=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.25)
        ax.grid(axis="x", visible=False)

    fig.suptitle(
        "ClinGen HI genes carry less predicted regulatory variance than background genes",
        fontsize=12.5,
        fontweight="bold",
        y=1.02,
    )
    save_manuscript_figure(fig, "4_4_vg_depletion_violin_publication")
    plt.close(fig)


def main() -> None:
    df_all, df_real, _, df_hi, df_bg = load_gene_tables()
    build_regional_boxplot(df_all)
    build_cluster_figure(df_real)
    build_af_overview(df_hi, df_bg)


if __name__ == "__main__":
    main()
