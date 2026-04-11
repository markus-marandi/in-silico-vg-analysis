import sys
from pathlib import Path
from typing import Optional, List, Tuple

import polars as pl
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

console = Console()

PROJECT_ROOT = Path(__file__).resolve().parent
EXPERIMENTS_DATA = PROJECT_ROOT / 'experiments_data'


def _latest_matching(pattern: str) -> Path:
    """Return newest match inside experiments data."""
    matches = sorted(EXPERIMENTS_DATA.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No files match pattern: {pattern}")
    return matches[-1]


VARIANT_PATHS = {
    'background': _latest_matching('dataset4_background/Background_Gnomad_variants_dedup_perm_21022026.parquet'),
    'background_null': _latest_matching('dataset5_NULL/Background_Synth_variants_downsampled_perm_21022026.parquet'),
    'clingen': _latest_matching('dataset3_ClinGen/ClinGen_HI_Gnomad_variants_dedup_21022026.parquet'),
    'clingen_null': _latest_matching('dataset5_NULL/ClinGen_HI_Synth_variants_downsampled_perm_21022026.parquet'),
}

GENE_PATHS = {
    'background': _latest_matching('dataset4_background/Background_Gnomad_genes_21022026.parquet'),
    'background_null': _latest_matching('dataset5_NULL/Background_Synth_genes_21022026.parquet'),
    'clingen': _latest_matching('dataset3_ClinGen/ClinGen_HI_Gnomad_genes_21022026.parquet'),
    'clingen_null': _latest_matching('dataset5_NULL/ClinGen_HI_Synth_genes_21022026.parquet'),
}

# ── Display / plotting config ──────────────────────────────────────────────────
# Vg spans multiple orders of magnitude; log scale is the canonical display mode.
# Toggle VG_LOG_SCALE = False to revert to linear in any notebook.
VG_LOG_SCALE: bool = True
# Floor value clipped before log-transform to avoid log(0) / -Inf artefacts.
VG_MIN_CLIP: float = 1e-12

SOURCE_PALETTE = {
    'background': '#4F46E5', 'background_null': '#1e1b4b',
    'clingen': '#10B981', 'clingen_null': '#064E3B',
}

DISPLAY_NAMES = {
    'background': 'Background',
    'background_null': 'Background\n(Synthetic)',
    'clingen': 'ClinGen HI',
    'clingen_null': 'ClinGen HI\n(Synthetic)',
}

# Simplified display names for tight plots
DISPLAY_NAMES_SHORT = {
    'background': 'Bg',
    'background_null': 'Bg-Synth',
    'clingen': 'HI',
    'clingen_null': 'HI-Synth',
}



POLYA_ASSAY_TITLE = 'polyA plus RNA-seq'


def load_variants_processed(
        dataset_name: str,
        columns: Optional[List[str]] = None,
        filter_polya: bool = True,
) -> pl.DataFrame:
    """Load variant data directly from parquet.

    args:
        dataset_name (str): key into VARIANT_PATHS
        columns (list[str] | None): subset of columns to return; core columns always included
        filter_polya (bool): when true, restrict to rows where 'Assay title' == POLYA_ASSAY_TITLE

    returns:
        pl.DataFrame: filtered variant dataframe
    """
    if dataset_name not in VARIANT_PATHS:
        raise ValueError(f"Unknown dataset: {dataset_name}")

    path = VARIANT_PATHS[dataset_name]
    df = pl.read_parquet(path)

    if filter_polya and 'Assay title' in df.columns:
        df = df.filter(pl.col('Assay title') == POLYA_ASSAY_TITLE)

    if columns:
        base = ['variant_id', 'gene_id', 'raw_score']
        target_cols = [c for c in (base + columns) if c in df.columns]
        df = df.select(target_cols)

    return df


def load_variant_pairs_matched(
        real_dataset: str,
        null_dataset: str,
        verbose: bool = True,
        filter_polya: bool = True,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Load matching real/null variant pairs with a progress bar and summary table.

    args:
        real_dataset (str): key for the real dataset in VARIANT_PATHS
        null_dataset (str): key for the null dataset in VARIANT_PATHS
        verbose (bool): print summary table after loading
        filter_polya (bool): passed through to load_variants_processed

    returns:
        tuple[pl.DataFrame, pl.DataFrame]: (real_df, null_df)
    """
    with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            console=console,
            transient=True
    ) as progress:
        task1 = progress.add_task(f"[bold blue]Loading {real_dataset}...", total=100)
        real_df = load_variants_processed(real_dataset, filter_polya=filter_polya)
        progress.update(task1, completed=100)

        task2 = progress.add_task(f"[bold magenta]Loading {null_dataset}...", total=100)
        null_df = load_variants_processed(null_dataset, filter_polya=filter_polya)
        progress.update(task2, completed=100)

    if verbose:
        _print_matching_summary(real_dataset, real_df, null_dataset, null_df)

    return real_df, null_df



def _print_matching_summary(real_name: str, real_df: pl.DataFrame, null_name: str, null_df: pl.DataFrame):
    """Prints a clean summary table of the loaded pair."""
    table = Table(title="[bold]Dataset Loading Summary", title_justify="left", border_style="bright_black")

    table.add_column("Metric", style="cyan")
    table.add_column(real_name.capitalize(), style=SOURCE_PALETTE.get(real_name, "white"))
    table.add_column(null_name.capitalize(), style=SOURCE_PALETTE.get(null_name, "white"))

    table.add_row("Total Variants", f"{real_df.height:,}", f"{null_df.height:,}")
    table.add_row("Unique Genes", f"{real_df['gene_id'].n_unique():,}", f"{null_df['gene_id'].n_unique():,}")

    # Calculate overlap
    common_genes = len(set(real_df['gene_id'].unique()) & set(null_df['gene_id'].unique()))
    table.add_row("Gene Overlap", f"{common_genes:,}", f"{common_genes:,}")

    console.print(table)
    console.print(f"[dim]Columns available: {', '.join(real_df.columns[:5])}... ({len(real_df.columns)} total)[/]\n")