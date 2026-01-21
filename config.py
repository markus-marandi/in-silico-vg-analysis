import sys
from pathlib import Path
from typing import Optional, List, Tuple

import polars as pl
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

# Initialize Rich console for beautiful logging
console = Console()

# --- Path Resolution ---
PROJECT_ROOT = Path(__file__).resolve().parent
EXPERIMENTS_DATA = PROJECT_ROOT / 'experiments_data'

# Import local utilities
from config_helpers import dedup_scores_by_variant, downsample_null_to_real

def _latest_matching(pattern: str) -> Path:
    """Return newest match inside experiments data."""
    matches = sorted(EXPERIMENTS_DATA.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No files match pattern: {pattern}")
    return matches[-1]

# --- Dataset Registries ---
VARIANT_PATHS = {
    'background': _latest_matching('dataset4_background/background_variants_20260119.parquet'),
    'background_null': _latest_matching('dataset5_NULL/*Background_NULL_variant_level_summary_1901.parquet'),
    'clingen': _latest_matching('dataset3_ClinGen/*variant_level_summary.parquet'),
    'clingen_null': _latest_matching('dataset5_NULL/*ClinGen_NULL_variant_level_summary_1901.parquet'),
}

GENE_PATHS = {
    'background': _latest_matching('dataset4_background/background_genes_20260102.parquet'),
    'background_null': _latest_matching('dataset5_NULL/dataset5_Background_NULL_gene_level_summary_1901.parquet'),
    'clingen': _latest_matching('dataset3_ClinGen/clingen_genes_20260102.parquet'),
    'clingen_null': _latest_matching('dataset5_NULL/dataset5_ClinGen_NULL_gene_level_summary_1901.parquet'),
}

SOURCE_PALETTE = {
    'background': '#4F46E5', 'background_null': '#1e1b4b',
    'clingen': '#10B981', 'clingen_null': '#064E3B',
}

# --- Core Loading Functions ---

def load_variants_processed(
    dataset_name: str,
    verbose: bool = True,
    columns: Optional[List[str]] = None,
) -> pl.DataFrame:
    """Load and deduplicate variant data, preserving schema flexibility."""
    if dataset_name not in VARIANT_PATHS:
        raise ValueError(f"Unknown dataset: {dataset_name}")
    
    path = VARIANT_PATHS[dataset_name]
    
    # Use our updated dedup function that preserves all columns
    df = dedup_scores_by_variant(
        path=path,
        label=dataset_name if verbose else None,
        verbose=False # We handle our own logging here
    )
    
    if columns:
        # Ensure core columns are always present
        base = ['variant_id', 'gene_id', 'raw_score']
        target_cols = [c for c in (base + columns) if c in df.columns]
        df = df.select(target_cols)
    
    return df

def load_variant_pairs_matched(
    real_dataset: str,
    null_dataset: str,
    downsample: bool = True,
    seed: int = 42,
    verbose: bool = True,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Load matching pairs with a visual progress bar and summary table."""
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
        transient=True
    ) as progress:
        
        # 1. Load Real Dataset
        task1 = progress.add_task(f"[bold blue]Loading {real_dataset}...", total=100)
        real_df = load_variants_processed(real_dataset, verbose=False)
        progress.update(task1, completed=100)

        # 2. Load Null Dataset
        task2 = progress.add_task(f"[bold magenta]Loading {null_dataset}...", total=100)
        null_df = load_variants_processed(null_dataset, verbose=False)
        progress.update(task2, completed=100)

        # 3. Downsampling
        if downsample:
            task3 = progress.add_task("[bold yellow]Downsampling null set...", total=100)
            gene_key = real_dataset.replace('_null', '')
            real_gene_counts = pl.read_parquet(GENE_PATHS[gene_key])
            
            null_final = downsample_null_to_real(
                null_df=null_df,
                real_counts_df=real_gene_counts,
                seed=seed
            )
            progress.update(task3, completed=100)
        else:
            null_final = null_df

    if verbose:
        _print_matching_summary(real_dataset, real_df, null_dataset, null_final)

    return real_df, null_final

# --- Reporting Utilities ---

def _print_matching_summary(real_name: str, real_df: pl.DataFrame, null_name: str, null_df: pl.DataFrame):
    """Prints a clean summary table of the loaded pair."""
    table = Table(title="[bold]Dataset Matching Summary", title_justify="left", border_style="bright_black")
    
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