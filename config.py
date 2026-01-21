from pathlib import Path

import polars as pl

# central hub for local project paths
PROJECT_ROOT = Path(__file__).resolve().parent
EXPERIMENTS_DATA = PROJECT_ROOT / 'experiments_data'

# import preprocessing utilities
from config_helpers import dedup_scores_by_variant, downsample_null_to_real


def _latest_matching(pattern: str) -> Path:
    """return newest match inside experiments data."""
    matches = sorted(EXPERIMENTS_DATA.glob(pattern))
    if not matches:
        raise FileNotFoundError(f'no files match pattern: {pattern}')
    return matches[-1]


CLINGEN_VAR_PATH = _latest_matching('dataset3_ClinGen/*variant_level_summary.parquet')
BACKGROUND_VAR_PATH = _latest_matching('dataset4_background/background_variants_20260119.parquet')
BACKGROUND_NULL_VAR_PATH = _latest_matching(
    'dataset5_NULL/*Background_NULL_variant_level_summary_1901.parquet'
)
CLINGEN_NULL_VAR_PATH = _latest_matching(
    'dataset5_NULL/*ClinGen_NULL_variant_level_summary_1901.parquet'
)

VARIANT_PATHS = {
    'background': BACKGROUND_VAR_PATH,
    'background_null': BACKGROUND_NULL_VAR_PATH,
    'clingen': CLINGEN_VAR_PATH,
    'clingen_null': CLINGEN_NULL_VAR_PATH,
}

CLINGEN_GENE_PATH = _latest_matching('dataset3_ClinGen/clingen_genes_20260102.parquet')
BACKGROUND_GENE_PATH = _latest_matching('dataset4_background/background_genes_20260102.parquet')
BACKGROUND_NULL_GENE_PATH = _latest_matching(
    'dataset5_NULL/dataset5_Background_NULL_gene_level_summary_1901.parquet'
)
CLINGEN_NULL_GENE_PATH = _latest_matching(
    'dataset5_NULL/dataset5_ClinGen_NULL_gene_level_summary_1901.parquet'
)

GENE_PATHS = {
    'background': BACKGROUND_GENE_PATH,
    'background_null': BACKGROUND_NULL_GENE_PATH,
    'clingen': CLINGEN_GENE_PATH,
    'clingen_null': CLINGEN_NULL_GENE_PATH,
}

# standard palette for dataset sources
SOURCE_PALETTE = {
    'background': '#4F46E5',       # Indigo (Vibrant Blue)
    'background_null': '#1e1b4b',  # Midnight (Deep Navy)
    'clingen': '#10B981',          # Emerald (Vibrant Green)
    'clingen_null': '#064E3B',     # Jungle (Deep Green)
}


def load_variants_processed(
    dataset_name: str,
    af_column: str | None = None,
    verbose: bool = True,
) -> pl.DataFrame:
    """load and deduplicate variant data.
    
    args:
        dataset_name (str): one of 'background', 'background_null', 'clingen', 'clingen_null'.
        af_column (str | None): af column to keep ('AF' for real, 'perm_AF' for null, None for none).
        verbose (bool): print loading stats.
    
    returns:
        pl.DataFrame: deduplicated variants with gene_id, variant_id, raw_score, and af column if specified.
    """
    if dataset_name not in VARIANT_PATHS:
        raise ValueError(f'unknown dataset: {dataset_name}. choose from {list(VARIANT_PATHS.keys())}')
    
    path = VARIANT_PATHS[dataset_name]
    columns = [af_column] if af_column else None
    
    return dedup_scores_by_variant(
        path=path,
        label=dataset_name if verbose else None,
        columns=columns,
        verbose=verbose,
    )


def load_variant_pairs_matched(
    real_dataset: str,
    null_dataset: str,
    downsample: bool = True,
    seed: int = 42,
    verbose: bool = True,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """load and process matching real and null variant pairs.
    
    args:
        real_dataset (str): real dataset name ('background' or 'clingen').
        null_dataset (str): null dataset name ('background_null' or 'clingen_null').
        downsample (bool): downsample null to match real counts per gene.
        seed (int): random seed for downsampling.
        verbose (bool): print processing stats.
    
    returns:
        tuple[pl.DataFrame, pl.DataFrame]: (real_dedup, null_sampled) dataframes.
    """
    # infer af columns
    real_af_col = 'AF'
    null_af_col = 'perm_AF'
    
    # load and deduplicate real variants
    if verbose:
        print(f'loading {real_dataset}...')
    real_dedup = load_variants_processed(
        dataset_name=real_dataset,
        af_column=real_af_col,
        verbose=verbose,
    )
    
    # load and deduplicate null variants
    if verbose:
        print(f'loading {null_dataset}...')
    null_dedup = load_variants_processed(
        dataset_name=null_dataset,
        af_column=null_af_col,
        verbose=verbose,
    )
    
    # downsample null to match real if requested
    if downsample:
        if verbose:
            print(f'downsampling {null_dataset} to match {real_dataset} counts...')
        
        # load gene counts
        gene_key = real_dataset.replace('_null', '')
        real_gene = pl.read_parquet(GENE_PATHS[gene_key])
        
        null_sampled = downsample_null_to_real(
            null_df=null_dedup,
            real_counts_df=real_gene,
            seed=seed,
        )
        
        if verbose:
            print(f'  {null_dataset} after downsampling: {null_sampled.shape}')
        
        return real_dedup, null_sampled
    else:
        return real_dedup, null_dedup