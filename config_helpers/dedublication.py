"""deduplication utilities for variant data."""
from pathlib import Path

import polars as pl


def dedup_scores_by_variant(
    path: str | Path,
    label: str | None = None,
    columns: list[str] | None = None, # Kept for compatibility, though redundant now
    verbose: bool = True,
) -> pl.DataFrame:
    """deduplicate scores by variant, keeping all original columns."""
    df = pl.read_parquet(path)
    
    if verbose and label:
        print(f'  loaded {label}: {df.shape}')
    
    # 1. Detect score column name
    schema_cols = df.columns
    if 'raw_score' in schema_cols:
        score_col = 'raw_score'
    elif 'score' in schema_cols:
        score_col = 'score'
    else:
        raise ValueError(f'missing raw_score/score column in {path}')
    
    # 2. Rename score column to 'raw_score' if needed (keeping all other columns)
    if score_col != 'raw_score':
        df = df.rename({score_col: 'raw_score'})
    
    # 3. Filter out null raw_scores
    df = df.filter(pl.col('raw_score').is_not_null())
    
    # 4. Deduplicate: keep the row with the maximum absolute raw_score
    # unique(subset=['variant_id']) will keep all original columns for the chosen row
    dedup = (
        df
        .with_columns(pl.col('raw_score').abs().alias('_abs_score'))
        .sort('_abs_score', descending=True)
        .unique(subset=['variant_id'], keep='first', maintain_order=False)
        .drop('_abs_score')
    )
    
    if verbose and label:
        print(f'  after dedup: {dedup.shape}')
        print(f'  columns kept: {len(dedup.columns)}')
    
    return dedup