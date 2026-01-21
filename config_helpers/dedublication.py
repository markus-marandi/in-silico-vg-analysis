"""deduplication utilities for variant data."""
from pathlib import Path

import polars as pl


def dedup_scores_by_variant(
    path: str | Path,
    label: str | None = None,
    columns: list[str] | None = None,
    verbose: bool = True,
) -> pl.DataFrame:
    """deduplicate scores by variant, keeping max absolute score.
    
    args:
        path (str | Path): parquet file path.
        label (str | None): dataset label for logging.
        columns (list[str] | None): additional columns to keep (e.g. ['AF'], ['perm_AF']).
        verbose (bool): print loading and dedup stats.
    
    returns:
        pl.DataFrame: deduplicated variant table with gene_id, variant_id, raw_score, and requested columns.
    """
    df = pl.read_parquet(path)
    
    if verbose and label:
        print(f'  loaded {label}: {df.shape}')
    
    # detect score column name
    schema_cols = df.columns
    if 'raw_score' in schema_cols:
        score_col = 'raw_score'
    elif 'score' in schema_cols:
        score_col = 'score'
    else:
        raise ValueError(f'missing raw_score/score column in {path}')
    
    # required columns
    required = ['variant_id', 'gene_id', score_col]
    
    # add optional columns (only if they exist in the file)
    if columns:
        for col in columns:
            if col in schema_cols:
                required.append(col)
            elif verbose and label:
                print(f'  warning: column {col} not found in {label}, skipping')
    
    # select needed columns and rename score column to raw_score
    df = df.select(required)
    if score_col != 'raw_score':
        df = df.rename({score_col: 'raw_score'})
    
    # filter out null raw_scores before deduplication
    df = df.filter(pl.col('raw_score').is_not_null())
    
    # deduplicate by keeping variant with max absolute score
    dedup = (
        df
        .with_columns(pl.col('raw_score').abs().alias('_abs_score'))
        .sort('_abs_score', descending=True)
        .unique(subset=['variant_id'], keep='first')
        .drop('_abs_score')
    )
    
    if verbose and label:
        print(f'  after dedup: {dedup.shape}')
    
    return dedup