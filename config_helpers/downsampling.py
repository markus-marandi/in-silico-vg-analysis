"""downsampling utilities for null variant data."""
import pandas as pd
import polars as pl


def downsample_null_to_real(
    null_df: pl.DataFrame,
    real_counts_df: pl.DataFrame,
    seed: int = 42,
) -> pl.DataFrame:
    """downsample null variants to match real variant counts per gene.
    
    args:
        null_df (pl.DataFrame): null variant dataframe with gene_id column.
        real_counts_df (pl.DataFrame): real gene counts with gene_id and n_variants columns.
        seed (int): random seed for reproducibility.
    
    returns:
        pl.DataFrame: downsampled null variants matching real counts per gene.
    """
    n_map = (
        real_counts_df
        .select(['gene_id', 'n_variants'])
        .to_pandas()
        .set_index('gene_id')['n_variants']
        .to_dict()
    )
    
    pdf = null_df.to_pandas()
    sampled_parts = []
    
    for gene, group in pdf.groupby('gene_id'):
        n_req = int(n_map.get(gene, 0))
        if n_req <= 0:
            continue
        
        avail = len(group)
        if avail == 0:
            continue
        
        if avail <= n_req:
            # take all available
            sampled_parts.append(group)
        else:
            # sample without replacement
            taken = group.sample(n=n_req, random_state=seed)
            sampled_parts.append(taken)
    
    if sampled_parts:
        sampled_pdf = pd.concat(sampled_parts, axis=0)
        return pl.from_pandas(sampled_pdf)
    else:
        return pl.DataFrame(schema=null_df.schema)
