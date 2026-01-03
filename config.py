from pathlib import Path

# central hub for local project paths
PROJECT_ROOT = Path(__file__).resolve().parent
EXPERIMENTS_DATA = PROJECT_ROOT / 'experiments_data'


def _latest_matching(pattern: str) -> Path:
    """return newest match inside experiments data."""
    matches = sorted(EXPERIMENTS_DATA.glob(pattern))
    if not matches:
        raise FileNotFoundError(f'no files match pattern: {pattern}')
    return matches[-1]


CLINGEN_VAR_PATH = _latest_matching('dataset3_ClinGen/*variant_level_summary.parquet')
BACKGROUND_VAR_PATH = _latest_matching('dataset4_background/background_variants_*.parquet')
BACKGROUND_NULL_VAR_PATH = _latest_matching('dataset5_NULL/*Background_NULL_variant_level_summary.parquet')
CLINGEN_NULL_VAR_PATH = _latest_matching('dataset5_NULL/*ClinGen_NULL_variant_level_summary.parquet')

VARIANT_PATHS = {
    'background': BACKGROUND_VAR_PATH,
    'background_null': BACKGROUND_NULL_VAR_PATH,
    'clingen': CLINGEN_VAR_PATH,
    'clingen_null': CLINGEN_NULL_VAR_PATH,
}

CLINGEN_GENE_PATH = _latest_matching('dataset3_ClinGen/clingen_genes_*.parquet')
BACKGROUND_GENE_PATH = _latest_matching('dataset4_background/background_genes_*.parquet')
BACKGROUND_NULL_GENE_PATH = _latest_matching('dataset5_NULL/dataset5_Background_NULL_gene_level_summary.parquet')
CLINGEN_NULL_GENE_PATH = _latest_matching('dataset5_NULL/dataset5_ClinGen_NULL_gene_level_summary.parquet')

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