from pathlib import Path

# central hub for cluster paths
BASE_DIR = Path('/cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background')
RESULTS_DIR = BASE_DIR / '03_results'
INPUTS_DIR = BASE_DIR / '01_inputs'

def get_latest_gene_file() -> Path:
    # return the most recent gene-level summary
    return sorted(RESULTS_DIR.glob('*_genes_*.parquet'))[-1]