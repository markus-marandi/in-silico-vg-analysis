# in-silico-vg-analysis

this repository is the analysis laboratory for data produced by the [in-silico-genetic-variance](https://github.com/markus-marandi/in-silico-genetic-variance) pipeline. it contains notebooks for quality control, benchmarking against ground truth, and spatial visualization of variant effects.

## project structure

the notebooks are organized by analysis stage to ensure a reproducible research flow:

```text
in-silico-vg-analysis/
├── config.py              # the 'golden link': central paths to the pipeline repo
├── 01_qc/                 # technical validation: distributions, nulls, and coverage
├── 02_benchmarking/       # validation: predicted vg vs ground truth (aneva-dot/gtex)
├── 03_visualisation/      # spatial analysis: tss decay, promoter vs body binning
├── 04_case_studies/       # deep dives: high-impact gene candidates
├── figures/               # high-resolution exports for thesis/publication
└── dashboard.py           # interactive streamlit dashboard for global exploration

```

## setup: the golden link

to maintain synchronization with the pipeline outputs on the hpc cluster, all notebooks import paths from `config.py`. edit this file to point to your latest experiment run.

```python
from pathlib import Path

# central hub for cluster paths
BASE_DIR = Path('/cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background')
RESULTS_DIR = BASE_DIR / '03_results'
INPUTS_DIR = BASE_DIR / '01_inputs'

def get_latest_gene_file() -> Path:
    # return the most recent gene-level summary
    return sorted(RESULTS_DIR.glob('*_genes_*.parquet'))[-1]

```

## analysis workflow

### 01_qc: technical sanity

before biological analysis, we verify the data integrity.

### 02_benchmarking: validation

comparing our predictions with observed genetic variance.

### 03_visualization: spatial effects

visualizing the spatial architecture of variant impacts.

## interactive exploration

for rapid exploratory data analysis (eda), we use `pygwalker` inside the notebooks. this allows for tableau-like interactions directly with polars dataframes.

for a global, binned overview across all parameters, run the local streamlit dashboard:

```bash
streamlit run dashboard.py

```
