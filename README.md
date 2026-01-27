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