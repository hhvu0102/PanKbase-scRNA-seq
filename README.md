# PanKbase single cell map generation
This directory contains Code and docs for PanKbase islet single cell map analysis.

## Content
1. `1_preprocessing`: <br>
This directory contains Jupyter notebooks with details about how to execute code, along with scripts (R/Python/Snakemake). <br>
- `scripts`: directory contains scripts to run analyses. <br>
- `0_collecting_data.ipynb`: Jupyter notebook with instructions on how data was downloaded. <br>
- `1_run_processing_pipeline.ipynb`: Jupyter notebook with instructions on how to run Nextflow pipeline to process scRNA-seq data. <br>
- `2_barcode_qc.ipynb`: Jupyter notebook with instructions on how to pull together multiple QC metrics, generate QC plots and choosing barcodes that likely correspond to cells. <br>
- `2_python_req.yml`: File with conda environment requirements to generate plots like in `3_barcode_qc.ipynb`.
- `3_run_modifed_cellbender.R.ipynb`: Jupyter notebook with instructions on how to set up pipeline to run Cellbender with modified settings. <br>
- `4_HTODemux.R.ipynb`: Jupyter notebook with instructions on how to analyze HTO data. <br>

2. `2_doublet_detection/`:
This directory contains Jupyter notebooks with details about how to execute code, along with scripts (R/Python/Snakemake). <br>
- `scripts`: directory contains scripts to run analyses <br>
- `1_run_doubletfinder.R.ipynb`: Jupyter notebook with instructions on how to run doublet detection pipeline <br>
- `2_doublet_rates_per_source.R.ipynb`: Jupyter notebook with instructions on how to obtain doublet rate per cluster. <br>

3. `3_integration/`:
This directory contains Jupyter notebooks with details about how to execute code, along with scripts (R/Python/Snakemake). <br>
- `scripts`: directory contains scripts to run analyses. <br>
- `1_getMetrics.ipynb`: Jupyter notebook with instructions on how to pull together multiple QC metrics. This script is different from `1_preprocessing/2_barcode_qc.ipynb` as it obtains metrics from CellBender modified setting runs, while at the earlier steps, metrics come from CellBender with default settings runs. <br>
- `2_integration.R.ipynb`: Jupyter notebook with instructions on how to integrate data across multiple tissue sources, how to remove outlier clusters, and annotate cell populations.

## Contact
Ha Vu <vthihong at umich.edu> and Stephen Parker <scjp at umich.edu>
