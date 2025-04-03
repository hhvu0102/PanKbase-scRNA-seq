# PanKbase single cell map generation
This directory contains Code and docs for PanKbase islet single cell map analysis.

## Content
0. `0_preprocessing`: <br>
This directory contains Jupyter notebooks with details about how to execute code, along with scripts (R/Python/Snakemake). <br>
	`scripts`: directory contains scripts to run analyses <br>
	`0_collecting_data.ipynb`: Jupyter notebook with instructions on how data was downloaded <br>
	`1_run_processing_pipeline.ipynb`: Jupyter notebook with instructions on how to run Nextflow pipeline to process scRNA-seq data <br>
	`2_fastQC_inspection.R.ipynb`: Jupyter notebook with instructions on how Per Tile Quality plots were analyzed <br>
	`3_barcode_qc.ipynb`: Jupyter notebook with instructions on how to pull together multiple QC metrics, generate QC plots and choosing barcodes that likely correspond to cells.
	`python_req.yml`: File with conda environment requirements to generate plots like in `3_barcode_qc.ipynb`.
