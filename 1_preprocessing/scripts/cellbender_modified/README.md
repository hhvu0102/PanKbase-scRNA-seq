# NextFlow pipeline for CellBender modified settings

## Dependencies
Singularity (v. 3) and NextFlow (>= v. 20.10.0). Containers with the software for each step are pulled from the Sylabs cloud library (https://cloud.sylabs.io/library).


## Configuration

When launching the pipeline, as shown in the `nextflow` command below, you'll also need to set the following:

1. The location of the results directory (e.g., `--results /path/to/results`)
2. A library config file which includes information about each library. Each column in this file is: 
* `sample`: often in the format `{library}-{genome}` to match with sample names from upstream steps
* `expected_cells`: number of cells expected in the dataset (see more at https://cellbender.readthedocs.io/en/latest/reference/index.html)
* `source`: (optional) tissue source
* `learning_rate`: learning rate
* `solo_out`: path to STARsolo ouput as obtained upstream
* `total_droplets_included`: the number of droplets from the rank-ordered UMI plot that will have their cell probabilities inferred as an output. 

## Running
Once you have all of the above information, you can run the pipeline as follows (in this case, indicating the path to the results on the command line):

```bash
sbatch --job-name=cellbender --mem=500M --time=24:00:00 --account=scjp99 --mail-user=vthihong@umich.edu --mail-type=END,FAIL --signal=B:TERM@60 --wrap="exec ~/tools/nextflow run -resume --library_info library_info.tsv --results /path/to/results/cellbender /path/to/main.nf"
```

