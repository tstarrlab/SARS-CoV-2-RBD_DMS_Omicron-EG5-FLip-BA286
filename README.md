# Deep mutational scanning of the Omicron EG.5, HK.3 (i.e. FLip), and BA.2.86 SARS-CoV-2 RBDs

Analysis of deep mutational scanning of barcoded codon variants of SARS-CoV-2 Omicron variants for ACE2-binding affinity and RBD expression.

Study and analysis by Tyler Starr, Ashley Taylor, and co-authors.

Interactive visualizations of the key deep mutational scanning data are available [here](https://tstarrlab.github.io/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ/).

## Summary of workflow and results
For a summary of the workflow and links to key results files, [click here](results/summary/summary.md).
Reading this summary is the best way to understand the analysis.

## Running the analysis

### Computing environment

We will generally use group conda environments that are built into the shared `starr-group1/software` directory on CHPC. If you need an updated conda environment, talk to Tyler about updating or installing a new version of our common environments. This repository uses the `ysd-dms1` environment, which is specified in the [run_CHPC_cluster.bash](./run_CHPC_cluster.bash) script with the absolute path to the prebuilt conda environment under the `--conda-prefix` tag.

You can activate this environment for interactive work from a CHPC bash prompt with:

    conda activate ysd-dms1
    
To activate this environment for an interactive Jupyter Notebook session launched via CHPC's [OnDemand](https://www.chpc.utah.edu/documentation/software/ondemand.php) portal, choose "Custom (Environment Setup below)" for the Jupyter Python version option, and enter the following into the Environment Setup for Custom Python text box:
	
	module load miniconda3/latest
	conda activate ysd-dms1
    
The `snakemake` pipeline will automatically output a YAML file detailing all packages and versions loaded into the computing environment with each execution of the pipeline. This file [environment.yml](./environment.yml) (and its version history) should retain a record of the exact software versions that were used over the course of pipeline execution.

Setting up the `conda` environment above loads everything to run all parts of the analysis **except** the `R` markdown notebooks.
For those, the pipeline currently uses the CHPC computing cluster module `R/4.1.1` as specified in `Snakefile`. All the `R` packages are listed at the beginning of their output in the [summary results](results/summary/summary.md). Similar to Jupyter Notebooks above, interactive RStudio sessions to work on individual scripts can be run launched via CHPC's [OnDemand](https://www.chpc.utah.edu/documentation/software/ondemand.php) portal.

Now you can run the entire analysis.


The analysis consists primarily of a series of Jupyter notebooks and R markdown in to the top-level directory along with some additional code in [Snakefile](Snakefile).
You can run the analysis by using [Snakemake](https://snakemake.readthedocs.io) to run [Snakefile](Snakefile), specifying the conda environment in `./env`, as in:

    snakemake --use-conda --conda-prefix /uufs/chpc.utah.edu/common/home/starr-group1/software/pkg/miniconda3/envs/ysd-dms1

However, you probably want to use the server to help with computationally intensive parts of the analysis.
To run using the cluster configuration for the CHPC server, simply run the bash script [run_CHPC_cluster.bash](run_CHPC_cluster.bash), which executes [Snakefile](Snakefile) in a way that takes advantage of the CHPC server resources.

You likely want to submit [run_CHPC_cluster.bash](run_CHPC_cluster.bash) itself to the cluster (since it takes a while to run) with:

    sbatch -t 0-2 run_CHPC_cluster.bash -A starr-np -p starr-shared-np

You can leave out the -A and -p flags and run the script `run_CHPC_cluster_shared.bash` if you want to run on general shared nodes instead of our starr-lab personal compute nodes.

### Configure `.git` to not track Jupyter notebook metadata
To simplify git tracking of Jupyter notebooks, we have added the filter described [here](https://stackoverflow.com/questions/28908319/how-to-clear-an-ipython-notebooks-output-in-all-cells-from-the-linux-terminal/58004619#58004619) to strip notebook metadata to [.gitattributes](.gitattributes) and [.gitconfig](.gitconfig).
The **first time** you check out this repo, run the following command to use this configuration (see [here](https://stackoverflow.com/a/18330114)):
```
   git config --local include.path ../.gitconfig
```
Then don't worry about it anymore.

### Configuring the analysis
The configuration for the analysis is specifed in [config.yaml](config.yaml).
This file defines key variables for the analysis, and should be relatively self-explanatory.

Whenever possible, you should modify the analysis by changing this configuration file instead of hard-coding important experiment-specific details into the Python/R scripts or Snakefile.

The input files pointed to by [config.yaml](config.yaml) are in the [./data/](data) subdirectory.
See the [./data/README.md](./data/README.md) file for details.


### Cluster configuration
There is a cluster configuration file [cluster.yaml](cluster.yaml) that configures [Snakefile](Snakefile) for the CHPC cluster, as recommended by the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
The [run_CHPC_cluster.bash](run_CHPC_cluster.bash) script uses this configuration to run [Snakefile](Snakefile). If you add new `snakemake` rules that require different computing requirements than the default given in this [cluster.yaml](cluster.yaml) file, add a new element to this YAML file.
If you are using a different cluster than CHPC, you may need to modify the cluster configuration file.

### Notebooks that perform the analysis
The Jupyter notebooks and R markdown dscripts that perform most of the analysis are in this top-level directory with the extension `*.ipynb` or `*.Rmd`.
These notebooks read the key configuration values from [config.yaml](config.yaml).

There is also a [./scripts/](scripts) subdirectory with related scripts, as needed.

The notebooks need to be run in the order described in [the workflow and results summary](results/summary/summary.md).
This will occur automatically if you run them via [Snakefile](Snakefile) as described above.

### Results
Results are placed in the [./results/](results) subdirectory.
Many of the files created in this subdirectory are not tracked in the `git` repo as they are very large.
However, key results files are tracked as well as a summary that shows the code and results.
Click [here](./results/summary/summary.md) to see that summary.

The large results files are tracked via [git-lfs](https://git-lfs.github.com/).
This requires `git-lfs` to be installed, which it is in our group `conda` environments.
The following commands were then run to activate `git-lfs`:

    git lfs install

You may need to run this if you are tracking these files and haven't installed `git-lfs` in your user account.
Then the large results files were added for tracking with:
```
git lfs track results/variants/codon_variant_table_*.csv
git lfs track results/counts/variant_counts.csv
```
