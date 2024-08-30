# Usage and configuration

Here is a rough overview:
1. Install [conda](https://docs.conda.io/en/latest/miniconda.html) (mamba or miniconda is fine).
2. Install snakemake with:
```bash
conda install -c conda-forge -c bioconda snakemake
```
3. Download checkm2 database (via `wget https://zenodo.org/api/files/fd3bc532-cd84-4907-b078-2e05a1e46803/checkm2_database.tar.gz`)
4. Download GTDB-Tk database (via `wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz`)
3. [Download the latest release from this repo](https://github.com/richardstoeckl/basecallNanopore/releases/latest) and cd into it
4. Edit the `config/config.yaml` to provide the paths to your results/logs directories, and the paths to the databases you downloaded, as well as any parameters you might want to change.
5. Edit the `config/sampleData.csv` file with the specific details for each assembly you want to check. Depending on what you enter here, the pipeline will automatically adjust what will be done.

---

# General configuration

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

## "Main" section

Here you should provide the paths to your intermediary/results/logs directories. The `interim` directory will contain larger intermediary files. The `results` directory will contain the final output of the pipeline. The `log`directory will be used to store the log files for each step.
Here you should also write the name of your sample data file (see [relevant section below](#sampleData-file-setup)).

## "Tools" section

Here you should give the paths to the databases needed for some of the tools.


# sampleData file setup

The setup of the samples is specified via comma-separated values files (`.csv`).
You can use the `config/sampleData.csv`file as a template.
