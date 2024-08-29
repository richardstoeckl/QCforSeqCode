# Snakemake workflow: `QCforSeqCode`

Author: richard.stoeckl@ur.de

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.10.0-brightgreen.svg)](https://snakemake.github.io)

## About
[Snakemake](https://snakemake.github.io) Pipeline to check the requirements for a prokaryotic assembly to be included in the [SeqCode](https://registry.seqco.de/) initiative.

The requirements are outlined in [APPENDIX I](https://registry.seqco.de/page/seqcode#data-quality-necessary-for-completion-of-seqcode-registryb) of the SeqCode.

## Usage

**[Check out the usage instructions in the snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=richardstoeckl/QCforSeqCode)**

But here is a rough overview:
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
5. Open a terminal in the main dir and start a dry-run of the pipeline with the following command. This will download and install all the dependencies for the pipeline (this step takes may take some time) and it will show you if you set up the paths correctly:

```bash
snakemake --sdm conda -n --cores
```
6. Run the pipeline with
```bash
snakemake --sdm conda --cores
```
---

## TODO and planned features
- add 16S rRNA gene truncation check
- add automatic switches for Kingdom specific modes of some tools
- automate checkm2 and gtdb-tk database downloads
- add checks if the config file and the sample file are correctly filled


## Notes on the Test data:
- `data/GCF_000007305.1_ASM730v1_genomic.fna` - This is the reference genome of Pyrococcus furiosus, which does fit the criteria of SeqCode. It was acquired from the [RefSeq database](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000007305.1/).
- `data/GCA_015662175.1_ASM1566217v1_genomic.fna` - This is the assembly of Thermococcus paralvinellae, which does not fit the criteria of SeqCode. It was acquired from [GenBank database](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_015662175.1/)
- `data/SRR8767914_subsampled.fastq.gz` is a [DNA-Seq of Pyrococcus furiosus DSM 3638](https://www.ncbi.nlm.nih.gov/sra/SRR8767914) dataset, that was subsampled for quicker testing via `zcat SRR8767914.fastq.gz | seqkit sample --rand-seed 42 -p 0.1 -o SRR8767914_subsampled.fastq.gz`.

```
Copyright Richard Stöckl 2024.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE or copy at 
https://www.boost.org/LICENSE_1_0.txt)
```