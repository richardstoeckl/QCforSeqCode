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
3. Download checkm2 database (via `wget https://zenodo.org/records/14897628/files/checkm2_database.tar.gz`)
4. Download GTDB-Tk database (via `wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz`)
3. [Download the latest release from this repo](https://github.com/richardstoeckl/QCforSeqCode/releases/latest) and cd into it
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

## Tools used in the pipeline and reasoning. Please cite these tools if you use this pipeline.
- **Taxonomy**
    - **[GTDB-Tk v2.4.1](https://github.com/Ecogenomics/GTDBTk/) - toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes.** *Used to get full genome taxonomic classification.*
    - **[Infernal v1.1.5](https://github.com/EddyRivasLab/infernal) - RNA secondary structure/sequence profiles for homology search and alignment.** *Used to find and extract rRNA genes in the genomes.*
    - **[DECIPHER v3.2.0](https://doi.org/doi:10.18129/B9.bioc.DECIPHER) - Tools for curating, analyzing, and manipulating biological sequences.** *Used to get 16S rRNA gene taxonomic classification by comparing to SILVA db.*
    - **[SILVA r138](https://www.arb-silva.de/) - rRNA database.** *Used as source of rRNA gene taxonomy*
- **Contamination and Completeness**
    - **[CheckM2 v1.1.0](https://github.com/chklovski/CheckM2/) - Assessing the quality of metagenome-derived genome bins using machine learning.** *Used to get completeness and contamination stats.* Unlike CheckM1 (one of the most popular tools for completeness and contamination prediction), CheckM2 has universally trained machine learning models it applies regardless of taxonomic lineage. This allows it to work better with organisms that have only few known representative genomes.
- **tRNA gene occurence**
    - **[tRNAscan-SE v2.0.12](https://github.com/UCSC-LoweLab/tRNAscan-SE) - An improved tool for transfer RNA detection.** *Used to find tRNA genes in the genomes.*
- **General stats, file manipulation, alignment, and reporting**
    - **[seqkit v2.10.0](https://github.com/shenwei356/seqkit) - ultrafast toolkit for FASTA/Q file manipulation.** *Used for quick and easy general stat gathering and sequence concatination.*
    - **[minimap2 v2.29](https://github.com/lh3/minimap2) - versatile pairwise aligner for genomic and spliced nucleotide sequences.** *Used to align sequencing reads to assembly to get coverage stats.*
    - **[samtools v1.21](https://github.com/samtools/samtools) - Tools for manipulating next-generation sequencing data** *Used to calculate coverage stats.*
    - **[tidyverse v2.0.0](https://github.com/tidyverse) - R packages for data science** *Used for general data manipulation for reporting*
    - **[fs v1.6.6](https://github.com/r-lib/fs/) - cross platform file operations** *Used for file manipulation for reporting*
    - **[tinytable v0.8.0](https://github.com/vincentarelbundock/tinytable) - Simple and Customizable Tables** *Used to generate the final report*


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
