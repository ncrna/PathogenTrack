[![Build Status](https://app.travis-ci.com/ncrna/PathogenTrack.svg?branch=master)](https://app.travis-ci.com/ncrna/PathogenTrack)
[![PYPI](https://img.shields.io/pypi/v/pathogentrack.svg)](https://pypi.org/project/pathogentrack/)
[![The MIT License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/ncrna/PathogenTrack/blob/master/LICENSE)
# PathogenTrack
PathogenTrack is an unsupervised computational software that uses `unmapped single-cell RNAseq reads` to characterize `intracellular pathogens` at the single-cell level. It is a python-based script that can be used to identify and quantify intracellular pathogenic `viruses` and `bacteria` reads at the single-cell level.
PathogenTrack has been tested on various scRNA-seq datasets derived from simulated and real datasets and performed robustly. The detailes are described in our paper *`Decoding Intracellular Pathogens of scRNA-seq experiments with PathogenTrack and SCKIT`*.

### System Requirements

PathogenTrack has been tested on Linux platform with CentOS 7 operation system. The RAM is 120 GB, with 40 computational threads.

## Installation

### PathogenTrack can be installed in two steps:

1 . Installing Miniconda on Linux Platform. For details, please refer to [Miniconda Installation](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent).
```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

2 . Installing PathogenTrack.
```sh
conda env create -f environment.yml
```
Users can install the dependencies manually. The dependencies and test versions are listed below.

Package|Version
--|:--:
python|3.6.10
biopython|1.78
star|2.7.5a
umi_tools|1.1.1
kraken2|2.1.1

## Databases Preparation

### 1. Prepare the Human genome database
Download the Human GRCh38 genome and genome annotation file, and then decompress them:
```sh
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
gzip -d Homo_sapiens.GRCh38.101.gtf.gz
```

Build STAR Index with the following command:
```sh
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ./ \
     --genomeFastaFiles ./Homo_sapiens.GRCh38.dna.toplevel.fa \
     --sjdbGTFfile ./Homo_sapiens.GRCh38.101.gtf \
     --sjdbOverhang 100
```

### 2. Prepare Kraken2 database

```sh
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
tar zxf minikraken_8GB_202003.tgz
```

## How to use PathogenTrack?
Before running this tutorial, you should run `cellranger` or `alevin` to get the single cells' gene expression matrix. Here, we take the simulated 10X sequencing data as an example:

First, we use cellranger to get scRNA-seq expression matrix and valid barcodes:
```sh
cellranger count --id cellranger_out --transcriptom /path/to/cellranger_database/
```

Then we run PathogenTrack to identify and quantify pathogen expression at the single-cell level:
```sh
conda activate PathogenTrack
python PathogenTrack.py count --project_id PathogenTrack_out --pattern CCCCCCCCCCCCCCCCNNNNNNNNNN \
                              --min_reads 10 --confidence 0.11 --star_index ~/database/STAR_index/ \
                              --kraken_db ~/database/minikraken_8GB_20200312/ --barcode barcodes.tsv \
                              --read1 simulation_S1_L001_R1_001.fastq.gz \
                              --read2 simulation_S1_L001_R2_001.fastq.gz 
```
**IMPORTANT**: The Read 1 in the example is made up of 16 bp CB and 10 bp UMI, so the --pattern is *CCCCCCCCCCCCCCCCNNNNNNNNNN* (16C and 10N). Users must adjust the pattern with their own Read 1 accordingly.

*Note:* It may take 4-6 hours to complete one sample, and it depends on the performance of computational resources and the size of the raw single-cell data.

### Please see [QUICK_START.md](https://github.com/ncrna/PathogenTrack/blob/master/doc/QUICK_START.md) for a complete tutorial.

## Questions

For questions and suggestions about the pipeline or the code, please contact [admin@ncrna.net](mailto:admin@ncrna.net) and [ty12260@rjh.com.cn](mailto:ty12260@rjh.com.cn). We will try our best to provide support, address new issues, and keep improving this software.
