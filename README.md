[![Build Status](https://app.travis-ci.com/ncrna/PathogenTrack.svg?branch=master)](https://app.travis-ci.com/ncrna/PathogenTrack)
[![The MIT License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/ncrna/PathogenTrack/blob/master/LICENSE)
[![PYPI](https://img.shields.io/pypi/v/pathogentrack.svg)](https://pypi.org/project/pathogentrack/)
[![Conda](https://anaconda.org/bioconda/pathogentrack/badges/installer/conda.svg)](https://anaconda.org/bioconda/pathogentrack)
[![Conda](https://anaconda.org/bioconda/pathogentrack/badges/version.svg)](https://anaconda.org/bioconda/pathogentrack)
[![Conda Downloads](https://anaconda.org/bioconda/pathogentrack/badges/downloads.svg)](https://anaconda.org/bioconda/pathogentrack)
[![Platform](https://img.shields.io/badge/platform-any-ec2eb4.svg)](https://github.com/ncrna/PathogenTrack)
[![check in Biotreasury](https://img.shields.io/badge/Biotreasury-collected-brightgreen)](https://biotreasury.rjmart.cn/#/tool?id=20721)

# PathogenTrack <img src="https://github.com/ncrna/PathogenTrack/blob/master/doc/PathogenTrack_logo.png" align="right" height=142 width=132/>
PathogenTrack is an unsupervised computational software that uses `unmapped single-cell RNAseq reads` to characterize `intracellular pathogens` at the single-cell level. It is a python-based script that can be used to identify and quantify intracellular pathogenic `viruses` and `bacteria` reads at the single-cell level.
PathogenTrack has been tested on various scRNA-seq datasets derived from simulated and real datasets and performed robustly. The detailes are described in our paper [*`PathogenTrack and Yeskit: tools for identifying intracellular pathogens from single-cell RNA-sequencing datasets as illustrated by application to COVID-19`*](https://journal.hep.com.cn/fmd/EN/10.1007/s11684-021-0915-9).

### System Requirements

PathogenTrack has been tested on Linux platform with CentOS 7 and Mac platform with macOS 11.6.1.

## Installation

### PathogenTrack can be installed in two steps:

1 . Installing Miniconda on Linux / MacOS Platform. For details, please refer to [Miniconda Installation](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent).
```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh    # For Linux users
bash Miniconda3-latest-Linux-x86_64.sh                                        # For Linux users
-----------------------------------------------------------------------------------------------
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh   # For MacOS users
bash Miniconda3-latest-MacOSX-x86_64.sh                                       # For MacOS users
```

2 . Installing PathogenTrack and the dependencies.
```sh
conda env create -f environment.yml
```
Users are strongly sugguested to install these software with conda. 
The dependencies and test versions are listed below.

Package|Version
--|:--:
python|3.6.10
biopython|1.78
fastp|0.12.4
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

## Run PathogenTrack
Before running PathogenTrack, you should run `cellranger` or `alevin` to get the single cells' gene expression matrix. Here, we take the simulated 10X sequencing data as an example:

First, we use cellranger to get scRNA-seq expression matrix and valid barcodes:
```sh
cellranger count --id cellranger_out --transcriptom /path/to/cellranger_database/
```
**Attention** 

Three files must be ready to run PathogenTrack: 1) the valid barcode.tsv file; 2) the raw scRNA-seq fastq file (xxx_R1.fastq.gz; xxx_R2.fastq.gz).

Then we run PathogenTrack to identify and quantify pathogen expression at the single-cell level:

(Users should change the '/path/to/' in the following command to the databases' real paths)
```sh
conda activate PathogenTrack
PathogenTrack count --project_id PathogenTrack_out \
                    --pattern CCCCCCCCCCCCCCCCNNNNNNNNNN \
                    --min_reads 10 --confidence 0.11 \
                    --star_index /path/to/STAR_index/ \
                    --kraken_db /path/to/minikraken_8GB_20200312/ \
                    --barcode barcodes.tsv \
                    --read1 test_S1_L001_R1_001.fastq.gz \
                    --read2 test_S1_L001_R2_001.fastq.gz 
```
**IMPORTANT**: The Read 1 in the example is made up of 16 bp CB and 10 bp UMI, so the --pattern is *CCCCCCCCCCCCCCCCNNNNNNNNNN* (16C and 10N). Users must adjust the pattern with their own Read 1 accordingly. (for 10X Genomics scRNA-seq Chemistry Version 2: 16 bp CB and 10 bp UMI; for Version 3: 16 bp CB and 12 bp UMI)

*Note:* It may take 4-6 hours to complete one sample, and it depends on the performance of computational resources and the size of the raw single-cell data.

### Please see [QUICK_START.md](https://github.com/ncrna/PathogenTrack/blob/master/doc/QUICK_START.md) for a complete tutorial.

## Questions
If you have any questions/problems with PathogenTrack, feel free to leave an issue! We will try our best to provide support, address new issues, and keep improving this software.

## Citation
Wei Zhang, Xiaoguang Xu, Ziyu Fu, Jian Chen, Saijuan Chen, Yun Tan. PathogenTrack and Yeskit: tools for identifying intracellular pathogens from single-cell RNA-sequencing datasets as illustrated by application to COVID-19. Front. Med., https://doi.org/10.1007/s11684-021-0915-9

The preprint version can be found [here](https://journal.hep.com.cn/fmd/EN/10.1007/s11684-021-0915-9).
