# PathogenTrack
PathogenTrack is an unsupervised computational software that uses unmapped single-cell RNAseq reads to characterize intracellular pathogens at the single-cell level. It is a python-based script that can be used to identify and quantify intracellular pathogenic `viruses` and `bacteria` reads at the single-cell level.
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
3. Get PathogenTrack
```sh
git clone git@github.com:rstatistics/PathogenTrack.git
```
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
STAR --runThreadN 16 --runMode genomeGenerate --limitGenomeGenerateRAM 168632691637 --genomeDir ./ \
     --genomeFastaFiles ./Homo_sapiens.GRCh38.dna.toplevel.fa --sjdbGTFfile ./Homo_sapiens.GRCh38.101.gtf \
     --sjdbOverhang 100
```

### 2. Prepare Kraken2 database

```sh
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
tar zxf minikraken_8GB_202003.tgz
```

## How to use PathogenTrack?
Before running this tutorial, you should run cellranger or other tools (e.g., umi_tools) to get the single cells' gene expression matrix. Here, we took cellranger results as an example:
### Example 1:
Let's take the simulated 10X sequencing data as an example:
```sh
conda activate PathogenTrack
python PathogenTrack.py count --project_id TMP --pattern CCCCCCCCCCCCCCCCNNNNNNNNNN --min_reads 10 --confidence 0.11 --star_index ~/database/STAR_index/ --kraken_db ~/database/minikraken_8GB_20200312/ --barcode barcodes.tsv --read1 simulation_S1_L001_R1_001.fastq.gz --read2 simulation_S1_L001_R2_001.fastq.gz 
```

### Example 2:
If we have extracted the CB (cell barcodes) and UMI for Read1 and added it to the read name of Read2, we can run it as follows:
```sh
conda activate PathogenTrack
python PathogenTrack.py --clean Input_R2.fastq_addCB_fastp.gz --barcode barcodes.tsv --thread 8 --starindex /db/human/STAR-index/ --krakendb /db/minikraken_8GB_20200312/ --taxondb taxons.db --output Input_matrix.tsv
```

**IMPORTANT**: The Read 1 in the example is made up of 16 bp CB and 12 bp UMI, so the --bcpattern is CCCCCCCCCCCCCCCCNNNNNNNNNNNN (16C and 12N). Users must adjust the bcpattern with their own Read 1 accordingly.

*Note:* It may take 4-6 hours to complete one sample, and it depends on the performance of computational resources and the size of the raw single-cell data.

## Questions

For questions and suggestions about the pipeline or the code, please contact [rstatistics@sjtu.edu.cn](mailto:rstatistics@sjtu.edu.cn) and [ty12260@rjh.com.cn](mailto:ty12260@rjh.com.cn). We will try our best to provide support, address new issues, and keep improving this software.
