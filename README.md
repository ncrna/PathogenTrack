# PathogenTrack
PathogenTrack is a python-based computational software based on **UMI-tools** and **Kraken2** developed to detect and identify pathogenic microorganisms from single-cell RNA-sequencing (scRNA-seq) raw data. We have tested PathogenTrack on various scRNA-seq datasets derived from human normal and tumor lung samples as described in our paper *'Detecting and studying pathogenic microorganisms invasion at the single-cell resolution using PathogenTrack'*.

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
3. Get PathogenTrack tool
```sh
git clone git@github.com:rstatistics/PathogenTrack.git
```
## Databases Preparation

### 1. Prepare the Human genome database
Download the Human GRCh38 genome and genome annotation file, and then decompress them:
```sh
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr_patch_hapl_scaff.gtf.gz
gzip -d Homo_sapiens.GRCh38.101.chr_patch_hapl_scaff.gtf.gz
```

Build STAR Index with the following command:
```sh
STAR --runThreadN 16 --runMode genomeGenerate --limitGenomeGenerateRAM 168632691637 --genomeDir ./ \
     --genomeFastaFiles ./Homo_sapiens.GRCh38.dna.toplevel.fa --sjdbGTFfile ./Homo_sapiens.GRCh38.101.chr_patch_hapl_scaff.gtf \
     --sjdbOverhang 100
```
*Note:* It was executed on a CentOS 7 system with 120GiB of memory, and cost 150GiB disk space in 10 hours.

### 2. Prepare Kraken2 database

```sh
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
tar zxf minikraken_8GB_202003.tgz
```
### 3. Prepare Taxonomy database
```sh
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-09-01.zip
unzip taxdmp_2020-09-01.zip
# only species-level was kept for further use
grep -wE 'forma|forma specialis|varietas|subspecies|strain|species' nodes.dmp | cut -f 1 > taxid.txt
grep 'scientific name' names.dmp | cut -f 1,3 > taxid2organism.txt
awk -F'\t' 'NR==FNR{a[$1]=$2; next}; {print $1"\t"a[$1]}' taxid2organism.txt taxid.txt > taxons.db
```

## How to use PathogenTrack?
Before running this tutorial, you should run cellranger or other tools (e.g., umi_tools) to get the single cells' gene expression matrix. Here, we took cellranger results as an example:
### Example 1:
Let's take the raw 10X sequencing data as an example:
```sh
conda activate PathogenTrack
python PathogenTrack.py --bcpattern CCCCCCCCCCCCCCCCNNNNNNNNNN --read1 Input_R1.fastq.gz --read2 Input_R2.fastq.gz --barcode barcodes.tsv --thread 8 --starindex /db/human/STAR-index/ --krakendb /db/minikraken_8GB_20200312/ --taxondb taxons.db --output Input_matrix.tsv
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
