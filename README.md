# PathogenTrack
PathogenTrack is a python-based computational software based on **UMI-tools** and **Kraken2** developped to detect and identify pathogenic microorganisms from single-cell RNA-sequencing (scRNA-seq) raw data. This tool was tested on various scRNA-seq datasets derived from human normal and tumor lung samples as described in our paper *'Detecting and studying pathogenic microorganisms invasion at the single-cell resolution using PathogenTrack'*.

## Installation

#### PathogenTrack can be installed in two steps:

1 . Installing Miniconda on Linux Platform. For details, please refer to [Miniconda Installation](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent).
```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

2 . Installing PathogenTrack.
```sh
conda env create -f environment.yml
```

## Databases Preparation

### 1. Prepare Human genome database

```sh
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr_patch_hapl_scaff.gtf.gz
gzip -d Homo_sapiens.GRCh38.101.chr_patch_hapl_scaff.gtf.gz
STAR --runThreadN 16 --runMode genomeGenerate --limitGenomeGenerateRAM 168632691637 --genomeDir ./ --genomeFastaFiles ./Homo_sapiens.GRCh38.dna.toplevel.fa --sjdbGTFfile ./Homo_sapiens.GRCh38.101.chr_patch_hapl_scaff.gtf --sjdbOverhang 100
```
This was performed on CentOS7, 120G RAM, 150G Disk, in 10 hours.

### 2. Prepare Kraken2 database

```sh
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
tar zxf minikraken_8GB_202003.tgz
```
### 3. Prepare Taxon database
```sh
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-09-01.zip
unzip taxdmp_2020-09-01.zip
# only species-level was kept for further use
grep -wE 'forma|forma specialis|varietas|subspecies|strain|species' nodes.dmp | cut -f 1 > taxid.txt
grep 'scientific name' names.dmp | cut -f 1,3 > taxid2organism.txt
awk -F'\t' 'NR==FNR{a[$1]=$2; next}; {print $1"\t"a[$1]}' taxid2organism.txt taxid.txt > taxons.db
```


## Tutorial

Before running this tutorial, you must run cellranger or other tools (eg., umi_tools) to quant the gene expression of single cells. We take cellranger as an example:

#### 1. Prepare the barcode file

```sh
sed 's/-.*//' barcodes.tsv > Input_barcodes.tsv
```

#### 2. Extract the barcodes and filter the reads

The next step is to extract the CB (cell barcodes) and UMI from Read 1 and add it to the Read 2 read name. We will also filter out reads that do not match one of the accepted cell barcode.

The most basic form of this is executed with:
```sh
umi_tools extract --bc-pattern CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
                  --stdin Input_1.fastq.gz \
                  --stdout /dev/null \
                  --read2-in Input_2.fastq.gz \
                  --read2-out Input_extracted_2.fq.gz \
                  --filter-cell-barcode \
                  --whitelist Input_barcodes.tsv
```
#### 3. Classify reads by taxon
Kraken is a taxonomic sequence classifier that assigns taxonomic labels to NGS sequences. Kraken examines the k-mers within a query sequence and uses the information within those k-mers to query a database. That database maps k-mers to the lowest common ancestor (LCA) of all genomes known to contain a given k-mer.
```sh
kraken2 --db minikraken_8GB_20200312 \
        --threads 40 \
        --report Input.kreport2 \
        --classified-out Input_extracted.fq \
        Input_extracted_2.fq.gz 1> Input.kraken2
```
```sh
awk '$1=="C"' Input.kraken2 > Input.kraken
```

#### 4. Deduplication and Quantification
```sh
python umi_count-v0.1.py -b Input_barcodes.tsv \
                         -i Input_extracted.fq \
                         -k Input.kraken \
                         -t taxons.db \
                         -o Input
```
