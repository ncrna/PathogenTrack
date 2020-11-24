# PathogenTrack
PathogenTrack is a python-based computational software based on **UMI-tools** and **Kraken2** developed to detect and identify pathogenic microorganisms from single-cell RNA-sequencing (scRNA-seq) raw data. We have tested PathogenTrack on various scRNA-seq datasets derived from human normal and tumor lung samples as described in our paper *'Detecting and studying pathogenic microorganisms invasion at the single-cell resolution using PathogenTrack'*.

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

## Tutorial

Before running this tutorial, you should run cellranger or other tools (e.g., umi_tools) to get the single cells' gene expression matrix. Here, we took cellranger results as an example:

#### 1. Modify the barcode file

If cellranger output files are in gzip format, please decompress them first. Then remove the suffix '-1' of the barcode file: 

```sh
sed -i 's/-.*//' barcodes.tsv > barcodes.tsv
```

#### 2. Extract the barcodes and filter the reads

The next step is to extract the CB (cell barcodes) and UMI from Read 1 and add it to the read name of Read 2. We also filter out reads of which the UMI does not exist in the accepted cell barcodes. 
**IMPORTANT**: The Read 1 here is made up of 16 bp CB and 12 bp UMI, so the --bc-pattern is CCCCCCCCCCCCCCCCNNNNNNNNNNNN (16C and 12N). Users must adjust the pattern with their own Read 1 accordingly.

The most basic form of this is executed with:
```sh
umi_tools extract --bc-pattern CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
                  --stdin Input_R1.fq.gz \
                  --stdout /dev/null \
                  --read2-in Input_R2.fq.gz \
                  --read2-out Input_extracted_R2.fq.gz \
                  --filter-cell-barcode \
                  --whitelist barcodes.tsv
```
#### 3. Filter out reads with low quality or low complexity

The barcoded reads with low quality and low complexity were filtered out using fastp with the following command:

```sh
fastp --thread 8 --low_complexity_filter -i Input_extracted_R2.fq.gz -o Input_R2.fp.fq.gz
```

#### 4. Filter out reads from the host reference genome

Then we aligned the quality-filtered Read 2 from step3 to the human reference genome using STAR under the following command:

```sh
STAR --genomeDir ./STAR-index --readFilesIn Input_R2.fp.fq.gz --readFilesCommand zcat --runThreadN 16 \
     --outFilterMismatchNmax 6 --outSAMtype None --outFilterMultimapNmax 20 --outFilterIntronMotifs RemoveNoncanonical \
     --quantMode - --outFileNamePrefix Input_ --outReadsUnmapped Fastx
```

We renamed the unmapped reads with 'Input_rmHost.fq':
```sh
mv Input_Unmapped.out.mate1 Input_rmHost.fq
```
#### 5. Classify unmapped reads by taxonomy
Kraken2 is an excellent taxonomic sequence classifier that assigns taxonomic labels to NGS sequences. Here we employed it to classify the taxonomy of the unmapped reads:
```sh
kraken2 --db minikraken_8GB_20200312 \
        --threads 40 \
        --report Input.kreport2 \
        --classified-out Input_rmHost_kraken.fq \
        Input_rmHost.fq 1> Input.kraken2
```
Only classified entries are kept for further use:
```sh
awk '$1=="C"' Input.kraken2 > Input.kraken
```

#### 6. Reads de-Duplication and Quantification

The script 'PathogenTrack.py' was designed for reads de-duplication and pathogen species' abundance quantification at the single-cell level. It output a matrix with rows represent pathogen species and columns represent cells.

```sh
python PathogenTrack.py -b barcodes.tsv \
                        -i Input_rmHost.fq \
                        -k Input.kraken \
                        -t taxons.db \
                        -o Input
```
