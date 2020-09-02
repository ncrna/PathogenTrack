# Pathogen-Track
Pathogen-Track is a python-based computational software based on **---** and **---** developped to detect and identify pathogenic microorganisms from single-cell RNA-sequencing (scRNA-seq) raw data. This tool was tested on various scRNA-seq datasets derived from human tumor samples as described in our paper *'Detecting and studying pathogenic microorganisms invasion at the single-cell resolution using Pathogen-Track'*.

## Installation

Before running Pathogen-Track, several dependencies must be installed :

1 . The first step is to install [**UMI-tools**](https://github.com/CGATOxford/UMI-tools). Umi_tools is dependent on python>=3.5, numpy, pandas, scipy, cython, pysam, future, regex and matplotlib, to install it you should start an ssh session and type :

```sh
conda install -c bioconda -c conda-forge umi_tools
```
or
```sh
pip install umi_tools
```


## Tutorial

### Pre-processing of the scRNA-seq data

Before running this tutorial, you must run cellranger or other tools to quant the gene expression of single cells. We take cellranger as an example:
you got an barcodes.tsv in the output, such as:

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
#### 3. 


#### 4. Extract the UMIs
UMIs are strings of random nucleotides attached to the start of reads. Before the reads are mapped the random nucleotides must be removed from the reads, but the sequence must be kept. The 'extract' command of UMI-Tools moves the UMI from the read to the read name.

Cell barcodes are short nucleotide sequences, very much like UMIs, except instead of identifying independent molecules, they identify independent cells. We generally observe more of them in an experiment than there were cells. This could be for several reasons, including sequencing or PCR errors and the sequencing of empty droplets or those containing contaminants. Thus we must identify which cell barcodes we wish to use downstream. UMI-Tools whitelist command is used to produce a list of CB to use downstream.

whitelist currently allows the common method of taking the top X most abundant barcodes. X can be estimated automatically from the data using the knee method (for more detail see this blog post). However, it is just an estimate and for this data we've been told that there were 100 cells, so we can just supply that number (see variations section for performing the estimation for data sets where cell number is unknown).
