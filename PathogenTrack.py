'''
This tool is only designed to work with single cell RNA-seq fastq reads.
'''

import sys
import re
import os
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b','--barcode',help='input barcode.tsv',required=True)
parser.add_argument('-i','--fastq',help='input fastq file',required=True)
parser.add_argument('-k','--kraken',help='input kraken file',required=True)
parser.add_argument('-t','--taxon',help='input taxon file',required=True)
parser.add_argument('-f','--fasta',help='output fasta file',required=True)
parser.add_argument('-o','--output',help='output file prefix',required=True)

args = parser.parse_args()
barcodeFile = args.barcode
fastqFile = args.fastq
krakenFile = args.kraken
taxonFile = args.taxon
fastaFile = args.fasta + ".fa"
outFile = args.output + "_matrix.txt"

# Read taxon file
taxons = {}
with open(taxonFile) as TX:
    for line in TX:
        taxid = line.strip().split("\t")[0]
        taxnm = line.strip().split("\t")[1]
        if taxid == "9606":
            continue
        taxons[taxid] = taxnm
TX.close()

# Define valid cell barcodes from cellranger barcodes.tsv file
cells = set()
with open(barcodeFile) as CB:
    for line in CB:
        cb = line.strip().split("-")[0]
        cells.add(cb)
CB.close()

Dict = {}
FA = open(fastaFile, "w")
#with gzip.open("NS.fq.gz", "rt") as FASTQ, open("NS.kraken2") as KK:
with open(fastqFile) as FASTQ, open(krakenFile) as KK:
    for record, row in zip(SeqIO.parse(FASTQ, "fastq"), KK):
        _, cb, umi = record.id.split("_")
        if cb not in cells:
            continue
        taxon = row.split("\t")[2]
        taxid = re.findall(r"taxid (\d+)", taxon)
        if taxid is None:
           taxid = taxon
        else:
           taxid = taxid[0]
        if taxid not in taxons: # only keep entries in taxon database
            continue
        FA.write('>' + record.id + ' ' + taxons[taxid] + '\n' + str(record.seq) + '\n')
        seq = cb + '_' + umi + '_' + record.seq
        #print("%s\t%s" % (organism, seq))
        organism = taxons[taxid]
        if organism in Dict:
            Dict[organism].add(str(seq))
        else:
            Dict[organism] = set()
            Dict[organism].add(str(seq))

FA.close()

genes = set()
gene_counts = {}

for gene in Dict:
    for seq in Dict[gene]:
        cell = seq.split("_")[0]
        count = 1
        genes.add(gene)
        if gene not in gene_counts:
            gene_counts[gene] = {}
        if cell not in gene_counts[gene]:
            gene_counts[gene][cell] = 1
        else:
            gene_counts[gene][cell] = gene_counts[gene][cell] + 1

with open(outFile, "w") as FO: 

    FO.write("gene" + '\t' + '\t'.join(sorted(cells)) + '\n')

    for gene in sorted(genes):
        counts = []
        for cell in sorted(cells):
            if cell in gene_counts[gene]:
                counts.append(gene_counts[gene][cell])
            else:
                counts.append(0)
        FO.write(gene + '\t' + '\t'.join(map(str, counts)) + '\n')

FO.close()
