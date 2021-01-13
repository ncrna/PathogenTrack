'''
This tool is only designed to work with single cell RNA-seq fastq reads.
'''

import sys
import re
import os
import gzip
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import subprocess

VERSION = 'v0.2'
USAGE = '''%(prog)s [options]'''

def checkstatus(cmd):
    ''' get command status '''
    status = 0
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        status = e.returncode
        output = e.output
    output = output.decode()
    if output[-1] == '\n':
        output = output[:-1]
    return (status, output)

def runCellRanger(args): ## run cellranger
    print()

def processCB(args): ## process cell barcode
    barcodeInput = args.barcode
    cells = set()
    with open(barcodeInput) as CB:
        for line in CB:
            cb = line.strip().split("-")[0]
            cells.add(cb)
    CB.close()

    barcodeFile = ''
    if args.read2.endswith('.gz'):
        barcodeFile = args.read2.replace('.gz', '_barcodes.tsv')
    else:
        barcodeFile = args.read2 + '_barcodes.tsv'
    
    CB = open(barcodeFile, "w")
    for cb in cells:
        CB.write(cb + '\n')
    CB.close()
    return(barcodeFile)

def addCB(barcodeFile, args): ## add cell barcode
    output = ''
    if args.read2.endswith('.gz'):
        output = args.read2.replace('.gz', '_addCB.gz')
    else:
        output = args.read2 + '_addCB.gz'
    cmd = 'umi_tools extract --bc-pattern ' + str(args.bcpattern)
    cmd += ' --stdin ' + args.read1 + ' --stdout /dev/null' + ' --read2-in ' + args.read2
    cmd += ' --read2-out ' + output + ' --filter-cell-barcode' + ' --whitelist ' + barcodeFile
    status, stdout = checkstatus(cmd)
    print("status is %s and output is %s" % (status, stdout))
    return(output)

def qc(raw, args): ## quality control
    output = raw.replace('.gz', '_fastp.gz')
    #cmd = 'fastp --thread ' + str(args.thread) + ' -i ' + raw + ' -o ' + output + ' --low_complexity_filter' + ' --complexity_threshold 50 --trim_poly_x' + ' && rm fastp.json && rm fastp.html'
    cmd = 'fastp --thread ' + str(args.thread) + ' -i ' + raw + ' -o ' + output + ' --low_complexity_filter' + ' && rm fastp.json && rm fastp.html'
    status, stdout = checkstatus(cmd)
    print("status is %s and output is %s" % (status, stdout))
    return(output)

def star(clean, args): ## do star alignment
    rmHostFq = 'rmHost.fq'
    cmd = 'STAR --genomeDir ' + args.starindex + ' --readFilesIn ' + clean + ' --readFilesCommand zcat '
    cmd += '--runThreadN ' + str(args.thread) + ' --outFilterMismatchNmax 6 ' + ' --outSAMtype None ' + ' --outFilterMultimapNmax 20 '
    cmd += '--outFilterIntronMotifs RemoveNoncanonical ' + ' --quantMode - ' + ' --outFileNamePrefix ' + '_' + ' --outReadsUnmapped Fastx && mv _Unmapped.out.mate1 rmHost.fq && rm _Log.final.out && rm _Log.out && rm _Log.progress.out && rm _SJ.out.tab'
    status, stdout = checkstatus(cmd)
    print("status is %s and output is %s" % (status, stdout))
    return(rmHostFq)

def kraken(rmHostFq, args): ## run kraken2
    krakenFq = 'rmHost_kraken.fq'
    krakenFile = 'kraken.out'
    cmd = 'kraken2 --db ' + args.krakendb + ' --threads ' + str(args.thread) + ' --classified-out ' + krakenFq + ' ' + rmHostFq
    cmd += ' | grep -w ^C > ' + krakenFile
    status, stdout = checkstatus(cmd)
    print("status is %s and output is %s" % (status, stdout))
    return(krakenFq, krakenFile)

def readCB(args):
    cells = set()
    with open(args.barcode) as CB:
        for line in CB:
            cb = line.strip().split("-")[0]
            cells.add(cb)
    CB.close()
    return cells

def dedup(cells, krakenFq, krakenFile, args): ## deduplication
    taxons = {}
    with open(args.taxondb) as TX:
        for line in TX:
            taxid = line.strip().split("\t")[0]
            taxnm = line.strip().split("\t")[1]
            if taxid == "9606":
                continue
            taxons[taxid] = taxnm
    TX.close()

    Dict = {}
    fastaFile = 'out.fa'

    #FA = open(fastaFile, "w")
    with open(krakenFq) as FASTQ, open(krakenFile) as KK:
        for record, row in zip(SeqIO.parse(FASTQ, "fastq"), KK):
            _, cb, umi = record.id.split("_")
            if cb not in cells:
                continue
            taxid = row.split("\t")[2]
            if taxid not in taxons: # only keep entries in taxon database
                continue
            #FA.write('>' + record.id + ' ' + taxons[taxid] + '\n' + str(record.seq) + '\n')
            seq = cb + '_' + umi + '_' + record.seq
            organism = taxons[taxid]
            if organism in Dict:
                Dict[organism].add(str(seq))
            else:
                Dict[organism] = set()
                Dict[organism].add(str(seq))
    #FA.close()
    return(Dict)

def quant(cells, Dict, args):
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

    with open(args.output, "w") as FO:
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

    return

def get_args():
    parser = argparse.ArgumentParser(usage=USAGE)
    parser.add_argument('--version', action='version', version=VERSION)
    parser.add_argument('--bcpattern', action='store', default='CCCCCCCCCCCCCCCCNNNNNNNNNN',
                        help='barcode pattern, please refer to umi_tools help manual')
    parser.add_argument('--read1', action='store', default='',
                        help='file name for Read1', required=False)
    parser.add_argument('--read2', action='store', default='',
                        help='file name for Read2', required=False)
    parser.add_argument('--clean', action='store', default=None,
                        help='clean single-end Read if Read1 and Read2 are not set', required=False)
    parser.add_argument('--barcode', action='store', default='',
                        help='file name for valid barcodes', required=True)
    parser.add_argument('--thread', action='store', default=8,
                        help='number of threads to use')
    parser.add_argument('--starindex', action='store', default='',
                        help='the path for star index', required=True)
    parser.add_argument('--krakendb', action='store', default='',
                        help='kraken2 database', required=True)
    parser.add_argument('--taxondb', action='store', default='',
                        help='taxon database', required=True)
    parser.add_argument('--output', action='store', default='',
                        help='output name for pathogen quantification matrix', required=True)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    if (args.clean == None):
        barcode = processCB(args)
        raw = addCB(barcode, args)
        print("raw file is %s" % raw)
        clean = qc(raw, args)
        print("clean is %s" % clean)
    else:
        clean = args.clean
    rmHostFq = star(clean, args)
    krakenFq, krakenFile = kraken(rmHostFq, args)
    cells = readCB(args)
    Dict = dedup(cells, krakenFq, krakenFile, args)
    quant(cells, Dict, args)
    return

if __name__ == '__main__':
    main()

