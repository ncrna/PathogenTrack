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

VERSION = '0.2.3'
USAGE = '''%(prog)s [options]'''

def check_status(cmd):
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

def check_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return(0)

def check_deps():
    status, _ = check_status('fastp --help')
    if (status != 0):
        sys.stderr.write('fastp is required for PathogenTrack!\n')

    status, _ = check_status('STAR --help')
    if (status != 0):
        sys.stderr.write('STAR is required for PathogenTrack!\n')

    status, _ = check_status('umi_tools --help')
    if (status != 0):
        sys.stderr.write('umi_tools is required for PathogenTrack!\n')

    status, _ = check_status('kraken2 --help')
    if (status != 0):
        sys.stderr.write('kraken2 is required for PathogenTrack!\n')

def trim(project_id, barcode, output): ## trim -1 at the end of barcode
    check_deps()
    check_dir(project_id)
    output = project_id + '/' + output
    fout = open(output, "w")
    with open(barcode, 'r') as fin:
        for line in fin:
            cb = line.strip().split("-")[0]
            fout.write(cb + '\n')
    fin.close()
    fout.close()

def extract(project_id, read1, read2, pattern, barcode, output, ignore_suffix): ## add cell barcode
    check_deps()
    check_dir(project_id)
    pattern = str(pattern)
    barcode = project_id + '/' + barcode
    output = project_id + '/' + output
    cmd = 'umi_tools extract --bc-pattern ' + str(pattern)
    cmd += ' --stdin ' + read1 + ' --stdout /dev/null' + ' --read2-in ' + read2
    cmd += ' --read2-out ' + output + ' --filter-cell-barcode' + ' --whitelist ' + barcode
    if (ignore_suffix):
        cmd += ' --ignore-read-pair-suffixes'
    status, stdout = check_status(cmd)
    print("status is %s and output is %s" % (status, stdout))
    return(0)

def filter(project_id, input_reads, thread, output, trim_poly_x, filter_low_complexity, complexity_threshold): ## filter low quality / complexity reads
    check_deps()
    check_dir(project_id)
    input_reads = project_id + '/' + input_reads
    output = project_id + '/' + output
    cmd = 'fastp --thread ' + str(thread) + ' -i ' + input_reads + ' -o ' + output
    if (trim_poly_x):
        cmd += ' --trim_poly_x'
    if (filter_low_complexity):
        cmd += ' --low_complexity_filter'
        if (complexity_threshold):
            cmd += ' --complexity_threshold ' + str(complexity_threshold)
    cmd += ' && rm fastp.json && rm fastp.html'
    status, stdout = check_status(cmd)
    print("status is %s and output is %s" % (status, stdout))
    return(0)

def align(project_id, star_index, input_reads, thread, output): ## do star alignment
    check_deps()
    check_dir(project_id)
    input_reads = project_id + '/' + input_reads
    output = project_id + '/' + output
    cmd = 'STAR --genomeDir ' + star_index + ' --readFilesIn ' + input_reads + ' --readFilesCommand gunzip -c '
    cmd += '--runThreadN ' + str(thread) + ' --outFilterMismatchNmax 6 ' + ' --outSAMtype None ' + ' --outFilterMultimapNmax 20 '
    cmd += '--outFilterIntronMotifs RemoveNoncanonical ' + ' --quantMode - ' + ' --outFileNamePrefix ' + project_id + '/ --outReadsUnmapped Fastx'
    cmd += ' && mv ' + project_id + '/Unmapped.out.mate1 ' + output + ' && rm ' + project_id + '/Log.final.out && rm ' + project_id + '/Log.out && rm ' + project_id + '/Log.progress.out && rm ' + project_id + '/SJ.out.tab'
    status, stdout = check_status(cmd)
    print("status is %s and output is %s" % (status, stdout))
    return(0)

def classify(project_id, kraken_db, input_reads, thread, confidence, out_reads, out_table, out_report): ## run kraken2
    check_deps()
    check_dir(project_id)
    input_reads = project_id + '/' + input_reads
    out_reads = project_id + '/' + out_reads
    out_table = project_id + '/' + out_table
    out_report = project_id + '/' + out_report
    
    cmd = 'kraken2 --db ' + kraken_db + ' --confidence ' + str(confidence) + ' --threads ' + str(thread) + ' --use-names --classified-out ' + out_reads + ' --report ' + out_report + ' '
    cmd += input_reads + ' | grep -w ^C > ' + out_table
    status, stdout = check_status(cmd)
    print("status is %s and output is %s" % (status, stdout))
    return(0)

def check_rank(rank):
    ranks = ['U', 'R', 'D', 'P', 'C', 'O', 'F', 'G', 'S']
    level = 0
    if rank in ranks:
        level = ranks.index(rank) # rank in comman order, eg 'F', 'G', 'S'
    elif rank[0] in ranks: # rank in sub-comman order, eg 'F1', 'G1', 'S1'
        level = ranks.index(rank[0]) + 0.5
    return(level)

def quant(project_id, barcode, input_reads, input_table, input_report, min_reads, output_reads, output_matrix): ## deduplication and quantification microbes at single-cell levels
    check_deps()
    check_dir(project_id)
    min_reads = int(min_reads)
    barcode = project_id + '/' + barcode
    input_reads = project_id + '/' + input_reads
    input_table = project_id + '/' + input_table
    input_report = project_id + '/' + input_report
    output_reads = project_id + '/' + output_reads
    output_matrix = project_id + '/' + output_matrix

    cells = set()
    with open(barcode) as fcb:
        for line in fcb:
            cb = line.strip().split("-")[0]
            cells.add(cb)
    fcb.close()

    ranks = ['U', 'R', 'D', 'P', 'C', 'O', 'F', 'G', 'S']
    Taxons = dict()
    level_taxid = '1' # initialize
    level_taxname = 'root'
    last_rank = 'R'
    current_rank = 'R'
    selected = False
    selected_taxid = ''
    selected_taxname = ''
    tmp_dict = dict() # store rows containing species/strain
    with open(input_report) as freport:
        for line in freport:
            total_count = int(line.strip().split("\t")[1])
            rank_count = int(line.strip().split("\t")[2])
            current_rank = line.strip().split("\t")[3]
            taxid = line.strip().split("\t")[4]
            taxname = line.strip().split("\t")[5].strip()
            if current_rank == 'U': # skip first line in report
                continue
            if check_rank(current_rank) < check_rank(last_rank) or check_rank(current_rank) == check_rank(last_rank) == 8: # changed organism
                last_rank = current_rank
                if selected == True:
                    for row in dict(sorted(tmp_dict.items(), key=lambda item: item[1], reverse=True)).keys():
                        if not (selected_taxid and selected_taxname):
                            selected_taxid = row.strip().split("\t")[4]
                            selected_taxname = row.strip().split("\t")[5].strip()
                            selected_reads = int(row.strip().split("\t")[2])
                        tmp_taxid = row.strip().split("\t")[4]
                        tmp_taxname = row.strip().split("\t")[5].strip()
                        if selected_reads >= min_reads:
                            Taxons[tmp_taxid] = selected_taxid + '|' + selected_taxname
                    selected = False
                    selected_taxid = ''
                    selected_taxname = ''
                    tmp_dict = dict()
                if check_rank(current_rank) >= 8:
                    if total_count == rank_count:
                        selected = True
                    tmp_dict[line] = rank_count
            else:
                last_rank = current_rank
                if check_rank(current_rank) >= 8:
                    if total_count == rank_count:
                        selected = True
                    tmp_dict[line] = rank_count
        if selected == True:
            for row in dict(sorted(tmp_dict.items(), key=lambda item: item[1], reverse=True)).keys():
                if not (selected_taxid and selected_taxname):
                    selected_taxid = row.strip().split("\t")[4]
                    selected_taxname = row.strip().split("\t")[5].strip()
                    selected_reads = int(row.strip().split("\t")[2])
                tmp_taxid = row.strip().split("\t")[4]
                tmp_taxname = row.strip().split("\t")[5].strip()

                if selected_reads >= min_reads:
                    Taxons[tmp_taxid] = selected_taxid + '|' + selected_taxname
    freport.close()

    Dict = {}
    Pool = {}
    fout = open(output_reads, "w")
    with open(input_reads) as ffastq, open(input_table) as fkraken:
        for record, row in zip(SeqIO.parse(ffastq, "fastq"), fkraken):
            cb, umi = record.id.split(" ")[0].split("_")[-2:]
            if cb not in cells:
                print("CB is: " + cb + " and UMI is: " + umi)
                continue
            taxonomy = row.split("\t")[2]
            taxname, taxid = taxonomy.split(" (taxid ")
            taxid = taxid.split(")")[0]
            if taxid == "9606" or taxid == "10090" or taxid == "28384":
                continue
            if taxid not in Taxons.keys():
                continue
            taxname = Taxons[taxid].split('|')[1]
            fout.write('>' + record.id + ' ' + taxname + '\n' + str(record.seq) + '\n')
            ids = cb + '_' + umi
            if ids not in Pool:
                Pool[ids] = {}
            if taxid not in Pool[ids]:
                Pool[ids][taxname] = 1
            else:
                Pool[ids][taxname] = Pool[ids][taxname] + 1
            if taxname in Dict:
                Dict[taxname].add(ids)
            else:
                Dict[taxname] = set()
                Dict[taxname].add(ids)
    fout.close()

    # find optimal organism with same barcode + umi
    organisms = {}
    for ids in Pool:
        organisms[ids] = sorted(Pool[ids].items(), reverse=True)[0][0]

    # store organisms in set
    organisms_set = {organisms[k] for k in organisms}

    genes = set()
    gene_counts = {}
    for gene in Dict:
        if gene not in organisms_set:
            continue
        for ids in Dict[gene]:
            cell = ids.split("_")[0]
            count = 1
            genes.add(gene)
            if gene not in gene_counts:
                gene_counts[gene] = {}
            if cell not in gene_counts[gene]:
                gene_counts[gene][cell] = 1
            else:
                gene_counts[gene][cell] = gene_counts[gene][cell] + 1

    with open(output_matrix, "w") as fmatrix:
        fmatrix.write("gene" + '\t' + '\t'.join(sorted(cells)) + '\n')
        for gene in sorted(genes):
            counts = []
            for cell in sorted(cells):
                if cell in gene_counts[gene]:
                    counts.append(gene_counts[gene][cell])
                else:
                    counts.append(0)
            if sorted(counts, reverse=True)[0] >= 2: # remove species/strains with maximum UMIs < 2
                positive_counts = str(sum(i > 0 for i in counts))
                fmatrix.write(gene + '\t' + '\t'.join(map(str, counts)) + '\n')
                print(gene + '\t' + positive_counts)
    fmatrix.close()

    return

def count(project_id, star_index, kraken_db, barcode, read1, read2, pattern, thread, confidence, min_reads, ignore_suffix, trim_poly_x, filter_low_complexity, complexity_threshold):
    check_deps()
    check_dir(project_id)
    trim(project_id, barcode, output='barcodes.tsv')
    extract(project_id, read1, read2, pattern, barcode='barcodes.tsv', output='reads_barcoded.fq.gz', ignore_suffix=ignore_suffix)
    filter(project_id, input_reads='reads_barcoded.fq.gz', thread=thread, output='reads_barcoded_clean.fq.gz', trim_poly_x=trim_poly_x, filter_low_complexity=filter_low_complexity, complexity_threshold=complexity_threshold)
    align(project_id, star_index, input_reads='reads_barcoded_clean.fq.gz', thread=thread, output='unmapped_reads.fq')
    classify(project_id, kraken_db, input_reads='unmapped_reads.fq', thread=thread, confidence=confidence, out_reads='classified_reads.fq', out_table='classified_table.txt', out_report='classified_report.txt')
    quant(project_id, barcode='barcodes.tsv', input_reads='classified_reads.fq', input_table='classified_table.txt', input_report='classified_report.txt', min_reads=min_reads, output_reads='microbes.fa', output_matrix='microbes.tsv')
    return(0)

def get_args():
    class Parser(argparse.ArgumentParser):
        def print_message(self, message, color=None):
            if message:
                sys.stderr.write('\033[91m' + message + '\033[0m\n')
        def error(self, message):
            #sys.stderr.write('\nerror: %s\n\n' % message)
            #sys.stderr.write('\033[31m', '\nerror: %s\n\n' % message)
            self.print_message('\nerror: %s\n' % message)
            self.print_help(sys.stderr)
            sys.exit(2)
    #parser = argparse.ArgumentParser(usage=USAGE)
    parser = Parser()
    # Global options
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    subparsers = parser.add_subparsers(dest='command')
    
    # Command: count
    parser_count = subparsers.add_parser('count',
                        help='one command to count microbes at single-cell levels')
    parser_count.add_argument('--project_id', default='./', action='store', 
                        help='project id to store output files', required=True)
    parser_count.add_argument('--star_index', action='store', default='',
                        help='star indices of the reference genome', required=True)
    parser_count.add_argument('--kraken_db', action='store', default='',
                        help='kraken2 database name', required=True)
    parser_count.add_argument('--barcode', action='store', default='barcodes.tsv',
                        help='file name for raw barcodes from cellranger, decompress it before use', required=True)
    parser_count.add_argument('--read1', default='', 
                        help='file name for Read1', action='store', required=True)
    parser_count.add_argument('--read2', default='', 
                        help='file name for Read2', action='store', required=True)
    parser_count.add_argument('--pattern', action='store', default='CCCCCCCCCCCCCCCCNNNNNNNNNN', 
                        help='barcode pattern, default="CCCCCCCCCCCCCCCCNNNNNNNNNN", please refer to umi_tools help manual', required=False)
    parser_count.add_argument('--thread', action='store', default=8, type=int, 
                        help='number of threads to use, default is 8', required=False)
    parser_count.add_argument('--confidence', action='store', default=0.11, type=float, 
                        help='confidence threshold to use [0-1], default is 0.11', required=False)
    parser_count.add_argument('--min_reads', action='store', default=10, 
                        help='minimum number of reads that supports the taxonomy, default is 10', required=False)
    parser_count.add_argument('--trim_poly_x', action='store_true', 
                        help='enable polyA/T/C/G/N trimming in 3\' ends', required=False)
    parser_count.add_argument('--filter_low_complexity', action='store_true', 
                        help='enable filter out low complexity reads, default is False', required=False)
    parser_count.add_argument('--complexity_threshold', action='store', default=30,
                        help='the threshold for low complexity filter [0-100], default is 30', required=False)
    parser_count.add_argument('--ignore_suffix', action='store_true',
                        help='ignore the read name suffixes, eg., /1 and /2 in the read are ignored, default is False', required=False)

    # Command: trim
    parser_trim = subparsers.add_parser('trim',
                        help='trim -1 at the end of barcode')
    parser_trim.add_argument('--project_id', default='./', action='store', 
                        help='project id to store output files', required=True)
    parser_trim.add_argument('--barcode', default='', action='store', 
                        help='file name for raw barcodes from cellranger, decompress it before use', required=True)
    parser_trim.add_argument('--output', default='barcodes.tsv', action='store', 
                        help='file to output barcodes to, default is barcodes.tsv', required=False)

    # Command: extract
    parser_extract = subparsers.add_parser('extract', 
                        help='extract and append barcode and UMI to the read2 name')
    parser_extract.add_argument('--project_id', default='./', action='store', 
                        help='project id to store output files', required=True)
    parser_extract.add_argument('--read1', default='', 
                        help='file name for Read1', action='store', required=True)
    parser_extract.add_argument('--read2', default='', 
                        help='file name for Read2', action='store', required=True)
    parser_extract.add_argument('--pattern', action='store', default='CCCCCCCCCCCCCCCCNNNNNNNNNN', 
                        help='barcode pattern, default="CCCCCCCCCCCCCCCCNNNNNNNNNN", please refer to umi_tools help manual', required=False)
    parser_extract.add_argument('--barcode', default='barcodes.tsv',
                        help='file name for valid barcodes', action='store', required=False)
    parser_extract.add_argument('--output', action='store', default='reads_barcoded.fq.gz',
                        help='file to output processed read to', required=False)
    parser_extract.add_argument('--ignore_suffix', action='store_true',
                        help='ignore the read name suffixes, eg., /1 and /2 in the read are ignored, default is False', required=False)
    # Command: filter
    parser_filter = subparsers.add_parser('filter',
                        help='filter out low quality / complexity reads')
    parser_filter.add_argument('--project_id', default='./', action='store', 
                        help='project id to store output files', required=True)
    parser_filter.add_argument('--input_reads', action='store', default='reads_barcoded.fq.gz',
                        help='input single-end fastq file produced by "extract" command', required=False)
    parser_filter.add_argument('--thread', action='store', default=8, type=int, 
                        help='number of threads to use, default is 8', required=False)
    parser_filter.add_argument('--output', action='store', default='reads_barcoded_clean.fq.gz',
                        help='file to output processed read to, default is reads_barcoded_clean.fq.gz', required=False)
    parser_filter.add_argument('--trim_poly_x', action='store_true', 
                        help='enable polyA/T/C/G/N trimming in 3\' ends', required=False)
    parser_filter.add_argument('--filter_low_complexity', action='store_true', 
                        help='enable filter out low complexity reads, default is False', required=False)
    parser_filter.add_argument('--complexity_threshold', action='store', default=30,
                        help='the threshold for low complexity filter [0-100], default is 30', required=False)
    # Command: align
    parser_align = subparsers.add_parser('align',
                        help='align clean reads to the reference genome')
    parser_align.add_argument('--project_id', default='./', action='store', 
                        help='project id to store output files', required=True)
    parser_align.add_argument('--star_index', action='store', default='',
                        help='star indices of the reference genome', required=True)
    parser_align.add_argument('--input_reads', action='store', default='reads_barcoded_clean.fq.gz',
                        help='input single-end fastq file', required=False)
    parser_align.add_argument('--thread', action='store', default=8, type=int, 
                        help='number of threads to use, default is 8', required=False)
    parser_align.add_argument('--output', action='store', default='unmapped_reads.fq',
                        help='file to output unmapped reads to', required=False)
    # Command: classify
    parser_classify = subparsers.add_parser('classify',
                        help='classify unmapped reads to taxons')
    parser_classify.add_argument('--project_id', default='./', action='store', 
                        help='project id to store output files', required=True)
    parser_classify.add_argument('--kraken_db', action='store', default='',
                        help='kraken2 database name', required=True)
    parser_classify.add_argument('--input_reads', action='store', default='unmapped_reads.fq',
                        help='input single-end fastq file produced by "align" command', required=False)
    parser_classify.add_argument('--thread', action='store', default=8, type=int, 
                        help='number of threads to use, default is 8', required=False)
    parser_classify.add_argument('--confidence', action='store', default=0.11, type=float, 
                        help='confidence threshold to use [0-1], default is 0.11', required=False)
    parser_classify.add_argument('--out_reads', action='store', default='classified_reads.fq',
                        help='file to output processed reads to', required=False)
    parser_classify.add_argument('--out_table', action='store', default='classified_table.txt',
                        help='file to output kraken results to', required=False)
    parser_classify.add_argument('--out_report', action='store', default='classified_report.txt',
                        help='file to output kraken report to', required=False)
    # Command: quant
    parser_quant = subparsers.add_parser('quant',
                        help='deduplication and quantification of microbes at single-cell levels')
    parser_quant.add_argument('--project_id', default='./', action='store', 
                        help='project id to store output files', required=True)
    parser_quant.add_argument('--barcode', action='store', default='barcodes.tsv',
                        help='file name for valid barcodes', required=False)
    parser_quant.add_argument('--input_reads', action='store', default='classified_reads.fq',
                        help='input fastq file produced by "classify" command', required=False)
    parser_quant.add_argument('--input_table', action='store', default='classified_table.txt',
                        help='input table file produced by "classify" command', required=False)
    parser_quant.add_argument('--input_report', action='store', default='classified_report.txt',
                        help='input report file produced by "classify" command', required=False)
    parser_quant.add_argument('--min_reads', action='store', default=10, 
                        help='minimum number of reads that supports the taxonomy, default is 10', required=False)
    parser_quant.add_argument('--output_reads', action='store', default='microbes.fa',
                        help='file to output classified reads to, default is microbes.fa', required=False)
    parser_quant.add_argument('--output_matrix', action='store', default='microbes.tsv',
                        help='file to output microbe quantification matrix to, default is microbes.tsv', required=False)

    return parser

def main():
    parser = get_args()
    args = parser.parse_args()
    if not args.command:
        parser.parse_args(["--help"])
        sys.exit(0)
    if args.command == "count":
        count(args.project_id, args.star_index, args.kraken_db, args.barcode, args.read1, args.read2, args.pattern, args.thread, args.confidence, args.min_reads, args.ignore_suffix, args.trim_poly_x, args.filter_low_complexity, args.complexity_threshold)
    elif args.command == "trim":
        trim(args.project_id, args.barcode, args.output)
    elif args.command == "extract":
        extract(args.project_id, args.read1, args.read2, args.pattern, args.barcode, args.output, args.ignore_suffix)
    elif args.command == "filter":
        filter(args.project_id, args.input_reads, args.thread, args.output, args.trim_poly_x, args.filter_low_complexity, args.complexity_threshold)
    elif args.command == "align":
        align(args.project_id, args.star_index, args.input_reads, args.thread, args.output)
    elif args.command == "classify":
        classify(args.project_id, args.kraken_db, args.input_reads, args.thread, args.confidence, args.out_reads, args.out_table, args.out_report)
    elif args.command == "quant":
        quant(args.project_id, args.barcode, args.input_reads, args.input_table, args.input_report, args.min_reads, args.output_reads, args.output_matrix)
    return

if __name__ == '__main__':
    main()

