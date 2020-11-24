File=$1
Thread=16

ulimit -SHn 4000

echo "$File starts at `date`"

umi_tools extract --bc-pattern CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin ${File}_R1.fastq.gz \
		  --stdout /dev/null \
		  --read2-in ${File}_R2.fastq.gz \
		  --read2-out ${File}_R2.fq.gz \
		  --filter-cell-barcode \
		  --whitelist ${File}/barcodes.tsv

fastp --thread $Thread --low_complexity_filter \
      -i ${File}_R2.fq.gz \
      -o ${File}_R2.fp.fq.gz

STAR --genomeDir ./STAR-index \
     --readFilesIn ${File}_R2.fp.fq.gz \
     --readFilesCommand zcat \
     --runThreadN $Thread \
     --outFilterMismatchNmax 6 \
     --outSAMtype None \
     --outFilterMultimapNmax 20 \
     --outFilterIntronMotifs RemoveNoncanonical \
     --quantMode - \
     --outFileNamePrefix ${File}_ \
     --outReadsUnmapped Fastx

mv ${File}_Unmapped.out.mate1 ${File}_rmHost.fq

kraken2 --db ./minikraken_8GB_20200312 \
        --threads $Thread \
	--report ${File}.kreport2 \
	--classified-out ${File}_rmHost_kraken.fq \
	${File}_rmHost.fq \
	1> ${File}.kraken2

awk '$1=="C"' ${File}.kraken2 > ${File}.kraken

python PathogenTrack.py -b ${File}/barcodes.tsv \
                        -i ${File}_rmHost_kraken.fq \
			-k ${File}.kraken \
			-t taxons.db \
			-f $File \
			-o $File

echo "$File completes at `date`"
