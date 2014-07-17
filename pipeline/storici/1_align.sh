#!/usr/bin/env bash

#BSUB -J align[1-13]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -P storici
#BSUB -n 6

<<DOC
Trim the UMI from the FASTQ, align trimmed reads using bowtie suppressing 
all reads that align more than once, then remove UMI duplicates from the
alignment.
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

unprocessed_fastq=$DATA/$sample.fq.gz
fastq=$DATA/$sample.umi.fq.gz

# trim the UMI
if [[ ! -f $fastq ]]; then
    umitools trim $unprocessed_fastq $UMI \
        | gzip -c \
        > $fastq
fi

results=$RESULT/$sample/alignment
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

for idx in ${!ALIGN_MODES[@]}; do

    align_mode=${ALIGN_MODES[$idx]}
    align_arg=${ALIGN_ARGS[$idx]}

    umibam=$results/$sample.UMIs_not_removed.align.$align_mode.bam
    bam=$results/$sample.align.$align_mode.bam
    stats=$results/$sample.align.$align_mode.alignment_stats.txt

    # align the reads
    if [[ ! -f $umibam ]]; then
        zcat $fastq \
            | bowtie $align_arg --sam $BOWTIEIDX -p 6 - \
            2> $stats \
            | samtools view -ShuF4 - \
            | samtools sort -o - $sample.temp -m 8G \
            > $umibam
        samtools index $umibam
    fi

    # process the UMIs
    umibed="$results/$sample.umi-report.align.$align_mode.bed.gz"
    if [[ ! -f $bam ]]; then
        umitools rmdup $umibam $bam | gzip -c > $umibed
        samtools index $bam
    fi
done

