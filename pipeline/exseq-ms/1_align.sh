#! /usr/bin/env bash

#BSUB -J align[1-10]
#BSUB -o align.%J.%I.out
#BSUB -e align.%J.%I.err
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 8

set -o nounset -o pipefail -o errexit -x

source $CONFIG

NUM_THREADS=8

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fastq="$FASTQDIR/$sample.fastq.gz"

results=$RESULT/$sample/alignments

if [[ ! -d $results ]]; then
    mkdir -p $results
fi

for idx in ${!ALIGN_MODES[@]}; do

    align_mode=${ALIGN_MODES[$idx]}
    align_arg=${ALIGN_ARGS[$idx]}

    bam=$results/$sample.align.$align_mode.bam
    stats=$results/$sample.stats.txt

    if [[ ! -f $bam ]]; then
        zcat $fastq \
            | bowtie $align_arg --sam $BOWTIEIDX -p $NUM_THREADS - \
            2> $stats \
            | samtools view -ShuF4 - \
            | samtools sort -o - $results/$sample.temp -m 8G \
            > $bam
        samtools index $bam
    fi

done
