#!/usr/bin/env bash
#BSUB -J align[1-56]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -R "span[hosts=1]"
#BSUB -q normal
#BSUB -P stivers
#BSUB -n 12

<<DOC
Trim the UMI from the FASTQ, align trimmed reads using bowtie suppressing 
all reads that align more than once, then remove UMI duplicates from the
alignment.
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/devel/modmap/pipeline/stivers/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

fastq=$DATA/$sample.fq.gz

results=$RESULT/$sample/alignment
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

threads=12

for bwt_idx in ${!BOWTIEIDXS[@]}; do

    bwt_index=${BOWTIEIDXS[$bwt_idx]}    
    index_type=${BOWTIEIDX_TYPES[$bwt_idx]}

    for aln_idx in ${!ALIGN_MODES[@]}; do

        align_mode=${ALIGN_MODES[$aln_idx]}
        align_arg=${ALIGN_ARGS[$aln_idx]}

        if [[ $align_mode == "all" && $index_type == "hg19" ]]; then
            continue
        fi

        bam=$results/$sample.$index_type.$align_mode.bam
        stats=$results/$sample.$index_type.$align_mode.alignment.txt

        # align the reads
        if [[ ! -f $bam ]]; then
            zcat $fastq \
                | bowtie $align_arg -p $threads --sam $bwt_index - \
                2> $stats \
                | samtools view -ShuF4 - \
                | samtools sort -o - $results/$sample.temp -m 8G \
                > $bam
            samtools index $bam
        fi
    done
done
