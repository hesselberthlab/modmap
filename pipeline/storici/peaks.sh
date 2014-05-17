#! /usr/bin/env bash

#BSUB -J peaks[1-10]
#BSUB -e peaks.%J.%I.err
#BSUB -o peaks.%J.%I.out
#BSUB -q short

<<DOC
call peaks using macs2.
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

results=$RESULT/$sample
peakresults=$results/peaks

# yeast genome size
genomesize=12e6

# narrowPeak autosql
ucscdir=/vol1/software/modules-sw/ucsc/build/v286
asfile=$ucscdir/kent/src/hg/lib/encode/narrowPeak.as

if [[ ! -f $peakresults ]]; then
    mkdir -p $peakresults
fi

strands=("all" "pos" "neg")

for strand in ${strands[@]}; do
    for align_mode in ${ALIGN_MODES[@]}; do

        exp_name=$peakresults/$sample.$strand
        peak=${exp_name}_peaks.bed
        narrowpeak=${exp_name}_peaks.narrowPeak
        bigbed=${exp_name}_peaks.bb

        # combined peaks with appropriate strand column
        bam=$results/alignment/$sample.align.$align_mode.bam

        if [[ ! -f $bam ]]; then
            echo "bam not found for $sample"
            exit 1
        fi

        macs2 callpeak -t $bam \
            -n $exp_name \
            --keep-dup all \
            --nomodel \
            -s 25 \
            --extsize 5 \
            --gsize $genomesize \
            --call-summits
        
        bedToBigBed -type=bed6+4 -as=$asfile \
            $narrowpeak $CHROM_SIZES $bigbed
    done
done
