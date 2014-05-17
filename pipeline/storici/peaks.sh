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

if [[ ! -f $peakresults ]]; then
    mkdir -p $peakresults
fi

strands=("all" "pos" "neg")

for strand in ${strands[@]}; do
    for align_mode in ${ALIGN_MODES[@]}; do

        exp_name=$peakresults/$sample.$strand
        peak=${exp_name}_peaks.bed
        narrowpeak=${exp_name}_peaks.narrowPeak
        summit=${exp_name}_summits.bed
        xls=${exp_name}_peaks.xls

        # some peaks were extending outside of genomic coords
        clipped_peak=${exp_name}_peaks.bed.clipped

        # combined peaks with appropriate strand column
        bam=$results/alignment/$sample.align.$align_mode.bam

        if [[ ! -f $bam ]]; then
            echo "bam not found for $sample"
            exit 1
        fi

        if [[ ! -f "$peak.gz" ]]; then
            macs2 callpeak -t $bam \
                -n $exp_name \
                --keep-dup all \
                --nomodel \
                -s 25 \
                --extsize 5 \
                --gsize $genomesize \
                --call-summits \
            cut -f1-4 $narrowpeak > $peak
            gzip -f $peak $narrowpeak
            # rm -f $xls $summit
        fi
    done
done
