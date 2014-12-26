#! /usr/bin/env bash

#BSUB -J peaks.summary.table
#BSUB -e peaks.summary.table.%J.err
#BSUB -o peaks.summary.table.%J.out
#BSUB -q normal
#BSUB -P storici

<<DOC
generate summary table
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG
cd $BIN

results=$RESULT/peaks_summary_table
if [[ ! -d $results ]]; then
    mkdir $results
fi

strands=("pos" "neg")
for strand in ${strands[@]}; do
    for align_mode in ${ALIGN_MODES[@]}; do

        samples=""
        names=""
        for sample in ${SAMPLES[@]}; do
            peakdir=$RESULT/$sample/peaks
            exp_name="$sample.align.$align_mode.strand.$strand"
            peakbed="$peakdir/$exp_name""_peaks.narrowPeak"
            samples="$samples $peakbed"
            names="$names $sample"
        done

        result="$results/peaks_summary_table.align.$align_mode.strand.$strand.tab"

        bedtools multiinter -i $samples -names $names \
            -g $CHROM_SIZES \
            > $result
    done
done
