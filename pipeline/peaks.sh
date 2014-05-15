#! /usr/bin/env bash

#BSUB -J peaks[1-10]
#BSUB -e peaks.%J.%I.err
#BSUB -o peaks.%J.%I.out
#BSUB -q short

<<DOC
call peaks using macs2.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/devel/modmap/pipeline/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
results=$RESULT/$sample
peakresults=$results/peaks

if [[ ! -f $peakresults ]]; then
    mkdir -p $peakresults
fi

strands=("all" "pos" "neg")

for strand in ${strands[@]}; do
    for align_mode in ${ALIGN_MODES[@]}; do

        peakbase=$peakresults.$strand
        peak=${peakbase}_peaks.bed
        narrowpeak=${peakbase}_peaks.narrowPeak
        summit=${peakbase}_summits.bed
        xls=${peakbase}_peaks.xls

        # some peaks were extending outside of genomic coords
        clipped_peak=${peakbase}_peaks.bed.clipped

        # combined peaks with appropriate strand column
        peak=$result/${sample}_peaks.bed.gz
   
        bam=$results/alignment/$sample.align.$align_mode.bam

        if [[ ! -f $bam ]]; then
            echo "bam not found for $sample"
            exit 1
        fi

        if [[ ! -f $peak.gz ]]; then
            macs2 callpeak -t $bam -n $peakbase --keep-dup auto \
                --nomodel -s 25 --extsize 5 --call-summits
            cut -f1-4 $narrowpeak > $peak
            bedClip $peak $CHROM_SIZES $clipped_peak
            mv $clipped_peak $peak
            gzip -f $peak $narrowpeak
            rm -f $xls $summit
        fi
    done
done
