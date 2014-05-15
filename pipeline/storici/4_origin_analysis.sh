#! /usr/bin/env bash

#BSUB -J ori.analysis[1-10]
#BSUB -e ori.analysis.%J.%I.err
#BSUB -o ori.anlaysis.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
origin analysis
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

# output directory
results=$RESULT/$sample/origin_analysis
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

# origins and timing
origin_bed=$DATA/$ASSEMBLY/oridb.confirmed.bed
timing_bedgraph=$DATA/$ASSEMBLY/yabuki.timing.bedgraph

# origin variables
timing_cutoffs=(20 25 30 35)
flank_sizes=(1000 2500 5000 10000)

# XXX need to run module out of bin directory
cd $BIN

aligndir=$RESULT/$sample/alignment

for align_mode in ${ALIGN_MODES[@]}; do
    BAM=$aligndir/$sample.align.$align_mode.bam
    result_tab="$results/origin_analysis.align.$align_mode.tab"

    # delete old results if exists - going to loop and append so need a
    # fresh empty file
    if [[ -f $result_tab ]]; then
        rm -f $result_tab
    fi

    for timing in ${timing_cutoffs[@]}; do
        for flank_size in ${flank_sizes[@]}; do
            python -m modmap.origin_analysis \
                $origin_bed $timing_bedgraph \
                $BAM $FASTA $CHROM_SIZES \
                --max-timing $timing \
                --flank-size $flank_size \
                --verbose \
                >> $result_tab
        done
    done
done
