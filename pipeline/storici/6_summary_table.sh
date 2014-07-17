#! /usr/bin/env bash

#BSUB -J summary.table
#BSUB -e summary.table.%J.err
#BSUB -o summary.table.%J.out
#BSUB -q normal
#BSUB -P storici

<<DOC
generate summary table
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG
results=$RESULT/summary_table
if [[ ! -d $results ]]; then
    mkdir $results
fi

strands=("pos" "neg")

# keep a running tally of the bedgraph files and sample IDs
bedgraphs=""
samplenames=""

for strand in ${strands[@]}; do
    for align_mode in ${!ALIGN_MODES[@]}; do

        result="$results/summary_table.align.$align_mode.strand.$strand.tab"

        for sample in $SAMPLES; do
            bedgraph_dir=$RESULT/$sample/bedgraphs
            bedgraph="$bedgraph_dir/$sample.align.$align_mode.strand.$strand.bg"
            bedgraphs="$bedgraphs $bedgraph"

            sampleid="$sample"
            samplenames="$samplenames $sampleid"
        done

        bedtools multiinter \
            -i $bedgraphs \
            -names $samplenames \
            -empty -g $CHROM_SIZES -header \
        > $result 
    done
done
