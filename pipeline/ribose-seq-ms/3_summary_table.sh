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
cd $BIN

results=$RESULT/summary_table
if [[ ! -d $results ]]; then
    mkdir $results
fi

strands=("pos" "neg")

for strand in ${strands[@]}; do

    for align_mode in ${ALIGN_MODES[@]}; do

        # keep a running tally of the bedgraph files and sample IDs
        args=""

        result="$results/summary_table.align.$align_mode.strand.$strand.tab"

        for sample in ${SAMPLES[@]}; do
            bedgraph_dir=$RESULT/$sample/bedgraphs
            bedgraph="$bedgraph_dir/$sample.align.$align_mode.strand.$strand.counts.bg"

            sampleid="$sample"

            args="$args -s $sampleid=$bedgraph "
        done

        if [[ $strand == "neg" ]]; then
            args="$args --revcomp"
        fi

        python -m modmap.summary_table \
            --fasta $FASTA $args --verbose \
            > $result
        
    done
done
