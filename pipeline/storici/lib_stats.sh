#! /usr/bin/env bash

#BSUB -J lib.stats 
#BSUB -e lib.stats.%J.err
#BSUB -o lib.stats.%J.out
#BSUB -q normal
#BSUB -P storici

<<DOC
generate library stats
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG
cd $BIN

results=$RESULT/library_stats
if [[ ! -d $results ]]; then
    mkdir $results
fi

# calc genome size for norm later
genomesize=$(awk '{N += $2} END {print N}' < $CHROM_SIZES)

for align_mode in ${ALIGN_MODES[@]}; do

    result="$results/library_stats.align.$align_mode.tab"

    for sample in ${SAMPLES[@]}; do
        bedgraph_dir=$RESULT/$sample/bedgraphs
        posbedgraph="$bedgraph_dir/$sample.align.$align_mode.strand.pos.counts.bg"
        negbedgraph="$bedgraph_dir/$sample.align.$align_mode.strand.neg.counts.bg"

        sampleid="$sample"

        cat $posbedgraph $negbedgraph \
            | awk '{SUM += $4} END {print SUM}'
    done

done
