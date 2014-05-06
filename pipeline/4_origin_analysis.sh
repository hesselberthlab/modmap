#!/usr/bin/env bash

#BSUB -J ori.analysis[1-10]
#BSUB -e ori.analysis.%J.%I.err
#BSUB -o ori.anlaysis.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
origin analysis
1. select origin flanks for analysis based on timing data
2. analyze signal in flanks and assign to leading or laggin
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/collab/storici-lab/bin/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

# count data
bedgraphs=$RESULT/$sample/bedgraphs
# output directory
results=$RESULT/$sample/origin_analysis

# origins and timing
oribed=$DATA/$ASSEMBLY/oridb.confirmed.bed
timingbg=$DATA/$ASSEMBLY/yabuki.timing.bedgraph

# XXX origin variables
max_timing=25
flank_size=10000

strands=("pos" "neg")
# XXX: upstream and downstream maybe more acc
sides=("left" "right")
directions=("leading" "lagging")

if [[ $ASSEMBLY == "sacCer2" ]]; then
    ignore_modes=("all" "only-mito" "no-mito" "only-2micron")
    ignore_args=("" "--only-chrom chrM"
                 "--ignore-chrom chrM"
                 "--only-chrom 2micron")
else
    ignore_modes=("all" "only-mito" "no-mito")
    ignore_args=("" "--only-chrom chrM"
                 "--ignore-chrom chrM")
fi

if [[ ! -d $results ]]; then
    mkdir -p $results
fi

# selecte origins for analysis
selectoritab="$results/selected.origins.tab"
selectoribed="$results/selected.origins.bed"
bedtools intersect -wb -a $oribed -b $timingbg \
    | bedtools groupby -g 4 -c 10 -o mean \
    | awk -v TIMING=$max_timing '$2 < TIMING' \
    | cut -f1 \
    > $selectoritab

for oriname in $(cat $selectoritab); do
    cat $oribed \
    | awk -v NAME=$oriname '$4 == NAME'
done > $selectoribed

# make flanks for origins - writing loop not worth it
left_origins="$results/origins.left.flank.bed"
right_origins="$results/origins.right.flank.bed"
bedtools flank -i $selectoribed -l $flank_size -r 0 -g $CHROM_SIZES > $left_origins
bedtools flank -i $selectoribed -r $flank_size -l 0 -g $CHROM_SIZES > $right_origins

# analysis is run for each alignment mode
for align_mode in ${ALIGN_MODES[@]}; do

    # generate signal for each strand and side
    for strand in ${strands[@]}; do
        signal_bg="$bedgraphs/$sample.align.$align_mode.strand.$strand.counts.bg"

        for side in ${sides[@]}; do
            origin_bed="$results/origins.$side.flank.bed"
            result_bg="$results/$sample.align.$align_mode.strand.$strand.side.$side.origin.counts.bg"

            bedtools intersect -a $signal_bg -b $origin_bed > $result_bg
        done
    done

    # analyze signal on leading and lagging directions
    for ignore_idx in ${!ignore_modes[@]}; do

        ignore_mode=${ignore_modes[$ignore_idx]}
        ignore_arg=${ignore_args[$ignore_idx]}

        for direction in ${directions[@]}; do

            if [[ $direction == "lagging" ]]; then
                # - lagging strand = -p $pos_left_bedgraph -n $neg_right_bedgraph
                pos_bg="$results/$sample.align.$align_mode.strand.pos.side.left.origin.counts.bg"
                neg_bg="$results/$sample.align.$align_mode.strand.neg.side.right.origin.counts.bg"
                result_tab="$results/lagging.align.$align_mode.ignore.$ignore_mode.origin.nuc_counts.tab"

            elif [[ $direction == "leading" ]]; then
                # - leading strand = -p $pos_right_bedgraph -n $neg_left_bedgraph
                pos_bg="$results/$sample.align.$align_mode.strand.pos.side.right.origin.counts.bg"
                neg_bg="$results/$sample.align.$align_mode.strand.neg.side.left.origin.counts.bg"
                result_tab="$results/leading.align.$align_mode.ignore.$ignore_mode.origin.nuc_counts.tab"
            fi

            # delete old results if exists
            if [[ -f $result_tab ]]; then
                rm -f $result_tab
            fi

            # offsets are with respect to base in question:
            # 0 = ribo, -1 = upstream, 1 = downstream
            python $BIN/nuc_frequencies.py \
                --offset-min -1 --offset-max 1 --region-size 1 \
                -p $pos_bg -n $neg_bg -f $FASTA \
                $ignore_arg \
                | awk -v ID=$sample -v DIR=$direction \
                    '{print $0, "\t", DIR, "\t", ID}' \
                > $result_tab
        done
        # combine results into 1 file
        combined="$results/combined.align.$align_mode.ignore.$ignore_mode.tab"
        cat $results/*align.$align_mode.ignore.$ignore_mode.*.tab >> $combined
    done
done
