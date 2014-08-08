#! /usr/bin/env bash

#BSUB -J dist.analysis
#BSUB -e dist.analysis.%J.err
#BSUB -o dist.analysis.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
Generate data for distribution analysis
DOC

set -o nounset -o pipefail -o errexit -x

# XXX debugging
ASSEMBLY="sacCer1"
DEBUG="-debug"
RESULT=$HOME/projects/collab/storici-lab/results/common$DEBUG/$ASSEMBLY
PIPELINE=$HOME/devel/modmap/pipeline/storici
CONFIG=$PIPELINE/config.sh
# XXX

source $CONFIG

if [[ $ASSEMBLY == "sacCer2" ]]; then
    ignore_modes=("nuclear" "only-mito" "only-2micron")
else
    ignore_modes=("nuclear" "only-mito")
fi

results=$RESULT/distribution_analysis
if [[ ! -d $results ]]; then
    mkdir -p $results
fi
strands=("pos" "neg")
output="$results/dist_data.tab"

if [[ -f $output ]]; then
    rm -f $output
fi

for sample in ${SAMPLES[@]}; do

    bgresults=$RESULT/$sample/bedgraphs

    for align_mode in ${ALIGN_MODES[@]}; do
        for strand in ${strands[@]}; do

            counts=$bgresults/$sample.align.$align_mode.strand.$strand.counts.bg

            if [[ ! -f $counts ]]; then
                echo "counts file not found: $counts"
                continue
            fi

            for ignore_mode in ${ignore_modes[@]}; do

                if [[ $ignore_mode == "only-mito" ]]; then
                    chrom="mito"
                    awk '$1 == "chrM"' < $counts \
                        | awk -v chrom=$chrom -v strand=$strand \
                              -v align=$align_mode -v sample=$sample \
                                    'BEGIN {OFS="\t"}
                                     {print $4, chrom, align,
                                     sample, strand}'

                elif [[ $ignore_mode == "only-2micron" ]]; then
                    chrom="2micron"
                    awk '$1 == "2micron"' < $counts \
                         | awk -v chrom=$chrom -v strand=$strand \
                              -v align=$align_mode -v sample=$sample \
                                    'BEGIN {OFS="\t"}
                                     {print $4, chrom, align,
                                     sample, strand}'
                elif [[ $ignore_mode == "nuclear" ]]; then
                    chrom="nuclear"
                    grep '^chr' < $counts \
                        | awk '$1 != "2micron" && $1 != "chrM"' \
                        | awk -v chrom=$chrom -v strand=$strand \
                              -v align=$align_mode -v sample=$sample \
                                    'BEGIN {OFS="\t"}
                                     {print $4, chrom, align,
                                     sample, strand}'
                fi

            done
        done
    done
done
