#!/usr/bin/env bash

#BSUB -J coverage[1-10]
#BSUB -e coverage.%J.%I.err
#BSUB -o coverage.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
Calculate coverage from bam files
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/devel/modmap/pipeline/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

results=$RESULT/$sample
bgresults=$RESULT/$sample/bedgraphs
if [[ ! -d $bgresults ]]; then
    mkdir -p $bgresults
fi

umi_types=("removed" "UMIs_not_removed")

for umi_type in ${umi_types[@]}; do

    for align_mode in ${ALIGN_MODES[@]}; do

        if [[ $umi_type == 'removed' ]]; then
            bam=$results/alignment/$sample.align.$align_mode.bam
            countsbg=$bgresults/$sample.align.$align_mode.strand.both.counts.bg
            countsposbg=$bgresults/$sample.align.$align_mode.strand.pos.counts.bg
            countsnegbg=$bgresults/$sample.align.$align_mode.strand.neg.counts.bg
        else
            bam=$results/alignment/$sample.$umi_type.align.$align_mode.bam
            countsbg=$bgresults/$sample.$umi_type.align.$align_mode.strand.both.counts.bg
            countsposbg=$bgresults/$sample.$umi_type.align.$align_mode.strand.pos.counts.bg
            countsnegbg=$bgresults/$sample.$umi_type.align.$align_mode.strand.neg.counts.bg
        fi

        if [[ ! -f $countsbg ]]; then
            bedtools genomecov -5 -bg -g $CHROM_SIZES -ibam $bam \
                > $countsbg
            bedtools genomecov -5 -bg -g $CHROM_SIZES -ibam $bam \
                -strand "+" > $countsposbg
            bedtools genomecov -5 -bg -g $CHROM_SIZES -ibam $bam \
                -strand "-" > $countsnegbg
        fi
    done
done

# create bigwigs
for bgfile in $(ls $bgresults/*.bg); do
    bwfile="$bgresults/$(basename $bgfile .bg).bw"
    bedGraphToBigWig $bgfile $CHROM_SIZES $bwfile
done

