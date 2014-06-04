#! /usr/bin/env bash

#BSUB -J exp.analysis[1-10]
#BSUB -e exp.analysis.%J.%I.err
#BSUB -o exp.anlaysis.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
origin analysis
DOC

set -o nounset -o pipefail -o errexit -x

# XXX testing
# CONFIG=$HOME/devel/modmap/pipeline/storici/config.sh
source $CONFIG
#ASSEMBLY=sacCer2
#RESULT=$HOME/projects/collab/storici-lab/results/common$DEBUG/$ASSEMBLY
#CHROM_SIZES=$HOME/ref/genomes/$ASSEMBLY/$ASSEMBLY.chrom.sizes

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

# output directory
results=$RESULT/$sample/expression_analysis
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

# origins and timing
region_types=('mrna' 'promoters')
region_filenames=("$DATA/$ASSEMBLY/exp.fpkm.bed"
                  "$DATA/$ASSEMBLY/promoters.250bp.fpkm.bed")

# XXX need to run module out of bin directory
cd $BIN

aligndir=$RESULT/$sample/alignments

for align_mode in ${ALIGN_MODES[@]}; do
    BAM=$aligndir/$sample.align.$align_mode.bam

    for region_idx in ${!region_types[@]}; do

        region_type=${region_types[$region_idx]}
        region_filename=${region_filenames[$region_idx]}

        result_tab="$results/exp_analysis.$region_type.align.$align_mode.tab"

        # delete old results if exists - going to loop and append so need a
        # fresh empty file
        if [[ -f $result_tab ]]; then
            rm -f $result_tab
        fi

        python -m modmap.signal_analysis \
            $BAM $region_filename $CHROM_SIZES \
            --verbose \
            >> $result_tab
    done
done
