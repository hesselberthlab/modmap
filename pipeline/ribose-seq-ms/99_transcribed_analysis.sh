#! /usr/bin/env bash

#BSUB -J txn.analysis[1-13]
#BSUB -e txn.analysis.%J.%I.err
#BSUB -o txn.anlaysis.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
transcription analysis
DOC

set -o nounset -o pipefail -o errexit -x

# XXX testing
# source debug-config.sh

source $CONFIG
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

# output directory
results=$RESULT/$sample/transcription_analysis
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

region_types=('nuc-genic' 'nuc-intergenic'
              'mito-genic' 'mito-intergenic')
              
region_filenames=("$RESULT/transcribed_regions/sgdGene.nuc.bed"
                  "$RESULT/transcribed_regions/sgdGene.complement.nuc.bed"
                  "$RESULT/transcribed_regions/sgdGene.mito.bed"
                  "$RESULT/transcribed_regions/sgdGene.complement.mito.bed")

# XXX need to run module out of bin directory
cd $BIN

aligndir=$RESULT/$sample/alignment

for align_mode in ${ALIGN_MODES[@]}; do

    BAM=$aligndir/$sample.align.$align_mode.bam

    result_tab="$results/txn_analysis.align.$align_mode.tab"
    if [[ -f $result_tab ]]; then
        rm -f $result_tab
    fi

    for region_idx in ${!region_types[@]}; do

        region_type=${region_types[$region_idx]}
        region_filename=${region_filenames[$region_idx]}

        python -m modmap.signal_analysis \
            $BAM $region_filename $CHROM_SIZES \
            --region-type $region_type \
            --normalize \
            --verbose \
            >> $result_tab
    done
done

