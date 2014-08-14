#! /usr/bin/env bash

#BSUB -J txn.analysis[1-13]
#BSUB -e txn.analysis.%J.%I.err
#BSUB -o txn.anlaysis.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
origin analysis
DOC

set -o nounset -o pipefail -o errexit -x

# XXX testing
CONFIG=$HOME/devel/modmap/pipeline/storici/config.sh
source $CONFIG
ASSEMBLY=sacCer2
RESULT=$HOME/projects/collab/storici-lab/results/common$DEBUG/$ASSEMBLY
CHROM_SIZES=$HOME/ref/genomes/$ASSEMBLY/$ASSEMBLY.chrom.sizes
LSB_JOBINDEX=1
# XXX

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

# output directory
results=$RESULT/$sample/transcription_analysis
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

region_types=('transcribed' 'not-transcribed')
region_filenames=("$DATA/$ASSEMBLY/sgdGene.bed"
                  "$DATA/$ASSEMBLY/sgdGene.complement.bed")

# make complement bed if required
bedfile=${region_filenames[0]}
complementbed=${region_filenames[1]}
if [[ ! -f $complementbed ]]; then
    bedtools complement -i $bedfile -g $CHROM_SIZES \
        > $complementbed 
fi

# XXX need to run module out of bin directory
cd $BIN

aligndir=$RESULT/$sample/alignment

for align_mode in ${ALIGN_MODES[@]}; do

    BAM=$aligndir/$sample.align.$align_mode.bam

    result_tab="$results/txn_analysis.align.$align_mode.tab"

    # delete old results if exists - going to loop and append so need a
    # fresh empty file
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

