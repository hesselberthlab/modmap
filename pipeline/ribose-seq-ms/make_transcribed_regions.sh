#! /usr/bin/env bash

#BSUB -J txn.regions
#BSUB -e txn.regions.%J.err
#BSUB -o txn.regions.%J.out
#BSUB -q normal
#BSUB -P storici

<<DOC
transcription analysis
DOC

set -o nounset -o pipefail -o errexit -x

# XXX
# source debug-config.sh

source $CONFIG

# output directory
results=$RESULT/transcribed_regions
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

genesbed=$DATA/$ASSEMBLY/sgdGene.bed

genes="$results/$(basename $genesbed .bed).bed"
comp_genes="$results/$(basename $genesbed .bed).complement.bed"

strands="+ -"

for strand in $strands; do

    awk -v strand=$strand '$6 == strand' \
        < $genesbed \
        | bedtools sort -i - \
        | bedtools merge -i - \
        | awk -v strand=$strand \
            'BEGIN {OFS="\t"} {print $0, ".", ".", strand}' \

done | bedtools sort -i - > $genes

bedtools complement -i $genes -g $CHROM_SIZES \
    | bedtools sort -i - \
    | awk 'BEGIN {OFS="\t"} {print $0, ".", ".", "."}' \
    > $comp_genes

# get mito regions
grep '^chrM' $genes \
    > "$results/$(basename $genes .bed).mito.bed"
grep '^chrM' $comp_genes \
    > "$results/$(basename $comp_genes .bed).mito.bed"

# get nuc regions
grep -v '^chrM' $genes \
    | grep -v '^2micron' \
    > "$results/$(basename $genes .bed).nuc.bed"
grep -v '^chrM' $comp_genes \
    | grep -v '^2micron' \
    > "$results/$(basename $comp_genes .bed).nuc.bed"

