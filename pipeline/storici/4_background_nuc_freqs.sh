#! /usr/bin/env bash

#BSUB -J background.nuc.freqs
#BSUB -e bkgd.nuc.freqs.%J.err
#BSUB -o bkgd.nuc.freqs.%J.out
#BSUB -q normal

<<DOC
Calculate nucleotide frequencies
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG

# need to be in BIN to run module
BIN=$HOME/devel/modmap
cd $BIN

output="$RESULT/$ASSEMBLY.genome.nuc.freqs.tab"
if [[ ! -f $output ]]; then
    python -m modmap.genome_nuc_freqs \
        $FASTA \
        --region-size-minimum 1 \
        --region-size-maximum 3 \
        --verbose \
        > $output
fi

output="$RESULT/$ASSEMBLY.chrM.nuc.freqs.tab"
if [[ ! -f $output ]]; then
    python -m modmap.genome_nuc_freqs \
        $FASTA \
        --region-size-minimum 1 \
        --region-size-maximum 3 \
        --only-chrom chrM \
        --verbose \
        > $output
fi

if [[ $ASSEMBLY == "sacCer2" ]]; then
    output="$RESULT/$ASSEMBLY.2micron.nuc.freqs.tab"
    if [[ ! -f $output ]]; then
        python -m modmap.genome_nuc_freqs \
            $FASTA \
            --region-size-minimum 1 \
            --region-size-maximum 3 \
            --only-chrom 2micron \
            --verbose \
            > $output
    fi    
fi
