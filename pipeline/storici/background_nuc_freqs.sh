#! /usr/bin/env bash

#BSUB -J background.nuc.freqs[1-3]
#BSUB -e bkgd.nuc.freqs.%J.%I.err
#BSUB -o bkgd.nuc.freqs.%J.%I.out
#BSUB -q normal

<<DOC
Calculate nucleotide frequencies
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG

# need to be in BIN to run module
BIN=$HOME/devel/modmap
cd $BIN

assembly=${ASSEMBLIES[$(($LSB_JOBINDEX - 1))]}

output="$RESULT/$assembly.genome.nuc.freqs.tab"
if [[ ! -f $output ]]; then
    python -m modmap.genome_nuc_freqs \
        $FASTA \
        --region-size-minimum 1 \
        --region-size-maximum 3 \
        --verbose \
        > $output
fi

output="$RESULT/$assembly.chrM.nuc.freqs.tab"
if [[ ! -f $output ]]; then
    python -m modmap.genome_nuc_freqs \
        $FASTA \
        --region-size-minimum 1 \
        --region-size-maximum 3 \
        --only-chrom chrM \
        --verbose \
        > $output
fi

if [[ $assembly == "sacCer2" ]]; then
    output="$RESULT/$assembly.2micron.nuc.freqs.tab"
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
