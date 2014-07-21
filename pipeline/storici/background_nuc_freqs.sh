#! /usr/bin/env bash

#BSUB -J background.nuc.freqs[1-3]
#BSUB -e bkgd.nuc.freqs.%J.%I.err
#BSUB -o bkgd.nuc.freqs.%J.%I.out
#BSUB -q normal

<<DOC
Calculate nucleotide frequencies
DOC

set -o nounset -o pipefail -o errexit -x

# need to be in BIN to run module
BIN=$HOME/devel/modmap
cd $BIN

ASSEMBLIES=(sacCer1 sacCer2 sacCer3)
assembly=${ASSEMBLIES[$(($LSB_JOBINDEX - 1))]}

output="$HOME/ref/genomes/$assembly/$assembly.genome.nuc.freqs.tab"
FASTA="$HOME/ref/genomes/$assembly/$assembly.fa"

if [[ ! -f $output ]]; then
    python -m modmap.genome_nuc_freqs \
        $FASTA \
        --region-size-minimum 1 \
        --region-size-maximum 3 \
        --verbose \
        > $output
fi

output="$HOME/ref/genomes/$assembly/$assembly.chrM.nuc.freqs.tab"
FASTA="$HOME/ref/genomes/$assembly/$assembly.chrM.fa"

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
    output="$HOME/ref/genomes/$assembly/$assembly.2micron.nuc.freqs.tab"
    FASTA="$HOME/ref/genomes/$assembly/$assembly.2micron.fa"

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
