#!/usr/bin/env bash
#BSUB -J master[1-3]
#BSUB -e master.%J.%I.err
#BSUB -o master.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/collab/storici-lab/bin/config.sh

ASSEMBLIES=("sacCer1" "sacCer2" "sacCer3")

for assembly in ${ASSEMBLIES[@]}; do

    export BOWTIEIDX=$HOME/ref/genomes/$assembly/$assembly
    export CHROM_SIZES=$HOME/ref/genomes/$assembly/$assembly.chrom.sizes
    export GTF=$HOME/ref/genomes/$assembly/sgdGene.$assembly.gtf
    export FASTA=$HOME/ref/genomes/$assembly/$assembly.fa
    export RESULT=$HOME/projects/collab/storici-lab/results/common/$assembly

done
