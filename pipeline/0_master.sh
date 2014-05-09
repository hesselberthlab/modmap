#!/usr/bin/env bash
#BSUB -J master
#BSUB -e master.%J.err
#BSUB -o master.%J.out
#BSUB -q normal
#BSUB -P storici

<<DOC
DOC

set -o nounset -o pipefail -o errexit -x

ASSEMBLIES=("sacCer1" "sacCer2" "sacCer3")
PIPLINE=$HOME/devel/modmap/pipeline

for assembly in ${ASSEMBLIES[@]}; do

    source $HOME/projects/collab/storici-lab/bin/config.sh

    export BOWTIEIDX=$HOME/ref/genomes/$assembly/$assembly
    export CHROM_SIZES=$HOME/ref/genomes/$assembly/$assembly.chrom.sizes
    export GTF=$HOME/ref/genomes/$assembly/sgdGene.$assembly.gtf
    export FASTA=$HOME/ref/genomes/$assembly/$assembly.fa
    export RESULT=$HOME/projects/collab/storici-lab/results/common/$assembly

    bsub -J align.dep \
        < $PIPELINE/1_align.sh

    bsub -J coverage.dep -w align.dep \
        < $PIPELINE/2_coverage.sh 

    bsub -J nuc.freqs.dep -w coverage.dep \
        < $PIPELINE/3_nuc_freqs.sh

    bsub -J origin.anal.dep -w nuc.freqs.dep \
        < $PIPELINE/4_origin_analysis.sh

    bsub -J plots.dep -w origin.anal.dep \
        < $PIPELINE/5_plots.sh
done
