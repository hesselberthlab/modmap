#!/usr/bin/env bash

#BSUB -J master
#BSUB -e master.%J.err
#BSUB -o master.%J.out
#BSUB -q normal
#BSUB -P stivers

<<DOC
master analysis loop for stivers mopmap pipeline
DOC

set -o nounset -o pipefail -o errexit -x

export PIPELINE=$HOME/devel/modmap/pipeline/stivers
export CONFIG=$PIPELINE/config.sh

source $CONFIG

# XXX DEBUG provides a directory extension ("common-debug"), or
# nothing ("common")
export RESULT=$HOME/projects/collab/stivers-lab/results/common$DEBUG

export BOWTIEIDX=$HOME/ref/genomes/$assembly/$assembly
export CHROM_SIZES=$HOME/ref/genomes/$assembly/$assembly.chrom.sizes
export GTF=$HOME/ref/genomes/$assembly/sgdGene.$assembly.gtf
export FASTA=$HOME/ref/genomes/$assembly/$assembly.fa

job_array="[1-$NUM_SAMPLES]"

bsub -J "align$job_array" \
    < $PIPELINE/1_align.sh

bsub -J "coverage$job_array" \
    -w "done('align[*]')" \
    < $PIPELINE/2_coverage.sh 

bsub -J "summarize" \
    -w "done('coverage[*]')" \
    < $PIPELINE/6_summary_table.sh

bsub -J "plots_$ASSEMBLY$job_array" \
    -w "done('origin_anal_$ASSEMBLY[*]') && \
        done('nuc_freqs_$ASSEMBLY[*]') && \
        done('exp_anal_$ASSEMBLY[*]')" \
    < $PIPELINE/5_plots.sh

