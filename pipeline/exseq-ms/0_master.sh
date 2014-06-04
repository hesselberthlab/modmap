#!/usr/bin/env bash

#BSUB -J master
#BSUB -e master.%J.err
#BSUB -o master.%J.out
#BSUB -q normal
#BSUB -P exseq

<<DOC
master analysis loop for exset manuscript pipeline
DOC

set -o nounset -o pipefail -o errexit -x

export ASSEMBLIES=("sacCer1" "sacCer2" "sacCer3")
export PIPELINE=$HOME/devel/modmap/pipeline/exseq-ms
export CONFIG=$PIPELINE/config.sh

for assembly in ${ASSEMBLIES[@]}; do

    source $CONFIG

    # reassign assembly-specific variables
    export ASSEMBLY=$assembly
    # XXX DEBUG provides a directory extension ("common-debug"), or
    # nothing ("common")
    export RESULT=$HOME/projects/exseq-ms/results/common$DEBUG/$assembly

    export BOWTIEIDX=$HOME/ref/genomes/$assembly/$assembly
    export CHROM_SIZES=$HOME/ref/genomes/$assembly/$assembly.chrom.sizes
    export GTF=$HOME/ref/genomes/$assembly/sgdGene.$assembly.gtf
    export FASTA=$HOME/ref/genomes/$assembly/$assembly.fa
    export GDARCHIVE="$RESULT/genomedata/exseq.ms.genomedata"
    export BKGD_FREQS="$HOME/ref/genomes/$assembly/$assembly.genome.nuc.freqs.tab"
    
    job_array="[1-$NUM_SAMPLES]"

    bsub -J "align_$ASSEMBLY$job_array" \
        < $PIPELINE/1_align.sh

    bsub -J "coverage_$ASSEMBLY$job_array" \
        -w "done('align_$ASSEMBLY[*]')" \
        < $PIPELINE/2_coverage.sh 

    bsub -J "nuc_freqs_$ASSEMBLY$job_array" \
        -w "done('align_$ASSEMBLY[*]')" \
        < $PIPELINE/3_nuc_freqs.sh

    bsub -J "exp_anal_$ASSEMBLY$job_array" \
        -w "done('align_$ASSEMBLY[*]')" \
        < $PIPELINE/expression_analysis.sh

    bsub -J "plots_$ASSEMBLY$job_array" \
        -w "done('nuc_freqs_$ASSEMBLY[*]') && \
            done('exp_anal_$ASSEMBLY[*]')" \
        < $PIPELINE/5_plots.sh

    #bsub -J "tracklines_$ASSEMBLY" \
    #    < $PIPELINE/tracklines.sh

done
