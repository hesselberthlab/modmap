#!/usr/bin/env bash

#BSUB -J master
#BSUB -e master.%J.err
#BSUB -o master.%J.out
#BSUB -q normal
#BSUB -P storici

<<DOC
master analysis loop for storici mopmap pipeline
DOC

set -o nounset -o pipefail -o errexit -x

export ASSEMBLIES=("sacCer1" "sacCer2" "sacCer3")
export PIPELINE=$HOME/devel/modmap/pipeline/storici
export CONFIG=$PIPELINE/config.sh

for assembly in ${ASSEMBLIES[@]}; do

    source $CONFIG

    # reassign assembly-specific variables
    export ASSEMBLY=$assembly
    # XXX DEBUG provides a directory extension ("common-debug"), or
    # nothing ("common")
    export RESULT=$HOME/projects/collab/storici-lab/results/common$DEBUG/$assembly

    export BOWTIEIDX=$HOME/ref/genomes/$assembly/$assembly
    export CHROM_SIZES=$HOME/ref/genomes/$assembly/$assembly.chrom.sizes
    export FASTA=$HOME/ref/genomes/$assembly/$assembly.fa
    
    job_array="[1-$NUM_SAMPLES]"

    # these are annotation jobs that don't wait, run them first
    bsub -J "txn_regions_$ASSEMBLY" \
        < $PIPELINE/make_transcribed_regions.sh

    bsub -J "bkgd_freqs_$ASSEMBLY" \
        < $PIPELINE/4_background_nuc_freqs.sh

    bsub -J "tracklines_$ASSEMBLY" \
        < $PIPELINE/99_tracklines.sh

    # job names look like: align_sacCer1[1-10]
    bsub -J "align_$ASSEMBLY$job_array" \
        < $PIPELINE/1_align.sh

    bsub -J "coverage_$ASSEMBLY$job_array" \
        -w "done('align_$ASSEMBLY[*]')" \
        < $PIPELINE/2_coverage.sh 

    bsub -J "peaks_$ASSEMBLY$job_array" \
        -w "done('align_$ASSEMBLY[*]') && \
            done('coverage_$ASSEMBLY[*]')" \
        < $PIPELINE/peaks.sh

    # get job IDs
    align_jobid=$(bjobs -w | grep "align_$ASSEMBLY" | cut -f1 -d' ' | uniq)
    coverage_jobid=$(bjobs -w | grep "coverage_$ASSEMBLY" | cut -f1 -d' ' | uniq)
    peaks_jobid=$(bjobs -w | grep "peaks_$ASSEMBLY" | cut -f1 -d' ' | uniq)

    bsub -J "peaks_summary_$ASSEMBLY" \
        -w "numdone($peaks_jobid, == $NUM_SAMPLES)" \
        < $PIPELINE/peaks_summary_table.sh

    # wait until *all* coverage jobs are complete
    bsub -J "summary_table_$ASSEMBLY" \
        -w "numdone($coverage_jobid, == $NUM_SAMPLES)" \
        < $PIPELINE/3_summary_table.sh

    bsub -J "geo_bundle_$ASSEMBLY"  \
        -w "numdone($align_jobid, == $NUM_SAMPLES) && \
            numdone($coverage_jobid, == $NUM_SAMPLES)" \
        < $PIPELINE/99_geo_bundle.sh

    bsub -J "rand_dist_$ASSEMBLY$job_array" \
        -w "numdone($coverage_jobid, == $NUM_SAMPLES)" \
        < $PIPELINE/99_random_dist.sh

    bsub -J "origin_anal_$ASSEMBLY$job_array" \
        -w "numdone($align_jobid, == $NUM_SAMPLES)" \
        < $PIPELINE/6_origin_analysis.sh

    bsub -J "nuc_freqs_$ASSEMBLY$job_array" \
        -w "numdone($align_jobid, == $NUM_SAMPLES) && \
            done('bkgd_freqs_$ASSEMBLY')" \
        < $PIPELINE/5_nuc_freqs.sh

    bsub -J "exp_anal_$ASSEMBLY$job_array" \
        -w "numdone($align_jobid, == $NUM_SAMPLES)" \
        < $PIPELINE/7_expression_analysis.sh

    bsub -J "txn_anal_$ASSEMBLY$job_array" \
        -w "numdone($align_jobid, == $NUM_SAMPLES)" \
        < $PIPELINE/99_transcribed_analysis.sh

    # get job IDs
    nuc_freqs_jobid=$(bjobs -w | grep "nuc_freqs_$ASSEMBLY" | cut -f1 -d' ' | uniq)
    origin_anal_jobid=$(bjobs -w | grep "origin_anal_$ASSEMBLY" | cut -f1 -d' ' | uniq)
    exp_anal_jobid=$(bjobs -w | grep "exp_anal_$ASSEMBLY" | cut -f1 -d' ' | uniq)
    txn_anal_jobid=$(bjobs -w | grep "txn_anal_$ASSEMBLY" | cut -f1 -d' ' | uniq)

    bsub -J "plots_$ASSEMBLY$job_array" \
        -w "numended($origin_anal_jobid, == $NUM_SAMPLES) && \
            numended($nuc_freqs_jobid, == $NUM_SAMPLES) && \
            numended($exp_anal_jobid, == $NUM_SAMPLES) && \
            numended($txn_anal_jobid, == $NUM_SAMPLES)" \
        < $PIPELINE/8_plots.sh

done
