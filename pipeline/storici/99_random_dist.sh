#! /usr/bin/env bash

#BSUB -J random.dist[1-13]
#BSUB -e random.dist.%J.%I.err
#BSUB -o random_dist.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
DOC

set -o nounset -o pipefail -o errexit -x

# XXX testing
#CONFIG=$HOME/devel/modmap/pipeline/storici/config.sh
#ASSEMBLY=sacCer2
#DEBUG="-debug"
#RESULT=$HOME/projects/collab/storici-lab/results/common$DEBUG/$ASSEMBLY
#CHROM_SIZES=$HOME/ref/genomes/$ASSEMBLY/$ASSEMBLY.chrom.sizes
#LSB_JOBINDEX=1
# XXX

source $CONFIG
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

# output directory
results=$RESULT/$sample/random_distribution
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

region_types=('nuc' 'mito')
region_args=('--ignore-chrom chrM' '--only-chrom chrM')
              
# XXX need to run module out of bin directory
cd $BIN

bedgraphdir=$RESULT/$sample/bedgraphs

interval_sizes=(500 1000 2500 5000 10000)

for align_mode in ${ALIGN_MODES[@]}; do

    bedgraph="$bedgraphdir/$sample.align.$align_mode.strand.all.counts.bg"

    for interval_size in ${interval_sizes[@]}; do
        result_tab="$results/random_dist.align.$align_mode.interval.$interval_size.tab"
        if [[ -f $result_tab ]]; then
            rm -f $result_tab
        fi

        for region_idx in ${!region_types[@]}; do

            region_type=${region_types[$region_idx]}
            region_arg=${region_args[$region_idx]}

            python -m modmap.random_dist \
                $bedgraph $CHROM_SIZES \
                --region-type $region_type \
                $region_arg \
                --interval-size $interval_size \
                --verbose \
                >> $result_tab
        done
    done
done
