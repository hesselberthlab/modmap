#!/usr/bin/env bash

#BSUB -J coverage[1-56]%12
#BSUB -e coverage.%J.%I.err
#BSUB -o coverage.%J.%I.out
#BSUB -q normal
#BSUB -P stivers

<<DOC
Calculate coverage from bam files
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/devel/modmap/pipeline/stivers/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

results=$RESULT/$sample
bgresults=$RESULT/$sample/bedgraphs

if [[ ! -d $bgresults ]]; then
    mkdir -p $bgresults
fi

strands=(pos neg both)
strand_args=("-strand +" "-strand -" "")

for bwt_idx in ${!BOWTIEIDXS[@]}; do

    bwt_index=${BOWTIEIDXS[$bwt_idx]}    
    index_type=${BOWTIEIDX_TYPES[$bwt_idx]}

    for aln_idx in ${!ALIGN_MODES[@]}; do

        align_mode=${ALIGN_MODES[$aln_idx]}

        bam=$results/alignment/$sample.$index_type.$align_mode.bam
        stats=$results/alignment/$sample.$index_type.$align_mode.alignment.txt
        
        if [[ ! -f $bam ]]; then
            echo ">> bam file not found: $bam" 1>&2
            continue
        fi

        num_align=$(cat $stats | grep Reported | cut -f2 -d ' ')
        scale_pm=$(echo "1 / ($num_align / 1000000)" | bc -l)

        for strand_idx in ${!strands[@]}; do

            strand=${strands[$strand_idx]}
            strand_arg=${strand_args[$strand_idx]}

            bedgraph=$bgresults/$sample.$index_type.$align_mode.strand.$strand.counts.bg.gz
            normbedgraph=$bgresults/$sample.$index_type.$align_mode.strand.$strand.counts.norm.bg.gz
            tab=$bgresults/$sample.$index_type.$align_mode.strand.$strand.counts.tab.gz
            normtab=$bgresults/$sample.$index_type.$align_mode.strand.$strand.counts.norm.tab.gz

            if [[ ! -f $bedgraph ]]; then
                bedtools genomecov -bg -g $CHROM_SIZES \
                    -ibam $bam $strand_arg \
                    | gzip -c \
                    > $bedgraph
            fi

            if [[ ! -f $normbedgraph ]]; then
                bedtools genomecov -bg -g $CHROM_SIZES \
                    -ibam $bam $strand_arg -scale $scale_pm \
                    | gzip -c \
                    > $normbedgraph
            fi

            if [[ ! -f $tab ]]; then
                bedtools genomecov -dz -g $CHROM_SIZES \
                    -ibam $bam $strand_arg \
                    | gzip -c \
                    > $tab
            fi

            if [[ ! -f $normtab ]]; then
                bedtools genomecov -dz -g $CHROM_SIZES \
                    -ibam $bam $strand_arg -scale $scale_pm \
                    | gzip -c \
                    > $normtab
            fi

        done
    done
done


