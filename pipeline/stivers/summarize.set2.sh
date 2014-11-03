#!/usr/bin/env bash

#BSUB -J tab
#BSUB -o %J.out
#BSUB -e %J.err

set -o nounset -o pipefail -o errexit -x

samples=(JS57 JS58 JS59 JS60 JS61 JS62)

PROJECT="$HOME/projects/collab/stivers-lab"
RESULT=$PROJECT/results/common
out=$RESULT/combined.set2.tab.gz

index_type="virus"
ALIGN_MODES=("uniq" "all")
strands=("pos" "neg")

for sample in ${samples[@]}; do

    samplenum=$(echo $sample | sed 's/JS//')
    if [[ $samplenum -gt 63 ]]; then
        treatment="udg"
    else
        treatment="none"
    fi

    if [[ $samplenum -eq 57 || $samplenum -eq 63 ]]; then
        expname="1:-V"
    elif [[ $samplenum -eq 58 || $samplenum -eq 64 ]]; then
        expname="2:-V RTX"
    elif [[ $samplenum -eq 59 || $samplenum -eq 65 ]]; then
        expname="3:+V"
    elif [[ $samplenum -eq 60 || $samplenum -eq 66 ]]; then
        expname="4:+V RTX"
    elif [[ $samplenum -eq 61 || $samplenum -eq 67 ]]; then
        expname="5:+V RAL"
    elif [[ $samplenum -eq 62 || $samplenum -eq 68 ]]; then
        expname="6:+V RAL RTX"
    fi

    bgresults=$RESULT/$sample/bedgraphs

    for aln_mode in ${ALIGN_MODES[@]}; do
        for strand in ${strands[@]}; do

            countstab=$bgresults/$sample.$index_type.$aln_mode.strand.$strand.counts.tab.gz
            normcountstab=$bgresults/$sample.$index_type.$aln_mode.strand.$strand.counts.norm.tab.gz

            zcat $countstab \
                | awk -v strand=$strand -v mode=$aln_mode \
                    -v sample=$sample -v treatment=$treatment \
                    -v hyb=$hyb -v expname="$expname" -v scale="raw" \
                    'BEGIN {OFS="\t"} {print $2, $3, sample, strand, mode,
                    treatment, hyb, expname, scale}'

             zcat $normcountstab \
                | awk -v strand=$strand -v mode=$aln_mode \
                    -v sample=$sample -v treatment=$treatment \
                    -v hyb=$hyb -v expname="$expname" -v scale="norm" \
                    'BEGIN {OFS="\t"} {print $2, $3, sample, strand, mode,
                    treatment, hyb, expname, scale}'
       done
    done
done | gzip -c > $out 

