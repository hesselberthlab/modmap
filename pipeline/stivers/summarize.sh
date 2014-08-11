#!/usr/bin/env bash

#BSUB -J tab
#BSUB -o %J.out
#BSUB -e %J.err

set -o nounset -o pipefail -o errexit -x

samples=(JS33 JS34 JS35 JS36 JS37 JS38 JS39 JS40
         JS41 JS42 JS43 JS44 JS45 JS46 JS47 JS48
         JS49 JS50 JS51 JS52 JS53 JS54 JS55 JS56)

unenriched=(JS33 JS34 JS35 JS36 JS37 JS38 JS39 JS40)
udg=(JS49 JS50 JS51 JS51 JS53 JS54 JS55 JS56)

PROJECT="$HOME/projects/collab/stivers-lab"
RESULT=$PROJECT/results/common
out=$RESULT/combined.tab.gz

index_type="virus"
ALIGN_MODES=("uniq" "all")
strands=("pos" "neg")

for sample in ${samples[@]}; do

    samplenum=$(echo $sample | sed 's/JS//')
    if [[ $samplenum -gt 40 ]]; then
        hyb="enriched"
    else
        hyb="unenriched"
    fi
    if [[ $samplenum -gt 48 ]]; then
        treatment="udg"
    else
        treatment="none"
    fi

    if [[ $samplenum -eq 33 || $samplenum -eq 41 || $samplenum -eq 49 ]]; then
        expname="1:-V"
    elif [[ $samplenum -eq 34 || $samplenum -eq 42 || $samplenum -eq 50 ]]; then
        expname="2:RTX -V"
    elif [[ $samplenum -eq 35 || $samplenum -eq 43 || $samplenum -eq 51 ]]; then
        expname="3:V-WT"
    elif [[ $samplenum -eq 36 || $samplenum -eq 44 || $samplenum -eq 52 ]]; then
        expname="4:RTX V-WT"
    elif [[ $samplenum -eq 37 || $samplenum -eq 45 || $samplenum -eq 53 ]]; then
        expname="5:V-D64E"
    elif [[ $samplenum -eq 38 || $samplenum -eq 46 || $samplenum -eq 54 ]]; then
        expname="6:RTX V-D64E"
    elif [[ $samplenum -eq 39 || $samplenum -eq 47 || $samplenum -eq 55 ]]; then
        expname="7:RAL V-WT"
    elif [[ $samplenum -eq 40 || $samplenum -eq 48 || $samplenum -eq 56 ]]; then
        expname="8:RTX RAL V-WT"
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

