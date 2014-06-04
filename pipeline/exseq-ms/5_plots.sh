#! /usr/bin/env bash

#BSUB -J plots[1-10]
#BSUB -e plots.%J.%I.err
#BSUB -o plots.%J.%I.out
#BSUB -q normal

<<DOC
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
results=$RESULT/$sample

plotdir=$results/plots
if [[ $ASSEMBLY == "sacCer2" ]]; then
    include_modes=("all" "only-mito" "no-mito" "only-2micron")
else
    include_modes=("all" "only-mito" "no-mito")
fi

# count thresholds
count_thresh="1 10 30"

for aln_idx in ${!ALIGN_MODES[@]}; do

    align_mode=${ALIGN_MODES[$aln_idx]}

    # --- expression plots --------------------------------------
    subplotdir="$plotdir/expression_analysis"
    if [[ ! -d $subplotdir ]]; then
        mkdir -p $subplotdir
    fi

    region_types=('mrna' 'promoters')
    for region_type in ${region_types[@]}; do
        counts="$results/expression_analysis/exp_analysis.$region_type.align.$align_mode.tab"
        sampleid="$sample.align-$align_mode.region-$region_type"
        Rscript $RSCRIPTS/exp.plots.R $counts "$sampleid" $subplotdir
    done

    for mincount in $count_thresh; do
        for ig_idx in ${!include_modes[@]}; do

            include_mode=${include_modes[$ig_idx]}

            # --- nuc_freq plots ------------------------------------
            subplotdir="$plotdir/nuc_freqs"
            if [[ ! -d $subplotdir ]]; then
                mkdir -p $subplotdir
            fi

            counts="$results/nuc_freqs/$sample.align.$align_mode.include.$include_mode.mincount.$mincount.nuc_freqs.tab"
            sampleid="$sample.align-$align_mode.subset-$include_mode.mincount-$mincount"
            Rscript $RSCRIPTS/nuc.freqs.R $counts "$sampleid" $subplotdir
        done
    done
done

# convert all PDFs to PNGs
for pdffile in $(find $plotdir -name '*.pdf' -print); do
    pngdirname=$(dirname $pdffile)
    pngfilename="$(basename $pdffile .pdf).png"
    pngfilepath="$pngdirname/$pngfilename"
    convert $pdffile $pngfilepath
done

