#! /usr/bin/env bash

#BSUB -J plots[1-13]
#BSUB -e plots.%J.%I.err
#BSUB -o plots.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
results=$RESULT/$sample

plotdir=$results/plots

if [[ $ASSEMBLY == "sacCer2" ]]; then
    ignore_modes=("all" "only-mito" "no-mito" "only-2micron")
    ignore_args=("" "--only-chrom chrM"
                 "--ignore-chrom chrM"
                 "--only-chrom 2micron")
else
    ignore_modes=("all" "only-mito" "no-mito")
    ignore_args=("" "--only-chrom chrM"
                 "--ignore-chrom chrM")
fi

# for nuc freq plots
offset_maxs=(100 50 15)

for aln_idx in ${!ALIGN_MODES[@]}; do

    align_mode=${ALIGN_MODES[$aln_idx]}

    # --- origin plots --------------------------------------
    subplotdir="$plotdir/origin_analysis"
    if [[ ! -d $subplotdir ]]; then
        mkdir -p $subplotdir
    fi

    counts="$results/origin_analysis/origin_analysis.align.$align_mode.tab"
    sampleid="$sample.align-$align_mode"
    Rscript $RSCRIPTS/origin.plots.R $counts "$sampleid" $subplotdir

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

    for ig_idx in ${!ignore_modes[@]}; do

        ignore_mode=${ignore_modes[$ig_idx]}

        # --- nuc_freq plots ------------------------------------
        subplotdir="$plotdir/nuc_freqs"
        if [[ ! -d $subplotdir ]]; then
            mkdir -p $subplotdir
        fi

        for offset_max in ${offset_maxs[@]}; do
            
            counts="$results/nuc_freqs/$sample.align.$align_mode.ignore.$ignore_mode.nuc_freqs.tab"
            sampleid="$sample.align-$align_mode.subset-$ignore_mode"

            Rscript $RSCRIPTS/nuc.freqs.R \
                -n "$sampleid" \
                -d $subplotdir \
                --offsetmax $offset_max \
                $counts
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

