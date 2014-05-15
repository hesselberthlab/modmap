#! /usr/bin/env bash

#BSUB -J tracklines
#BSUB -o tracklines.%J.out
#BSUB -e tracklines.%J.err

<<DOC
Generate tracklines for UCSC
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/collab/storici-lab/bin/config.sh

# this grabs colors for each sample
num_colors=${#SAMPLES[@]}
color_str=$(python -c "import colorbrewer; print colorbrewer.Paired[$num_colors]" \
    | sed 's/ //g; s/(/"/g;s/)/"/g; s/","/ /g; s/"//g; s/\[//g; s/\]//g')
colors=(`echo ${color_str}`)

urlbase="http://amc-sandbox.ucdenver.edu"
strands=("both" "pos" "neg")
umi_types=("removed" "UMIs_not_removed")

tracklinefile=$RESULT/tracklines.txt

if [[ -f $tracklinefile ]]; then
    rm -f $tracklinefile
fi

for sample_idx in ${!SAMPLES[@]}; do

    sample=${SAMPLES[$sample_idx]}
    webpath="~jhessel/projects/storici/results/common/$ASSEMBLY/$sample/bedgraphs"
    color="color=${colors[$sample_idx]}"
    height="maxHeightPixels=50:35:35"

    # keep these loops in order so that align modes are group together for
    # visual comparison
    for umi_type in ${umi_types[@]}; do
        if [[ $umi_type == 'UMIs_not_removed' ]]; then
            vis="visibility=hide"
        else
            vis="visibility=full"
        fi
        for align_mode in ${ALIGN_MODES[@]}; do
            for strand in ${strands[@]}; do

                if [[ $umi_type == 'removed' ]]; then
                    countsbw="$urlbase/$webpath/$sample.align.$align_mode.strand.$strand.counts.bw"
                    name="name='$sample $strand $align_mode umi.removed'"
                    description="description='sample=$sample strand=$strand \
                                align.mode=$align_mode umi=removed'"
                else
                    countsbw="$urlbase/$webpath/$sample.$umi_type.align.$align_mode.strand.$strand.counts.bw"
                    name="name='$sample $strand $align_mode umi.not.removed'"
                    description="description='sample=$sample strand=$strand \
                                align.mode=$align_mode umi=not removed'"
                fi

                url="bigDataUrl=$countsbw"
                trackbase="track type=bigWig"
                trackline="$trackbase $url $name $description $color $vis $height"

                # print the line
                echo $trackline >> $tracklinefile

            done
        done
    done
done
