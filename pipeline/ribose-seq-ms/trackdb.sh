#! /usr/bin/env bash

#BSUB -J trackdb
#BSUB -o trackdb.%J.out
#BSUB -e trackdb.%J.err

<<DOC
Generate tracklines for UCSC
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG

urlbase="http://amc-sandbox.ucdenver.edu"

trackdbroot=$RESULT/tracklines.txt
webroot="$urlbase/~jhessel/projects/storici/results/common$DEBUG/hub"

hubfile="$webroot/hub.txt"
echo "hub Storici Data" >> $hubfile
echo "shortLabel Ribo mapping" >> $hubfile
echo "longLabel Storici lab Ribonucleotide mapping" >> $hubfile
echo "genomesFile genomes.txt" >> $hubfile
echo "email jay.hesselberth@gmail.com" >> $hubfile

genomesfile="$webroot/genomes.txt"
for assembly in $ASSEMBLIES; do
    echo "genome $assembly" >> $genomesfile
    echo "trackdb $assembly/trackDb.txt"
done

strands=("both" "pos" "neg")
for assembly in $ASSEMBLIES; do

    for sample_idx in ${!SAMPLES[@]}; do

        # provided by config.sh
        sample=${SAMPLES[$sample_idx]}
        color=${COLORS[$sample_idx]}
        descrip=${DESCRIPS[$sample_idx]}

        webroot="$urlbase/~jhessel/projects/storici/results/common$DEBUG/$ASSEMBLY/$sample"

        bigwigdir="$webroot/bedgraphs"
        bigbeddir="$webroot/peaks"

        for align_mode in ${ALIGN_MODES[@]}; do
            for strand in ${strands[@]}; do

                track_color="color=$color"

                #bigbed tracks
                bigbed="$bigbeddir/${sample}.${strand}.align.${align_mode}_peaks.bb"
                bigbed_name="name='${sample} ${strand} ${align_mode} peaks"
                bigbed_descrip="description='$descrip PEAKS sample=$sample \
                                strand=$strand \
                                align.mode=$align_mode '"

                bigbed_track="track type=bigBed"
                bigbed_url="bigDataUrl=$bigbed"
                bigbed_vis="visibility=squish"
                bigbed_trackline="$bigbed_track $bigbed_url $bigbed_name \
                                  $bigbed_descrip $track_color $bigbed_vis"

                # print the line
                echo $bigbed_trackline >> $tracklinefile

                # bigWig tracks
                bigwig="$bigwigdir/${sample}.align.${align_mode}.strand.${strand}.counts.bw"
                bigwig_name="name='$sample $strand $align_mode'"
                bigwig_descrip="description='$descrip sample=$sample strand=$strand \
                                align.mode=$align_mode'"
                bigwig_height="maxHeightPixels=50:35:35"

                bigwig_track="track type=bigWig"
                bigwig_url="bigDataUrl=$bigwig"
                bigwig_vis="visibility=full"
                bigwig_trackline="$bigwig_track $bigwig_url $bigwig_name \
                                  $bigwig_descrip $track_color \
                                  $bigwig_vis $bigwig_height"

                # print the line
                echo $bigwig_trackline >> $tracklinefile

            done
        done
    done
done
