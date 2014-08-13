#!/usr/bin/env bash

#BSUB -J geo.bundle[1-13]
#BSUB -e geo.bundle.%J.%I.err
#BSUB -o geo.bundle.%J.%I.out
#BSUB -P storici

source $CONFIG

bundle_dir=$RESULT/geo_submission
if [[ ! -d $bundle_dir ]]; then
    mkdir -p $bundle_dir
    mkdir -p $bundle_dir/fastq
    mkdir -p $bundle_dir/bedgraphs
fi

bundle_file=$bundle_dir/storici-geo-submission.tar.gz

# running list of tarfiles
for sample in ${SAMPLES[@]}; do
    fastq=$DATA/$sample.fq.gz
    posbg=$bgresults/$sample.align.all.strand.pos.counts.bg
    negbg=$bgresults/$sample.align.all.strand.neg.counts.bg

    cp $fastq $bundle_dir/fastq
    cp $posbg $bundle_dir/bedgraphs
    cp $negbg $bundle_dir/bedgraphs
done

# gzip bedgraphs
for bgfile in $bundle_dir/bedgraphs; do
    gzip $bgfile
done

tar zvcf $bundle_file $bundle_dir
