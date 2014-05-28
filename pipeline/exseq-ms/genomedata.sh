#! /usr/bin/env bash

#BSUB -J gd.create
#BSUB -e log/gd.create.%J.err
#BSUB -o log/gd.create.%J.out

set -o nounset -o pipefail -o errexit -x
source config.sh

GEOBEDGRAPHS="$PROJECT/ncbi-geo/sacCer1/bedgraphs"
MISCBEDGRAPHS="$PROJECT/stat-model/bedgraphs"

GDARCHIVE="$RESULT/genomedata/exseq.ms.genomedata"

genomedata-load \
    -t cpd.rep1="$GEOBEDGRAPHS/CPD.rep1.bg.gz" \
    -t cpd.rep2="$GEOBEDGRAPHS/CPD.rep2.bg.gz" \
    -t 64.rep1="$GEOBEDGRAPHS/6_4.rep1.bg.gz" \
    -t 64.rep1="$GEOBEDGRAPHS/6_4.rep2.bg.gz" \
    -t udg.afu.bg234.predig="$GEOBEDGRAPHS/Afu_BG234_predig.bg.gz" \
    -t udg.afu.y7092.predig="$GEOBEDGRAPHS/Afu_Y7092_predig.bg.gz" \
    -t udg.postdig="$GEOBEDGRAPHS/UDG_postdig.bg.gz" \
    -t udg.predig="$GEOBEDGRAPHS/UDG_predig.bg.gz" \
    -t udg.ung1.5fu.postdig="$GEOBEDGRAPHS/ung1_5fu_predig.bg.gz" \
    -t udg.ung1.postdig="$GEOBEDGRAPHS/ung1_predig.bg.gz" \
    -t rep.timing.yabuki="$MISCBEDGRAPHS/yabuki.timing.bedgraph.gz" \
    -t rep.timing.raghu="$MISCBEDGRAPHS/raghu.timing.bedgraph.gz" \
    -t dnasei="$MISCBEDGRAPHS/dnasei.bedgraph.gz" \
    -t mrna.levels="$MISCBEDGRAPHS/regev.combined.logscores.bedgraph.gz" \
    -s $FASTA \
    $GDARCHIVE \
    --verbose

