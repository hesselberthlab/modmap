PROJECTID=stivers-lab

SAMPLES=(JS1 JS2 JS3 JS4 JS5 JS6 JS7 JS8
         JS9 JS10 JS11 JS12 JS13 JS14 JS15 JS16 
         JS17 JS18 JS19 JS20 JS21 JS22
         JS23 JS24
         JS25 JS26 JS27 JS28
         JS29 JS30 JS31 JS32
         JS33 JS34 JS35 JS36 JS37 JS38 JS39 JS40
         JS41 JS42 JS43 JS44 JS45 JS46 JS47 JS48
         JS49 JS50 JS51 JS52 JS53 JS54 JS55 JS56
         JS57 JS58 JS59 JS60 JS61 JS62)

PROJECT="$HOME/projects/collab/stivers-lab"
DOC=$PROJECT/doc
RESULT=$PROJECT/results/common
DATA=$PROJECT/data/common

ALIGN_MODES=("uniq" "all")
ALIGN_ARGS=("-m 1" "--all")

BOWTIEIDXS=($HOME/ref/genomes/hg19/hg19
            $DATA/index/NL4-3-integrant
            $DATA/index/NL4-3_replication_competent_no_GFP
            $DATA/index/UGI)
BOWTIEIDX_TYPES=("hg19" "virus" "plasmid" "UGI")

CHROM_SIZES=$HOME/ref/genomes/hg19/hg19.chrom.sizes
CHROM_SIZES_VIRUS=$DATA/index/NL4-3-integrant.chrom.sizes
CHROM_SIZES_PLASMID=$DATA/index/NL4-3_replication_competent_no_GFP.chrom.sizes
CHROM_SIZES_UGI=$DATA/index/UGI.chrom.sizes

METADATA=$DOC/metadata.tsv
