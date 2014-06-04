# this sets the RESULT directory to "common-debug"
DEBUG="-debug"

# 10 yeast samples
SAMPLES=(6_4.rep1 Afu_BG234_predig CPD.rep1
         UDG_predig ung1_predig
         6_4.rep2  Afu_Y7092_predig  CPD.rep2
         UDG_postdig  ung1_5fu_predig)
NUM_SAMPLES=${#SAMPLES[@]}

# XXX update this with project specific location
DATA=$HOME/projects/collab/storici-lab/data/common

PROJECT=$HOME/devel/lab-projects/excision-seq-ms
FASTQDIR=$PROJECT/ncbi-geo/fastq
LOG=$PROJECT/log

# from modmap pipeline
BIN=$HOME/devel/modmap
RSCRIPTS=$HOME/devel/modmap/R

ALIGN_MODES=("uniq" "all")
ALIGN_ARGS=("-m 1" "--all")
