PROJECTID=storici-lab

# this sets the RESULT directory to "common-debug"
DEBUG="-debug"

# XXX SAMPLES and DESCRIPS and COLORS must be in same order

SAMPLES=(FS1 FS2 FS3 FS4 FS5 FS6 FS7 FS8 FS9 FS10 FS14 FS15 FS16)
NUM_SAMPLES=${#SAMPLES[@]}
DESCRIPS=("YFP17, rnh1 rnh201 / RE"
          "YFP17, rnh1 rnh201 / Fragmentase"
          "YFP17, rnh1 rnh201 ung1"
          "E134, rnh201, rep1"
          "E134, pol3-5DV rnh201"
          "E134, pol2-4 rnh201"
          "E134, rnh201, rep2"
          "YFP17, rnh201"
          "YFP17, pol2-M644G rnh201"
          "E134, rnh1 rnh201"
          "YFP17, rnh1 rnh201",
          "E134, rnh201",
          "YFP17, rnh1 rnh201 / Taq")

COLORS=("228,26,28",
        "228,26,28",
        "55,126,184",
        "77,175,74",
        "152,78,163",
        "152,78,163",
        "77,175,74",
        "166,86,40"
        "152,78,163",
        "255,127,0",
        "228,26,28",
        "228,26,28",
        "228,26,28"
        )

DATA=$HOME/projects/collab/storici-lab/data/common
BIN=$HOME/devel/modmap
RSCRIPTS=$HOME/devel/modmap/R
UMI=NNNNNNNN
ALIGN_MODES=("uniq" "all")
ALIGN_ARGS=("-m 1" "--all")
