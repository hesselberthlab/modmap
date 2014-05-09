PROJECTID=storici-lab

SAMPLES=(FS1 FS2 FS3 FS4 FS5 FS6 FS7 FS8 FS9 FS10)
NUM_SAMPLES=${#SAMPLES[@]}
# SAMPLES and DESCRIPS must be in same order
DESCRIPS=("YFP17, rnh1 rnh201 / RE"
          "YFP17, rnh1 rnh201 / Fragmentase"
          "YFP17, rnh1 rnh201 ung1"
          "E134, rnh201"
          "E134, pol3-5DV rnh201"
          "E134, pol2-4 rnh201"
          "E134, rnh201"
          "YFP17, rnh201"
          "YFP17, pol2-M644G rnh201"
          "E134, rnh1 rnh201")

#ASSEMBLY=sacCer2
#BOWTIEIDX=$HOME/ref/genomes/$ASSEMBLY/$ASSEMBLY
#RESULT=$HOME/projects/collab/storici-lab/results/common/$ASSEMBLY
#FASTA=$HOME/ref/genomes/$ASSEMBLY/$ASSEMBLY.fa
#CHROM_SIZES=$HOME/ref/genomes/$ASSEMBLY/$ASSEMBLY.chrom.sizes
# GTF=$HOME/ref/genomes/$ASSEMBLY/sgdGene.$ASSEMBLY.gtf

DATA=$HOME/projects/collab/storici-lab/data/common
BIN=$HOME/devel/modmap/modmap
RSCRIPTS=$HOME/devel/modmap/R
UMI=NNNNNNNN
ALIGN_MODES=("uniq" "all")
ALIGN_ARGS=("-m 1" "--all")
EXPPOS=$DATA/$ASSEMBLY/regev.exp.pos.bg
EXPNEG=$DATA/$ASSEMBLY/regev.exp.neg.bg
METADATA=$DATA/metadata.tsv
