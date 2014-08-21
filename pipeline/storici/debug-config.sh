# config file for debuggin purposes.
# source this from any script outside of the master to run on a single
# example (sacCer2, first sample (FS1))
CONFIG=$HOME/devel/modmap/pipeline/storici/config.sh
ASSEMBLY=sacCer2
DEBUG="-debug"
RESULT=$HOME/projects/collab/storici-lab/results/common$DEBUG/$ASSEMBLY
CHROM_SIZES=$HOME/ref/genomes/$ASSEMBLY/$ASSEMBLY.chrom.sizes
LSB_JOBINDEX=1

