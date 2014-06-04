#! /usr/bin/env bash

#BSUB -J nuc.freqs[1-10]
#BSUB -o nuc.freqs.%J.%I.out
#BSUB -e nuc.freqs.%J.%I.err

set -o nounset -o pipefail -o errexit -x

source $CONFIG

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

if [[ $ASSEMBLY == "sacCer2" ]]; then
    include_modes=("all" "only-mito" "no-mito" "only-2micron")
    include_args=("" "--only-chrom chrM"
                 "--ignore-chrom chrM"
                 "--only-chrom 2micron")
else
    include_modes=("all" "only-mito" "no-mito")
    include_args=("" "--only-chrom chrM"
                 "--ignore-chrom chrM")
fi

# mono, di and trinucleotides
sizes="1 2 3"

# count thresholds
count_thresh="1 10 30"

alignments=$RESULT/$sample/alignments
results=$RESULT/$sample/nuc_freqs

if [[ ! -d $results ]]; then
    mkdir -p $results
fi

# need to be in BIN to run module
cd $BIN

for idx in ${!ALIGN_MODES[@]}; do
    align_mode=${ALIGN_MODES[$idx]}
    align_arg=${ALIGN_ARGS[$idx]}

    BAM=$alignments/$sample.align.$align_mode.bam

    for inc_idx in ${!include_modes[@]}; do

        include_mode=${include_modes[$inc_idx]}
        include_arg=${include_args[$inc_idx]}

        for mincount in $count_thresh; do

            output="$results/$sample.include.$include_mode.mincount.$mincount.nuc_freqs.tab"

            if [[ -f $output ]]; then
                rm -f $output
            fi

            # signals need to be reverse complemented because the sequence is
            # the reverse complement of the captured strand
            for size in $sizes; do
                python -m modmap.nuc_frequencies \
                    $BAM $FASTA \
                    --region-size $size \
                    $include_arg \
                    --background-freq-table $BKGD_FREQS \
                    --verbose \
                    >> $output
            done
        done
    done
done

