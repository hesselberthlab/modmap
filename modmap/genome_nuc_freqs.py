#! /usr/bin/env python

''' modmap.genome_nuc_freqs: calculate background nucleotide frequencies
and report a table
'''

import sys
from collections import Counter, defaultdict
from pyfasta import Fasta, complement

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

def genome_nuc_freqs(fasta_filename, region_size_min, region_size_max,
                     verbose):

    header = ('#region.size','nuc','count','freq')
    print '\t'.join(header)

    nuc_counts = calc_nuc_counts(fasta_filename, region_size_min,
                                 region_size_max, verbose)

    for region_size in nuc_counts:

        nuc_sum = float(sum(nuc_counts[region_size].values()))

        for nuc, count in sorted(nuc_counts[region_size].items()):
            nuc_freq = count / nuc_sum
            print '%d\t%s\t%s\t%s' % (region_size, nuc, str(count), str(nuc_freq))

def calc_nuc_counts(fasta_filename, region_size_min,
                    region_size_max, verbose):
    ''' calculate nuc frequencies for normalization.

        Returns: dict of nucleotide frequencies.
    '''

    nuc_counts = defaultdict(Counter)

    fasta = Fasta(fasta_filename)

    for chrom, seq in fasta.items():

        for idx, pos in enumerate(seq):

            for region_size in range(region_size_min,
                                     region_size_max + 1):

                nucs = seq[idx:idx+region_size]

                if len(nucs) < region_size: continue

                nuc_counts[region_size][nucs] += 1

    return nuc_counts

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... FASTA_FILENAME"
    version = "%%prog %s" % __version__
    description = ("pre-calculate nucleotide frequences for genome")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("--region-size-minimum", action="store", type='int',
        default=1,
        help="minimum region size "
        " (default: %default)")

    group.add_option("--region-size-maximum", action="store", type='int',
        default=1,
        help="maximum region size "
        " (default: %default)")

    group.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if len(args) != 1:
        parser.error("specify FASTA")

    return options, args

def main(args=sys.argv[1:]):

    options, args = parse_options(args)

    kwargs = {'region_size_min':options.region_size_minimum,
              'region_size_max':options.region_size_maximum,
              'verbose':options.verbose}

    fasta_filename = args[0]

    return genome_nuc_freqs(fasta_filename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
