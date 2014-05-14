#! /usr/bin/env python

''' modmap.nuc_frequencies: calculate nucleotide frequencies at modified base
sites.
'''

import sys
import pdb

from operator import itemgetter
from itertools import izip
from collections import Counter, defaultdict

from pybedtools import BedTool
from pyfasta import Fasta, complement

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

# Copyright 2013,2014 Jay R. Hesselberth

def nuc_frequencies(bam_filename, fasta_filename, 
                    revcomp_strand, min_counts, 
                    offset_min, offset_max, region_size,
                    ignore_chroms, only_chroms, verbose):

    pos_signal_bedtool = BedTool.genomecov(ibam=bam_filename, five=True,
                                           strand='+')
    neg_signal_bedtool = BedTool.genomecov(ibam=bam_filename, five=True,
                                           strand='-')

    nuc_counts = calc_nuc_counts(pos_signal_bedtool, neg_signal_bedtool,
                                 seq_fasta, 
                                 revcomp_strand, min_counts, 
                                 offset_min, offset_max, region_size,
                                 ignore_chroms, only_chroms, verbose)

    print_report(nuc_counts, verbose)

def calc_nuc_counts(pos_signal_bedtool, neg_signal_bedtool, fasta_filename, 
                    revcomp_strand, min_counts, 
                    offset_min, offset_max, region_size,
                    ignore_chroms, only_chroms, verbose):

    ''' main routine for calculating nuc_counts) '''

    if verbose:
        msg =  ">> analyzing sequences ...\n"
        msg += ">> ignore:%s only:%s\n" % \
            (str(ignore_chroms), str(only_chroms))
        msg += ">> offset range: %d to %d\n" % (offset_min, offset_max)
        msg += ">> region size: %d\n" % (region_size)
        msg += ">> revcomp strand: %s\n" % str(revcomp_strand)
        print >>sys.stderr, msg

    seq_fasta = Fasta(fasta_filename)

    nuc_counts = defaultdict(Counter)

    bedtools = (pos_signal_bedtool, neg_signal_bedtool)
    strands = ('+', '-')

    # total number of sites examined
    total_sites = 0

    for bedtool, strand in izip(bgfiles, strands):

        for row in bedtool:

            # skip data based on specified chromosomes
            if row.chrom in ignore_chroms:
                continue
            if only_chroms and row.chrom not in only_chroms:
                continue

            # skip data if counts are too low
            if row.count < min_counts: continue

            # sites in bedgraph examined - must come after all checks
            # above
            total_sites += 1

            for offset in range(offset_min, offset_max + 1):

                # upstream offsets are negative values
                if strand == '+':
                    start = row.start + offset
                elif strand == '-':
                    start = row.start - offset

                if region_size == 1:
                    # half open at the position of interest
                    end = start + region_size
                else:
                    # make sure that the 3' most position in a region
                    # is the base of interest
                    if strand == '+':
                        end = start + 1 # include position with + 1
                        start = end - region_size
                    else:
                        # negative strand
                        end = start + region_size

                # XXX: does this ever happen?
                if start < 0: continue

                nucs = seq_fasta[row.chrom][start:end]

                #  1. libs where the captured strand is sequenced
                #     are the correct polarity as-is (i.e. Excision-seq
                #     libs)
                #  2. libs where the *copy* of the captured strand
                #     is sequenced should be revcomplemented (i.e.
                #     circularization-based libs)

                if (strand == '+' and revcomp_strand) or \
                   (strand == '-' and not revcomp_strand):
                    nucs = complement(nucs[::-1])

                if len(nucs.strip()) != region_size: continue

                nuc_counts[offset][nucs] += row.count

    return nuc_counts

def print_report(nuc_counts, verbose):
    ''' print report of nuc_counts with frequency calculations '''

    if verbose:
        debug_report(nuc_counts)

    # report the results
    header = ('#nuc','offset','region.size','count','freq','total.sites') 
    print '\t'.join(header)

    for offset, counts in sorted(nuc_counts.items()):
        sum_counts = sum(counts.values())

        for nuc, count in sorted(counts.items(), key=itemgetter(1), \
                                 reverse=True):
            freq = float(count) / float(sum_counts)
            vals = map(str, [nuc, offset, region_size, count, freq,
                             total_sites])
            print '\t'.join(vals)

def debug_report(nuc_counts):
    ''' debug counts at specific positions. currently calculates the sum
    of frequencies a the 3' most position '''

    for offset, counts in sorted(nuc_counts.items()):

        sum_counts = sum(counts.values())
        sums = Counter()

        for nuc, count in sorted(counts.items(), key=itemgetter(1), \
                                 reverse=True):
            
            freq = float(count) / float(sum_counts)
            sums[nuc[-1]] += freq

        for nuc, nucsum in sums.items():
            vals = map(str, ['#>>debug', nuc, offset, nucsum])
            print >>sys.stderr, '\t'.join(vals)

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... BAM_FILENAME FASTA_FILENAME"
    version = "%%prog %s" % __version__
    description = ("calculate nucleotide frequences at modified base "
                   "sites")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("--minimum-counts", action="store", type='int',
        metavar="MINIMUM_COUNTS", default=0,
        help="minimum number of counts to be included in analysis"
        " (default: %default)")

    group.add_option("--revcomp-strand", action="store_true", 
        metavar="CAPTURED STRAND", default=False,
        help="read is reverse comp of captured strand"
        " (default: %default)")

    group.add_option("--offset-min", action="store", type='int', 
        metavar="MIN_OFFSET", default=-15,
        help="min offset from current position to examine. 5' positions are "
        " more negative. (default: %default)")

    group.add_option("--offset-max", action="store", type='int', 
        metavar="MAX_OFFSET", default=15,
        help="max offset from current position to examine. 5' positions are "
        " more negative. (default: %default)")

    group.add_option("--region-size", action="store", type='int', 
        metavar="REGION_SIZE", default=1,
        help="size of region to examine"
        " (default: %default)")

    group.add_option("--ignore-chrom", action="append",
        metavar="CHROM", default=[],
        help="list of chroms to ignore"
        " (default: %default)")

    group.add_option("--only-chrom", action="append",
        metavar="CHROM", default=[],
        help="list of chroms to include"
        " (default: %default)")

    group.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if len(args) != 2:
        parser.error("specify BAM and FASTA")

    if options.region_size <= 0:
        parser.error("region_size must be > 0")
       
    if len(options.ignore_chrom) != 0 and len(options.only_chrom) != 0:
        parser.error("--ignore-chrom and --only-chrom are mutually exclusive ")

    return options, args

def main(args=sys.argv[1:]):

    options, args = parse_options(args)

    ignore_chroms = tuple(options.ignore_chrom)
    only_chroms = tuple(options.only_chrom)

    kwargs = {'revcomp_strand':options.revcomp_strand,
              'min_counts':options.minimum_counts,
              'offset_min':options.offset_min,
              'offset_max':options.offset_max,
              'region_size':options.region_size,
              'verbose':options.verbose,
              'ignore_chroms':ignore_chroms,
              'only_chroms':only_chroms}

    bam_filename = args[0]
    fasta_filename = args[1]

    return nuc_frequencies(bam_filename, fasta_filename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
