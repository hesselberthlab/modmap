#! /usr/bin/env python

''' modmap.random_dist '''

import ipdb
import sys

from collections import Counter, defaultdict

from pybedtools import BedTool

from .common import load_coverage

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

# Copyright 2014 Jay R. Hesselberth

def random_dist(bedgraph_filename, chrom_size_filename, interval_size,
                only_chroms, ignore_chroms, verbose):

    bedtool = BedTool(bedgraph_filename)

    result = interval_counts(bedtool, interval_size, chrom_size_filename,
                             only_chroms, ignore_chroms, verbose)

    ipdb.set_trace()

def interval_counts(bedtool, interval_size, chrom_size_filename,
                    only_chroms, ignore_chroms, verbose):

    result = defaultdict()

    # make windows for analysis
    windows = BedTool().window_maker(w=interval_size,
                                     g=chrom_size_filename)

    # collapse per inteval (comma delim counts, or '0')
    mapresult = windows.map(bedtool, o='collapse', c=4, null=0)

    for idx, row in enumerate(mapresult):

        if (only_chroms and row.chrom not in only_chroms) or \
           (ignore_chroms and row.chrom in ignore_chroms):
            continue

        nums = row.name.split(',')
        counts = Counter(nums)

        # find number of non-zero counts
        total_counts = sum([i for i in counts.values() if i > 0])

        total_size = int(row.end - row.start)
        num_zeros = total_size - total_counts

        counts[0] = num_zeros

        result[idx] = counts

    return result

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... BEDGRAPH_FILENAME CHROM_SIZE_FILENAME"
    version = "%%prog %s" % __version__
    description = ("")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("--interval-size", type='int',
        metavar="SIZE", default=10000,
        help="size of interval for counts"
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
        parser.error("specify BEDGRAPH and CHROM_SIZES")

    if len(options.ignore_chrom) != 0 and len(options.only_chrom) != 0:
        parser.error("--ignore-chrom and --only-chrom are mutually exclusive ")

    return options, args

def main(args=sys.argv[1:]):

    options, args = parse_options(args)

    ignore_chroms = tuple(options.ignore_chrom)
    only_chroms = tuple(options.only_chrom)

    kwargs = {'interval_size':options.interval_size,
              'ignore_chroms':ignore_chroms,
              'only_chroms':only_chroms,
              'verbose':options.verbose}

    bedgraph_filename = args[0]
    chrom_size_filename = args[1]

    return random_dist(bedgraph_filename, chrom_size_filename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
