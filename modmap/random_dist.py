#! /usr/bin/env python

''' modmap.random_dist '''

import ipdb
import sys

from collections import Counter, defaultdict

from toolshed import reader
from pybedtools import BedTool

from rpy2.robjects.packages import importr
stats = importr('stats')
# poisson distribution 
rpois = stats.rpois

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

# Copyright 2014 Jay R. Hesselberth

def random_dist(bedgraph_filename, chrom_size_filename, interval_size,
                region_type, only_chroms, ignore_chroms, verbose):

    print >>sys.stderr, ">> only chroms: %s" % str(only_chroms)
    print >>sys.stderr, ">> ignore chroms: %s" % str(ignore_chroms)
    print >>sys.stderr, ">> interval size: %s" % str(interval_size) 

    bedtool = BedTool(bedgraph_filename)

    counts = interval_counts(bedtool, interval_size, chrom_size_filename,
                             only_chroms, ignore_chroms, verbose)

    genome_size = calc_genome_size(chrom_size_filename,
                                   only_chroms, ignore_chroms, verbose)

    plambda = calc_lambda(bedtool, genome_size,
                          only_chroms, ignore_chroms, verbose)

    rand_counts = random_counts(genome_size, interval_size, plambda,
                                verbose)

    write_table(counts, data_type='obs', region_type=region_type,
                interval_size=interval_size, verbose=verbose)
    write_table(rand_counts, data_type='rand', region_type=region_type,
                interval_size=interval_size, verbose=verbose)

def write_table(counts, data_type, region_type, interval_size, verbose):

    header = ('#obs.counts','num.intervals','data.type',
              'interval.size', 'region.type')
    print '\t'.join(map(str, header))

    num_interval = Counter() # number of intervals with this count

    for interval_num in counts:
        for site_count in counts[interval_num].keys():
            num_interval[site_count] += 1

    for obs_counts, num_intervals in sorted(num_interval.items()):

        fields = (obs_counts, num_intervals, data_type, interval_size, 
                  region_type)
        print '\t'.join(map(str, fields))

def interval_counts(bedtool, interval_size, chrom_size_filename,
                    only_chroms, ignore_chroms, verbose):

    result = defaultdict()

    # make windows for analysis
    windows = BedTool().window_maker(w=interval_size,
                                     g=chrom_size_filename).sort()

    # collapse per inteval (comma delim counts, or '0')
    mapresult = windows.map(bedtool, o='collapse', c=4, null=0)

    for idx, row in enumerate(mapresult):

        if (only_chroms and row.chrom not in only_chroms) or \
           (ignore_chroms and row.chrom in ignore_chroms):
            continue

        print row
        nums = [int(i) for i in row.name.split(',')]
        counts = Counter(nums)

        # find number of non-zero counts
        total_counts = sum([i for i in counts.values() if i > 0])

        total_size = int(row.end - row.start)
        num_zeros = total_size - total_counts

        # change the 0 counts to the calculated number 
        counts[0] = num_zeros

        result[idx] = counts

    return result

def random_counts(genome_size, interval_size, plambda, verbose):

    result = defaultdict()

    rand_obs = rpois(genome_size, plambda)

    interval_count = 1
    for idx in range(0, int(genome_size), interval_size):

        obs = rand_obs[idx:idx+interval_size]

        result[interval_count] = Counter(obs)

        interval_count += 1

    return result

def calc_lambda(bedtool, genome_size, only_chroms,
                ignore_chroms, verbose):
    ''' generate counts from poisson from observed lambda and genome size
    '''
    # find lambda
    total_counts = 0.0
    for row in bedtool:
        if (only_chroms and row.chrom not in only_chroms) or \
           (ignore_chroms and row.chrom in ignore_chroms):
            continue

        total_counts += float(row['name'])

    plambda = total_counts / genome_size

    if verbose:
        print >>sys.stderr, ">> lambda calc: %s" % str(plambda)

    return plambda

def calc_genome_size(chrom_size_filename, only_chroms,
                     ignore_chroms, verbose):

    genome_size = 0.0

    for row in reader(chrom_size_filename, header=['chrom','size']):
        if (only_chroms and row['chrom'] not in only_chroms) or \
           (ignore_chroms and row['chrom'] in ignore_chroms):
            continue

        genome_size += float(row['size'])

    if verbose:
        print >>sys.stderr, ">> genome size: %s" % str(genome_size)

    return genome_size

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

    group.add_option("--region-type", type='str',
        metavar="CHROM", default='generic',
        help="label for region type in output"
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

    ignore_chroms = set(options.ignore_chrom)
    only_chroms = set(options.only_chrom)

    kwargs = {'interval_size':options.interval_size,
              'region_type':options.region_type,
              'ignore_chroms':ignore_chroms,
              'only_chroms':only_chroms,
              'verbose':options.verbose}

    bedgraph_filename = args[0]
    chrom_size_filename = args[1]

    return random_dist(bedgraph_filename, chrom_size_filename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
