#! /usr/bin/env python

''' modmap.origin_analysis: XXX
'''

import sys
import pdb

from operator import itemgetter
from collections import defaultdict

from pybedtools import BedTool

from .common import load_coverage
from .nuc_frequencies import calc_nuc_counts

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

SIDES = ('left','right')
STRANDS = ('pos','neg')

def origin_analysis(origin_bed, timing_bedgraph, bam_filename,
                    fasta_filename, chrom_sizes_filename,
                    max_timing, flank_size, verbose):

    origin_bedtool = select_origins(origin_bed, timing_bedgraph, max_timing,
                                  verbose)

    pos_signal_bedtool = load_coverage(bam_filename, strand='pos',
                                       verbose=verbose) 
    neg_signal_bedtool = load_coverage(bam_filename, strand='neg',
                                       verbose=verbose) 

    ori_signals = calc_origin_signals(origin_bedtool,
                                      pos_signal_bedtool,
                                      neg_signal_bedtool,
                                      chrom_sizes_filename,
                                      flank_size, verbose)

    origin_result = calc_origin_nuc_counts(ori_signals,
                                           pos_signal_bedtool,
                                           neg_signal_bedtool,
                                           fasta_filename,
                                           verbose)

    print_report(origin_result, max_timing, flank_size, verbose)

def print_report(origin_result, max_timing, flank_size, verbose):
    ''' print report of nuc_counts with frequency calculations '''

    header = ('#nuc','offset','count','freq','total.sites',
              'max.timing','flank.size','strand',) 
    print '\t'.join(header)

    for strand, tup in origin_result.items():
        total_sites, nuc_counts = tup

        for offset, counts in sorted(nuc_counts.items()):
            sum_counts = sum(counts.values())

            for nuc, count in sorted(counts.items(), key=itemgetter(1), \
                                     reverse=True):

                freq = float(count) / float(sum_counts)

                vals = map(str, [nuc, offset, count, freq, total_sites,
                                 max_timing, flank_size, strand])

                print '\t'.join(vals)

def calc_origin_nuc_counts(ori_signals, pos_signal_bedtool,
                           neg_signal_bedtool, fasta_filename, verbose):

    ''' calculate nucleotide frequencies within origins.

    leading strand = -p $neg_right_bedgraph -n $pos_left_bedgraph
    lagging strand = -p $neg_left_bedgraph -n $pos_right_bedgraph
    '''

    # args for nuc_frequencies()
    kwargs = {'offset_min':-1,
              'offset_max':1,
              'region_size':1,
              'revcomp_strand':True,
              'min_counts':1,
              'ignore_chroms':[],
              'only_chroms':[],
              'verbose':verbose}

    result = defaultdict(dict)

    # keys of the dicts refer to the pos and neg signals
    # for calc_nuc_counts()

    # XXX these are static to avoid ambiguous assigment later
    strand_args = {}
    strand_args['leading'] = {'pos':('neg','right'),
                              'neg':('pos','left')}
    strand_args['lagging'] = {'pos':('neg','left'),
                              'neg':('pos','right')}

    if verbose:
        print >>sys.stderr, ">> calculating nuc frequencies ..."

    for ori_strand in strand_args.keys():

        # index 0 = strand; index 1 = side
        pos_strand_arg = strand_args[ori_strand]['pos'][0]
        neg_strand_arg = strand_args[ori_strand]['neg'][0]
        pos_side_arg = strand_args[ori_strand]['pos'][1]
        neg_side_arg = strand_args[ori_strand]['neg'][1]

        pos_signal_bedtool = ori_signals[pos_strand_arg][pos_side_arg]
        neg_signal_bedtool = ori_signals[neg_strand_arg][neg_side_arg]

        total_sites, nuc_counts = calc_nuc_counts(pos_signal_bedtool,
                                 neg_signal_bedtool, 
                                 fasta_filename,
                                 **kwargs)

        result[ori_strand] = (total_sites, nuc_counts)

    return result

def calc_origin_signals(origin_bedtool, pos_signal_bedtool,
                        neg_signal_bedtool, chrom_sizes_filename,
                        flank_size, verbose):

    ''' calcualte signals within flanks up and downstream of selected
    origins '''

    # left = upstream; right = downstream
    flanks = {}
    flanks['left'] = origin_bedtool.flank(l=flank_size, r=0,
                                            g=chrom_sizes_filename)
    flanks['right'] = origin_bedtool.flank(l=0, r=flank_size,
                                            g=chrom_sizes_filename)

    # get positive and negative signals within flanks
    signals = {}
    signals['pos'] = pos_signal_bedtool
    signals['neg' ]= neg_signal_bedtool

    result = defaultdict(dict)

    for strand in STRANDS:
        for side in SIDES:

            signal_bedtool = signals[strand]
            flank_bedtool = flanks[side]

            result[strand][side] = signal_bedtool.intersect(flank_bedtool)

    return result

def select_origins(origin_bed, timing_bedgraph, max_timing, verbose):
    ''' select origins based on max_timing.
    
        returns: BedTool
        '''

    origin_bedtool = BedTool(origin_bed)
    timing_bedtool = BedTool(timing_bedgraph)

    # intersect and group to calculate timing near origins
    int_bedtool = origin_bedtool.intersect(timing_bedtool, wb=True)
    group_bedtool = int_bedtool.groupby(g=4, c=10, o='mean')
  
    # XXX hack to get a filehandle for looping. have to do this because
    # the groupby result is only 2 columns and doesn't appear to allow
    # iteration by the std BedTool iter() method
    group_data = open(group_bedtool.TEMPFILES[-1])

    select_origins = set() 
    for row in group_data:
        origin, timing = row.strip().split('\t')
        timing = float(timing)

        if timing <= max_timing:
            select_origins.add(origin)
   
    if verbose:
        print >>sys.stderr, ">> %d origins selected" % len(select_origins)

    # now reconstruct a BedTool with selected origins
    select_intervals = []
    for row in origin_bedtool:
        if row.name in select_origins:
            select_intervals.append(row)

    result = BedTool(select_intervals)

    return result 

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... ORIGIN_BED TIMING_BEDGRAPH BAM_FILENAME FASTA_FILENAME CHROM_SIZE"
    version = "%%prog %s" % __version__
    description = ("origin analysis")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("--max-timing", action="store", type='float',
        default=20, help="maximum timing in Trep (default: %default)")

    group.add_option("--flank-size", action="store", type='int',
        default=500, help="origin flank size in bp (default: %default)")

    group.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if len(args) != 5:
        parser.error("specify 5 required files")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    kwargs = {'max_timing':options.max_timing,
              'flank_size':options.flank_size,
              'verbose':options.verbose}

    origin_bed = args[0]
    timing_bedgraph = args[1]
    bam_filename = args[2]
    fasta_filename = args[3]
    chrom_sizes_filename = args[4]

    return origin_analysis(origin_bed, timing_bedgraph,
                           bam_filename, fasta_filename,
                           chrom_sizes_filename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 

