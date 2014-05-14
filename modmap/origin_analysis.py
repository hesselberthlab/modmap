#! /usr/bin/env python

''' modmap.origin_analysis: XXX
'''

import sys
import pdb

from collections import defaultdict

from pybedtools import BedTool

from .nuc_frequencies import calc_nuc_counts

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

SIDES = ('left','right')
STRANDS = ('pos','neg')

def origin_anlaysis(timing_bedgraph, origin_bed, max_timing, flank_size,
                    bam_filename, fasta_filename, chrom_size_filename,
                    verbose):

    selected_ori = select_origins(origin_bed, timing_bedgraph, max_timing,
                                  verbose)

    # XXX duplicative in modmap.nuc_frequencies
    pos_signal_bedtool = BedTool.genomecov(ibam=bam_filename, five=True,
                                           strand='+')
    neg_signal_bedtool = BedTool.genomecov(ibam=bam_filename, five=True,
                                           strand='-')

    ori_signals = calc_origin_signals(selected_ori, pos_signal_bedtool,
                                      neg_signal_bedtool, flank_size,
                                      verbose)

def nuc_counts(ori_signals, verbose):
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

    for ori_strand in strand_args.keys():

        # index 0 = strand; index 1 = side
        pos_strand_arg = strand_args[ori_strand]['pos'][0]
        neg_strand_arg = strand_args[ori_strand]['neg'][0]
        pos_side_arg = strand_args[ori_strand]['pos'][1]
        neg_side_arg = strand_args[ori_strand]['neg'][1]

        pos_signal = ori_signals[pos_strand_arg][pos_side_arg]
        neg_signal = ori_signals[neg_strand_arg][neg_side_arg]

        counts = calc_nuc_counts(pos_signal, neg_signal, 
                                 fasta_filename, **kwargs)

        result[ori_strand] = counts

    return result

def calc_origin_signals(selected_ori, pos_signal_bedtool,
                        neg_signal_bedtool, flank_size, verbose):
    ''' calcualte signals within flanks up and downstream of selected
    origins '''

    # left = upstream; right = downstream
    flanks = {}
    flanks['left'] = selected_ori.flank(l=flank_size, r=0,
                                            g=chrom_sizes)
    flanks['right'] = selected_ori.flank(l=0, r=flank_size,
                                            g=chrom_sizes)

    # get positive and negative signals within flanks
    signals = {}
    siganls['pos'] = signal_pos_bedtool
    signals['neg' ]= signal_neg_bedtool

    result = defaultdict(dict)

    for strand in STRANDS:
        for side in SIDES:

            signal_bt = signals[strand]
            flank_bt = flanks[side]

            result[strand][side] = signal_bt.intersect(flank_bt)

    return result

def select_origins(origin_bed, timing_bedgraph, max_timing, verbose):
    ''' select origins based on max_timing '''

    origin_bt = BedTool(origin_bed)
    timing_bt = BedTool(timing_bedgraph)

    # intersect and group to calculate timing near origins
    int_result = origin_bt.intersect(timing_bt, wb=True)
    group_result = int_result.groupby(g=5, c=10, o=mean)
   
    select_origins = []
    for row in group_result:
        
        # XXX
        pdb.set_trace()

        if timing <= max_timing:
            select_origins.append(origin)
   
    if verbose:
        print >>sys.stderr, ">> %d origins selected" % len(select_origins)

    return select_origins

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... BAM_FILENAME FASTA_FILENAME"
    version = "%%prog %s" % __version__
    description = ("origin analysis")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if len(args) != 2:
        parser.error("specify BAM and FASTA")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    kwargs = {'verbose':options.verbose}

    bam_filename = args[0]
    fasta_filename = args[1]

    return nuc_frequencies(bam_filename, fasta_filename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 

