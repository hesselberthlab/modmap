#! /usr/bin/env python

''' modmap.signal_analysis: XXX
'''

import sys
import pdb

from itertools import izip
from collections import defaultdict

from pybedtools import BedTool

from .common import load_coverage, STRANDS

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

def signal_analysis(bam_filename, region_bed_filename, chrom_sizes_filename,
                    pos_strand_label, neg_strand_label, signal_colnum,
                    signal_operation, verbose):

    pos_signal_bedtool = load_coverage(bam_filename, strand='pos',
                                       verbose=verbose) 
    neg_signal_bedtool = load_coverage(bam_filename, strand='neg',
                                       verbose=verbose) 

    # region_bed_filename has FPKM as score
    signals = calc_signals(region_bed_filename,
                           pos_signal_bedtool,
                           neg_signal_bedtool,
                           signal_operation,
                           signal_colnum,
                           verbose)

    header_fields = ('#region.name','region.score','region.strand',
                     'signal.pos','signal.neg')
    print '\t'.join(header_feilds)

    for fields in signals:
        region_name, region_score, region_strand, signal_pos, signal_neg = fields
        print '\t'.join(map(str, fields))

def calc_signals(region_bed_filename, pos_signal_bedtool,
                 neg_signal_bedtool, signal_operation, signal_colnum,
                 verbose):

    ''' generator to calculate signals from BED regions mapped onto positive and
    negative strand data.'''

    region_bedtool = BedTool(region_bed_filename)

    # strand flags are irrelevant, treating pos and neg as separate files
    signal_pos_bedtool = _map_signals(region_bedtool, pos_signal_bedtool,
                                      signal_operation, signal_colnum)

    signal_neg_bedtool = _map_signals(region_bedtool, neg_signal_bedtool,
                                      signal_operation, signal_colnum)

    zipped = izip(region_bedtool, signal_pos_bedtool, signal_neg_bedtool)

    for region_row, signal_pos_row, signal_neg_row in zipped:

        # XXX: need less ugly accessor, maybe using BED and bedGraph specs
        region_name = region_row[3]
        region_score = region_row[4]
        region_strand = region_row[5]

        signal_pos = signal_pos_row[6]
        signal_neg = signal_neg_row[6]

        yield (region_name, region_score, region_strand, signal_pos, signal_neg)

def _map_signals(region_bedtool, signal_bedtool, operation, signal_colnum):
    ''' aux function to improve readability '''
    map_bedtool = region_bedtool.map(signal_bedtool, o=operation,
                                     c=signal_colnum, null=0)
    return map_bedtool

def print_report(signal_result, verbose):
    ''' DOC '''
    pass

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = ("%prog [OPTION]... BAM_FILENAME REGION_BED "
            "CHROM_SIZE")
    version = "%%prog %s" % __version__
    description = ("signal analysis")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("--signal-colnum", action="store", type='int',
        default=4, help="column num for signal (default: %default)")

    group.add_option("--signal-operation", action="store", type='str',
        default='sum', help="operation for bedtool map (default: %default)")

    group.add_option("--pos-strand-label", action="store", type='str',
        default='pos', help="alternate pos strand label (default: %default)")

    group.add_option("--neg-strand-label", action="store", type='str',
        default='neg', help="alternate neg strand label (default: %default)")

    group.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if len(args) != 3:
        parser.error("specify 3 required files")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    kwargs = {'pos_strand_label':options.pos_strand_label,
              'neg_strand_label':options.neg_strand_label,
              'signal_colnum':options.signal_colnum,
              'signal_operation':options.signal_operation,
              'verbose':options.verbose}

    bam_filename = args[0]
    region_bedfilename = args[1]
    chrom_sizes_filename = args[2]

    return signal_analysis(bam_filename,
                           region_bedfilename,
                           chrom_sizes_filename,
                           **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 

