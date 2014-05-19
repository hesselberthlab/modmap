#! /usr/bin/env python

''' modmap.signal_analysis: XXX
'''

import sys
import pdb

from operator import itemgetter
from itertool import izip
from collections import defaultdict

from pybedtools import BedTool

from .common import load_coverage, STRANDS

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

def expression_analysis(bam_filename, region_bed_filename,
                        chrom_sizes_filename,
                        pos_strand_label, neg_strand_label,
                        verbose):

    pos_signal_bedtool = load_coverage(bam_filename, strand='pos',
                                       verbose=verbose) 
    neg_signal_bedtool = load_coverage(bam_filename, strand='neg',
                                       verbose=verbose) 

    # XXX: make this a var
    operation = 'sum'
    # region_bed_filename has FPKM as score
    signals = calc_signals(region_bed_filename,
                           pos_signal_bedtool,
                           neg_signal_bedtool,
                           operation,
                           verbose)
                                   
def calc_signals(region_bed_filename, pos_signal_bedtool,
                 neg_signal_bedtool, operation, verbose):
    ''' DOC '''
    result = defaultdict(dict)

    # signal is bedGraph format
    signal_colnum = 4

    region_bedtool = BedTool(region_bed_filename)

    # strand flags are irrelevant, treating pos and neg as separate files
    signal_pos_bedtool = region_bedtool.map(pos_signal_bedtool,
                                            o=operation, c=signal_colnum)

    signal_neg_bedtool = region_bedtool.map(neg_signal_bedtool,
                                            o=operation, c=signal_colnum)

    zipped = izip(region_bedtool, signal_pos_bedtool, signal_neg_bedtool)

    for row in zipped:
        pdb.set_trace()

    return result

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

