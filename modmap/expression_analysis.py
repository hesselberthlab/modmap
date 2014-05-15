#! /usr/bin/env python

''' modmap.expression_analysis: XXX
'''

import sys
import pdb

from operator import itemgetter
from collections import defaultdict

from pybedtools import BedTool

from .common import load_coverage

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

STRANDS = ('pos','neg')

def expression_analysis(bam_filename, mrna_bed_filename,
                        pos_exp_bedgraph, neg_exp_bedgraph,
                        chrom_sizes_filename,
                        verbose):

    pos_signal_bedtool = load_coverage(bam_filename, strand='pos',
                                       verbose=verbose) 
    neg_signal_bedtool = load_coverage(bam_filename, strand='neg',
                                       verbose=verbose) 

    exp_signals = calc_exp_signals(mrna_bed_filename,
                                   pos_signal_bedtool, neg_signal_bedtool,
                                   pos_exp_bedgraph, neg_exp_bedgraph)

def calc_exp_signals(mrna_bed_filename, pos_signal_bedtool,
                     neg_signal_bedtool, pos_exp_bedgraph,
                     neg_signal_bedgraph, verbose):
    pass

def calc_promoter_signals(mrna_bed_filename, pos_signal_bedtool,
                          neg_signal_bedtool, pos_exp_bedgraph,
                          neg_signal_bedgraph, verbose):
    pass

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = ("%prog [OPTION]... BAM_FILENAME MRNA_BED POS_EXP_BEDGRAPH "
            "NEG_EXP_BEDGRAPH CHROM_SIZE")
    version = "%%prog %s" % __version__
    description = ("expression analysis")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if len(args) != 5:
        parser.error("specify 5 required files")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    kwargs = {'verbose':options.verbose}

    bam_filename = args[0]
    mrna_bedfile = args[1]

    return expression_analysis(bam_filename,
                               mrna_bedfile,
                               pos_exp_bedgraph,
                               neg_exp_bedgraph)

if __name__ == '__main__':
    sys.exit(main()) 

