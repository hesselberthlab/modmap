#! /usr/bin/env python

''' modmap.signal_analysis: XXX
'''

import sys
import ipdb

from itertools import izip
from collections import defaultdict

from pybedtools import BedTool

from .common import load_coverage, STRANDS

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

def signal_analysis(bam_filename, region_bed_filename,
                    chrom_sizes_filename,
                    pos_strand_label, neg_strand_label,
                    signal_colnum, region_type, 
                    normalize, verbose):

    # region_bed_filename has FPKM as score
    signals = calc_signals(bam_filename,
                           region_bed_filename,
                           signal_colnum,
                           region_type,
                           normalize,
                           verbose)

    header_fields = ('#region.name','region.score','region.strand',
                     'region.type', 'signal.strand','operation',
                     'signal', 'signal.type')
    print '\t'.join(header_fields)

    for fields in signals:
        print '\t'.join(map(str, fields))

def calc_signals(bam_filename, region_bed_filename, signal_colnum,
                 region_type, normalize, verbose):

    ''' generator to calculate signals from BED regions mapped onto positive and
    negative strand data.'''

    region_bedtool = BedTool(region_bed_filename)

    # bedtools.map operations
    operations = ('sum','count')

    signal_type = 'raw'
    if normalize:
        signal_type = 'norm'

    for signal_strand in STRANDS:

        signal_bedtool = load_coverage(bam_filename, strand=signal_strand,
                                       verbose=verbose)
        for oper in operations:

            map_bedtool = region_bedtool.map(signal_bedtool, o=oper,
                                             c=signal_colnum, null=0)

            for region_row, signal_row in izip(region_bedtool, map_bedtool):
   
                try:
                    region_name = region_row[3]
                    region_score = region_row[4]
                    region_strand = region_row[5]

                except IndexError:
                    region_name = '%s-%s-%d-%d' % (region_type,
                                                   region_row.chrom,
                                                   region_row.start,
                                                   region_row.end)
                    region_score = 0
                    region_strand = 'none'

                if region_strand == '+':
                    region_strand = 'pos'
                elif region_strand == '-':
                    region_strand = 'neg'

                # last field is the calculated signal
                signal = float(signal_row[-1])

                if normalize and signal != 0:
                    region_size = float(region_row.end - region_row.start)
                    signal = signal / region_size

                result = (region_name, region_score, region_strand,
                          region_type, signal_strand, oper, signal, signal_type)

                yield result

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = ("%prog [OPTION]... BAM_FILENAME REGION_BED "
            "CHROM_SIZE")
    version = "%%prog %s" % __version__
    description = ("signal analysis")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("--signal-colnum", action="store", type='str',
        default=4, help="column num for signal (default: %default)")

    group.add_option("--pos-strand-label", action="store", type='str',
        default='pos', help="alternate pos strand label (default: %default)")

    group.add_option("--neg-strand-label", action="store", type='str',
        default='neg', help="alternate neg strand label (default: %default)")

    group.add_option("--region-type", action="store", type='str',
        default='generic', help="region type (default: %default)")

    group.add_option("--normalize", action="store_true",
        default=False, help="normalize signal to region length (default: %default)")

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
              'region_type':options.region_type,
              'normalize':options.normalize,
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

