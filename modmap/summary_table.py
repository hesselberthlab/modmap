#! /usr/bin/env python

''' modmap.summary_table: XXX
'''

import sys
import ipdb

from itertools import izip
from collections import Counter, defaultdict

from pybedtools import BedTool, Interval
from toolshed import reader
from pyfaidx import Fasta

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

# Copyright 2014 Jay R. Hesselberth

def summary_table(sample_names, bedgraph_filenames, fasta_filename, revcomp, verbose):

    signal_bedtools = generate_bedtools(bedgraph_filenames, verbose)

    intersection = calc_intersection(signal_bedtools, verbose)

    intersection_bedtool = generate_intersection_bedtool(intersection,
                                                         revcomp,
                                                         verbose)

    fasta = generate_fasta(intersection_bedtool, fasta_filename,
                           revcomp, verbose)
 
    maps = generate_maps(intersection_bedtool, signal_bedtools, verbose)

    strand = 'pos'
    if revcomp:
        strand = 'neg'

    headers = ['#chrom','start','end','strand','seq']
    headers.extend(sample_names)
    print '\t'.join(headers)

    for els in izip(intersection_bedtool, fasta, *maps):

        region = els[0]
        seq = str(els[1])

        signals = []
        for interval in els[2:]:
            signals.append(interval.fields[-1])

        fields = [region.chrom, region.start, region.end, strand, seq]
        fields.extend(signals)

        print '\t'.join(map(str, fields))

def generate_fasta(intersection_bedtool, fasta_filename, revcomp, verbose):

    if verbose:
        print >>sys.stderr, ">> generating fasta of positions ..."

    # -s: force strandedness
    fasta_seqs = intersection_bedtool.sequence(fi=fasta_filename, s=True)

    fasta = Fasta(fasta_seqs.seqfn)

    return fasta

def generate_maps(intersection_bedtool, signal_bedtools, verbose):
    
    if verbose:
        print >>sys.stderr, ">> generating signal maps ..."

    map_bedtools = []
    for signal_bedtool in signal_bedtools:

        signal_map = intersection_bedtool.map(signal_bedtool, c=4,
                                              o='sum', null='0')
        map_bedtools.append(signal_map)

    return map_bedtools

def generate_intersection_bedtool(intersection_bedtool, revcomp, verbose):

    intervals = []
    for row in intersection_bedtool:

        chrom, start, stop = row.fields[:3]
        name = score = '.'

        if revcomp:
            strand = '-'
        else:
            strand = '+'

        intervals.append(Interval(chrom, int(start), int(stop), name, score, strand))

    bedtool = BedTool(intervals)

    return bedtool

def generate_bedtools(bedfilenames, verbose):

    result = []
    for bedfilename in bedfilenames:
        bedtool = BedTool(bedfilename)
        result.append(bedtool)

    return result

def calc_intersection(bedtools, verbose):

    intersect_tool = BedTool()

    if verbose:
        print >>sys.stderr, ">> generating intersection ... "

    result = intersect_tool.multi_intersect(i=[bt.fn for bt in bedtools])

    return result

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog --fasta FASTA [OPTION]... -s sample=bedgraph "
    version = "%%prog %s" % __version__
    description = ("")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Required")

    group.add_option("--fasta", action="store", type='str',
        metavar="Fasta file", default=None,
        help="fasta filename"
        " (default: %default)")

    group.add_option("-s", "--signal", action="append", 
        metavar="mapping of sample=bedgraph file", default=[],
        help="-s sample=bedgraph "
        " (default: %default)")

    parser.add_option_group(group)

    group = OptionGroup(parser, "Variables")

    group.add_option("--revcomp", action="store_true", 
        metavar="CAPTURED STRAND", default=False,
        help="read is reverse comp of captured strand"
        " (default: %default)")

    group.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if len(options.signal) < 2:
        parser.error("specify more than 1 sample=BEDGRAPH mapping")
  
    if not options.fasta:
        parser.error("specify FASTA file")

    return options, args

def main(args=sys.argv[1:]):

    options, args = parse_options(args)

    samples = {}
    for opt in options.signal:
        sample_name, bedgraph_filename = opt.split('=')
        samples[sample_name] = bedgraph_filename

    sample_names = samples.keys()
    bedgraph_filenames = samples.values()

    kwargs = {'fasta_filename':options.fasta,
              'revcomp':options.revcomp,
              'verbose':options.verbose}
    
    return summary_table(sample_names, bedgraph_filenames, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
