#! /usr/bin/env python

''' modmap.origin_analysis: XXX
'''

import sys
import pdb

from operator import itemgetter
from collections import defaultdict

from pybedtools import BedTool

from .common import load_coverage, STRANDS
from .nuc_frequencies import calc_nuc_counts

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

# Copyright 2013,2014 Jay R. Hesselberth

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
                
                              __________                 
                             /<---5'    \
             left           /            \          right
    5' --------------------/              \--------------------- 3'

    3' --------------------\              /--------------------- 5'
                            \            /
                             \_____5'-->/

    Leading signals are the captured negative strand (i.e. pos
    strand sequences) upstream (left) ofthe origin and captured 
    positive strand (i.e. neg strand sequences) downstream
    (right) of origin.
    '''

    # args for modmap.nuc_frequencies(). must be total of 8 args
    kwargs = {'offset_min':-15,
              'offset_max':15,
              'region_size':1,
              'revcomp_strand':True,
              'min_counts':1,
              'ignore_chroms':[],
              'only_chroms':[],
              'verbose':verbose}

    # keys are ori_strand, values are (sites, nuc_counts)
    result = defaultdict(dict)

    ori_strands = ('leading', 'lagging')

    # XXX refer to drawing on board for these updated settings
    # `leading` and `lagging` represent the **template** strand
    for ori_strand in ori_strands:

        if ori_strand == 'lagging':
            pos_signal_bedtool = ori_signals['pos']['right']
            neg_signal_bedtool = ori_signals['neg']['left']        
       
        elif ori_strand == 'leading':
            pos_signal_bedtool = ori_signals['pos']['left']
            neg_signal_bedtool = ori_signals['neg']['right']        
        
        if verbose:
            print >>sys.stderr, ">> nuc.freq on %s strand ..." % ori_strand

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
    origins.
    
    result is a dict where keys are strands, values are dicts where keys
    are sides and values are bedtools containg intercted bedGraph data,
    e.g.:
        
        result[pos][left] = bedtool

    `pos` and `neg` refer to the raw coverage data, must be reverse comped
    later
    '''

    result = defaultdict(dict)

    origin_sides = ('left','right')

    for strand in STRANDS:

        if strand == 'pos':
            signal_bedtool = pos_signal_bedtool

        elif strand == 'neg':
            signal_bedtool = neg_signal_bedtool

        for side in origin_sides:

            # left = upstream; right = downstream
            if side == 'left':
                flank_bedtool = origin_bedtool.flank(l=flank_size, r=0,
                                                     g=chrom_sizes_filename)
            elif side == 'right':
                flank_bedtool = origin_bedtool.flank(l=0, r=flank_size,
                                                     g=chrom_sizes_filename)

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
    with open(group_bedtool.TEMPFILES[-1]) as group_data:

        select_origins = set() 
        for row in group_data:
            origin, timing = row.strip().split('\t')
            timing = float(timing)

            if timing <= max_timing:
                select_origins.add(origin)

    if verbose:
        print >>sys.stderr, ">> %d origins selected:" % len(select_origins)

        # these are printed in the data table as comments
        print "# origins selected (%d total):" % len(select_origins)
        for oriname in select_origins:
            print "# %s" % oriname

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

