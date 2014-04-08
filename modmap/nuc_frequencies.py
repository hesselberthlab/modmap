#! /usr/bin/env python

''' nuc_frequencies: calculate nucleotide frequencies at modified base
sites.

Part of the modmap analysis pipeline. Assumes generation of pos and neg
strand bedgraph counts via e.g.:

$ bedtools genomecov -5 -strand + > pos.bg 
$ bedtools genomecov -5 -strand - > neg.bg
'''

import sys
import pdb

from operator import itemgetter
from itertools import izip
from collections import Counter, defaultdict

from pyfasta import Fasta, complement

from toolshed import reader

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

# Copyright 2013,2014 Jay R. Hesselberth

def nuc_frequencies(posbedgraph, negbedgraph, fastafilename, 
                    revcomp_strand, offset_min, offset_max, region_size,
                    ignore_chroms, only_chroms, verbose):

    if verbose:
        print >>sys.stderr, ">> analyzing sequences ..."
        print >>sys.stderr, ">> ignore:%s only:%s" % \
            (str(ignore_chroms), str(only_chroms))
        print >>sys.stderr, ">> offset range: %d to %d" % \
            (offset_min, offset_max)
        print >>sys.stderr, ">> region size: %d" % \
            (region_size)
        print >>sys.stderr, ">> revcomp strand: %s" % \
            str(revcomp_strand)
     
    seqs = Fasta(fastafilename)

    nuc_counts = defaultdict(Counter)

    bgfiles = (posbedgraph, negbedgraph)
    strands = ('+', '-')

    for bgfile, strand in izip(bgfiles, strands):

        # filenames can be None
        if not bgfile: continue

        for row in reader(bgfile, header=BedGraphLine):

            if row.chrom in ignore_chroms:
                continue
            if only_chroms and row.chrom not in only_chroms:
                continue
          
            for offset in range(offset_min, offset_max + 1):

                # -15 ------------- 0 ------------- +15
                #     TAGCTAGGCTAGCGACGTAGCGACTAGCA
                # (+) -----------------------------
                #                   >>>>>>>>>>>>>>>
                #      <<<<<<<<<<<<<<
                # (-) -----------------------------

                if strand == '+':
                    start = row.start + offset
                elif strand == '-':
                    # this will be revcomped later
                    start = row.start - offset

                end = start + region_size

                if start < 0: continue

                nucs = seqs[row.chrom][start:end]

                #  1. libs where the captured strand is sequenced
                #     are the correct polarity as-is (i.e. Excision-seq
                #     libs)
                #  2. libs where the *copy* of the captured strand
                #     is sequenced should be revcomplemented (i.e.
                #     circularization-based libs)

                if (strand == '+' and revcomp_strand) or \
                   (strand == '-' and not revcomp_strand):
                    nucs = complement(nucs[::-1])

                if len(nucs.strip()) != region_size: continue

                nuc_counts[offset][nucs] += row.count

    # report the results
    header = ('#nuc','offset','region.size','count','freq') 
    print '\t'.join(header)

    for offset, counts in sorted(nuc_counts.items()):
        sum_counts = sum(counts.values())

        for nuc, count in sorted(counts.items(), key=itemgetter(1), \
                                 reverse=True):
            freq = float(count) / float(sum_counts)
            vals = map(str, [nuc, offset, region_size, count, freq])
            print '\t'.join(vals)

class BedGraphLine:
    ''' class for bedgraph row conversion '''
    def __init__(self, inlist):
        self.chrom = inlist[0]
        self.start = int(inlist[1])
        self.end = int(inlist[2])
        self.count = int(inlist[3])
    
def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... -p POS_BEDGRAPH -n NEG_BEDGRAPH -f FASTA"
    version = "%%prog %s" % __version__
    description = ("calculate nucleotide frequences at modified base "
                   "sites")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Required")
    group.add_option("-p", "--pos-bedgraph", action="store", type='string', 
        metavar="POS_BEDGRAPH", default=None,
        help="bedgraph data, positive strand"
        " (default: %default)")
    group.add_option("-n", "--neg-bedgraph", action="store", type='string', 
        metavar="NEG_BEDGRAPH", default=None,
        help="bedgraph data, negative strand"
        " (default: %default)")
    group.add_option("-f", "--fasta", action="store", type='string', 
        metavar="FASTA", default=None,
        help="FASTA file"
        " (default: %default)")

    parser.add_option_group(group)

    group = OptionGroup(parser, "Variables")
    group.add_option("--revcomp-strand", action="store_true", 
        metavar="CAPTURED STRAND", default=False,
        help="read is reverse comp of captured strand"
        " (default: %default)")

    group.add_option("--offset-min", action="store", type='int', 
        metavar="MIN_OFFSET", default=-15,
        help="min offset from current position to examine. 5' positions are "
        " more negative. (default: %default)")

    group.add_option("--offset-max", action="store", type='int', 
        metavar="MAX_OFFSET", default=15,
        help="max offset from current position to examine. 5' positions are "
        " more negative. (default: %default)")

    group.add_option("--region-size", action="store", type='int', 
        metavar="REGION_SIZE", default=1,
        help="size of region to examine"
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

    if not options.fasta:
        parser.error("specify fasta")
  
    if options.region_size <= 0:
        parser.error("region_size must be > 0")
       
    if len(options.ignore_chrom) != 0 and len(options.only_chrom) != 0:
        parser.error("--ignore-chrom and --only-chrom are mutually exclusive ")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    ignore_chroms = tuple(options.ignore_chrom)
    only_chroms = tuple(options.only_chrom)

    kwargs = {'revcomp_strand':options.revcomp_strand,
              'offset_min':options.offset_min,
              'offset_max':options.offset_max,
              'region_size':options.region_size,
              'verbose':options.verbose,
              'ignore_chroms':ignore_chroms,
              'only_chroms':only_chroms}
    pos_bedgraph = options.pos_bedgraph
    neg_bedgraph = options.neg_bedgraph
    fastafilename = options.fasta
    return nuc_frequencies(pos_bedgraph, neg_bedgraph, fastafilename, 
                           **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
