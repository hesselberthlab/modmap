#! /usr/bin/env python

''' origin_analysis: XXX
'''

import sys
import pdb

from pybedtools import BedTool

def origin_anlaysis(timing_bedgraph, origin_bed, max_timing, flank_size,
                    signal_pos_bedgraph, signal_neg_bedgraph,
                    verbose):

    SIDES = ['left','right']
    STRANDS = ['pos','neg']

    selected_ori = select_origins(origin_bed, timing_bedgraph, max_timing,
                                  verbose)

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

