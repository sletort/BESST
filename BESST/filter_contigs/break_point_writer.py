#!/usr/bin/env python
# encoding: utf-8
"""
break_point_writer.py

Created by Måns Magnusson on 2013-01-31.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os


##
# This class only outputs base pairs which seem to be a breakpoint
# of a structural variation, by looking at the number of reads that
# are clipped.
#
class BreakPointWriter(object):
    """Constructor.
    
     @param bam_file The .bam from which the stats are computed.
     @param num_clipped_threshold The number of clipped reads required to
                                  for a base pair to be a breakpoint.
     @param output_file The file to write the base pairs to.
    """
    def __init__(self, bam_file, num_clipped_threshold, output_file):
        super(BreakPointWriter, self).__init__()
        self.bam_file = bam_file
        self.num_clipped_threshold = num_clipped_threshold
        self.output_file = output_file
    ##
    # Writes base pairs which are above the threshold, and depending
    # on whether clipped reads start or end at the base pair 'start'
    # or 'end' is printed.
    #
    # @param tid Contig id.
    # @param pos Position of the base pair.
    # @param bp_stats Stats of the base pair.
    # 
    def __call__(self, tid, pos, bp_stats=None):
        event = None
        if bp_stats == None:
            self.output_file.write(">{0}\t{1}\n".format(self.bam_file.getrname(tid), pos))

        elif bp_stats.num_starting >= self.num_clipped_threshold:
            event = ("start", pos)
        elif bp_stats.num_ending >= self.num_clipped_threshold:
            event = ("end", pos)
        elif bp_stats.span_coverage:
            event = ('span_cov', pos, bp_stats.span_coverage)

        if event:
            self.output_file.write("{0}\t{1}\t{2}\t{3}\n".format(self.bam_file.getrname(tid), event[ 0 ], event[ 1 ], event[ 2 ]))

if __name__ == '__main__':
    main()
