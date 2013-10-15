from collections import defaultdict
from BESST.filter_contigs import get_splits

##
# This class holds statistics of reads mapped to a reference
# genome. The idea is to feed it aligned reads and for each
# read it is fed, the statistics are updated. You may then at
# any point get the statistics you want.
#
class BasePairTracker:
    ##
    # Constructor.
    #
    # @param handle_bp Takes care of a BasePairStats object that have no more 
    #                        stats to update. Takes a tid and BasePairStats object.
    #
    def __init__(self, handle_bp):
        # Stats for each active base pair in their respective chromosome
        self.active_bps = defaultdict(lambda: defaultdict(BasePairStats))

        # Handles base pairs that become inactive
        self.handle_bp = handle_bp

    ##
    # Determines whether we should visit a read or not. In
    # this case we ignore reads that have not been mapped.
    # 
    # @param read The read to visit.
    # 
    # @return True if the read is mapped, false otherwise.
    #
    def accept(self, read):
        if read.rname != read.mrnm or read.is_unmapped or read.mpos <= read.pos:
            return False
        else:
            return True

    ##
    # Updates the statistics by being fed a new aligned read.
    #
    # @param read An aligned read.
    #
    def visit(self, read):
        #self.update_stats_for_soft_clipped_reads(read)
        self.update_stats_for_span_coverage(read)

        start_pos = read.positions[0]
        if read.cigar[0][0] != 0:
            start_pos = read.positions[0] - read.cigar[0][1]
        self.output_bp(read.tid, start_pos)

    ##
    # Output trailing base pairs.
    #
    def finish(self):
        for chrom, contig_active_bps in self.active_bps.iteritems():
            if len(contig_active_bps) > 0:
                max_index = max(contig_active_bps.iterkeys())
                self.output_bp(chrom, max_index + 1)

    ##
    # Update the split read statistics from the given read.
    #
    # Assumes that both ends of the read are aligned.
    #
    # @param read An aligned read.
    #
    def update_stats_for_soft_clipped_reads(self, read):
        contig_active_bps = self.active_bps[ read.tid ]
        for position, is_start in get_splits.get_split_positions(read, 1):
            if is_start:
                contig_active_bps[ position ].num_starting += 1
            else:
                contig_active_bps[ position ].num_ending += 1


    def update_stats_for_span_coverage(self, read):
        contig_active_bps = self.active_bps[ read.tid ]
        #for i in range(min(read.pos, read.mpos), max(read.pos, read.mpos), 100):
        start = int(round(read.pos + 100, -2))
        end = int(round(read.mpos, -2))

        for position in range(start, end, 100):
            #print position
            contig_active_bps[ position ].span_coverage += 1
    ##
    # Feed base pairs with index less than end to the
    # result handler and then remove the base pairs
    # from the active list.
    #
    # @param chr Chromosome of the position.
    # @param end Base pairs with index less than this will 
    #            become inactive.
    #
    def output_bp(self, chrom, end):
        contig_active_bps = self.active_bps[ chrom ]
        inactive_bps = (i for i in contig_active_bps.keys() if i < end)
        for i in sorted(inactive_bps):
            self.handle_bp(chrom, i, contig_active_bps[ i ])
            del contig_active_bps[ i ]

##
# This is a simple class for holding different sufficient
# statistics of base pairs in the reference genome.
#
class BasePairStats:
    def __init__(self):
        # Number of reads ending at the given position, split reads
        # are considered ending in multiple positions.
        self.num_ending = 0.0

        # Number of reads starting at the given position
        self.num_starting = 0.0

        self.span_coverage = 0

##
# For each base pair with soft clipped reads, call the handle_bp
# function.
#
# @param bam_file An opened .bam that is sorted.
# @param handle_bp A callback function that takes a pysam tid, contig position
#                  and a BasePairStats object.
#
def find_soft_clipped_base_pairs(bam_file, handle_bp):
    base_pair_tracker = BasePairTracker(handle_bp)
    for read in bam_file:
        if base_pair_tracker.accept(read):
            base_pair_tracker.visit(read)

    base_pair_tracker.finish()
