'''
Created on Oct 15, 2013

@author: ksahlin
'''


import argparse
import sys

import pysam

from BESST.filter_contigs.break_point_writer import BreakPointWriter
from BESST.filter_contigs.base_pair_tracker import find_soft_clipped_base_pairs
from BESST.util.bamio import open_bam_file

##
# Finds putative structural variations by looking at break
# points.
#
# @param reads An iterator over reads, i.e. a bam_file.
# @param output_file The file where the break points will be written.
def predict_sv(reads, output_file):
    writer = BreakPointWriter(reads, 5, output_file)
    find_soft_clipped_base_pairs(reads, writer)



##
# Main
#
def main(infile, outfile):
    bam_file = open_bam_file(infile)
    predict_sv(bam_file, outfile)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Find regions that seem to contain a misassembly.")
    arg_parser.add_argument("alignment_file", type=str, help="The aligned reads to infer contig quality from.")
    arg_parser.add_argument("output_file", type=argparse.FileType("w"), help="Output file containing predicted SV regions.")
    args = arg_parser.parse_args()
    main(args.alignment_file, args.output_file)

