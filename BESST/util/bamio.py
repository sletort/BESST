import os
import pysam

##
# Opens a .bam or .sam file and returns the file
# object.
#
# @param bam_file_path Path to the .bam or .sam.
# 
# @return File object for .bam or .sam.
#
def open_bam_file(bam_file_path):
    bam_file_name, bam_file_ext = os.path.splitext(bam_file_path)
    if bam_file_ext == ".bam":
        return pysam.Samfile(bam_file_path, 'rb')
    elif bam_file_ext == ".sam":
        return pysam.Samfile(bam_file_path, 'r')
    else:
        return IOError("open_bam_file: File must be either .bam or .sam.")
