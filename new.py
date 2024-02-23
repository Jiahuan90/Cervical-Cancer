#!/usr/bin/env python
# encoding: utf-8

"""
barcodecollapsepe.py
Created by Gabriel Pratt
reads in a .bam file where the first 9 nt of the read name are the barcode
and merge reads mapped to the same position that have the same barcode
"""


from __future__ import print_function


from collections import Counter
import itertools
from optparse import OptionParser  # TODO replace with argparse
import sys
import pysam


help_message = """
barcodecollapse_pe reads in a .bam file where the first 9 nt of the read name
are the barcode and merge reads mapped to the same position that have the same
barcode
"""


def stranded_read_start(read):
    if read.is_reverse:
        return read.reference_end
    else:
        return read.reference_start

def stranded_read_end(read):
    if read.is_reverse:
        return read.reference_start
    else:
        return read.reference_end
		
def output_metrics(metrics_file, total_count, removed_count):
    with open(metrics_file, 'w') as metrics:
        metrics.write("\t".join(["randomer",
                                 "total_count",
                                 "removed_count"])
                      + "\n")
        for barcode in total_count.keys():
            metrics.write("\t".join(map(str, [barcode,
                                              total_count[barcode],
                                              removed_count[barcode]]))
                          + "\n")


def barcode_collapse(in_bam, out_bam):
    number_of_unmapped_mate_pairs = 0
    different_chroms = 0
    removed_count = Counter()
    total_count = Counter()
    result_dict = {}
    with pysam.Samfile(in_bam, 'r') as samfile_read1:
        #samfile_read1 = itertools.islice(samfile1, 0, None, 2)
        for read1 in samfile_read1:
            randomer = read1.qname.split(":")[0]
            start = stranded_read_start(read1)
            stop = stranded_read_end(read1)
            # read1.is_read1
            strand = "-" if read1.is_reverse else "+"
            unique_location = (read1.rname, start, stop, strand, randomer)
            # increment appropriate counter
            total_count[randomer] += 1
            if unique_location in result_dict:
                removed_count[randomer] += 1
                continue
            result_dict[(read1.rname, start, stop, strand, randomer)] = read1
        # ouput barcode collapsed reads
        with pysam.Samfile(out_bam, 'wb', template=samfile_read1) as out_bam:
            for key, (read1) in result_dict.items():
                out_bam.write(read1)
        return total_count, removed_count

def main():
    description = """Paired End randomer aware duplciate removal algorithm."""
    usage  = """
Assumes paired end reads are adjacent in output file (ie needs unsorted bams)
Also assumes no multimappers in the bam file (otherwise behavior is undefined)
"""
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-b", "--bam",
                      dest="bam",
                      help="bam file to barcode collapse")
    parser.add_option("-o", "--out_file",
                      dest="out_file")
    parser.add_option("-m", "--metrics_file",
                      dest="metrics_file")
    (options, args) = parser.parse_args()
    if not (options.bam.endswith(".bam")):
        raise TypeError("%s, not bam file" % options.bam)
    total_count, removed_count = barcode_collapse(options.bam, options.out_file)
    output_metrics(options.metrics_file, total_count, removed_count)
    sys.exit(0)
	
if __name__ == "__main__":
    main()