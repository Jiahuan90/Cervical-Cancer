__author__ = 'gpratt'


"""

Converts randomer + barcoded fastq files into something that can be barcode collapsed and mapped

"""

from collections import Counter, defaultdict
from itertools import izip
import gzip
import os
from optparse import OptionParser

def reformat_read(name_2, seq_2, plus_2, quality_2,
                  RANDOMER_LENGTH=2):
    """ reformats read to have correct barcode attached
        name - read name
        seq - read sequence
        plus - +
        quality - read quality sequence this is a poor mans datastructure, designed for speed
        barcodes, dictionary of barcodes to search for

        returns str - barcode barcode found, str - randomer identified, str - reformateed read
    """

    randomer = seq_2[:RANDOMER_LENGTH]

    #name_1 = name_1[0] + randomer + ":" + name_1[1:]
    name_2 = name_2[0] + randomer + ":" + name_2[1:]
    seq_2 = seq_2[RANDOMER_LENGTH:]
    quality_2 = quality_2[RANDOMER_LENGTH:]


    #result_1 = name_1 + seq_1 + plus_1 + quality_1
    result_2 = name_2 + seq_2 + plus_2 + quality_2

    return randomer, result_2

if __name__ == "__main__":
    usage = """ takes raw fastq files and demultiplex inline randomer + adapter sequences  """
    parser = OptionParser(usage)
    #parser.add_option("--fastq_1", dest="fastq_1", help="fastq file to barcode seperate")
    parser.add_option("--fastq_2", dest="fastq_2", help="fastq file to barcode seperate")

    #parser.add_option("--out_file_1", dest="out_file_1")
    parser.add_option("--out_file_2", dest="out_file_2")
    parser.add_option("--length", type=int, dest="length", help="Number of randomers on the front of the second read", default=3)


    (options,args) = parser.parse_args()

    #if a gziped file then we reassign open to be the gzip open and continue using that for the rest of the
    #program
    my_open = gzip.open
    #creates different barcode files to assign reads to

    RANDOMER_LENGTH = options.length
    out_file = []

    out_file= [gzip.open(options.out_file_2, 'w'),]


    #reads through initial file parses everything out
    with my_open(options.fastq_2) as fastq_file_2:
        while True:
            try:
                #name_1 = fastq_file_1.next()
                #seq_1 = fastq_file_1.next()
                #fastq_file_1.next() #got to consume the read
                #plus = "+\n" #sometimes the descriptor is here, don't want it
                #quality_1 = fastq_file_1.next()

                name_2 = fastq_file_2.next()
                seq_2 = fastq_file_2.next()
                fastq_file_2.next() #got to consume the read
                plus = "+\n" #sometimes the descriptor is here, don't want it
                quality_2 = fastq_file_2.next()
                #if name_1.split()[0] != name_2.split()[0]:
                    #print name_1, name_2
                    #raise Exception("Read 1 is not same name as Read 2")

                randomer, result_2 = reformat_read(name_2, seq_2, plus, quality_2, RANDOMER_LENGTH)
                #out_file[0].write(result_1)
                out_file[0].write(result_2)
            except StopIteration:
                break

    #cleans up at the end
    for fn in out_file:
        fn.close()

