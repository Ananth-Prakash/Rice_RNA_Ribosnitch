import sys
import math
import pysam
import HTSeq
import utils

import countReads2M

# Python script to compute read stop counts
# Modified from RNcode pipeline by Fei
# https://gitlab.com/feideng/RNAcode

input_file = sys.argv[1]
output_count_file = sys.argv[2]

almnt_file = utils.openFile(input_file)
references = utils.getReference(input_file)
paired_end = True
stop_counts, coverage = countReads2M.calculateStopCounts(input_file, references, paired_end)
countReads2M.writeStopCounts(stop_counts, coverage, references, output_count_file)

