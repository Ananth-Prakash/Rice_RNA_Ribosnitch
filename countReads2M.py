import sys
import math
import pysam
import HTSeq
import utils

def getAlignmentPos(input_file, paired_end, output_file):
  """
  Write, for each meaningful alignment, its first and last alignment position into file.

  Args:
    input_file: string, input alignment file in SAM or BAM format.
    paired_end: bool, True if input_file contains paired_end reads.
    output_file: a file name to which the data are written

  Notes:
    One line for each alignment, with two columns representing the first and last aligned position.
    For paired-end reads, the two segments are considered together and one-line will be output. 
  """
  almnt_file = utils.openFile(input_file)
  fout = open(output_file, "w")
  if not paired_end:
    for almnt in almnt_file:
      if not almnt.aligned:
        continue
      start = almnt.iv.start + 1
      fout.write("%d\t%d\n" % (start, almnt.iv.end))
  else:
    for first, second in HTSeq.pair_SAM_alignments(almnt_file, bundle=False):
      #if len(bundle) != 1:
      # skip multiple alignments
      # continue
      #first, second = bundle[0]
      if first.aligned and second.aligned and first.iv.chrom == second.iv.chrom:
        start = min(first.iv.start, second.iv.start) + 1
        end = max(first.iv.end, second.iv.end)
        fout.write("%d\t%d\n" % (start, end))

  fout.close()   


def calculateStopCounts(input_file, references, paired_end):
    """
    Calculate stop counts and coverage.
    
    Args:
        input_file: string, input alignment file in SAM or BAM format.
        references: dict, Dictionary for reference sequences, mapping from reference_name to the length of reference sequence.
        paired_end: bool, True if input_file contains paired_end reads.
    
    Returns:
        stop_counts: GenomicArray for stop counts.
        coverage: GenomicArray for coverage.
    """
    almnt_file = utils.openFile(input_file)
    stop_counts = HTSeq.GenomicArray(references, stranded=False, typecode="i")
    coverage = HTSeq.GenomicArray(references, stranded=False, typecode="i")
    if not paired_end:
        for almnt in almnt_file:
            updateStopCountsForSingleEnd(almnt, stop_counts, coverage, False)
    else:
        for bundle in HTSeq.pair_SAM_alignments(almnt_file, bundle=True):
            if not bundle: # null case we pass only one or more
               continue 
            #len(bundle) != 1:
            #    continue
            #print (bundle[0])
            
            first, second = bundle[0] ## allow multiple Jitu
            if first and second: ## must not have aligned Jitu
                if first.aligned and second.aligned and first.iv.chrom == second.iv.chrom:
                   # both ends map to the same reference sequence
                   updateStopCountsForPairedEnd(first, second, stop_counts, coverage)
            #else:
            #    updateStopCountsForSingleEnd(first, stop_counts, coverage, True)
            #    updateStopCountsForSingleEnd(second, stop_counts, coverage, True)
    return stop_counts, coverage


def updateStopCountsForSingleEnd(almnt, stop_counts, coverage, paired_end):
    """
    Update stop counts and coverage for a single segment (eithr single-end read or one end of a paired-end read).
    
    Args:
        almnt: a HTSeq.Alignment object.
        stop_counts: a HTSeq.GenomicArray object for stop counts.
        coverage: a HTSeq.GenomicArray for coverage.
        paired_end: bool, True is almnt is part of a paired-end read.
    """
    if not almnt.aligned:
        return
    if (not paired_end) or (almnt.pe_which == "first"):
        # Update stop counts only when the read is from 1) single-end read; 2) first segment of a paired-end read
        modification_site = HTSeq.GenomicPosition(almnt.iv.chrom, almnt.iv.start, almnt.iv.strand)
        stop_counts[modification_site] += 1
    coverage[almnt.iv] += 1


def updateStopCountsForPairedEnd(first, second, stop_counts, coverage):
    """
    Update stop counts and coverage for paired-end read (assume that both ends are mappend to the same reference).
    
    Args:
        first: a HTSeq.Alignment object, pointing to the first end of a paired-end read.
        second: a HTSeq.Alignment object, pointing to the first end of a paired-end read.
        stop_counts: a HTSeq.GenomicArray object for stop counts.
        coverage: a HTSeq.GenomicArray object for coverage.
    """
    start = min(first.iv.start, second.iv.start)
    end = max(first.iv.end, second.iv.end)
    modification_site = HTSeq.GenomicPosition(first.iv.chrom, start, first.iv.strand)
    stop_counts[modification_site] += 1
    new_iv = HTSeq.GenomicInterval(first.iv.chrom, start, end, first.iv.strand)
    coverage[new_iv] += 1


## ***
def writeStopCounts(stop_counts, coverage, references, output_file):
    """
    Write stop counts and coverage to output file. Modified as Hongjing want: Jitender 071217
    ## added as she Want full length counts
    Args:
        stop_counts: GenomicArray for stop counts.
        coverage: GenomicArray for coverage.
        references: dict, Dictionary for reference sequences, mapping from reference_name to the length of reference sequence.
        output_file: string, output file to write stop counts and coverage.
    """
    fp = open(output_file, "w")
    for chrome, counts_vector in stop_counts.chrom_vectors.items():
        fp.write("%s\n" % chrome)
        chrom_size = references[chrome]
        for index in range(0, chrom_size):
            fp.write("%d\t" % stop_counts[HTSeq.GenomicPosition(chrome, index, ".")])
        #fp.write("0")
        fp.write("\n")
        for index in range(0, chrom_size):
            fp.write("%d\t" % coverage[HTSeq.GenomicPosition(chrome, index, ".")])
        #fp.write("0")    
        fp.write("\n\n")
    fp.close()

