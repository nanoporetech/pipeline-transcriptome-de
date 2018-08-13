#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import pandas as pd
from collections import OrderedDict, defaultdict

import pysam
import tqdm

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Count reads mapping to each reference in a BAM file.""")
parser.add_argument(
    '-a', metavar='min_aqual', type=int, help="Minimum mapping quality (0).", default=0)
parser.add_argument(
    '-t', metavar='tsv_file', type=str, help="Save results in tsv format in this file (bam_count_reads.tsv).", default="bam_count_reads.tsv")
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not print progress bar (False).", default=False)
parser.add_argument('bam', help='Input file')


def _count_reads(alignment_file, in_format='BAM', min_aln_qual=0, verbose=False):
    """Count reads mapping to references in a BAM file.
    """
    counts = defaultdict(int)
    if in_format == 'BAM':
        mode = "rb"
    elif in_format == 'SAM':
        mode = "r"
    else:
        raise Exception("Invalid format: {}".format(in_format))

    aln_iter = pysam.AlignmentFile(alignment_file, mode)

    if verbose and in_format == "BAM":
        try:
            total_reads = aln_iter.mapped + aln_iter.unmapped
        except:
            total_reads = None
        sys.stdout.write(
            "Gathering read statistics from file: {}\n".format(alignment_file))
        if in_format == "BAM":
            aln_iter = tqdm.tqdm(aln_iter, total=total_reads)

    for segment in aln_iter:
        if segment.is_unmapped:
            continue
        if segment.mapping_quality >= min_aln_qual:
            counts[segment.reference_name] += 1

    return counts


if __name__ == '__main__':
    args = parser.parse_args()

    counts = _count_reads(args.bam, min_aln_qual=args.a, verbose=not args.Q)
    df = pd.DataFrame(OrderedDict([('Reference', list(counts.keys())), ('Count', list(counts.values()))]))
    df = df.sort_values(['Count', 'Reference'], ascending=False)
    df.to_csv(args.t, sep='\t', index=False)
