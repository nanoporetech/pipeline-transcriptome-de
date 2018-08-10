#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from os import path
from functools import reduce

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Merge tab separated files on a given field using pandas.""")
parser.add_argument(
    '-j', metavar='join', type=str, help="Join type (outer).", default="outer")
parser.add_argument(
    '-f', metavar='field', type=str, help="Join on this field (Reference).", default="Reference")
parser.add_argument(
    '-o', metavar='out_tsv', type=str, help="Output tsv (merge_tsvs.tsv).", default="merge_tsvs.tsv")
parser.add_argument(
    '-z', action="store_true", help="Fill NA values with zero.", default=False)
parser.add_argument(
    'tsvs', metavar='input_tsvs', nargs='*', type=str, help="Input tab separated files.")


if __name__ == '__main__':
    args = parser.parse_args()

    dfs = [pd.read_csv(x, sep="\t").rename(columns={'Count': path.basename(x).split('.tsv')[0]}) for x in args.tsvs]

    df_merged = reduce(lambda left, right: pd.merge(left, right, on=args.f, how=args.j), dfs)
    if args.z:
        df_merged = df_merged.fillna(0)

    df_merged.to_csv(args.o, sep="\t", index=False)
