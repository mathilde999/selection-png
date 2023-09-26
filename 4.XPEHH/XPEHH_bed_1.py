#!/bin/python
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='This script will output all the score in the 99th percentile')
parser.add_argument('xpehh', help='The path to the xpehh.out.norm file.')
parser.add_argument("--chr", help="chromosome number", type=int)
args = parser.parse_args()
chrom = "chr"+str(args.chr)

df = pd.read_csv(args.xpehh, delimiter='\t', usecols=['pos', 'normxpehh'])

top = df[df['normxpehh'].le(df['normxpehh'].quantile(0.01))]
top['CHR'] = chrom

data = [top['CHR'], top['pos'], top['normxpehh']]

headers =['CHR', 'pos', 'NORM']

bed = pd.concat(data, axis=1, keys=headers)

bed['BP'] = bed['pos'] - 1

cols = bed.columns.tolist()
cols = ['CHR', 'BP', 'pos', 'NORM']
bed = bed[cols]

bed.to_csv('top_snp_XPEHH_1.bed', header=None, index=None, sep='\t')
