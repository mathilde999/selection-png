# !/usr/bin/python
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description='Will invert the fourth column of a bed file. Can be use to get the top XPEHH score for the reference pop')

parser.add_argument("-b", "--bed", help="bed file, each row is a snp CHR\tBP\tScore")

args = parser.parse_args()

df = pd.read_csv(args.bed, delimiter='\t',names=['CHR', 'START', 'END', 'NORM'])
df["NORM"] = -(df["NORM"])
df = df[["CHR", "START", "END", "NORM"]]
df.to_csv("inverted.xpehh.bed", index=None,header=None, sep='\t')
