#!/usr/bin/env python3

"""
This script is meant to convert output vcf from bcftools call to bed files.
It is a modified version from the bamCaller.py script on msmc-tools
https://github.com/stschiff/msmc-tools/blob/master/bamCaller.py
You need to store "utils.py" from the same source in the same directory as this script
"""

import sys
import re
import utils
import argparse
import gzip

parser = argparse.ArgumentParser()

parser.add_argument("mask")
args = parser.parse_args()

sites_parser = None

lastPos = 0
line_cnt = 0
chr_ = ""

line_list=[]
for line in sys.stdin:
    line_list.append(line)
    pass
last_line = line
if last_line[0]=='#':#this if block is to check if the inputed vcf is empty. In that cas it will output an empty bed file
    with gzip.open(args.mask, 'wb') as f:
        f.write(''.encode())
else :
    for line in line_list:
        if line[0] == '#':
            continue
        fields = line.strip().split('\t')
        if chr_ == "":
            chr_ = fields[0]
            mask = utils.MaskGenerator(args.mask, chr_)
        else:
            assert fields[0] == chr_, "found multiple chromosomes"
        pos = int(fields[1])
        refAllele = fields[3]
        altAllele = fields[4]
        info = fields[7]
        genotypes = fields[9]
        line_cnt += 1

        if sites_parser is not None:
            while not sites_parser.end and sites_parser.pos < pos:
                sites_parser.tick()

        if re.match("^[ACTGactg]$", refAllele) and re.match("^[ACTGactg\.]$", altAllele):
            dp_match = re.search("DP=(\d+)", info)
            mq_match = re.search("MQ=(\d+)", info)
            fq_match = re.search("FQ=([\d-]+)", info)
            if not (dp_match and mq_match and fq_match):
                continue
            dp = int(dp_match.group(1))
            mq = int(mq_match.group(1))
            fq = int(fq_match.group(1))
            mask.addCalledPosition(pos)


