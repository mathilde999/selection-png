#!/usr/bin/python
"""
This file to compute fisher score from XPEHH and 5.PBS. Will have to be updated to add more scores
"""

import argparse,os
import FisherScoreClass

parser = argparse.ArgumentParser(description='file format:\n\t'
                                             'Should have CHR   START END  SCORE   SNP.\n\t'
                                             'SNP column is optional.\n\t'
                                             , formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--xpehh',  help='the path for whole genome 4.XPEHH output files. It is expecting to have normalized')
parser.add_argument('--pbs',  help='the path for whole genome 5.PBS output files. ')
parser.add_argument('--out',  help='the path for the output file')

args = parser.parse_args()


if not os.path.exists(args.xpehh):
    print("The files could not be found. Please check: ",args.xpehh)
if not os.path.exists(args.pbs):
    print("The files could not be found. Please check: ",args.pbs)


scoring = FisherScoreClass.FisherScore()
scoring.wrapper(args)

