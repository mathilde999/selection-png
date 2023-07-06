#!/bin/python

import argparse

import Misc
import PBSClass

parser = argparse.ArgumentParser(description='This is the scikit version of PBS. Preferable over previous version')
parser.add_argument('vcf', help='the path for the file which has all the information.',
                    type=lambda x: Misc.args_valid_file(parser, x))
parser.add_argument('--popfile',
                    help='The path of population file. where first columns in indiviuals name present in vcf'
                         ' and the second columns in the populatio name i.e. inds1   pop1', required=True,
                    type=lambda x: Misc.args_valid_file(parser, x))
parser.add_argument('--pbspop', help="Put the name of the populations on which you'll do the pbs. pop1:pop2:pop3. "
                                     "sepearated by under score", required=True)
parser.add_argument('--dropna', help='Lines with nan will be removed', action="store_true")
args = parser.parse_args()

args = parser.parse_args()
pbs = PBSClass.PBS_Sickit.wrapper(vcf=args.vcf, popfile=args.popfile, pbspopstring=args.pbspop, dropna=args.dropna)
print(pbs.to_csv(index=False, sep='\t', na_rep='nan'), end='')
