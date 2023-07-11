#!/usr/bin/python
"""
This file to calculate Beta score (Genotype to Phenotype using GWAS) from a vcf file
"""

####Parsing the commands
import argparse
import sys

import NealeLabClass
import Misc
parser = argparse.ArgumentParser(
    description="Given a bed file with the snps it will extract the corresponding pval and beta of European GWAS score"
                "coming from UK Bio Bank (UKBB) Neale's lab files")
parser.add_argument('phenofile', help='the path for the Neales lab phenotype file for UKBB',
                    type=lambda x: Misc.args_valid_file(parser, x))
parser.add_argument('--bed', help='The bed file where a list of snps are written to be extract. remember the snp '
                                      'position in the end or last column and start should be snp-1 bp', required=True,
                    type = lambda x: Misc.args_valid_file(parser, x))
parser.add_argument('--beta_percent', help='In case you want to calculate the percentile score of the extracted beta '
                                           'values', action="store_true")
parser.add_argument('--plot', help='In case you want to plot histogram of beta values as well as the extracted beta '
                                   'values', action="store_true")
parser.add_argument('--tabixpath', help='The path of tabix. default is tabix ', default='tabix')
args = parser.parse_args()
if args.plot:
    if not args.beta_percent:
        print("--plot would not work unless --beta_percent is used")
        sys.exit(1)
betascore = NealeLabClass.UKBBExtract().wrapper(phenofile=args.phenofile, bedfile=args.bedfile,
                                             beta_percent=args.beta_percent, plot=args.plot, tabix=args.tabixpath)
betascore.to_csv(sys.stdout, sep='\t', na_rep='nan')
