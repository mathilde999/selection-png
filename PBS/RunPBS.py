#!/bin/python

import argparse
import os

from Common import InfoReaderClass
from Selection.PBS import PBSClass

parser = argparse.ArgumentParser(description='')
parser = argparse.ArgumentParser(description='infofilestructure:\n\t'
                                             'version=plink #either plink or vcf. Sorry this is not done for binary ped format because I am lazy\n\t'
                                             'input=it can be multiple files in a folder (to define a folder end it with "/". it will find all the vcf or vcf.gz in the folder), if you want to do it for single file put only the full path and prefix. Same for plink format\n\t'
                                             'popinfofile=should be a simple tab delimited file where first column is for inds and second column is for populations\n\t'
                                             'A,B,C= Is the three population that you will work with. A is the target population, B is the reference population and C is the outgroup\n\t'
                                             'plink=is the file path of plink. can be module load as well\n\t'
                                             'vcftools=is the file path of vcftools. can be module load as well\n\t'
                                             'python=is the file path of python. can be module load as well\n\t'
                                             'Fst=to define which Fst method to use. You can use "WC","H","R" for Weir and Crick, Hudson and Reynolds. Default is Hudson',formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('info',  help='the path for the info file which has all the information. ')
parser.add_argument('--dont_remove_missing',help='remember the default code will remove those position which have missing data. This will forecfully stop that. not recommended',action="store_true")
parser.add_argument('--cal',  help='This is mainly to calculate the PBS. This is internal dont have to used by yourself.',action="store_true")
args = parser.parse_args()

if args.dont_remove_missing:
    geno = 1
else:
    geno = 0
if args.cal:
    infos = args.info.split()
    # print infos
    pbs = PBSClass.PBS(infos)
    if args.dont_remove_missing:
        pbs.PBScalculationAbsence(infos[0], infos[1], infos[2], infos[3], infos[4], infos[5], infos[6])
    else:
        pbs.PBScalculation(infos[0], infos[1], infos[2], infos[3], infos[4], infos[5], infos[6])
else:
    infos = InfoReaderClass.Readinginfofiles(args.info)
    infodict = infos.reader()
    # print infodict
    pbs = PBSClass.PBS(infodict)
    pbs.wrapper(os.path.realpath(__file__), geno)
