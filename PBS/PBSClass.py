#!/usr/bin/python

"""
This file is specially meant to create clases for selscan
"""

import sys

import allel
import numpy
import pandas

import Misc


class PBS_Sickit():
    """
    This will also calculate PBS but using Scikit. Much more robust
    """

    @classmethod
    def wrapper(cls, vcf, popfile, pbspopstring, dropna=False):
        """
        This the wrapper for tyhe sickit pbs. Will produce a out.pbs.gz with required info
        :param vcf:  the path of the vcf file
        :param popfile: the pop file path
        :param pbspopstring: the name of the populations on which you'll do the pbs. pop1_pop2_pop3. sepearated by under
            score
        :param dropna:  Lines with nan will be removed
        :return: will return the pbs
        """

        pbspop = pbspopstring.split(':')
        popdict = cls.check(vcf=vcf, popfile=popfile, pbspop=pbspop)
        pbs = cls.sickitpbs(vcf=vcf, subpopdict=popdict, pbspop=pbspop)
        if dropna:
            pbs = pbs.dropna()
        return pbs

    @classmethod
    def check(cls, vcf, popfile, pbspop):
        """
        To check everything is ok or not
        :param vcf: the vcf file path
        :param popfile: the pop file path
        :param pbspop: the populations in a list format. [pop1,pop2,pop3]
        :return: will return pop dict for sicikt. {pop1:[0,1,2..],pop2:[33,53..]}
        """
        tobeprint, absentfiles = Misc.file_existence_checker(['popfile', 'vcffile'], [popfile, vcf])
        if len(absentfiles) > 0:
            print(tobeprint)
            sys.exit(1)
        popdict = cls.popdict4sickitallel(vcffile=vcf, popfile=popfile)
        if len(pbspop) != 3:
            print("The number of populations has to be 3 for pbs. Please check", pbspop)
        for pop in pbspop:
            if pop not in popdict:
                print(
                    'The population could not be found in popfile. Might be non ind left in the vcf file Please check:',
                    pop)
                sys.exit(1)
        return popdict

    @classmethod
    def popdict4sickitallel(cls, vcffile, popfile):
        """
        To create a popdict worthy of sickit allel class
        :param vcffile: The vcf file path
        :param popfile: the pop file path
        :return: will return pop dict for sicikt. {pop1:[0,1,2..],pop2:[33,53..]}
        """
        pops = pandas.DataFrame(allel.read_vcf_headers(vcffile)[-1])
        pops['pop'] = Misc.adding_pop(pops, popfile)
        subpopdict = {}
        for pop in list(pops['pop'].unique()):
            subpopdict[pop] = list(pops[pops['pop'] == pop].index)
        return subpopdict

    @classmethod
    def sickitpbs(cls, vcf, subpopdict, pbspop):
        """
        The main def for the pbs calculation. It will calcualte the pbs from a given vcf file
        :param vcf: the path of the vcf file
        :param subpopdict:pop dict for sicikt. {pop1:[0,1,2..],pop2:[33,53..]}
        :param pbspop: the populations in a list format. [pop1,pop2,pop3]
        :return:
        """

        callsets = allel.iter_vcf_chunks(vcf, fields=[
            'variants/CHROM', 'variants/POS', 'variants/ID', 'calldata/GT'], chunk_length=int(1e6))[-1]
        df = pandas.DataFrame()
        for callset in callsets:
            geno = allel.GenotypeArray(
                callset[0]['calldata/GT'])
            pbs = cls.pbs_cal(geno=geno, popdict=subpopdict, pbspop=pbspop)
            pbs = numpy.around(pbs, 2)
            temp = pandas.DataFrame({'CHR': callset[0]['variants/CHROM'], 'BP': callset[0]['variants/POS'],
                                     'SNPid': callset[0]['variants/ID'], 'PBS': pbs})
            df = pandas.concat([df, temp], ignore_index=True)
        df = df.reset_index(drop=True)
        return df

    @classmethod
    def pbs_cal(cls, geno, popdict, pbspop):
        """
        instead of using default hudson fst for pbs, this will calculate pbs using Weir and Cockerham
        :param geno: where information about the genotype of the vcf file is written in allel.GenotypeArray format
        :param popdict:  pop dict for sicikt. {pop1:[0,1,2..],pop2:[33,53..], ..}
        :param pbspop: the populations in a list format. [pop1,pop2,pop3]
        :return: will return a numpy array with pbs values for every snp present in the vcf file
        """
        t_ab = cls.fst2t(cls.cal_fst_weir(geno=geno, popdict=popdict, pop1=pbspop[0], pop2=pbspop[1]))
        t_ac = cls.fst2t(cls.cal_fst_weir(geno=geno, popdict=popdict, pop1=pbspop[0], pop2=pbspop[2]))
        t_bc = cls.fst2t(cls.cal_fst_weir(geno=geno, popdict=popdict, pop1=pbspop[1], pop2=pbspop[2]))
        pbs = (t_ab + t_ac - t_bc) / 2
        return pbs

    @classmethod
    def fst2t(cls, fst):
        """
        population divergence time from fst
        :param fst: fst values in numpy array
        :return: will return the divergence time in numpy array format
        """
        return -numpy.log(1 - fst)

    @classmethod
    def cal_fst_weir(cls, geno, popdict, pop1, pop2):
        """
        fst weir cockerham calculation
        :param geno: where information about the genotype of the vcf file is written in allel.GenotypeArray format
        :param popdict:  pop dict for sicikt. {pop1:[0,1,2..],pop2:[33,53..], ..}
        :param pop1: the name of the first pop. should be present in popdict key
        :param pop2: the name of the first pop. should be present in popdict key
        """
        a, b, c = allel.weir_cockerham_fst(geno, [popdict[pop1], popdict[pop2]])
        fst = (numpy.sum(a, axis=1) / (numpy.sum(a, axis=1) + numpy.sum(b, axis=1) + numpy.sum(c, axis=1)))
        return fst
