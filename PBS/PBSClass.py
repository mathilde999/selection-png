#!/usr/bin/python

"""
This file is specially meant to create clases for selscan
"""

import os
import pandas
from Selection import Manhattanplot
from Common import Misc
import allel
import sys
class PBS():
    """
    to understand the all the infos and then submit the command in sbatch system to calculate PBS
    """

    def __init__(self, infodict):
        self.infodict = infodict

    def wrapper(self, ownself, geno):
        import argparse
        from SubmittingJobs import SubmittingJobsClass
        self.infodictcleaner()
        modules, pop = self.checker()
        # print self.vcf
        ##creating all the folders
        if not os.path.exists("temp"):
            os.makedirs("temp")
        if not os.path.exists("logfiles"):
            os.makedirs("logfiles")
        if not os.path.exists("shfiles"):
            os.makedirs("shfiles")
        pop['family'] = pop['inds']
        pop = pop[['inds', 'family', 'pop']]
        pop.to_csv('temp/plink.pop', sep='\t', index=False, header=False)

        ##creating commands for per file
        if (self.infodict['version'] == 'vcf'):
            for file in self.vcf:
                filename = Misc.gettingfilename(file)
                commands = "\n".join(self.shscriptperfile(file, modules, ownself, geno))
                args = argparse.Namespace(command=commands, name=filename, amd=True, log=True, sh=True)
                submission = SubmittingJobsClass.sbatchoneliner(args)
                submitcommand = submission.argsparser()
                os.system(submitcommand)

        else:
            for file in self.map:
                filename = Misc.gettingfilename(file)
                commands = "\n".join(self.shscriptperfile(file, modules, ownself, geno))
                args = argparse.Namespace(command=commands, name=filename, log=True, sh=True)
                submission = SubmittingJobsClass.sbatchoneliner(args)
                submitcommand = submission.argsparser()
                os.system(submitcommand)

    def infodictcleaner(self):
        """
        To cleanup infodict from white spaces and other stuffs
        :return: returned modified and better version of infodict
        """
        self.infodict['version'] = self.infodict['version'].split()[0].lower()
        self.infodict['input'] = self.infodict['input'].split()[0]
        self.infodict['popinfofile'] = self.infodict['popinfofile'].split()[0]
        self.infodict['A'] = self.infodict['A'].split()[0]
        self.infodict['B'] = self.infodict['B'].split()[0]
        self.infodict['C'] = self.infodict['C'].split()[0]
        self.infodict['Fst'] = self.infodict['Fst'].split()[0].upper()[0]

    def checker(self):
        """
        To check if all the files are correct and everything is ok or not
        :return: it will return all the modules and pop dataframe. It will also automatically added input file
        """
        import glob, sys
        modules = []
        ####input problems
        if (self.infodict['version'] == 'vcf'):
            self.infodict['vcftools'] = Misc.module_checker(self.infodict['vcftools'], "vcftools", modules)
            if self.infodict['input'][-1] == "/":
                if len(glob.glob(self.infodict['input'] + "*.vcf")) == 0:
                    self.vcf = glob.glob(self.infodict['input'] + "*.vcf.gz")
                    self.vcfinputversion = 'gzvcf'
                else:
                    self.vcf = glob.glob(self.infodict['input'] + "*.vcf")
                    self.vcfinputversion = 'vcf'
            else:
                if os.path.exists(self.infodict['input'] + ".vcf"):
                    self.vcf = [self.infodict['input'] + ".vcf"]
                    self.vcfinputversion = 'vcf'
                elif os.path.exists(self.infodict['input'] + ".vcf.gz"):
                    self.vcf = [self.infodict['input'] + ".vcf.gz"]
                    self.vcfinputversion = 'gzvcf'
                else:
                    print("Could not find the input file", self.infodict['input'] + ".vcf or",
                          self.infodict['input'] + ".vcf.gz")
                    sys.exit(1)
        elif (self.infodict['version'] == 'plink'):
            if self.infodict['input'][-1] == "/":
                self.map = glob.glob(self.infodict['input'] + "*.map")
            else:
                if os.path.exists(self.infodict['input'] + ".map"):
                    self.map = [self.infodict['input'] + ".map"]
                else:
                    print("Could not find the input file", self.infodict['input'] + ".map")
                    sys.exit(1)
        else:
            print("This is only for vcf and plink version.", self.infodict['version'], " version does not exist yet")
            sys.exit(1)
        ###module problems continues
        self.infodict['plink'] = Misc.module_checker(self.infodict['plink'], "plink", modules)
        self.infodict['python'] = Misc.module_checker(self.infodict['python'], "python", modules)
        ###population informations
        pop = pandas.read_table(self.infodict['popinfofile'], names=['inds', 'pop'])
        popnames = pandas.Series.unique(pop['pop'])
        if self.infodict['A'] not in popnames:
            print("could not find the A population in pop info file", self.infodict['A'])
            sys.exit(1)
        if self.infodict['B'] not in popnames:
            print("could not find the B population in pop info file", self.infodict['B'])
            sys.exit(1)
        if self.infodict['C'] not in popnames:
            print("could not find the C population in pop info file", self.infodict['C'])
            sys.exit(1)
        if self.infodict['Fst'] != 'W' and self.infodict['Fst'] != 'H' and self.infodict['Fst'] != 'R':
            print(
                "Fst method can only be Fst=WC or Fst=H or Fst=R as it does not match with anything we put it as Fst=Hudson")
            self.infodict['Fst'] = 'H'

        return modules, pop

    def shscriptperfile(self, file, modules, ownself, geno):
        """
        This wil return commands per file, which then can be saved and run
        :param file: the input file
        :param modules: the list of modules
        :param ownself: the path of RunPBS and is need to be use for different stuff
        :param geno: to tell to remove missing information position. 0 means no missing information. and 1 mean it will keep all the missing informations
        :return: the total command that needs to be run
        """
        commands = []
        if len(modules) > 0:
            commands = commands + modules
        filename = Misc.gettingfilename(file)
        if self.infodict['version'] == 'vcf':
            if (self.vcfinputversion == "gzvcf"):
                filename = filename[:-7]
            else:
                filename = filename[:-4]
            commands.append(" ".join([self.infodict['vcftools'], "--" + self.vcfinputversion, file,
                                      "--keep temp/plink.pop --plink --remove-indels --non-ref-ac 1 --out temp/" + filename]))
            file = "temp/" + filename + ".map"
            removemap = ["rm -f " + "temp/" + filename + ".map", "rm -f " + "temp/" + filename + ".ped"]
            # removemap = []
        else:
            filename = Misc.gettingfilename(file)
            filename = filename[:-4]
            removemap = []
        commands.append(" ".join([self.infodict['plink'], "--file", file[:-4], "--geno", str(geno),
                                  "--freq --within temp/plink.pop --out temp/" + filename]))
        if geno == 0:
            commands.append(" ".join(
                [self.infodict['python'], ownself, '"'"temp/" + filename + ".frq.strat", file, self.infodict['A'],
                 self.infodict['B'], self.infodict['C'], self.infodict['Fst'], filename, '"', "--cal"]))
        else:
            commands.append(" ".join(
                [self.infodict['python'], ownself, '"'"temp/" + filename + ".frq.strat", file, self.infodict['A'],
                 self.infodict['B'], self.infodict['C'], self.infodict['Fst'], filename, '"', "--cal",
                 "--dont_remove_missing"]))
        removemap.append("rm -f " + "temp/" + filename + ".frq.strat")

        return commands + removemap

    def PBScalculation(self, freq, map, A, B, C, method, out):
        """
        To calculate PBS and write it an output file
        :param freq: the plink output of --freq --within input.pop which have freq by pop
        :param map: the map file from plink as freq file do not have the position
        :param A: A is the target population
        :param B: B is the reference population
        :param C: C is the outgroup
        :param method: Fst can be calculate in 3 mentods. Weir an Crick, Hudson and Raynolds. This is to define which Fst method to use
        :param out: The output file prefix
        :return: It will return nothing as it will already save the file
        """
        import math
        popA = ["A", "A"]
        popB = ["B", "B"]
        popC = ["C", "C"]
        with open(freq, "r") as f, open(out + '.pbs', "w") as o, open(map, 'r') as pos:
            o.write(
                Misc.joinginglistbyspecificstring(['CHR', 'BP', 'SNPid', 'Fa', 'Fb', 'Fc', 'Tab', 'Tac', 'Tbc', 'PBS'],
                                                  '\t') + '\n')
            for line in f:
                line = line.split()
                if line[2] == A:
                    popA = line
                if line[2] == B:
                    popB = line
                if line[2] == C:
                    popC = line
                if popA[1] == popB[1] and popA[1] == popC[1] and (line[2] == A or line[2] == B or line[2] == C):
                    posline = pos.readline().split()[3]
                    popA[5], popB[5], popC[5], popA[7], popB[7], popC[7] = float(popA[5]), float(popB[5]), float(
                        popC[5]), float(popA[7]), float(popB[7]), float(popC[7])
                    if popA[7] > 0:
                        fa = popA[5]
                    else:
                        fa = float('Nan')
                        popA[5]=float('Nan')
                    if popB[7] > 0:
                        fb = popB[5]
                    else:
                        fb = float('Nan')
                        popB[5] = float('Nan')
                    if popC[7] > 0:
                        fc = popC[5]
                    else:
                        fc = float('Nan')
                        popC[5] = float('Nan')

                    ab = round(self.fstcal(popA[7], popB[7], popA[5], popB[5], method), 4)
                    ac = round(self.fstcal(popA[7], popC[7], popA[5], popC[5], method), 4)
                    bc = round(self.fstcal(popC[7], popB[7], popC[5], popB[5], method), 4)

                    if math.isnan(ab) or math.isnan(ac) or math.isnan(bc):
                        o.write("\t".join([str(line[0]), posline, str(line[1]), str(fa), str(fb), str(fc), str(ab), str(ac), str(bc),str(float('NaN'))]) + "\n")
                    elif (ab == 1 or ac == 1) or bc == 1:
                        o.write("\t".join(
                            [str(line[0]), posline, str(line[1]), str(fa), str(fb), str(fc), str(ab), str(ac), str(bc),
                             str(float('NaN'))]) + "\n")
                    else:

                        tab = -math.log(1 - ab)
                        tac = -math.log(1 - ac)
                        tbc = -math.log(1 - bc)
                        pbs = (tab + tac - tbc) / 2
                        o.write("\t".join(
                            [str(line[0]), posline, str(line[1]), str(fa), str(fb), str(fc), str(ab), str(ac), str(bc),
                             str(pbs)]) + "\n")
                    ###in case no snpid given we need to clear it
                    if popA[1] == '.':
                        popA[1] = 'A'
                        popB[1] = 'B'
                        popC[1] = 'C'

    def fstcal(self, nA, nB, pA, pB, fstmethod):
        """
        With given information (number of individuals and frqeuency), it will calculate the Fst pos by pos
        :param nA: Number of individual from first population
        :param nB: Number of individual from second population
        :param pA: Frequency of SNP from first population
        :param pB: Frequency of SNP from second population
        :param fstmethod: It has 3 methods for Fst: Weir an Crick, Hudson or Reynolds
        :return: Will return the Fst value
        """
        import math
        if math.isnan(nA) or math.isnan(nB) or math.isnan(pA) or math.isnan(pB):
            return float('NaN')

        ###fst WC###
        def fstWC(nA, nB, pA, pB):
            nA = nA
            nB = nB
            if (nA > 0) and (nB > 0):
                num = 2 * ((nA * nB) / (nA + nB)) * (1 / (nA + nB - 2)) * ((nA * pA * (1 - pA)) + (nB * pB * (1 - pB)))
                den = (((nA * nB) / (nA + nB)) * ((pA - pB) ** 2)) + (2 * (nA * nB) / (nA + nB) - 1) * (
                        1 / (nA + nB - 2)) * (nA * pA * (1 - pA) + nB * pB * (1 - pB))
                if den == 0:
                    return float('NaN')
                else:
                    return 1 - (num / den)

            else:
                return float('NaN')

        ###fst Hudson###
        def fstH(nA, nB, pA, pB):
            nA = nA
            nB = nB
            X = (pA - pB) ** 2
            Y = (pA * (1 - pA)) / (nA - 1)
            W = (pB * (1 - pB)) / (nB - 1)
            Z = pA * (1 - pB)
            K = pB * (1 - pA)
            num = X - Y - W
            den = Z + K
            if den == 0:
                return float('NaN')
            else:
                return (num / den)

        ###fst Reynolds###
        def fstR(nA, nB, pA, pB):
            alfa1 = 1 - ((pA ** 2) + ((1 - pA) ** 2))
            alfa2 = 1 - ((pB ** 2) + ((1 - pB) ** 2))
            Al = (0.5 * (((pA - pB) ** 2) + (((1 - pA) - (1 - pB)) ** 2))) - (
                    ((nA + nB) * (nA * alfa1 + nB * alfa2)) / ((4 * nA * nB) * (nA + nB - 1)))
            AlBl = (0.5 * (((pA - pB) ** 2) + (((1 - pA) - (1 - pB)) ** 2))) + (
                    ((4 * nA * nB - nA - nB) * (nA * alfa1 + nB * alfa2)) / ((4 * nA * nB) * (nA + nB - 1)))
            if AlBl == 0:
                return float('NaN')
            else:
                return (Al / AlBl)

        #####
        if fstmethod == 'W':
            return fstWC(nA, nB, pA, pB)
        elif fstmethod == 'R':
            return fstR(nA, nB, pA, pB)
        else:
            return fstH(nA, nB, pA, pB)

    def biggest_count(self, freq):
        """
        To calcualte the biggest number of individuals present of avaliable in the file
        :param freq: the plink output of --freq --within input.pop which have freq by pop
        :return: will return the biggest number of haplotypes (2xindividuals) per population
        """
        totalnumbers = {}
        with open(freq, 'r') as f:
            f.readline()
            for line in f:
                line = line.split()
                line[7] = int(line[7])
                if line[2] in totalnumbers:
                    if totalnumbers[line[2]] < line[7]:
                        totalnumbers[line[2]] = line[7]
                else:
                    totalnumbers[line[2]] = line[7]
        return totalnumbers

    def PBScalculationAbsence(self, freq, map, A, B, C, method, out):
        """
        To calculate PBS and write it an output file. this is specfically if some of the place we have absence. Thus it will also the portion of data that is avaliable under CountA,CountB,CountCq
        :param freq: the plink output of --freq --within input.pop which have freq by pop
        :param map: the map file from plink as freq file do not have the position
        :param A: A is the target population
        :param B: B is the reference population
        :param C: C is the outgroup
        :param method: Fst can be calculate in 3 mentods. Weir an Crick, Hudson and Raynolds. This is to define which Fst method to use
        :param out: The output file prefix
        :return: It will return nothing as it will already save the file
        """
        import math
        popA = ["A", "A"]
        popB = ["B", "B"]
        popC = ["C", "C"]
        # max_count = self.biggest_count(freq)
        with open(freq, "r") as f, open(out + '.pbs', "w") as o, open(map, 'r') as pos:
            o.write(Misc.joinginglistbyspecificstring(
                ['CHR', 'BP', 'SNPid', 'CountA', 'CountB', 'CountC', 'Fa', 'Fb', 'Fc', 'Tab', 'Tac', 'Tbc', 'PBS'],
                '\t') + '\n')
            for line in f:
                line = line.split()
                if line[2] == A:
                    popA = line
                if line[2] == B:
                    popB = line
                if line[2] == C:
                    popC = line
                if popA[1] == popB[1] and popA[1] == popC[1] and (line[2] == A or line[2] == B or line[2] == C):

                    temp=pos.readline()
                    try:
                        posline = temp.split()[3]

                        popA[5], popB[5], popC[5], popA[7], popB[7], popC[7] = float(popA[5]), float(popB[5]), float(
                            popC[5]), float(popA[7]), float(popB[7]), float(popC[7])
                        if popA[7] > 0:
                            fa = popA[5]
                        else:
                            fa = float('Nan')
                            popA[5] = float('Nan')
                        if popB[7] > 0:
                            fb = popB[5]
                        else:
                            fb = float('Nan')
                            popB[5] = float('Nan')
                        if popC[7] > 0:
                            fc = popC[5]
                        else:
                            fc = float('Nan')
                            popC[5] = float('Nan')
                        ab = round(self.fstcal(popA[7], popB[7], popA[5], popB[5], method),4)
                        ac = round(self.fstcal(popA[7], popC[7], popA[5], popC[5], method),4)
                        bc = round(self.fstcal(popC[7], popB[7], popC[5], popB[5], method),4)

                        count_a = int(popA[7])
                        count_b = int(popB[7])
                        count_c = int(popC[7])
                        ##freq calculation
                        if popA[7] > 0:
                            fa = popA[5]
                        else:
                            fa = float('Nan')
                        if popB[7] > 0:
                            fb = popB[5]
                        else:
                            fb = float('Nan')
                        if popC[7] > 0:
                            fc = popC[5]
                        else:
                            fc = float('Nan')

                        if math.isnan(ab) or math.isnan(ac) or math.isnan(bc):
                            o.write(Misc.joinginglistbyspecificstring(
                                [line[0], posline, line[1], count_a, count_b, count_c, fa, fb, fc, ab, ac, bc,
                                 float('NaN')], '\t') + "\n")

                        elif (ab == 1 or ac == 1) or bc == 1:
                            o.write(Misc.joinginglistbyspecificstring(
                                [line[0], posline, line[1], count_a, count_b, count_c, fa, fb, fc, ab, ac, bc,
                                 float('NaN')], '\t') + "\n")
                        else:
                            tab = -math.log(1 - ab)
                            tac = -math.log(1 - ac)
                            tbc = -math.log(1 - bc)
                            pbs = (tab + tac - tbc) / 2
                            # tab = round(-math.log(1 - ab),4)
                            # tac = round(-math.log(1 - ac),4)
                            # tbc = round(-math.log(1 - bc),4)
                            # pbs = round((tab + tac - tbc) / 2,4)
                            o.write(Misc.joinginglistbyspecificstring(
                                [line[0], posline, line[1], count_a, count_b, count_c, fa, fb, fc, ab, ac, bc, pbs],'\t') + "\n")

                        ###in case no snpid given we need to clear it
                        if popA[1] == '.':
                            popA[1] = 'A'
                            popB[1] = 'B'
                            popC[1] = 'C'
                    except:
                        pass


class PBS_Sickit():
    """
    This will also calculate PBS but using Sickit. Much more robust
    """

    @classmethod
    def popdict4sickitallel(cls,vcffile, popfile):
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
    def check(cls,vcf,popfile,pbspop):
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
        if len(pbspop)!=3:
            print ("The number of populations has to be 3 for pbs. Please check",pbspop)
        for pop in pbspop:
            if pop not in popdict:
                print(
                    'The population could not be found in popfile. Might be non ind left in the vcf file Please check:',
                    pop)
                sys.exit(1)
        return popdict
    @classmethod
    def sickitpbs(cls,vcf, subpopdict, pbspop):
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
            ac = allel.GenotypeArray(
                callset[0]['calldata/GT']).count_alleles_subpops(subpopdict)
            Tac1 = ac[pbspop[0]].sum(1)
            Tac2 = ac[pbspop[1]].sum(1)
            Tac3 = ac[pbspop[2]].sum(1)
            sickitpbs = allel.pbs(ac[pbspop[0]], ac[pbspop[1]],
                                  ac[pbspop[2]], 1, normed=False)

            temp = pandas.DataFrame({'CHR': callset[0]['variants/CHROM'], 'BP': callset[0]
            ['variants/POS'], 'SNPid': callset[0]['variants/ID'], 'CountA': Tac1, 'CountB': Tac2, 'CountC': Tac3,
                                     'PBS': sickitpbs})
            df = pandas.concat([df, temp], ignore_index=True)
        df = df.reset_index(drop=True)
        return df
    @classmethod
    def wrapper(cls,vcf,popfile,pbspopstring,name=None):
        """
        This the wrapper for tyhe sickit pbs. Will produce a out.pbs.gz with required info
        :param vcf:  the path of the vcf file
        :param popfile: the pop file path
        :param pbspopstring: the name of the populations on which you'll do the pbs. pop1_pop2_pop3. sepearated by under score
        :param name:  output prefix for the file. Default is the vcf file name
        :return: wil lsave the pbs in a csv format
        """

        pbspop=pbspopstring.split('_')
        popdict=cls.check(vcf=vcf,popfile=popfile,pbspop=pbspop)
        pbs=cls.sickitpbs(vcf=vcf,subpopdict=popdict,pbspop=pbspop)
        if name is None:
            name = Misc.filenamewithoutextension_checking_zipped(vcf) + '.pbs.gz'
        pbs.to_csv(name,index=False)


class PBSPlot():
    def plotwrapper(self, data, filename):
        """
        The main and only def for this class
        :param args: args should have file information which is required to do the manhattan plot
        :return: will call the manhattan plot class to do the manhattan plots.
        """
        plotting = Manhattanplot.Manhattan()
        if ('CHR' in data) & ('BP' in data):
            if 'PBS' in data:
                plottingdata = data[['CHR', 'BP', 'PBS']]
                plottingdata = plottingdata.rename(columns={'PBS': 'VALUE'})
                print("les valeurs sont renommee BIS")
                return plotting.wrapper(plottingdata, filename, False)
            else:
                print("Could not find PBS column in the data. Please check.")
                sys.exit(1)
        else:
            print("Please check the format of input files. If it has header like: CHR   BP")
            sys.exit()


class PBSHighlightedPlot():
    def plotwrapper(self,data, highlight, filename):
        """
        The main and only def for this class
        :param args: args should have file information which is required to do the manhattan plot
        :return: will call the manhattan plot class to do the manhattan plots.
        """
        highlightdata = highlight
        filename = filename + "_highlight"
        plotting = Manhattanplot.Manhattan_Seaborn_Highlight()
        if ('CHR' in data) & ('BP' in data):
            if 'PBS' in data:
                plottingdata = data[['CHR', 'BP', 'PBS']]
                value = 'PBS'
                return plotting.wrapper(plottingdata, filename, value, highlightdata)
            else:
                print("Could not find PBS column in the data. Please check.")
                sys.exit(1)
        else:
            print("Please check the format of input files. If it has header like: CHR   BP")
            sys.exit()
