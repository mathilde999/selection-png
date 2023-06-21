#!/usr/bin/python

"""
This file is specially meant to create Manhattan plot. The input file format is same "CHR   POS VALUE SNPs". SNPs columns are not mandatory
"""
import sys
import pandas as pd
import math
from Selection import Manhattanplot
from Common import Misc

class FisherScore():

    def wrapper(self, args):
        """
        :param args.pbs: should have column CHR BP PBS
        :param args.xpehh: should have column CHR BP NORM
        :return: bed files with the snp with the 99% highest Fisher Score.
        """
        PBS = pd.read_csv(args.pbs, sep='\t')
        XPEHH = pd.read_csv(args.xpehh, sep='\t')
        filename = args.out
        for file in [PBS, XPEHH]:
            if ('CHR' in file) & ('START' in file) & ('END' in file):
                if 'PBS' in file:
                    file_pbs = file[['CHR', 'START', 'END', 'PBS']]
                elif 'NORM' in file:
                    file_xpehh = file[['CHR', 'START', 'END', 'NORM']]
                else:
                    print("Could not find PBS or NORM column in the data. Please check.")
                    sys.exit(1)
            else:
                print("Please check the format of input files. If it has header like: 'CHR', 'START', 'END'")
                sys.exit()

        XPEHH_percentile = self.PercentileRank(file_xpehh)
        PBS_percentile = self.PercentileRank(file_pbs)
        ALL = self.Score(XPEHH_percentile, PBS_percentile)
        ALL = self.pval(ALL)
        all_bed = self.create_bed(ALL)
        all_bed.to_csv(filename + '.all_SNP.bed', index=False, sep='\t')
        top = self.keep_top_snp(ALL)
        top.to_csv(filename + '.topSNP_all_info', index=False, sep='\t')
        bed = self.create_bed(top)
        bed.to_csv(filename + '.topSNP.bed', index=False,header=False, sep='\t')



    def PercentileRank(self, data):
        """
        lambda row: -math.log10(row.iloc[5])-math.log10(row.iloc[10])
        Remove site with NaN Score and set negative Score to 0
        """
        df = data.copy(deep=True)
        SCORE = list(df)[3]
        df = df[df[SCORE].notna()]
        df[SCORE] = df[SCORE].clip(lower=0)

        """
        Creation of the list of rank depending of the score.
        We keep only unique Score before assigning descending rank to each unique Score value
        """
        Score_rank = df[[SCORE]]
        Score_rank = Score_rank.drop_duplicates()
        Score_rank['Rank'] = Score_rank.rank(ascending=False)

        """
        Merging the rank column to global dataframe depending of Score value (when same Score value, same rank and same 
        percentile rank)
        Compute Percentile Rank = Percentage of SNP with rank lower than this rank in its frequency distribution)
        """
        Score_all = pd.merge(df, Score_rank, on=SCORE)
        Score_all['Percentile_Rank'] = Score_all.Rank.rank(pct=True)
        return Score_all

    def Score(self, XPEHH_percentile, PBS_percentile):
        """
        give the Fisher score to each snp depending of their XPEHH and PBS percentile rank
        """
        ALL = XPEHH_percentile.merge(PBS_percentile, how='inner', on=['CHR', 'START', 'END'])
        ALL['Fisher_Score'] = ALL.apply(
            lambda row: -math.log10(row.Percentile_Rank_x) - math.log10(row.Percentile_Rank_y), axis=1)

        return ALL

    def pval(self, ALL):
        """
        !!! Need to be updated - chi-square ?!!
        Computing empirical p-value for the fisher score : the genomic rank of the ith statistic divided by the total number
        of unique values obtained for this statistics in the entire genome (values exactly equal get the same rank qnd
        same p-value. (Deschamps et al.,2016)
        """
        """
        Creation of the list of rank depending of the Fisher score.
        We keep only unique Fisher score before assigning descending rank to each unique Fisher value
        """
        ALL_rank = ALL[['Fisher_Score']]
        ALL_rank = ALL_rank.drop_duplicates()
        ALL_rank['Rank_Fisher'] = ALL_rank.rank(ascending=False)

        """
        Merging the rank column to global dataframe depending of Fisher Score value (when same PBS value, same rank and same percentile 
        rank)
        """
        ALL = pd.merge(ALL, ALL_rank, on='Fisher_Score')
        fisher_count = ALL.Fisher_Score.nunique()
        ALL['Empirical_p'] = ALL.apply(lambda row: row.Rank_Fisher / fisher_count, axis=1)
        return ALL

    def keep_top_snp(self, df):
        """
        keep only the top  snp with highest PBS score
        """
        top = df[df['Fisher_Score'].ge(df['Fisher_Score'].quantile(0.99))]
        return top

    def create_bed(self, df):
        cols = ['CHR', 'START', 'END', 'Fisher_Score']
        fisher = df[cols]
        fisher = fisher.sort_values(['CHR', 'START'], ascending=[True, True])
        return fisher

class ManhattanPlot():
    """
    To plot the out put of selscan
    """

    def plotwrapper(self, args):
        """
        The main and only def for this class
        :param args: args should have file,raw and abs. information which is required to do the manhattan plot
        :return: will call the manhattan plot class to do the manhattan plots.
        """
        data = pd.read_csv(args.file, sep='\t')
        filename = Misc.gettingfilename(args.file)
        plotting = Manhattanplot.Manhattan()
        abs = args.abs
        if ('CHR' in data) & ('BP' in data) & ('Fisher_Score' in data):
            plottingdata = data[['CHR', 'BP', 'Fisher_Score']]
            plottingdata = plottingdata.rename(columns={'Fisher_Score': 'VALUE'})
            return plotting.wrapper(plottingdata, filename, abs)
        else:
            print("Please check the format of input files. If it has header like: CHR   BP  Fisher_Score")
            sys.exit()

class ManhattanPlotHighlighted():
    """
    To plot the out put of selscan with regions to highlight from a bed file
    """

    def plotwrapper(self, args):
        """
        The main and only def for this class
        :param args: args should have file,raw and abs. information which is required to do the manhattan plot
        :return: will call the manhattan plot class to do the manhattan plots.
        """
        data = pd.read_csv(args.file, sep='\t')
        highlightdata = args.highlight
        filename = Misc.gettingfilename(args.file)
        filename = filename + "_highlight"
        plotting = Manhattanplot.Manhattan_Seaborn()
        if ('CHR' in data) & ('BP' in data) & ('Fisher_Score' in data):
            plottingdata = data[['CHR', 'BP', 'Fisher_Score']]
            value = 'Fisher_Score'
            return plotting.wrapper(plottingdata, filename, value, highlightdata)
        else:
            print("Please check the format of input files. If it has header like: CHR   BP  Fisher_Score")
            sys.exit()




