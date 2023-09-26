#!/usr/bin/python
"""
This file specifically work with Neales lab GWAS score 10.UKBB UK bio bank files
"""
import numpy
import pandas
import Misc


class UKBBExtract():
    @classmethod
    def wrapper(cls, bedfile, phenofile, plot=False, beta_percent=False, tabix='tabix'):
        """
        This is the main wrapper of UKBBExtract class. Given input files it will return the gwas score of those snps
        and if needed can plot also
        Args:
            bedfile: The bed file where a list of snps are written to be extract. remember the snp position in the end
            or last column and start should be snp-1 bp
            phenofile: the path for the Neales lab phenotype file for 10.UKBB
            plot: In case you want to plot histogram of beta values as well as the extracted beta values. default false
            beta_percent: In case you want to calculate the percentile score of the extracted beta values. default false
            tabix: The path of tabix. default is tabix

        Returns: Will return a pandas dataframe with snps in the row and in columns snp position, pval_EUR, beta_EUR,
        beta_percentile (in case beta_percent=True)
        """
        cls.check_files(bedfile=bedfile, phenofile=phenofile)
        target = cls.extracting_snps(bedfile=bedfile, phenofile=phenofile, tabix=tabix)
        if beta_percent:
            beta = cls.reading_beta(phenofile=phenofile)
            target['beta_percentile'] = cls.percentile_target_beta(target=target, beta=beta)
            if plot:
                cls.plot_hist(target=target, beta=beta)
        return target

    @classmethod
    def check_files(cls, bedfile, phenofile):
        """
        Check all the input files exist or not.
        Args:
            bedfile: the path of the bed file
            phenofile: the path of the phenofile

        Returns: If not it will raise an error (for now print and stop the code) . if ok returns nothing

        """
        import sys
        tobeprint, absentfiles = Misc.file_existence_checker(['bedfile', 'phenofile', 'tbifile'],
                                                             [bedfile, phenofile, phenofile + '.tbi'])
        if len(absentfiles) > 0:
            print(tobeprint)
            sys.exit(1)

    @classmethod
    def extracting_snps(cls, bedfile, phenofile, tabix='tabix'):
        """
        Main commad to extract the snps exist in the bed file and extract it from phenofile so we know about the
        pval_EUR and beta_EUR of that snp
        Args:
            bedfile: the path of the bed file
            phenofile: the path of the phenofile
            tabix: the path of tabix

        Returns: will return pandas dataframe with snps in the row and in columns snp position, pval_EUR, beta_EUR,

        """
        command = Misc.joinginglistbyspecificstring([tabix, '-R', bedfile, phenofile])
        target = pandas.DataFrame([line.split('\t') for line in Misc.systemcode2output(command)])
        target.columns = pandas.read_csv(phenofile, nrows=10, compression='gzip', sep='\t').columns
        snp = target.chr.astype(str) + ":" + target.pos.astype(str)
        snp.name = 'snp'
        target.index = snp
        target = target.loc[:, ['pval_EUR', 'beta_EUR']]
        target = target.replace('NA', numpy.nan).astype(float)
        return target

    @classmethod
    def reading_beta(cls, phenofile):
        """
        in case beta_percent is true it will read the whole file and only extract the beta_EUR
        Args:
            phenofile: the path of the phenofile

        Returns: will return a dataframe of beta values only

        """
        beta = pandas.read_csv(phenofile, compression='gzip', sep='\t', usecols=['beta_EUR'])
        beta = beta.dropna()
        phenotype = Misc.filenamewithoutextension(phenofile).split(".")[0]
        beta.columns = [phenotype]
        return beta

    @classmethod
    def percentile_target_beta(cls, target, beta):
        """
        calculating the target beta score's percentile values to understand how deviated they are from normal
        Args:
            target: pandas dataframe with snps in the row and in columns snp position, pval_EUR, beta_EUR coming from
            extracting_snps()
            beta: the whole beta pandas dataframe coming from reading_beta()

        Returns:

        """
        from scipy import stats
        percentile = [stats.percentileofscore(beta, beta_value) for beta_value in target.beta_EUR.values]
        percentile = pandas.Series(percentile, index=target.index)
        return percentile

    @classmethod
    def plot_hist(cls, target, beta):
        """
        In case plot=True it will first plot all the values in a histogram to see the distribution of values and
        then it will only plot a narrow values where the target snps and or 1,99percentile are used for x limit (thus
        a zoomed plot). this plot will also have the information of target beta score.
        Args:
            target: pandas dataframe with snps in the row and in columns snp position, pval_EUR, beta_EUR coming from
            extracting_snps(), beta_percentile
            beta: the whole beta pandas dataframe coming from reading_beta()

        Returns: will not return anuthing but will plot two plots. one normal histogram (phenotype file name with .pdf
        extension) and another a zoomed plot (phenotype_file_name_zoomed.pdf)

        """
        from matplotlib import pyplot
        ax = beta.plot.hist(bins=1000)
        ax.figure.savefig(beta.columns[0] + '.pdf')
        min_x = min(list(target.dropna().beta_EUR.values) + [beta.quantile(0.01).values[0]])
        max_x = max(list(target.dropna().beta_EUR.values) + [beta.quantile(0.99).values[0]])
        pyplot.xlim(beta.quantile(0.01)[0], beta.quantile(0.99)[0])
        [pyplot.axvline(x=per, color='k', linestyle='dashed') for per in target.beta_EUR.values]
        ax.figure.savefig(beta.columns[0] + '_zoomed' + '.pdf')
