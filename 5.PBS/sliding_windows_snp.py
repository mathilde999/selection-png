import argparse
import pandas as pd
from statistics import mean
import numpy as np
parser = argparse.ArgumentParser(
    description='Will split the given bed file in windows of W number of SNP whit a step of S SNP\n This will output '
                'two files: bed file with the mean pos of each windows and the associated mean PBS (use to plot the '
                'ManhattanPlot\n bed file of the windows with the PBS score in the 99percentile')

parser.add_argument("-b", "--bed", help="bed file, each row is a snp CHR\tBP\tPBS")

parser.add_argument("-w", "--window", help="Number of SNP in each window", type=int)

parser.add_argument("-s", "--step", help="Step in number of SNP", type=int)

args = parser.parse_args()

PBS = pd.read_csv(args.bed, delimiter="\t", usecols=["CHR", "BP", "PBS"])

# sort for CHR and BP
PBS = PBS.sort_values(['CHR', 'BP'], ascending=[True, True])

# remove row where PBS = nan, inf or <0
PBS.replace([np.inf, -np.inf], np.nan, inplace=True)
PBS = PBS[(PBS["PBS"] >= 0) & (PBS["PBS"].notna())]


# create windows including X SNP with a step of Y SNP
# PBS for each of the windows is average of the PBS score of each SNP that were used to do the windows
window_size = args.window
step = args.step
end = window_size - 1
CHR = []
BP_START = []
BP_END = []
MEAN = []

for chr in PBS.iloc[:, 0].unique():
    data_chr = PBS[PBS.iloc[:, 0] == chr]
    for first_snp_i in range(0, len(data_chr.iloc[:, 0]) - end, step):
        list_window_PBS = data_chr.iloc[first_snp_i:first_snp_i + window_size, 2].tolist()
        PBS_window = round(mean(list_window_PBS), 6)
        CHR.append(chr)
        BP_START.append(data_chr.iloc[first_snp_i, 1])
        BP_END.append(data_chr.iloc[first_snp_i + window_size - 1, 1])
        MEAN.append(PBS_window)

dict = {'CHR': CHR, 'BP_START': BP_START, 'BP_END': BP_END, 'PBS': MEAN}
windows_df = pd.DataFrame(dict)
windows_df = windows_df.sort_values(['CHR', 'BP_START'], ascending=[True, True])

windows_df.to_csv('windows_'+str(window_size)+'SNP_step'+str(step)+'.bed', header=['CHR', 'START', 'END', 'PBS'], index=None, sep='\t')

# get top regions
top = windows_df[windows_df.iloc[:, 3].ge(windows_df.iloc[:, 3].quantile(0.99))].reset_index(drop=True)
# sort for CHR and BP
top = top.sort_values(['CHR', 'BP_START'], ascending=[True, True])
top.to_csv('sliding_windows_snp_top_regions_PBS.bed', header=None, index=None, sep='\t')

# get the mean pos of each windows to plot them
windows_df['BP'] = windows_df[['BP_START', 'BP_END']].mean(axis=1)
df = windows_df[['CHR', 'BP', 'PBS']]
# sort for CHR and BP
df = df.sort_values(['CHR', 'BP'], ascending=[True, True])
df.to_csv('mean_pos_PBS.bed', index=None, sep='\t')
