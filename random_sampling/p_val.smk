import pybedtools
import os
import pandas as pd

workdir: os.getcwd()
include: os.getcwd() + '/config.sk'
import scipy.stats
import numpy as np

"""
Functions
"""

def random_region(selection_score_df, score_name):
    # Convert start and end columns to integers
    selection_score_df[['start', 'end']] = selection_score_df[['start', 'end']].astype(int)
    print("this is the selection df \n", selection_score_df)
    # Sort the DataFrame by chromosomal position
    selection_score_df = selection_score_df.sort_values(by=['chrom', 'start'])
    print("this is the selection df sorted \n", selection_score_df)
    list_all_chr = []
    for chr in range(1,23):
        selection_score_df_chr = selection_score_df[selection_score_df["chrom"] == chr]
        print("this is the selection df per chr \n", selection_score_df)
        #list where we store the kept regions
        selected_regions = []
        # Iterate over the sorted BED file
        for index, row in selection_score_df_chr.iterrows():
            # If it's the first region, add it to the selected regions list
            if len(selected_regions) == 0:
                selected_regions.append([row['chrom'], row['start'], row['end'], row[score_name]])
            else:
                # Get the end position of the previous region
                prev_end = int(selected_regions[-1][2])
                # Get the start position of the current region
                current_start = int(row['start'])
                # Check if the current region is separated by 5MB from the previous region
                if current_start - prev_end > 5000000:
                    selected_regions.append([row['chrom'], row['start'], row['end'], row[score_name]])
        list_all_chr = list_all_chr + selected_regions
    df = pd.DataFrame(list_all_chr,columns=['chrom', 'start', 'end', score_name])
    return df


def get_pval(score_value, random_df, score_name):
    if score_name == "xpehh":
        score_name = score_name
    else:
        score_name = 'log10(' + score_name + ')'
    mean_normal_distrib = random_df[score_name].mean()
    sd_normal_distrib = random_df[score_name].std()
    observed_z_score = (score_value - mean_normal_distrib) / sd_normal_distrib
    if score_name == "xpehh":  #both negative and positive score are signe of selection
        p_value = scipy.stats.norm.sf(abs(observed_z_score))
    else:
        p_value = scipy.stats.norm.sf(observed_z_score)
    return observed_z_score, p_value


def read_selection_score(file_name,score_name):
    if score_name == "xpehh":
        selection_df = pd.read_csv(file_name,sep="\t",usecols=["id", "pos", "normxpehh"])
        chr=selection_df.iloc[0]['id'].split(':')[0]
        selection_df['id']="chr"+str(chr)
        selection_df.columns = ['chrom', 'end', 'xpehh']
        selection_df["start"] = selection_df["end"] -1
        selection_df = selection_df[selection_df['xpehh'].notnull()]
        selection_df = selection_df[["chrom", "start", "end", "xpehh"]]
    elif score_name == "pbs_mean":
        selection_df = pd.read_csv(file_name,sep="\t")
        selection_df.columns = ["chrom", "start", "end", "pbs_mean"]
        selection_df = selection_df[selection_df['pbs_mean'].notnull()]
    elif score_name =="fisher":
        selection_df = pd.read_csv(file_name,sep="\t")
        selection_df.columns = ["chrom", "start", "end", "fisher"]
        selection_df = selection_df[selection_df['fisher'].notnull()]
    return selection_df


"""
actual snakemake steps
"""

rule all:
    input:
        "out_" + score_name + "/top_" + score_name + "_all_region_pval.out",
        "out_" + score_name + "/top_" + score_name + "_pval.out"

# The rule to greate a temp file that will get the information for how many snp we have a selection score info for.
# We exclude the windows that have less than the inputed number of snp's
rule random_windows_to_keep:
    input:
        selection=selection_file
    params:
        score=score_name
    output:
        random_regions="out_" + score_name + "/output_kept_random_region.out"
    resources:
        mem_mb=100000
    run:
        score = params.score
        score_df = read_selection_score(input.selection,score)
        windows_to_keep = random_region(score_df,score)
        if score == 'xpehh':
            pass
        else:
            windows_to_keep['log10(' + score + ')'] = np.log10(windows_to_keep[score])
        windows_to_keep.to_csv(output.random_regions,sep="\t",index=False)

rule p_val_all:
    input:
        random_regions="out_" + score_name + "/output_kept_random_region.out",
        selection=selection_file
    output:
        out="out_" + score_name + "/output_enrich_" + score_name + "_all_pval.gz"
    params:
        score=score_name
    resources:
        mem_mb=200000
    run:
        score = params.score
        #open selection file depending on selection score used
        selection_df = read_selection_score(input.selection,score)
        random_regions_df = pd.read_csv(input.random_regions,sep="\t")
        if score == 'xpehh':
            #get pvalue from normal distribution of the random regions
            g = lambda x: pd.Series(get_pval(x[score],random_regions_df,score))
        else:
            selection_df['log10(' + score + ')'] = np.log10(selection_df[score])
            #get pvalue from normal distribution of the random regions
            g = lambda x: pd.Series(get_pval(x['log10(' + score + ')'],random_regions_df,score))
        selection_df[['Zscore', 'pval']] = selection_df.apply(g,axis=1)
        selection_df.to_csv(output.out,sep="\t",index=False,compression='gzip')


rule extract_p_val_top_region:
    input:
        pval_all = "out_" + score_name + "/output_enrich_" + score_name + "_all_pval.gz",
        top_regions = regions_bed
    output:
        top_region="out_" + score_name + "/top_" + score_name + "_all_region_pval.out",
        top_snp = "out_" + score_name + "/top_" + score_name + "_pval.out"
    params:
        score=score_name
    resources:
        mem_mb=10000
    run:
        score = params.score
        #open selection file depending on selection score used
        pval_df = pd.read_csv(input.pval_all,sep="\t")
        top_df = pd.read_csv(input.top_regions,sep="\t", header=None)
        top_bed = pybedtools.BedTool.from_dataframe(top_df)
        pval_bed = pybedtools.BedTool.from_dataframe(pval_df)
        intersection = pval_bed.intersect(top_bed, wb=True)  #this output all the sites in selection_bed included in regions_bed
        intersection_df = intersection.to_dataframe()
        if score == 'xpehh':
            columns_name = ['chrom', 'start', 'end', score,'Zscore', 'pval']
        else:
            columns_name = ['chrom', 'start', 'end', score,  'log10(' + score + ')', 'Zscore', 'pval']
        intersection_df.columns = columns_name + ['chrom_region', 'start_region', 'end_region']
        intersection_df['top_region'] = intersection_df['chrom_region'].apply(str) + '_' + intersection_df['start_region'].apply(str) + '_'+ intersection_df['end_region'].apply(str)
        intersection_df = intersection_df [columns_name + ['top_region']]
        intersection_df.to_csv(output.top_region,sep="\t",index=False)
        #get the info for the top windows (PBS and Fisher) or the top SNP (xpehh) for each top regions
        max_df = intersection_df.loc[intersection_df.groupby('top_region')[score].apply(lambda x: x.abs().idxmax())]
        max_df.to_csv(output.top_snp,sep="\t",index=False)

