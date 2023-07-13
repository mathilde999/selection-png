import pybedtools
import os
import pandas as pd
workdir: os.getcwd()
include: os.getcwd() + '/config.sk'
import scipy.stats
import statsmodels.api as sm
from statsmodels.graphics.gofplots import ProbPlot
import matplotlib.pyplot as plt
import glob
"""
Functions
"""

#find intesection regions on the introgression map
def intersection_region(name_region, region_dict, selection_df):
    chr=region_dict[name_region][0]
    start=region_dict[name_region][1]
    end=region_dict[name_region][2]
    region = " ".join([chr,start, end])
    region_bed = pybedtools.BedTool(region,from_string=True)
    selection_bed = pybedtools.BedTool.from_dataframe(selection_df)
    intersection = selection_bed.intersect(region_bed,u=True) #this output all the sites in selection_bed included in regions_bed
    return intersection

def random_region(selection_score_df, score_name):
    # Convert start and end columns to integers
    selection_score_df[['start', 'end']] = selection_score_df[['start', 'end']].astype(int)
    # Sort the DataFrame by chromosomal position
    selection_score_df = selection_score_df.sort_values(by=['chrom', 'start'])
    list_all_chr=[]
    for chr in range(1,23):
        selection_score_df_chr = selection_score_df[selection_score_df["chrom"] == chr]
        #list where we store the kept regions
        selected_regions=[]
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
        list_all_chr = list_all_chr +selected_regions
    df = pd.DataFrame(list_all_chr,columns=['chrom', 'start', 'end', score_name])
    return df

def read_selection_score(file_name,score_name):
    if score_name == "xpehh":
        selection_df = pd.read_csv(file_name,sep="\t",usecols=["CHR", "BP", "NORM"])
        selection_df["chrom"] = selection_df["CHR"]
        selection_df["start"] = selection_df["BP"] -1
        selection_df["end"] = selection_df["BP"]
        selection_df["xpehh"] = selection_df["NORM"]
        selection_df = selection_df[selection_df['xpehh'].notnull()]
        selection_df = selection_df[["chrom", "start", "end", "xpehh"]]
    elif score_name == "pbs":
        selection_df = pd.read_csv(file_name,sep="\t",usecols=["CHR", "BP", "PBS"])
        selection_df["chrom"] = selection_df["CHR"]
        selection_df["start"] = selection_df["BP"] - 1
        selection_df["end"] = selection_df["BP"]
        selection_df["pbs"] = selection_df["PBS"]
        selection_df = selection_df[selection_df['pbs'].notnull()]
        selection_df = selection_df[["chrom", "start", "end", "pbs"]]
    elif score_name == "pbs_mean":
        selection_df = pd.read_csv(file_name,sep="\t",names=["chrom", "start", "end", "pbs_mean"])
        selection_df = selection_df[selection_df['pbs_mean'].notnull()]
    elif score_name =="fisher":
        selection_df = pd.read_csv(file_name,sep="\t",names=["chrom", "start", "end", "fisher"])
        selection_df = selection_df[selection_df['fisher'].notnull()]
    return selection_df

"""
Get input ready
"""
bed = open(regions_bed)

####create the dict with all the info respectiv to the regions
regions_dic = {}
names_regions = []
for line in bed:
    CHR, START, END = line.split()
    name_region = "_".join([CHR,START,END])
    chrom=CHR.replace("chr","")
    names_regions.append(name_region)
    length = int(END) - int(START)
    regions_dic[name_region] = [CHR, START, END, length]

regions_bed = pd.DataFrame.from_dict(regions_dic, orient='index', columns=['chrom', 'start', 'end', 'length_region'] )


"""
actual snakemake steps
"""
rule all:
    input:
        expand("out_" + score_name + "/output_enrich_" + score_name + "_{name}.pval", name=names_regions)

rule intersection:
    input:
        score = selection_file
    resources:
        mem_mb = 10000
    params:
        score = score_name
    output:
        out="out_"+score_name+"/{name}_"+score_name+"_interesect.out"
    run:
        CHR = regions_dic[wildcards.name][0]
        name=wildcards.name
        score=params.score
        #open selection file depending on selection score used
        score_df = read_selection_score(input.score,score)
        #find intersection regions on the selection score file
        intersection_bed = intersection_region(wildcards.name, regions_dic, score_df)
        intersection_df = intersection_bed.to_dataframe()
        intersection_df.columns = ['chrom', 'start', 'end', score]
        intersection_df.to_csv(output.out, index=False, sep="\t")

# The rule to greate a temp file that will get the information for how many snp we have a selection score info for.
# We exclude the windows that have less than the inputed number of snp's
rule random_windows_to_keep:
    input:
        selection= selection_file
    params:
        score=score_name
    output:
        random_regions="out_"+score_name+"/output_kept_random_region.out"
    resources:
        mem_mb=100000
    run:
        score = params.score
        score_df = read_selection_score(input.selection,score)
        windows_to_keep=random_region(score_df, score)
        windows_to_keep.to_csv(output.random_regions,sep="\t",index=False)

rule p_val_enrich:
    input:
        random_regions = "out_"+score_name+"/output_kept_random_region.out",
        intersect="out_"+score_name+"/{name}_"+score_name+"_interesect.out"
    output:
        out = "out_"+score_name+"/output_enrich_"+score_name+"_{name}.pval"
    params:
        score= score_name
    resources:
        mem_mb=100000
    run:
        score = params.score
        random_regions_df = pd.read_csv(input.random_regions,sep="\t",usecols=[score])
        list_random_mean_score = random_regions_df[score].values.tolist()
        intersection_df = pd.read_csv(input.intersect,sep="\t",usecols=['chrom', 'start', 'end', score])
        if score == "xpehh" and intersection_df[score].mean() < 0:
            max_score_region_interest = intersection_df[score].min()
        else:
            max_score_region_interest = intersection_df[score].max()
        #get pvalue from normal distribution of the random regions
        mean_normal_distrib = random_regions_df[score].mean()
        sd_normal_distrib = random_regions_df[score].std()
        observed_z_score = (max_score_region_interest-mean_normal_distrib)/sd_normal_distrib
        p_value = scipy.stats.norm.sf(abs(observed_z_score))
        line_a = "P-val computed from the normal distribution of the random regions is : "+str(p_value) + "\n"
        line_b = "Zscore computed from the normal distribution of the random regions is : " + str(observed_z_score) + "\n"
        line_c = "Max value in the region " + str(wildcards.name) + " :" + str(max_score_region_interest) + "\n"
        file1 = open(output.out,"w")
        file1.writelines(line_a + line_b)
        file1.close()