import pybedtools
import os
import pandas as pd
workdir: os.getcwd()
include: os.getcwd() + '/config.sk'
import scipy.stats

"""
Functions
"""

def get_chr_end(CHR, autosome):
    len_chr_df = autosome.loc[autosome['chrom'] == CHR]
    len_chr= len_chr_df.iloc[0]['end']
    return len_chr

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

def random_region(selection_score_bed,autosomes_bed, len_region, minimal_number):
    #create non-overlapping windows of the lenght of the target region + 1 MB for all the genome
    windows_bed = autosomes_bed.window_maker(autosomes_bed,w=len_region, s=len_region+1000000)
    #convert to dataframe
    windows_df = windows_bed.to_dataframe()
    windows_df.columns=['chrom', 'start', 'tmp_end']
    #for every windows, we will change the end pos by start + len region
    #All of this is to have only windows of the lenght of the target region separated by 1MB
    windows_df["end"]=windows_df["start"]+len_region
    windows_df=windows_df[['chrom', 'start', 'end']]
    #convert back to bed format
    windows_bed = pybedtools.BedTool.from_dataframe(windows_df)
    #for each windows output the number of snp with selection scores
    intersection = windows_bed.intersect(selection_score_bed,c=True)
    intersection_df = intersection.to_dataframe()
    intersection_df.columns = ['chrom', 'start', 'end', 'count_sites']
    intersection_df = intersection_df[intersection_df['count_sites'] >= minimal_number]
    intersection_df = intersection_df[['chrom', 'start', 'end']]
    return intersection_df

def random_intersection(random_windows_df, selection_bed, score_name):
    list_of_list=[]
    for i in range(len(random_windows_df.index)):
        chr = random_windows_df.iloc[i]["chrom"]
        start = random_windows_df.iloc[i]["start"]
        end = random_windows_df.iloc[i]["end"]
        line = " ".join([str(chr), str(start), str(end)])
        random_region_bed = pybedtools.BedTool(line,from_string=True)
        # intersect with selection scores
        intersection = selection_bed.intersect(random_region_bed, u=True)
        intersection_df = intersection.to_dataframe()
        intersection_df.columns = ['chrom', 'start', 'end', score_name]
        #print("this is the intesection ", intersection_df.head())
        mean_score_region = intersection_df[score_name].mean()
        list_mean_score = [chr, start, end, mean_score_region]
        list_of_list.append(list_mean_score)
        #print("this is the output list ", list_of_list)
    return list_of_list

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


bed = open(autosome)
####create the dict with all the info respectiv to the chromosomes
chr_dic = {}
CHR_list = []
for line in bed:
    chrom, start, end = line.split()
    CHR_list.append(chrom)
    chr_dic[chrom] = end

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
        #open selection file depending of selection score used
        score_df = read_selection_score(input.score,score)
        #find intersection regions on the selection score file
        intersection_bed = intersection_region(wildcards.name, regions_dic, score_df)
        #from introgression to datafram so we can count the number of match
        intersection_df = intersection_bed.to_dataframe()
        intersection_df.columns = ['chrom', 'start', 'end', score]
        intersection_df.to_csv(output.out, index=False, sep="\t")

# The rule to greate a temp file that will get the information for how many snp we have a selection score info for.
# We exclude the windows that have less than the inputed number of snp's
rule random_windows_to_keep:
    input:
        selection= selection_file,
        autosomes=autosome
    params:
        score=score_name,
        minimal_line=minimal_line_inter
    output:
        random_regions="out_"+score_name+"/output_kept_random_region_{name}.out"
    resources:
        mem_mb=100000
    run:
        score = params.score
        score_df = read_selection_score(input.selection,score)
        score_bed = pybedtools.BedTool.from_dataframe(score_df)
        #load autosome bed file
        autosomes_df = pd.read_csv(input.autosomes,sep="\t",names=["chrom", "tmp_start", "tmp_end"])
        autosomes_df["start"] = autosomes_df["tmp_start"] + 10000
        autosomes_df["end"] = autosomes_df["tmp_end"] - 10000
        autosomes_df=autosomes_df[["chrom", "start", "end"]]
        #convert to bedtool object
        autosomes_bed = pybedtools.BedTool.from_dataframe(autosomes_df)
        len_region = regions_dic[wildcards.name][3]
        windows_to_keep=random_region(score_bed, autosomes_bed, len_region, params.minimal_line)
        random_region_list = random_intersection(windows_to_keep,score_bed,score)
        #print("this is the list of list again ", random_region_list)
        df = pd.DataFrame(random_region_list,columns=['CHR', 'start', 'end', score])
        #print("this is the df ",df.head())
        df.to_csv(output.random_regions,sep="\t",index=False)

rule p_val_enrich:
    input:
        random_regions = "out_"+score_name+"/output_kept_random_region_{name}.out",
        intersect="out_"+score_name+"/{name}_"+score_name+"_interesect.out"
    output:
        out = "out_"+score_name+"/output_enrich_"+score_name+"_{name}.pval"
    params:
        score= score_name
    resources:
        mem_mb=100000
    run:
        with open(input.random_regions) as myfile:
            sentence= "We didn't find random regions that have at least 5 sites"
            if sentence in myfile.read():
                f = open(output.out,"w")
                f.write(sentence)
                f.close()
            else:
                score = params.score
                random_regions_df = pd.read_csv(input.random_regions,sep="\t",usecols=[score])
                list_random_mean_score = random_regions_df[score].values.tolist()
                intersection_df = pd.read_csv(input.intersect,sep="\t",usecols=['chrom', 'start', 'end', score])
                mean_score_region_interest = intersection_df[score].mean()
                #get pvalue from normal distribution of the random regions
                mean_normal_distrib = random_regions_df[score].mean()
                sd_normal_distrib = random_regions_df[score].std()
                observed_z_score = (mean_score_region_interest-mean_normal_distrib)/sd_normal_distrib
                p_values = scipy.stats.norm.sf(abs(observed_z_score))
                line_d = "P-val computed from the normal distribution of the random regions is : "+str(p_values)
                file1 = open(output.out,"w")
                file1.writelines(line_d)
                file1.close()