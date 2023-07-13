import os
##input
##input info extraction
workdir: os.getcwd()
include: os.getcwd() + '/config.sk'
import pybedtools
import glob
import pandas as pd
from statsmodels.distributions.empirical_distribution import ECDF


def get_chr_lenght(chrom):
    return chr_dic[chrom]

def random_sampling(windows_df, number_to_select, blood_site_bed):
    #select X random windows in the windows that have at least one site associated to a phenotype
    sample_windows_df = windows_df.sample(n=number_to_select,replace=False)
    #convert to bed
    sample_windows_bed = pybedtools.BedTool.from_dataframe(sample_windows_df)
    #for each windows output the number of site that are significantly associated to a blodd pheno
    intersection = sample_windows_bed.intersect(blood_site_bed,c=True)
    intersection_df = intersection.to_dataframe()
    intersection_df.columns = ['chrom', 'start', 'end', 'count_sign_sites']
    #we don't actually care about the exact number of sites associated with blood pheno but only if there are at least one
    intersection_df['blood_assoc'] = (intersection_df['count_sign_sites'] >= 1).astype(int)
    #we return the number of windows with at least one site associated with blood
    number_site = intersection_df['blood_assoc'].sum()
    return number_site


"""
Get input ready
"""
bed = open(autosomes_file)
####create the dict with all the info respectiv to the chromosomes
chr_dic = {}
CHR = []
for line in bed:
    chrom, start, end = line.split()
    CHR.append(chrom)
    chr_dic[chrom] = end

dict_test = {}
list_info_regions = []

for i in range(len(n_region_sign_list)):
    n_region_sign_all = n_region_sign_list[i]
    n_region_sign_blood = n_region_sign_blood_list[i]
    info_list = [n_region_sign_blood, n_region_sign_all]
    key = str(n_region_sign_blood) + "_out_of_" + str(n_region_sign_all)
    dict_test[key]=info_list
    list_info_regions.append(key)



"""
We build the distribution of the count of how many windows include at least one significant association with a
blood phenotype in the UKbiobank by X random non-sliding windows of 200bp (100bp upstream and 100bp downstream) including
at least 1 site with a significant asssociation with any phenotype in the UKbiobank.
"""


rule all:
    input:
         expand("out/{info_region}_random_regions.out", info_region = list_info_regions)

rule split_chr:
    input:
        pheno_all=bed_pheno,
    output:
        pheno_all_chr= "pheno_pval_chr/pheno_pval_chr{chr}.tsv.bgz"
    params:
        end_chr = lambda wildcards: get_chr_lenght(wildcards.chr)
    resources:
        mem_mb = 100000
    shell:
        """
        tabix -h {input.pheno_all} {wildcards.chr}:0-{params.end_chr}|bgzip > {output.pheno_all_chr}
        """

rule clean_pheno_file:
    input:
        pheno_all = "pheno_pval_chr/pheno_pval_chr{chr}.tsv.bgz"
    output:
        pheno_all_sign= "pheno_sum/pheno_sign_sum_chr{chr}.bed.gz"
    resources:
        mem_mb = 100000
    run:
        # Read column names from file
        cols = list(pd.read_csv(input.pheno_all,nrows=1,sep='\t', compression='gzip'))
        all_pheno_list = cols[4:]

        #read and filter the file, it is huge so we will do it in chunks
        #what is the maximum size of chunck we can do at a time
        pheno_sample = pd.read_csv(input.pheno_all,nrows=10,sep='\t', compression='gzip')
        pheno_sample_size = pheno_sample.memory_usage(index=True).sum()#how much RAM for 10 rows, in byte
        # define a chunksize that would occupy a maximum of resources.mem_mb
        # we divide by 10 because we have selected 10 lines in our df_sample
        # we then get the integer part of the result
        max_ram_byte = (resources.mem_mb * 1000000) - 1000000000
        my_chunk = (max_ram_byte/ pheno_sample_size) / 10
        my_chunk = int(my_chunk // 1)
        # create the iterator
        iter_csv = pd.read_csv(input.pheno_all,iterator = True, chunksize = my_chunk, sep='\t',
            usecols =[i for i in cols if i != "REF" and i != "ALT"], compression='gzip')
        print("creations list of df done")
        #we apply filters on every chunk
        list_filtered_pheno_chunck=[]
        for chunk in iter_csv:
            #adding the header for each chunk so I can use dropna
            chunk.columns = ["#CHROM","POS"] + all_pheno_list
            #remove the site for which we only have nan at every pheno
            chunk = chunk.dropna(subset=all_pheno_list)
            #replace pval < -11.29 to 1 and the pval > -11.29 to 0
            chunk[all_pheno_list] = (chunk[all_pheno_list] < -11.29).astype(int)
            #I will keep only the site that are at least significantly associated to 1 pheno
            chunk = chunk[(chunk[all_pheno_list] == 1).any(axis=1)]
            list_filtered_pheno_chunck.append(chunk)

        #concat the chunks
        pheno_filtered = pd.concat(list_filtered_pheno_chunck)

        #make it bed format
        pheno_filtered["start"] = pheno_filtered["POS"] -1
        pheno_filtered["end"] = pheno_filtered["POS"]
        pheno_filtered["CHR"] = pheno_filtered["#CHROM"]
        col_to_use = ["CHR", "start", "end"] + all_pheno_list
        pheno_filtered = pheno_filtered[col_to_use]
        pheno_filtered.to_csv(output.pheno_all_sign,sep="\t",index=False, compression='gzip')


"""
Define the windows that have at least one site that is significantly associated to any pheno 
"""
rule create_non_overlapping_windows:
    input:
        pheno_sign = "pheno_sum/pheno_sign_sum_chr{chr}.bed.gz" #bed file with info for site that are associated significiantlyt to at least one pheno
    params:
        wind_size=size_windows_bp
    output:
        windows_sign = "pheno_sum/windows_"+str(size_windows_bp)+"_sign_any_pheno_chr{chr}.bed.gz" #windows of 2000bp that have at least one site that is associated to at least 1 pheno
    resources:
        mem_mb = 100000
    run:
        #create a bed file with the chromosome lenght info
        chr = wildcards.chr
        start = 0
        end = chr_dic[chr]
        region = " ".join([str(chr), str(start), str(end)])
        autosomes_bed = pybedtools.BedTool(region,from_string=True)
        #create windows of size X from the config file
        wind_size = params.wind_size
        windows_bed = autosomes_bed.window_maker(autosomes_bed,w=wind_size)
        #load pheno bed
        pheno_df = pd.read_csv(input.pheno_sign,sep='\t')
        pheno_bed = pybedtools.BedTool.from_dataframe(pheno_df)
        #for each windows output the number of site that have at least one selection score
        intersection = windows_bed.intersect(pheno_bed,c=True)
        intersection_df = intersection.to_dataframe()
        intersection_df.columns = ['chrom', 'start', 'end', 'count_sign_sites']
        #keep only the windows that minimum one site that is significantly associated to a pheno
        intersection_df = intersection_df[intersection_df['count_sign_sites'] >= 1]
        intersection_df.to_csv(output.windows_sign,sep="\t",index=False,compression='gzip')

"""
Get a bed of the site that have a significant association with at least one blood pheno
"""
rule blood_sign_site:
    input:
        pheno_sign = "pheno_sum/pheno_sign_sum_chr{chr}.bed.gz", #bed file with info for site that are significantly associated with at least one pheno
        blood_pheno = list_blood_pheno
    output:
        blood_sign = "pheno_sum/blood_sign_sum_chr{chr}.bed.gz" #bed file with info for site that are associated significantly associated with at least one blood pheno
    resources:
        mem_mb = 100000
    run:
        #from file with blood pheno name to list of the blood pheno name
        blood_pheno_list=[]
        blood_file = open(input.blood_pheno)
        for pheno in blood_file:
            pheno =pheno.replace("\n","")
            blood_pheno_list.append(pheno)
        col_to_use= ["CHR", "start", "end"] + blood_pheno_list
        pheno_blood_bed = pd.read_csv(input.pheno_sign, compression='gzip', sep="\t", usecols=col_to_use)
        pheno_blood_bed = pheno_blood_bed[(pheno_blood_bed[blood_pheno_list] == 1).any(axis=1)]
        pheno_bed = pheno_blood_bed[["CHR", "start", "end"]]
        pheno_bed.to_csv(output.blood_sign,sep="\t", index=False,compression='gzip')

rule concat_blood_sites:
    input:
        blood_sign_site = expand("pheno_sum/blood_sign_sum_chr{chr}.bed.gz", chr=CHR)
    params:
        blood_sign_site = "pheno_sum/blood_sign_sum_chr"
    output:
        out="pheno_sum/blood_sign_sum_all_chr.bed.gz"
    resources:
        mem_mb=100000
    run:
        blood_sites_all_list = [i for i in glob.glob(f"{params.blood_sign_site}*{'.bed.gz'}")]
        blood_sites_all = pd.concat([pd.read_csv(f,delimiter='\t') for f in blood_sites_all_list])
        blood_sites_all = blood_sites_all.sort_values(by=['CHR', 'start'])
        blood_sites_all.to_csv(output.out,sep="\t",index=False,compression={'method': 'gzip'})

rule concat_random_windows:
    input:
        windows_sign = expand("pheno_sum/windows_"+str(size_windows_bp)+"_sign_any_pheno_chr{chr}.bed.gz", chr=CHR)
    params:
        windows_sign = "pheno_sum/windows_"+str(size_windows_bp)+"_sign_any_pheno_chr"
    output:
        out= "pheno_sum/windows_"+str(size_windows_bp)+"_sign_any_pheno_all_chr.bed.gz"
    resources:
        mem_mb=100000
    run:
        windows_sign_all_list = [i for i in glob.glob(f"{params.windows_sign}*{'.bed.gz'}")]
        windows_sign_all = pd.concat([pd.read_csv(f,delimiter='\t') for f in windows_sign_all_list])
        windows_sign_all = windows_sign_all.sort_values(by=['chrom', 'start'])
        windows_sign_all.to_csv(output.out,sep="\t",index=False,compression={'method': 'gzip'})

rule get_pval:
    input:
        windows_sign = "pheno_sum/windows_"+str(size_windows_bp)+"_sign_any_pheno_all_chr.bed.gz",
        blood_sign_site="pheno_sum/blood_sign_sum_all_chr.bed.gz"
    params:
        number_random_sampling= n_resampling
    output:
        out= "out/{info_region}_random_regions.out"
    resources:
        mem_mb=100000
    run:
        number_random_windows = dict_test[wildcards.info_region][1]
        number_assoc_blood = dict_test[wildcards.info_region][0]
        windows=pd.read_csv(input.windows_sign, sep="\t", usecols=["chrom", "start", "end"])
        blood_df = pd.read_csv(input.blood_sign_site,sep="\t",usecols=["CHR", "start", "end"])
        blood_bed = pybedtools.BedTool.from_dataframe(blood_df)

        list_sum=[]
        for i in range(params.number_random_sampling):
            random_sum=random_sampling(windows,number_random_windows, blood_bed)
            list_sum.append(random_sum)

        # Get p-value of the nuber of blood assoc based on the random (null) distribution
        ecdf_score = ECDF(list_sum)
        pval_score_region = ecdf_score(number_assoc_blood)
        line_a = "Get p-value of the observed site sign assoc with blood (" + str(number_assoc_blood) + " for " +str(number_random_windows)+\
                 " regions with significant assoc based on the random (null) distribution"
        line_b = 'P(random_distrib < '+ str(number_assoc_blood) +' regions with at least one site sign assoc to blood pheno) = ' + str(pval_score_region) + "\n"
        lines = line_a + line_b

        file1 = open(output.out,"w")
        file1.writelines(lines)
        file1.close()

