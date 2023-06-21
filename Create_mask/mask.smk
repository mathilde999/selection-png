import os
import pandas as pd
##input
##input info extraction
workdir: os.getcwd()
include: os.getcwd()+'/config.sk'

# print (modules)
shell.prefix('module load '+ modules)

#create the name of the regions using the windows bed
def create_region_pos(row):
    return("chr"+str(row["chr"])+":"+str(row["beg"])+"-"+str(row["end"]))

def create_region_name(row):
    return("chr"+str(row["chr"])+"_"+str(row["beg"])+"_"+str(row["end"]))

def return_region_pos(wildcards):
    return region_dict[wildcards.region_name]


#read the bed files with the windows coordinate
windows_df = pd.read_csv(windows_bed, sep="\t", header = None, names=["chr", "beg", "end"])

#create a list to stor the region pos
region_pos_list = windows_df.apply (lambda row: create_region_pos(row), axis=1).tolist()

#create a list to stor the region name
region_name_list = windows_df.apply (lambda row: create_region_name(row), axis=1).tolist()

#create a dictionary with name as key and pos as value
zip_iterator = zip(region_name_list, region_pos_list)
region_dict = dict(zip_iterator)


rule all:
    input:
        "mask_all_sort.bed.gz"

rule mpileup:
    input:
        bam = bam_list,
        ref =ref_file
    output:
        bed = temp("temp/out_mask_{region_name}.bed.gz")
    params:
        region = return_region_pos,
        min_depth = mean_coverage/2,
        max_depth = mean_coverage*2,
        dir_script = convert_vcf
    shell:
        """
        bcftools mpileup -q 20 -Q 20 -C 50 -r {params.region} -f {input.ref} -b {input.bam} |\
        bcftools call -c -V indels|bcftools view --max-alleles 2 -i 'MIN(DP)>{params.min_depth} & MAX(DP)<{params.max_depth}'|\
        vcftools --vcf - --max-missing 0.95 -c --recode --recode-INFO-all|\
        python {params.dir_script} {output.bed}
   """
rule concat_bed:
    input:
         expand("temp/out_mask_{region_name}.bed.gz", region_name=region_name_list)
    output:
        no_sort = temp("temp/mask_all_no_sort.bed")
    shell:
         """
         zcat {input} >> {output.no_sort}
         """
rule sort_bed:
    input:
         "temp/mask_all_no_sort.bed"
    output:
        gz = "mask_all_sort.bed.gz"
    params:
        bed = "mask_all_sort.bed"
    shell:
        """
        sort -k1,1 -k2,2n {input} > {params.bed}
        gzip {params.bed}
        """



