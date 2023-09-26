from collections import defaultdict
import os
workdir: os.getcwd()
include: os.getcwd() + '/config.sk'

"""
Get input ready
"""
pop_to_extract = POP_TO_EXTRACT
dir_tree = TREES
bed = open(sites_bed)

"""
create the dictionary of sites list per chr
"""
chr_to_snp = defaultdict(list)
CHR=[]
for line in bed:
    chr, snp = line.split()
    CHR.append(str(chr))
    chr_to_snp[str(chr)].append(str(snp))

"""
Number of time clues is run per sites
"""
def createList(test_number):
    return [str(item) for item in range(1, int(number_of_test)+1)]

TEST = createList(number_of_test)

"""
get the sites list in function of chr
"""
def getTargetFiles(chrom_list, pos_dict, test_list, POP):
    file_list=[]
    for CHR in chrom_list:
        for SNP in pos_dict[CHR]:
            for TEST in test_list:
                file = "clues_output/{pop}/test_{test}/chr{chr}/CLUES_output_{pop}_chr{chr}_{snp}.txt".format(pop=POP, test=TEST, chr=CHR, snp=SNP)
                file_list.append(file)
    return file_list

"""
actual snakemake steps
"""
rule all:
    input:
        expand("{filename}", filename=getTargetFiles(CHR,chr_to_snp, TEST, POP_TO_EXTRACT))

#extract branch
rule create_clues_input:
    input:
        TREES + "/extracted_trees/" + POP_TO_EXTRACT + "/" + POP_TO_EXTRACT + "_chr{chr}.mut.gz",
        TREES + "/extracted_trees/" + POP_TO_EXTRACT + "/" + POP_TO_EXTRACT + "_chr{chr}.anc.gz",
        coal=TREES + "/pop_size/" + POP_TO_EXTRACT + "/" + POP_TO_EXTRACT + "_popsize.coal"
    output:
        "clues_input/" + POP_TO_EXTRACT + "/test_{test}/chr{chr}/" + POP_TO_EXTRACT + "_chr{chr}_{snp}.timeb",
        multiext("clues_input/" + POP_TO_EXTRACT  + "/test_{test}/chr{chr}/" + POP_TO_EXTRACT + "_chr{chr}_{snp}",".mut.gz",".anc.gz",".dist.gz")
    params:
        input=TREES + "/extracted_trees/" + POP_TO_EXTRACT + "/" + POP_TO_EXTRACT + "_chr{chr}",
        output="clues_input/" + POP_TO_EXTRACT + "/test_{test}/chr{chr}/" + POP_TO_EXTRACT + "_chr{chr}_{snp}",
        mut="clues_input/" + POP_TO_EXTRACT + "/test_{test}/chr{chr}/" + POP_TO_EXTRACT + "_chr{chr}_{snp}.mut",
        anc="clues_input/" + POP_TO_EXTRACT +  "/test_{test}/chr{chr}/" + POP_TO_EXTRACT + "_chr{chr}_{snp}.anc",
        dist="clues_input/" + POP_TO_EXTRACT + "/test_{test}/chr{chr}/" + POP_TO_EXTRACT + "_chr{chr}_{snp}.dist",
        relate = Relate_Dir
    resources:
        mem_mb = 4000
    shell:
        """
        {params.relate}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
                        -i {params.input} \
                        -o {params.output}\
                        -m 1.25e-8 \
                        --coal {input.coal} \
                        --format b \
                        --num_samples 200 \
                        --first_bp {wildcards.snp} \
                        --last_bp {wildcards.snp}
gzip {params.mut}
gzip {params.anc}
gzip {params.dist}
        """

rule clues:
    input:
        timeb="clues_input/" + POP_TO_EXTRACT + "/test_{test}/chr{chr}/" + POP_TO_EXTRACT + "_chr{chr}_{snp}.timeb",
        coal=TREES + "/pop_size/" + POP_TO_EXTRACT + "/" + POP_TO_EXTRACT + "_popsize.coal"
    output:
        "clues_output/" + POP_TO_EXTRACT + "/test_{test}/chr{chr}/CLUES_output_" + POP_TO_EXTRACT + "_chr{chr}_{snp}.txt"
    params:
        input = "clues_input/" + POP_TO_EXTRACT + "/test_{test}/chr{chr}/" + POP_TO_EXTRACT + "_chr{chr}_{snp}",
        clues = Clues_Dir
    resources:
        mem_mb = 10000
    shell:
        """
        python3 {params.clues}/inference.py \
        --times {params.input} \
        --coal {input.coal} \
        --thin 5 \
        --tCutoff 1000 \
        --burnin 100  >> {output}
        """
