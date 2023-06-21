import os
##input
##input info extraction
workdir: os.getcwd()
include: os.getcwd()+'/config.sk'

chrall = [str(i) for i in range(1,23)]

# print (modules)
shell.prefix('module load '+ modules)

rule all:
    input:
        "phased/phased_autosomes.vcf.gz.tbi"

rule breakbychr:
    input:
        vcf_file
    output:
        vcf=temp("chr/chr{chr}.vcf.gz"),
        tabix=temp("chr/chr{chr}.vcf.gz.tbi")
    resources:
        mem_mb = 100000, #means 10000 per cpu
        threads = 10
    shell:
        """
        bcftools view -r chr{wildcards.chr} --threads {resources.threads} -O z -o {output.vcf} {input}
        tabix {output.vcf}
        """

rule chr_filter:
    input:
        vcf="chr/chr{chr}.vcf.gz",
        ID=kept_ID
    output:
        vcf="chr/filtered_chr{chr}.vcf.gz"
    params:
        vcftools=vcftoolsfilter,
        bed=positiv_mask
    resources:
        mem_mb = 100000
    shell:
        """
        vcftools --gzvcf {input.vcf} --keep  {input.ID} {params.vcftools} --bed {params.bed} --recode --stdout|bgzip -c > {output.vcf}
        """
rule tabix:
    input:
        vcf="chr/filtered_chr{chr}.vcf.gz"
    output:
        tabix= "chr/filtered_chr{chr}.vcf.gz.tbi"
    resources:
        mem_mb = 100000
    shell:
        """
        tabix {input.vcf}
        """
rule format_map:
    input:
        map_file
    output:
        temp("maps/genetic_map_hg38_chr{chr}.txt.gz")
    resources:
        mem_mb = 4000
    shell:
        """
        zcat {input}|cut -d " " -f2,3,4|bgzip -c >{output}
        """

rule phasing:
    input:
        vcf="chr/filtered_chr{chr}.vcf.gz",
	tabix= "chr/filtered_chr{chr}.vcf.gz.tbi",
        map="maps/genetic_map_hg38_chr{chr}.txt.gz"
    output:
        vcf="phased/chr{chr}_phased.vcf.gz",
        tabix="phased/chr{chr}_phased.vcf.gz.tbi"
    resources:
        mem_mb = 100000, #means 10000 per cpu
        threads = 10
    shell:
        """
        shapeit4 --input {input.vcf} --map {input.map} --output {output.vcf} --region chr{wildcards.chr} --thread {resources.threads} --sequencing
        tabix {output.vcf}
        """

rule concat:
    input:
        expand("phased/chr{chr}_phased.vcf.gz",chr=chrall)
    output:
        vcf="phased/phased_autosomes.vcf.gz",
        tabix="phased/phased_autosomes.vcf.gz.tbi"
    resources:
        mem_mb = 100000, #means 10000 per cpu
        threads = 10
    shell:
        """
        bcftools concat --threads {resources.threads} -o {output.vcf} -O z {input}
        tabix {output.vcf}
        """


