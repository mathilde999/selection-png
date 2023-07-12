# selection-png
Scripts developed and used in [André et al., biorxiv 2022](https://doi.org/10.1101/2022.12.15.520226)

### Kinship analysis
#### Convert vcf to plink bed format
We used [plink v1.9](https://www.cog-genomics.org/plink/) ([Purcell et al., 2007](https://doi.org/10.1086/519795))
`plink --vcf all_pop_biallelic_PASS.vcf.gz --make-bed --double-id --recode --out all_pop_biallelic_PASS`
#### Relationship inference
We used [KING v2.2.4](https://www.kingrelatedness.com/manual.shtml) integrated relationship inference 
([Manichaikul et al.,2010](https://academic.oup.com/bioinformatics/article/26/22/2867/228512))to detect first and 
second-degree relatives pairs among our dataset\
\
`king -b all_pop_biallelic_PASS.bed --related --degree 2 --prefix all_relatedness`

### creation of a genomic mask 
#### create coverage mask 
We included in the positive coverage mask sites with a minimum base quality of 20, an alignment minimum mapping quality 
of 20 and downgrading mapping quality for reads containing excessive mismatches with a coefficient of 50. From these
sites, we excluded indels, sites with more than two alleles, sites with a maximum missing rate of 5%. We also masked
sites whose depth of coverage summed across all samples was higher or lower than the sum of the average depth across the
dataset by a factor of 2-fold.\
\
`snakemake -s mask.smk --profile [profile_file]`
\
We are using a [snakemake v7.24.0](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) workflow system 
in several of our pipeline for this step and other analysis 

#### create positive mask with coverage mask and mappability mask 
We create a positive mask using [bedtools v2.29.2](https://bedtools.readthedocs.io/en/latest/index.html) including the 
site present in both the coverage mask we created (mask_all_sort.bed) and the [mappability mask](https://share.eva.mpg.de/index.php/s/ygfMbzwxneoTPZj) use with [MSMC](https://github.com/stschiff/msmc-tools/tree/master)  liftover to GR38 \
\
`bedtools intersect -a mask_all_sort.bed -b gr38.mask.bed > mask_coverage_and_map.bed`
#### filtering out the variant without PASS flag from the variant calling
Keep only the sites that pass the PASS filter from the variant calling unfiltered vcf files.\
\
`vcftools --gzvcf raw_chr"$0".vcf.gz --remove-filtered-all --recode --stdout |bgzip > sites_PASS_chr"$0".vcf.gz`\
Output the site without PASS flag from the variant calling unfiltered vcf files.\
\
`bcftools isec -C -o sites_not_PASS_chr"$0".bed -Oz raw_chr"$0".vcf.gz sites_PASS_chr"$0".vcf.gz`\
\
Concat the not PASS sites\
\
`seq 22|awk '{print "zcat sites_not_PASS_chr"$0".bed.gz|cut -f1,2" >> sites_not_PASS_all.bed}' > concat.sh`\
\
`./concat.sh`\
\
Format of the output file\
\
`awk '{print $1 "\t" ($2 - 1) "\t" $2}' sites_not_PASS_all.bed|bedtools merge -i - > sites_not_PASS_all_good_format.bed`
#### remove the not PASS sites from the global positive mask 
`bedtools subtract -a mask_coverage_and_map.bed -b sites_not_PASS_all_good_format.bed > full_positive_mask.bed`
#### create the mask fasta file 
split the positive mask in chr \
\
`for chr in 'bedextract --list-chr full_positive_mask.bed'; do bedextract $chr full_positive_mask.bed > full_positive_mask_$chr.bed; done`
\
Generating the postive mask in a fasta format \
\
`./make_fa_mask.sh`

### LD pruning  
#### convert vcf file to ped and map file
We kept biallelic SNPs with a maximal missing rate of 5% filtered for our genomic mask and with a minor allele frequency > 5% . We excluded related individuals and two PNG samples with low call rate\
`vcftools --gzvcf combined.vcf.gz --max-alleles 2 --min-alleles 2 --remove-indels --max-missing 0.95 --remove-filtered-all --bed full_positive_mask.bed --maf 0.05 --keep all_unrelated.ID --plink --out all_biSNP_PASS_maf_unrelated`
#### LD pruning
We pruned variants in high linkage disequilibrium, we used [plink v1.9](https://www.cog-genomics.org/plink/) ([Purcell et al., 2007](https://doi.org/10.1086/519795))
with the default parameters of 50 variants count window shifting from five variants and a variance inflation factor (VIF) 
threshold of 2\
\
`plink --file all_biSNP_PASS_maf_unrelated --indep 50 5 2 --out all_biSNP_PASS_maf_unrelated_VIFLDPruned --noweb`
`plink --file all_biSNP_PASS_maf_unrelated_VIFLDPruned --extract all_biSNP_PASS_maf_unrelated_VIFLDPruned.prune.in  --out all_biSNP_PASS_maf_unrelated_VIFLDPruned --noweb --recode`

### PCA
We performed Principal Component Analysis (PCA) using the smartpca program from the [EIGENSOFT v.7.2.0](https://github.com/DReichLab/EIG/tree/master) package ([Patterson et al., 2006](https://doi.org/10.1371/journal.pgen.0020190)) 
on LD pruned SNPs.\
\
`smartpca -p par.file`

### Admixture
We ran [ADMIXTURE v1.3](https://dalexander.github.io/admixture/index.html) ([Alexander et al., 2009](10.1101/gr.094052.109)) 
for the components K=2 to K=10 on the same dataset as the PCA. For each component, ADMIXTURE computes the cross-validation 
error using k-fold cross-validation procedure. We set the k parameter to 100. In order to select the model with the most 
likely number of components, we generated the cross-validation error 10 times for each component.\
\
`for i in {1..10};do ./admixture.sh ${i}; done`


### phasing
Before phasing we filtered our vcf file. We kept biallelic SNPs with a maximal missing rate of 5% filtered for our genomic mask. We excluded related individuals
and two PNG samples with low call rate. We phased our sequences statistically without reference with the software [shapeit4 v4.2.2](https://odelaneau.github.io/shapeit4/) ([Delaneau, et al., 2019](https://www.nature.com/articles/s41467-019-13225-y)). 
We used the [genetic map](https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/) provided by [Eagle v2](https://alkesgroup.broadinstitute.org/Eagle/#x1-30002) ([Loh et al., 2016](https://www.nature.com/articles/ng.3679)). 
This phasing was performed per chromosome using the sequencing option that adjusts the default parameter to what the shapeit
authors consider the most adapted to sequenced data.\
\
`snakemake -s shapeit4.smk --profile [profile_file]`
\
That snakemake pipeline will create two directories in your working directory. Directory chr/ included unphased vcf files
per chr filtered for the filters indicated in the config.sk.The other directory created is phased/ that will include all the phased vcf files per chr
### XP-EHH scan 
To study signals of positive selection in the two target populations, we used [selscan V2.0](https://github.com/szpiech/selscan) ([Sabeti, et al., 2007](https://www.nature.com/articles/nature06250))
to compute the cross-extended haplotype homozygosity (XP-EHH) on the phased dataset.
#### create map files
`python Selscan/CreatingMapFiles.py --files_path chr/ --vcf`
#### run selscan
`selscan --xpehh --vcf High_chr{chr}.vcf.gz --vcf-ref Low_chr{chr}.vcf.gz --map High_Low_chr{chr}_GeneticMap.map --out High_chr{chr} --threads 4`
#### normalization
`norm --xpehh --files *.xpehh.out`
### top region score (target pop)
We kept the SNPs whose XPEHH score was in the 99th percentiles.\
\
`python XPEHH_bed_99.py` 
### merge top SNPs in genomic regions
We then merge top SNPs together when they were closer than 10kb with the closest windows  using [bedtools v2.29.2](https://bedtools.readthedocs.io/en/latest/index.html).
The XP-EHH score of the regions will be the top XP-EHH of the regions uses to create the bigger regions \
\
`sort -k1,1 -k2,2n top_snp_XPEHH.bed > sort_top_snp_XPEHH.bed`
`bedtools merge -i sort_top_snp_XPEHH.bed -d 10000 -c 4 -o max > merged_top_snp_XPEHH.bed`
####  get p-value for the 10 top regions with random sampling approach 
`snakemake -s p_val.smk -kp --jobs 100 --profile [profile_file]--groups intersection=group1 --group-component group1=10`
### PBS scan
#### Running PBS
`python RunPBS.py PBS.info`
#### PBS score for sliding windows
We create sliding windows of X snp with Y snps step and give to each of these regions a PBS that is the average 
of the PBS of the X SNP that the window is made. We keep only the windows whose PBS score in the 99th percentile: output 
"sliding_windows_snp_top_regions_PBS.bed". We make a file "mean_pos_PBS.bed" with the mean pos of the region (mean of 
the start and end bp) and the PBS of the region (average of the PBS of all the snp use to create the windows) 
We used windows of 20 SNPs sliding by 5 SNPs.\
\
`python sliding_windows_snp.py -b High_Low_EUR_AFR.pbs -w 20 -s 5`\
\
We then merge top regions together when they were closer than 10kb with the closest windows  using [bedtools v2.29.2](https://bedtools.readthedocs.io/en/latest/index.html).
The PBS score of the regions will be the top PBS of the regions uses to create the bigger regions \
\
`bedtools merge -i sliding_windows_snp_top_regions_PBS.bed -d 10000 -c 4 -o max > windows_and_PBS_top_windows_final_merged.bed`
####  get p-value for the 10 top regions with random sampling approach 
`snakemake -s p_val.smk -kp --jobs 100 --profile [profile_file]--groups intersection=group1 --group-component group1=10`
### Fisher score
We are computing Fisher for the same windows used to generate PBS. We then need to generate XP-EHH score for these windows too.
Make XP-EHH output as a bed file: \
\
`cat xpehh_all.out.norm|awk '{print $1 "\t" $2 "\t" $2 "\t" $4}'|tail -n+2 >xpehh.bed`\
\
Extract the PBS windows coordinates: \
\
`cut -f1,2,3 windows_20SNP_step5.bed > windows.bed`\
\
make a file that list each snp associated to each windows and their xpehh score: \
\
`bedtools intersect -a windows.bed -b xpehh.bed -wa -wb > windows_xpehh.bed`\
\
average the XPEHH score of all the snp in each region (top score in the region): \
\
`bedtools groupby -i windows_xpehh.bed -g 1-3 -c 7,7  -o max,count >windows_XPEHH_final.bed`\
\
keep only CHR  START   END     NORM + add header: \
\
`cut -f1,2,3,4 windows_XPEHH_final.bed > XPEHH_final.bed`\
\
create bed file with fisher score for each regions: \
keep only regions in the 99 percentile for the Fisher score: \
\
`python Fish.py --xpehh XPEHH_final.bed --pbs windows_20SNP_step5.bed --out FisherScore_windows`\
\
merge the regions together, score of the regions is the top score. Extract the top 10 regions from it: \
\
`bedtools merge -i FisherScore_windows.topSNP.bed -d 10000 -c 4 -o max > top_regions_Fisher_merged.bed`
#### get p-value for the 10 top regions with random sampling approach 
`snakemake -s p_val.smk -kp --jobs 100 --profile [profile_file]--groups intersection=group1 --group-component group1=10`
### Relate analysis
We used [Relate v1.8](https://myersgroup.github.io/relate/) ([Speidel et al., 2019](https://www.nature.com/articles/s41588-019-0484-x)) \
The relate snakemake pipeline is ran with the following command \
\
`snakemake -s Make_trees.smk -kp --profile [profile_file] `\
\
and will create several folders in the working directory: 
#### Create relate input file
The snakemake pipeline will create an input_relate folder including the *sample, *haps, *annot and *dist file used by Relate.
It will also includes the file poplabel.txt that is the population files correctly ordered for the further relate analysis 
#### Make trees
Will create recombination trees for all the samples included in the vcf files you have inputed.We used the default 
parameters of 1.25e-8 for the mutation rate and 3000 for the effective population size of haplotypes. The output are *anc and *mut files in the trees_all/ directory \
#### Extracting trees
Trees will be extracted for every population you included in config.sk. The extracted trees (*anc and *mut files) are 
found in the extracted_trees/ directory 
#### Effective population size
Effective population size will be estimated for each of the populations you included in the config.sk file. We used the 
default recombination rate of 1.25e-8 and we specified the first time interval from 0 to 10² years ago and 
then the periods between 10² and 10<sup>7</sup>  is split into successive bins with a step of 10<sup>0.1</sup>  years. We used Relate default 
parameters of 28 years per generation, the number of cycles and the fractions of trees to be dropped. Outputs (*coal) are found in the pop_size/ folder 

### Clues analysis 
We used [clues v1.0](https://github.com/standard-aaron/clues/) ([Stern et al., 2019)](https://journals.plos.org/plosgenetics/article/metrics?id=10.1371/journal.pgen.1008384)) \
`snakemake -s Clues.smk --profile [profile_file] --groups create_clues_input=group1 clues=group2 --group-components group1=20 group2=10`
#### Extract local trees with Relate
The first rule of the pipeline extract the local trees for each SNPs in the list of interest for the population indicated 
in the config.sk file. The Relate rule of the snakemake pipeline will sampled the branch length for those trees. The branch 
length of the focal SNP tree is resampled 200 times (“--num_samples 200” option of the “SampleBranchLengths.sh” script ), and each iteration result
is recorded to estimate the uncertainty in the age estimate of each node. To take the uncertainty in the branch length 
estimate into account, this can performed this several times (number_of_test in the config.sk file). For each SNP, we thus 
compute sampled branches X times. They output files *timeb can be found in the clues_input directory.
#### Running clues for each focal trees
Clues is run for each of the sampled branches output are found in the clues_output directory.

### Association genotype phenotype
We used [GEMMA v0.98.4](https://github.com/genetics-statistics/GEMMA/tree/master) ([Zhou & Stephens, 2012](https://www.nature.com/articles/ng.2310))
#### convert vcf files to bed files
Using [plink v1.9](https://www.cog-genomics.org/plink/).These files will be used to computed the relatedness matrices\
\
`plink --vcf unrelated_filtered_masked.vcf.gz --keep papuans.ID --output-missing-genotype 0 --double-id --make-bed --out unrelated_masked_with_pheno`
#### Extract the candidate SNPs and convert to bed file format.
These file will be used to compute association only for the top SNP\
\
`plink --vcf unrelated_filtered_masked.vcf.gz --keep papuans.ID  --extract range TopSNP.bed --output-missing-genotype 0 --double-id --set-missing-var-ids @:# --make-bed --out Papuan_unrelated_masked_with_pheno_Top_SNP`
#### add the pheno info to the fam files from the top SNP and the all site file
`python create_fam.py Papuan_unrelated_masked_with_pheno.fam residuals_pheno.cov`
`python create_fam.py Papuan_unrelated_masked_with_pheno_Top_SNP.fam residuals_pheno.cov`
#### compute the relatedness matrices
`gemma-0.98.4-linux-static-AMD64 -bfile Papuan_unrelated_masked_with_pheno -gk 1 -o Papuan_with_pheno_centered`
### gemma association
To run for every phenotype column\
\
`gemma-0.98.4-linux-static-AMD64 -bfile Papuan_unrelated_masked_with_pheno_Top_SNP -k output/Papuan_with_pheno_centered.cXX.txt -notsnp -miss 1 -lmm 4 -n [pheno_column] -outdir results_centered/ -o [pheno_name] `

### Enrichment for blood phenotypes
`snakemake -s enrich.smk --profile [profile_file]`
