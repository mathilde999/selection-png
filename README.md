# selection-png  
Scripts developed and used in [André et al., biorxiv 2022](https://doi.org/10.1101/2022.12.15.520226)\  
André, M., Brucato, N., Hudjasov, G., Pankratov, V., Yermakovich, D., …, Mondal, M., Ricaut, F.-X. (2022, December 15). Positive selection in the genomes of two Papua New Guinean populations at distinct altitude levels (p. 2022.12.15.520226). bioRxiv. https://doi.org/10.1101/2022.12.15.520226  
  
## 1. Creation of a genomic mask   
#### 1.1. Create coverage mask   
We included in the positive coverage mask sites with a minimum base quality of 20, an alignment minimum mapping quality   
of 20 and downgrading mapping quality for reads containing excessive mismatches with a coefficient of 50. From these  
sites, we excluded indels, sites with more than two alleles, sites with a maximum missing rate of 5%. We also masked  
sites whose depth of coverage summed across all samples was higher or lower than the sum of the average depth across the  
dataset by a factor of 2-fold.\  
\  
```shell
snakemake -s Create_mask/mask.smk --profile [profile_file]
```
#### 1.2. Create positive mask with coverage mask and mappability mask   
We create a positive mask using [bedtools v2.29.2](https://bedtools.readthedocs.io/en/latest/index.html) including the   
site present in both the coverage mask we created `mask_all_sort.bed` and `gr38.mask.bed`, the [mappability mask](https://share.eva.mpg.de/index.php/s/ygfMbzwxneoTPZj) used with [MSMC](https://github.com/stschiff/msmc-tools/tree/master)  liftover to GR38 \  
\  
```shell
bedtools intersect -a mask_all_sort.bed -b gr38.mask.bed > mask_coverage_and_map.bed
```
#### 1.3. Filtering out the variant without PASS flag from the variant calling  
Keep only the sites that pass the PASS filter from the variant calling unfiltered vcf files `raw_chr"$0".vcf.gz` \  
\  
`vcftools --gzvcf raw_chr"$0".vcf.gz --remove-filtered-all --recode --stdout |bgzip > sites_PASS_chr"$0".vcf.gz`\  
\  
Output the site without PASS flag from the variant calling unfiltered vcf files `raw_chr"$0".vcf.gz` \  
\  
`bcftools isec -C -o sites_not_PASS_chr"$0".bed -Oz raw_chr"$0".vcf.gz sites_PASS_chr"$0".vcf.gz`\  
\  
Concat bed files with the coordinate of the not PASS sites `sites_not_PASS_chr"$0".bed.gz`\  
\  
`seq 22|awk '{print "zcat sites_not_PASS_chr"$0".bed.gz|cut -f1,2" >> sites_not_PASS_all.bed}' > concat.sh`\  
\  
`./concat.sh`\  
Format of the output file `sites_not_PASS_all.bed`\  
\  
`awk '{print $1 "\t" ($2 - 1) "\t" $2}' sites_not_PASS_all.bed|bedtools merge -i - > sites_not_PASS_all_good_format.bed`  
#### 1.4. Remove the not PASS sites from the global positive mask  
`bedtools subtract -a mask_coverage_and_map.bed -b sites_not_PASS_all_good_format.bed > full_positive_mask.bed`  
#### 1.5. Create the mask fasta file   
split the positive mask `full_positive_mask.bed` per chr \  
\  
`for chr in 'bedextract --list-chr full_positive_mask.bed'; do bedextract $chr full_positive_mask.bed > full_positive_mask_$chr.bed; done`\  
\  
Convert the positive mask to fasta format using the [ancestral genome](https://ftp.ensembl.org/pub/release-93/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38/)   
from the Ensembl 93 release `homo_sapiens_ancestor_${chr}.fa.gz` \  
\  
`Create_mask/make_fa_mask.sh`  
## 2. Phasing  
We used [shapeit4 v4.2.2](https://odelaneau.github.io/shapeit4/) \  
\  
`snakemake -s Phasing/shapeit4.smk --profile [profile_file]`  
\  
That snakemake pipeline will create two directories in your working directory. Directory `chr/` includes unphased vcf files  
per chr filtered for the filters indicated in `Phasing/config.sk`.The other directory created is `phased/` that includes all the phased vcf files per chr  
## 3. XP-EHH scan   
We used [selscan V2.0](https://github.com/szpiech/selscan) with highlanders as the target population and lowlanders as the reference population.\  
Selscan can use phased vcf files separated by chromosome with no missing data. It will also need one population per file.\  
It will also need map files of those vcf files.  
#### 3.1. Create genetic map files for selscan from vcf files  
`python Selscan/CreatingMapFiles.py --files_path chr/ --vcf`  
#### 3.2. Run selscan  
`selscan --xpehh --vcf High_chr{chr}.vcf.gz --vcf-ref Low_chr{chr}.vcf.gz --map High_Low_chr{chr}_GeneticMap.map --out High_chr{chr} --threads 4`  
#### 3.3. Normalization  
`norm --xpehh --files *.xpehh.out`  
#### 3.4. Top region score  
Input is the concatenated files `xpehh_all.out.norm`   
This script will output the SNP with xpehhnorm in the 99th percentile `top_snp_XPEHH.bed`\  
\  
`python Selscan\XPEHH_bed_99.py`   
\  
If you are interested in the top SNP for the reference population (lowlanders in our case)\  
This script will output the SNP with xpehhnorm in the 1th percentile `top_snp_XPEHH.bed`\  
`python Selscan\XPEHH_bed_1.py`   
  
#### 3.5. Merge top SNPs in genomic regions  
We then merge top SNPs together when they were closer than 10kb with the closest windows  using [bedtools v2.29.2](https://bedtools.readthedocs.io/en/latest/index.html).  
The XP-EHH score of the regions will be the top XP-EHH of the top SNPs uses to create the bigger regions \  
\  
`sort -k1,1 -k2,2n top_snp_XPEHH.bed > sort_top_snp_XPEHH.bed`  
`bedtools merge -i sort_top_snp_XPEHH.bed -d 10000 -c 4 -o max > merged_top_snp_XPEHH.bed`  
#### 3.6. Get p-value for the 10 top regions with random sampling approach   
`snakemake -s random_sampling/p_val.smk -kp --jobs 100 --profile [profile_file]`  
  
## 4. PBS scan  
#### 4.1. Running PBS  
`python PBS/RunPBS.py --popfile HL_LL_YRI.pop --pbspop Mt_Wilhelm:Daru:YRI --dropna filtered_{chr}.vcf.gz`  
#### 4.2. PBS score for sliding windows  
We create sliding windows of X snp with Y snps step and give to each of these regions a PBS that is the average   
of the PBS of the X SNP that the window is made with `windows_20SNP_step5.bed`. We keep only the windows whose PBS score in the 99th percentile:   
`sliding_windows_snp_top_regions_PBS.bed`. It also output `mean_pos_PBS.bed` with the mean pos of the region (mean of   
the start and end bp) and the PBS of the region (average of the PBS of all the snp use to create the windows)   
We used windows of 20 SNPs sliding by 5 SNPs \  
\  
`python PBS/sliding_windows_snp.py -b High_Low_EUR_AFR.pbs -w 20 -s 5`  
\  
We then merge top windows `sliding_windows_snp_top_regions_PBS.bed` together when they were closer than 10kb with the closest windows  using [bedtools v2.29.2](https://bedtools.readthedocs.io/en/latest/index.html).  
The PBS score of the regions will be the top PBS of the regions uses to create the bigger regions \  
\  
`bedtools merge -i sliding_windows_snp_top_regions_PBS.bed -d 10000 -c 4 -o max > windows_and_PBS_top_windows_final_merged.bed`  
#### 4.3. Get p-value for the 10 top regions with random sampling approach   
`snakemake -s random_sampling/p_val.smk -kp --jobs 100 --profile [profile_file]`  
## 5. Fisher score  
We computed Fisher score  for the same windows used to generate PBS `windows_20SNP_step5.bed`.   
  
#### 5.1. Generate XP-EHH score for the sliding windows   
Make XP-EHH output as a bed file: \  
\  
`cat xpehh_all.out.norm|awk '{print $1 "\t" $2 "\t" $2 "\t" $4}'|tail -n+2 >xpehh.bed`\  
Extract the PBS windows coordinates: \  
\  
`cut -f1,2,3 windows_20SNP_step5.bed > windows.bed`\  
make a file that list each snp associated to each windows and their xpehh score: \  
\  
```shell
bedtools intersect -a windows.bed -b xpehh.bed -wa -wb > windows_xpehh.bed
```
average the XPEHH score of all the snp in each region (top score in the region): \  
\  
```shell
bedtools groupby -i windows_xpehh.bed -g 1-3 -c 7,7  -o max,count >windows_XPEHH_final.bed
```  
keep only CHR  START   END     NORM + add header: \  
\  
```shell
cut -f1,2,3,4 windows_XPEHH_final.bed > XPEHH_final.bed
```\  
#### 5.2. Create bed file with fisher score for each sliding windows: \  
keep only sliding windodws in the 99 percentile for the Fisher score `FisherScore_windows.topSNP.bed`: \  
\  
```shell
python FisherScore/Fish.py --xpehh XPEHH_final.bed --pbs windows_20SNP_step5.bed --out FisherScore_windows
```
\  
merge the top regions `FisherScore_windows.topSNP.bed` together, score of the regions is the top score. Extract the top   
10 regions from it: \  
\  
```shell
bedtools merge -i FisherScore_windows.topSNP.bed -d 10000 -c 4 -o max > top_regions_Fisher_merged.bed
```
#### 5.3. Get p-value for the 10 top regions with random sampling approach
```shell
snakemake -srandom_sampling/p_val.smk -kp --jobs 100 --profile [profile_file]
``` 
## 6. Relate analysis  
We used [Relate v1.8](https://myersgroup.github.io/relate/) \  
The relate snakemake pipeline is ran with the following command \  
\
```shell
snakemake -s Relate/Make_trees.smk -kp --profile [profile_file] 
```
\  
and will create several folders in the working directory:   
#### 6.1. Create relate input file  
The snakemake pipeline will create an input_relate folder including the `*sample`, `*haps`, `*annot` and `*dist` file used by Relate.  
It  also includes the file `poplabel.txt` that is the population files correctly ordered for the further relate analysis   
#### 6.2. Make trees  
Will create recombination trees for all the samples included in the vcf files. We used the default   
parameters of 1.25e-8 for the mutation rate and 3000 for the effective population size of haplotypes. The output are `*anc` and `*mut` files in the `trees_all/` directory   
#### 6.3. Extracting trees  
Trees will be extracted for every population you included in `config.sk`. The extracted trees (`*anc` and `*mut` files) are   
found in the `extracted_trees/` directory   
#### 6.4. Effective population size  
Effective population size will be estimated for each of the populations you included in the `config.sk` file. We used the   
default recombination rate of 1.25e-8 and we specified the first time interval from 0 to 10<sup>2</sup> years ago and   
then the periods between 10<sup>2</sup> and 10<sup>7</sup> is split into successive bins with a step of 10<sup>0.1</sup> years. We used Relate default   
parameters of 28 years per generation, the number of cycles and the fractions of trees to be dropped.\  
Outputs (`*coal`) are found in the `pop_size/` folder   
## 7. Clues analysis   
We used [clues v1.0](https://github.com/standard-aaron/clues/) \  
```shell
snakemake -s Clues/Clues.smk --profile [profile_file] --groups create_clues_input=group1 clues=group2 --group-components group1=20 group2=10
```
#### 7.1. Extract local trees with Relate  
The first rule `create_clues_input` of the snakemake script extracts the local trees for each SNPs in the list of interest   
for the population indicated in the `config.sk` file. It samples the branch length for those trees. The branch length of   
the focal SNP tree is resampled 200 times (`--num_samples 200` option of the `SampleBranchLengths.sh` script ), and each   
iteration result is recorded to estimate the uncertainty in the age estimate of each node. To take the uncertainty in the   
branch length estimate into account, this is performed this several times (`number_of_test` in the `config.sk` file).   
They output files `*timeb` is found in the `clues_input/` output directory   
#### 7.2. Running clues for each focal trees  
Clues is run for each of the sampled branches output are found in the `clues_output/`directory.  
  
## 8. Association genotype phenotype  
We used [GEMMA v0.98.4](https://github.com/genetics-statistics/GEMMA/tree/master)  
#### 8.1. Convert vcf files to bed files  
Using [plink v1.9](https://www.cog-genomics.org/plink/).These files will be used to computed the relatedness matrice `output/Papuan_with_pheno_centered.cXX.txt`\  
\  
```shell
plink --vcf unrelated_filtered_masked.vcf.gz --keep papuans.ID --output-missing-genotype 0 --double-id --make-bed --out unrelated_masked_with_pheno
``` 
#### 8.2. Extract the candidate SNPs and convert to bed file format.  
These file will be used to compute association only for the top SNP\  
\  
```shell
plink --vcf unrelated_filtered_masked.vcf.gz --keep papuans.ID  --extract range TopSNP.bed --output-missing-genotype 0 --double-id --set-missing-var-ids @:# --make-bed --out Papuan_unrelated_masked_with_pheno_Top_SNP
```
#### 8.3. Add the pheno info to the fam files from the top SNP and the all site file  
```shell
python create_fam.py Papuan_unrelated_masked_with_pheno.fam residuals_pheno.cov  
python create_fam.py Papuan_unrelated_masked_with_pheno_Top_SNP.fam residuals_pheno.cov
```
#### 8.4. Compute the relatedness matrices 
```shell
gemma-0.98.4-linux-static-AMD64 -bfile Papuan_unrelated_masked_with_pheno -gk 1 -o Papuan_with_pheno_centered
``` 
#### 8.5. Gemma association  
To run for every phenotype column\  
\  
```shell
gemma-0.98.4-linux-static-AMD64 -bfile Papuan_unrelated_masked_with_pheno_Top_SNP -k output/Papuan_with_pheno_centered.cXX.txt -notsnp -miss 1 -lmm 4 -n [pheno_column] -outdir results_centered/ -o [pheno_name]
```
## 9. UKBiobank Summary Statistics analysis
We first downloaded the summary statistics [files](http://www.nealelab.is/uk-biobank) of UK bio bank (UKBB) from 
publicly available. Files can be downloaded from https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/. Both 
data and index should be downloaded in the same folder to use our codes. 
#### 9.1. Download files 
```shell
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/*.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix/*.tbi
```
#### 9.2. Closest snp present in UKBB
To find the closest snp present in the UKBB file set use this command
```shell
bash UKBB/UKBB_ClosestSNP.sh <snp.bed> <ukbb.tsv.bgz> > <closest.bed> 
```
You need two files to run the code.  
##### input:   
First the targeted SNP, written in bed format. Remember UKBB is in gr37, thus your bed file should be also in gr37. Can
be multiple snp in a same file. Should be just in different lines.
Second you need any of the tsv file that you have downloaded. As all the files have exactly same snps inside, it does 
not matter which of the file you use for this command.   
##### output: 
It will produce a bed file, where the closest snp will be written. For now, if there is no SNP in UKBB within 1kb 
upstream or downstream, the code will ignore that snp.   
#### 9.3. Extract UKBB values for target snps 
```shell
python UKBB/Run_ExtractScore.py --bed <closest.bed> <ukbb.tsv.bgz> 
```
This will extract UKBB results of pval_EUR and beta_EUR for the targeted snp present in closest.bed created in previous
step. 
As these have to be done for all the UKBB phenotypes, you can use this shell script to extract it for all the 
phenotypes automatically:
```shell
sh <ExtractScore.sh> <closest_snp.bed> <UKBB_folder>
```
##### input: 
It has two input.  
First the closest_snp.bed that we have created in the previous step (or any bed file where all the snp is present in 
UKBB).  
Second is the downloaded folder where all the UKBB summary statistics files are kept together.  
##### output: 
It will create scores name folder. In which it will dump all the snps in the folder (chr_position.tsv). Inside every 
file, you will have all the phenotypes with pval_EUR and beta_EUR together, which then can be sorted using an Excel 
sheet to get the significant values. 
#### 9.5. Extract UKBB values for all snps
In case you want to extract UKBB values for all SNPs (which we have used for random sampling), you can use the code. 
```shell
sh <ExtractScore_AllSNP.sh> <UKBB_folder>
```
##### input: 
You just haave to give the folder path where all UKBB summary files are stored.
##### output:
It will create Phenotypes.pval.gz files which will be tabix indexed. Meaning you can extract any snp using: 
```shell
tabix -h Phenotypes.pval.gz 1:6623020-6623020
```
Remember all snp and all phenotype is reading and writing all of very big files. We will suggest rather than put all the
files is same folder, break it smaller folders and then run this code independently to make it faster.  

## 10. Enrichment for blood phenotypes
```shell
snakemake -s Enrichment_UKBB/enrich.smk --profile [profile_file]`
```
