# selection-png    
Scripts developed and used in [André et al., biorxiv 2022](https://doi.org/10.1101/2022.12.15.520226)
André, M., Brucato, N., Hudjasov, G., Pankratov, V., Yermakovich, D., …, Mondal, M., Ricaut, F.-X. (2022, December 15). Positive selection in the genomes of two Papua New Guinean populations at distinct altitude levels (p. 2022.12.15.520226). bioRxiv. https://doi.org/10.1101/2022.12.15.520226    

We used vcf of 1000 Genomes 30x on GRCh38 as example data. The used vcf file can be downloaded [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/) 
similarly for the [phased vcf](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/) and the cram files.

You can use conda to install the dependencies necessary to run the following analysis and listed in `png_selection.yml`. 
To install conda please visit anaconda webpage. 
```shell  
conda install --file png_selection.yml
```
## 1. Creation of a genomic mask   
#### 1.1. Create coverage mask     
We included in the positive coverage mask sites with a minimum base quality of 20, an alignment minimum mapping quality 
of 20 and downgrading mapping quality for reads containing excessive mismatches with a coefficient of 50. From these sites,
we excluded indels, sites with more than two alleles, sites with a maximum missing rate of 5%. We also masked sites whose 
depth of coverage summed across all samples was higher or lower than the sum of the average depth across the dataset by a factor of 2-fold.
```shell  
snakemake -s 1.GenomicMask/mask.smk
```  
#### 1.2. Create positive mask with coverage mask and mappability mask     
We create a positive mask using [bedtools v2.29.2](https://bedtools.readthedocs.io/en/latest/index.html) including the site
present in both the coverage mask we created in step 1.1. `1.GenomicMask/mask_all_sort_chr22.bed` and `gr38.mask.bed`, the [mappability mask](https://share.eva.mpg.de/index.php/s/ygfMbzwxneoTPZj) 
used with [MSMC](https://github.com/stschiff/msmc-tools/tree/master)  liftover to GR38 
```shell  
bedtools intersect -a 1.GenomicMask/mask_all_sort_chr22.bed -b gr38.mask.bed > 1.GenomicMask/mask_coverage_and_map_chr22.bed
```  
#### 1.3. Filtering out the variant without PASS flag from the variant calling    
Keep only the sites that pass the PASS filter from the variant calling unfiltered vcf files `example_dataset/example_raw_chr22.vcf.gz`   
```shell 
vcftools --gzvcf example_dataset/example_raw_chr22.vcf.gz --remove-filtered-all --recode --stdout |bgzip > 1.GenomicMask/sites_PASS_chr22.vcf.gz    
```    
Output the site without PASS flag from the variant calling unfiltered vcf files `example_dataset/example_raw_chr22.vcf.gz`
```shell
bcftools isec -C -o 1.GenomicMask/sites_not_PASS_chr22.bed -Oz example_dataset/example_raw_chr22.vcf.gz 1.GenomicMask/sites_PASS_chr22.vcf.gz   
```
Format of the output file `1.GenomicMask//sites_not_PASS_chr22.bed`
```shell
awk '{print $1 "\t" ($2 - 1) "\t" $2}' 1.GenomicMask/sites_not_PASS_chr22.bed|bedtools merge -i - > 1.GenomicMask/sites_not_PASS_chr22_format.bed
```
#### 1.4. Remove the not PASS sites from the global positive mask
```shell
bedtools subtract -a 1.GenomicMask/mask_coverage_and_map_chr22.bed -b 1.GenomicMask/sites_not_PASS_chr22_format.bed > 1.GenomicMask/full_positive_mask_chr22.bed
```  
#### 1.5. Create the mask fasta file      
Convert the positive mask to fasta format using the [ancestral genome](https://ftp.ensembl.org/pub/release-93/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38/)     
from the Ensembl 93 release `homo_sapiens_ancestor_22.fa.gz`. This will output the final positive mask `1.GenomicMask/full_positive_mask_chr22.fa` used in the further analyses.
```shell  
1.GenomicMask/make_fa_mask.sh
``` 
## 2. PCA and ADMIXTURE
#### 2.1. filtering
```shell
vcftools --gzvcf example_dataset/example_raw_chr22.vcf.gz --max-alleles 2 --min-alleles 2 --remove-indels --max-missing 0.95 --remove-filtered-all --bed 1.GenomicMask/full_positive_mask_chr22.bed --maf 0.05 --plink --out 2.PCA_Admixture/example_chr22
```
#### 2.2. LDpruning
We prune variants in high linkage disequilibrium using [plink v1.9](https://www.cog-genomics.org/plink/) and the default parameters of 50 variants count window 
shifting from five variants and a variance inflation factor (VIF) threshold of 2. 
```shell
plink --file 2.PCA_Admixture/example_chr22 --indep 50 5 2 --out 2.PCA_Admixture/example_chr22 --noweb #output example_chr22.prune.in
plink --file 2.PCA_Admixture/example_chr22 --extract 2.PCA_Admixture/example_chr22.prune.in  --out 2.PCA_Admixture/example_pruned_chr22 --noweb --recode #output example_pruned_chr22.ped and example_pruned_chr22.map
```
#### 2.3. Run PCA
We performed the Principal Component Analysis (PCA) on the unrelated LD pruned dataset `all_biSNP_PASS_maf_unrelated_VIFLDPruned.ped`  
using the smartpca program from the [EIGENSOFT v.7.2.0 package](https://github.com/DReichLab/EIG)
```shell
smartpca -p 2.PCA_Admixture/par.file # Will output example_pruned_chr22.evec and example_pruned_chr22.evac
```
#### 2.4. Convert ped to bed format
Admixture required a bed format
```shell
plink --file 2.PCA_Admixture/example_pruned_chr22 --make-bed --out 2.PCA_Admixture/example_pruned_chr22 # Will output example_pruned_chr22.bed
```
#### 2.4. Run ADMIXTURE
We ran [ADMIXTURE v1.3](https://dalexander.github.io/admixture/)  for the components K=2 to K=10 on the unrelated LD pruned dataset `2.PCA_Admixture/example_pruned.bed`
For each component, ADMIXTURE computes the cross-validation error using k-fold cross-validation procedure. We set the k parameter to 100.
In order to select the model with the most likely number of components, in thi example, we generated the cross-validation error 3 times for each component. 
```shell
./submit.sh 2.PCA_Admixture/example_pruned_chr22
```
This script will create 10 `trial` directory, each of which will have P and Q admixture output files for component k=2 to k=10

## 3. Phasing    
We used [shapeit4 v4.2.2](https://odelaneau.github.io/shapeit4/). Input parameter are found in `Phasing/config.sk`
We used the genetic map (`genetic_map_hg38_chr22.txt.gz`) from [Eagle](https://alkesgroup.broadinstitute.org/Eagle/downloads/tables) for GRCh38 to shapeit format 
```shell
snakemake -s 3.Phasing/shapeit4.smk
```    
That snakemake pipeline will create two directories in your working directory. Directory `chr/` includes unphased vcf files    
per chr filtered for the filters indicated in `3.Phasing/config.sk`.The other directory created is `phased/` that includes all the phased vcf files per chr

## 4. XP-EHH scan     
We used [selscan V2.0](https://github.com/szpiech/selscan) with highlanders as the target population and lowlanders as the reference population.
For this example run we used all the 99 CEU individuals from 1000G as the target population and the 103 CHB individuals from 1000G as the reference population. We used the already phased [phased vcf](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
for which we kept the biallelil SNP present with a max missing rate of 0.95. 
Selscan v2.0.0 is not available on Anaconda. You must install it following the instruction found in the selscan gitHub repository.
#### 4.1. Create genetic map files for selscan from vcf files
We used the genetic map (`genetic_map_hg38_chr22.txt.gz`) from [Eagle](https://alkesgroup.broadinstitute.org/Eagle/downloads/tables) for GRCh38. You need to convert it to plink .map format first
```shell
tail -n+2 genetic_map_hg38_chr22.txt|awk '{print $1 , $1":"$2,$4, $2}'> 4.XPEHH/genetic_map_hg38_chr22_plink_format.txt
vcftools --gzvcf 4.XPEHH/CEU_CHB_YRI_phased_chr22_biallel_snp.vcf.gz --plink --out 4.XPEHH/chr22_input # Will output chr22_input.map
Rscript 4.XPEHH/GeneticMapIntrapolate.r 4.XPEHH/chr22_input.map 4.XPEHH/CEU_CHB_YRI_chr22_GeneticMap.map  4.XPEHH/genetic_map_hg38_chr22_plink_format.txt # Will CEU_CHB_YRI_chr22_GeneticMap.map
```
#### 4.2. Run selscan 
Selscan uses phased vcf files separated by chromosome with no missing data. It will also need one population per file.
```shell
selscan --xpehh --vcf 4.XPEHH/CEU_phased_chr22_biallel_snp.vcf.gz --vcf-ref 4.XPEHH/CHB_phased_chr22_biallel_snp.vcf.gz --map 4.XPEHH/CEU_CHB_YRI_chr22_GeneticMap.map --out 4.XPEHH/CEU_CHB_chr22 # Output is 4.XPEHH/CEU_CHB_chr22.xpehh.out
``` 
#### 4.3. Normalization 
```shell
norm --xpehh --files 4.XPEHH/CEU_CHB_chr22.xpehh.out # Output is 4.XPEHH/CEU_CHB_chr22.xpehh.out.norm 
```   
#### 4.4. Top region score    
This script will output the SNP with xpehhnorm in the 99th percentile 
```shell
python 4.XPEHH\XPEHH_bed_99.py 4.XPEHH/CEU_CHB_chr22.xpehh.out.norm # Output is 4.XPEHH/CEU_top_snp_XPEHH_99.bed   
```  
If you are interested in the top SNP for the reference population. This script will output the SNP with xpehhnorm in the 1st percentile.    
```shell
python 4.XPEHH\XPEHH_bed_1.py 4.XPEHH/CEU_CHB_chr22.xpehh.out.norm # Output is 4.XPEHH/CHB_top_snp_XPEHH_1.bed   
```
#### 4.5. Merge top SNPs in genomic regions  
We then merge top SNPs together when they were closer than 10kb with the 
closest windows  using [bedtools v2.29.2](https://bedtools.readthedocs.io/en/latest/index.html)    
The XP-EHH score of the regions will be the top XP-EHH of the top SNPs uses to create the bigger regions  
```shell
bedtools merge -i 4.XPEHH/CEU_top_snp_XPEHH_99.bed -d 10000 -c 4 -o max > 4.XPEHH/CEU_merged_top_snp_XPEHH_99.bed
bedtools merge -i 4.XPEHH/CHB_top_snp_XPEHH_1.bed -d 10000 -c 4 -o min > 4.XPEHH/CHB_merged_top_snp_XPEHH_1.bed
```
#### 4.6. Get p-value for the 10 top regions with random sampling approach     
```shell
snakemake -s random_sampling/p_val.smk -kp --jobs 100 --profile [profile_file]
```  
## 5. PBS scan  
#### 5.1. Running PBS
```shell
python 5.PBS/RunPBS.py --popfile example_dataset/example.pop --pbspop CEU:CHB:YRI --dropna 5.PBS/CEU_CHB_YRI_chr22_biallel_snp.vcf.gz > 5.PBS/CEU_example_chr22.pbs
``` 
#### 5.2. PBS score for sliding windows    
We create sliding windows of X snp with Y snps step and give to each of these regions a PBS that is the average     
of the PBS of the X SNP that the window is made with (see `CEU_windows_20SNP_step5.bed`). We keep only the windows whose PBS score in 
the 99th percentile:`CEU_sliding_windows_snp_top_regions_PBS.bed`. It also output `CEU_mean_pos_PBS.bed` with the mean pos of the region (mean of     
the start and end bp) and the PBS of the region (average of the PBS of all the snp use to create the windows)     
We used windows of 20 SNPs sliding by 5 SNPs 
```shell
python 5.PBS/sliding_windows_snp.py -b 5.PBS/CEU_example_chr22.pbs -w 20 -s 5 
```
We then merge top windows `CEU_sliding_windows_snp_top_regions_PBS.bed` together when they were closer than 10kb with the closest 
windows  using [bedtools v2.29.2](https://bedtools.readthedocs.io/en/latest/index.html).    
The PBS score of the regions will be the top PBS of the regions uses to create the bigger regions 
```shell
bedtools merge -i 5.PBS/CEU_sliding_windows_snp_top_regions_PBS.bed -d 10000 -c 4 -o max > 5.PBS/CEU_windows_and_PBS_top_windows_final_merged.bed
```  
#### 5.3. Get p-value for the 10 top regions with random sampling approach     
```shell
snakemake -s random_sampling/p_val.smk -kp --jobs 100 --profile [profile_file]
```
## 6. Fisher score    
We computed Fisher score for the same windows used to generate PBS `5.PBS/CEU_windows_20SNP_step5.bed` or `5.PBS/CHB_windows_20SNP_step5.bed` depending
if we want to compute teh reference for the target population CEU or the reference population CHB.
    
#### 6.1. Generate XP-EHH score for the sliding windows 
Make XP-EHH output as a bed file:   
```shell    
cat 4.XPEHH/CEU_CHB_chr22.xpehh.out.norm |awk  '{print "chr22\t" $2-1 "\t" $2 "\t" $9}'|tail -n+2 > 6.FisherScore/CEU_CHB_chr22.xpehh.bed
```   
If we are interested to get the XPEHH score from the reference pop (CHB), we need to invert the value of XP-EH
```shell    
python get_inverted_xpehh.py -b 6.FisherScore/CEU_CHB_chr22.xpehh.bed # output inverted.xpehh.bed
```   

Extract the PBS windows coordinates: 
```shell
cut -f1,2,3 5.PBS/CEU_windows_20SNP_step5.bed|tail -n+2 > 6.FisherScore/CEU_windows.bed
cut -f1,2,3 5.PBS/CHB_windows_20SNP_step5.bed|tail -n+2 > 6.FisherScore/CHB_windows.bed
```
make a file that list each snp associated to each windows and their xpehh score:
```shell  
bedtools intersect -a 6.FisherScore/CEU_windows.bed -b 6.FisherScore/CEU_CHB_chr22.xpehh.bed -wa -wb > 6.FisherScore/CEU_windows_xpehh.bed
bedtools intersect -a 6.FisherScore/CHB_windows.bed -b 6.FisherScore/inverted.xpehh.bed -wa -wb > 6.FisherScore/CHB_windows_xpehh.bed
```  
Get the top XPEHH score of all the snp in each windows (top score in the region):
```shell  
bedtools groupby -i 6.FisherScore/CEU_windows_xpehh.bed -g 1-3 -c 7  -o max |awk 'BEGIN {print "CHR\tSTART\tEND\tNORM"} {print $0}'> 6.FisherScore/CEU_windows_XPEHH_final.bed
bedtools groupby -i 6.FisherScore/CHB_windows_xpehh.bed -g 1-3 -c 7  -o max |awk 'BEGIN {print "CHR\tSTART\tEND\tNORM"} {print $0}'> 6.FisherScore/CHB_windows_XPEHH_final.bed
```

#### 6.2. Create bed file with fisher score for each sliding windows:
keep only sliding windows in the 99 percentile for the Fisher score: 
```shell  
python 6.FisherScore/SubmitFish.py --xpehh 6.FisherScore/CEU_windows_XPEHH_final.bed --pbs 5.PBS/CEU_windows_20SNP_step5.bed --out 6.FisherScore/CEU_FisherScore_windows
python 6.FisherScore/SubmitFish.py --xpehh 6.FisherScore/CHB_windows_XPEHH_final.bed --pbs 5.PBS/CHB_windows_20SNP_step5.bed --out 6.FisherScore/CHB_FisherScore_windows
```  
merge the top regions together, score of the regions is the top score. Extract the top 10 regions from it: 
```shell  
bedtools merge -i 6.FisherScore/CEU_FisherScore_windows.topSNP.bed -d 10000 -c 4 -o max > 6.FisherScore/CEU_top_regions_Fisher_merged.bed
bedtools merge -i 6.FisherScore/CHB_FisherScore_windows.topSNP.bed -d 10000 -c 4 -o max > 6.FisherScore/CHB_top_regions_Fisher_merged.bed
```  
#### 6.3. Get p-value for the 10 top regions with random sampling approach  
```shell  
snakemake -s random_sampling/p_val.smk -kp --jobs 100 --profile [profile_file]
``` 
## 7. Relate analysis    
We used [Relate v1.8](https://myersgroup.github.io/relate/)  
The relate snakemake pipeline is ran with the following command. You must install Relate following relate gitHub directory.
```shell  
snakemake -s Relate/Make_trees.smk
```
and will create several folders in the working directory:     
#### 7.1. Create relate input file    
The snakemake pipeline will create an input_relate folder including the `*sample`, `*haps`, `*annot` and `*dist` file used by Relate.    
It  also includes the file `poplabel.txt` that is the population files correctly ordered for the further relate analysis     
#### 7.2. Make trees    
Will create recombination trees for all the samples included in the vcf files. We used the default     
parameters of 1.25e-8 for the mutation rate and 3000 for the effective population size of haplotypes. The output are `*anc` and `*mut` files in the `trees_all/` directory     
#### 7.3. Extracting trees    
Trees will be extracted for every population you included in `config.sk`. The extracted trees (`*anc` and `*mut` files) are     
found in the `extracted_trees/` directory     
#### 7.4. Effective population size    
Effective population size will be estimated for each of the populations you included in the `config.sk` file. We used the     
default recombination rate of 1.25e-8 and we specified the first time interval from 0 to 10<sup>2</sup> years ago and     
then the periods between 10<sup>2</sup> and 10<sup>7</sup> is split into successive bins with a step of 10<sup>0.1</sup> years. We used Relate default     
parameters of 28 years per generation, the number of cycles and the fractions of trees to be dropped.
Outputs (`*coal`) are found in the `pop_size/` folder     

## 8. Clues analysis
#### 8.1. Create a bed files with all the SNPs that are part of the genomic regions of interest
We have created a bed file `regions_of_interest.bed` including two random regions of the chr22 that we will consider as the genormic region of interest for which we want to run CLUES.
We are creating another bed files including all the SNPs that are included in these regions of interest. 
``shell
bcftools query -f '%CHROM\t%POS0\t%END\n' CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz -R 8.CLUES/regions_of_interest.bed  -o 8.CLUES/SNPs_of_interest.bed
``
#### 8.2. generate the input file for clues snakemake 
This will filter out the sites from `SNPs_of_interest.bed` that are not present in the population extracted tree
``shell
python bed_sites_to_kept_list.py -b 8.CLUES/SNPs_of_interest.bed -t 7.RELATE/extracted_trees/CEU -o 8.CLUES/SNPs_of_interest
``
#### 8.3 Running CLUES
We used [clues v1.0](https://github.com/standard-aaron/clues/)    
```shell  
snakemake -s Clues/Clues.smk
```
The first rule `create_clues_input` of the snakemake script extracts the local trees for each SNPs in the list of interest     
for the population indicated in the `config.sk` file. It samples the branch length for those trees using a RELATE script. The branch length of     
the focal SNP tree is resampled 200 times (`--num_samples 200` option of the `SampleBranchLengths.sh` script ), and each     
iteration result is recorded to estimate the uncertainty in the age estimate of each node. To take the uncertainty in the     
branch length estimate into account, this is performed this several times (`number_of_test` in the `config.sk` file).     
They output files `*timeb` is found in the `clues_input/` output directory.
Clues is run for each of the sampled branches output are found in the `clues_output/`directory.    
    
## 9. Association genotype phenotype 
We used [GEMMA v0.98.4](https://github.com/genetics-statistics/GEMMA/tree/master). None phenotype data are available for 1000 genomes individuals.
For the following example run we fill use the example datasets provided on GEMMA gitHub page. When working with your own vcf files and phenotype data, you must first: 
- Convert the main vcf file to plink bed file format. (`mouse_hs1940.bed`)
- Extract the candidate SNPs from the main vcf and convert to plink bed file format. (`target_SNP.bed`) 
- Add the pheno info to the fam files from the top SNP and the all site file. (`mouse_hs1940.fam`)

#### 9.1. Compute the centered relatedness matrice   
```shell  
gemma -bfile 9.GEMMA/mouse_hs1940 -gk 1 -o 9.GEMMA/mouse_hs1940_centered #will output 9.GEMMA/mouse_hs1940_centered.cXX.txt
``` 
#### 9.5. Gemma association    
These file will be used to compute association only for the target SNP. 
To run for every phenotype column of interest
```shell  
gemma -bfile 9.GEMMA/target_SNP -k 9.GEMMA/mouse_hs1940_centered.cXX.txt -notsnp -miss 1 -lmm 4 -n 1 -outdir results_centered/ -o pheno_1 # will output pheno_1.assoc.txt
```  
## 10. UKBiobank Summary Statistics analysis  
We first downloaded the summary statistics [files](http://www.nealelab.is/uk-biobank) of UK bio bank (UKBB) from   
publicly available. Files can be downloaded from https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/. Both   
data and index should be downloaded in the same folder to use our codes.   
#### 10.1. Download files   
```shell  
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/*.tsv.bgz wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix/*.tbi
```  
#### 10.2. Closest snp present in UKBB  
To find the closest snp present in the UKBB file set use this command  
```shell  
bash  UKBB_ClosestSNP.sh snp.bed biomarkers-30600-both_sexes-irnt.tsv.bgz closest.bed # Will output closest.bed, Present.bed and notpresent.snps
```  
You need two files to run the code.    
##### input:     
First the targeted SNP, written in bed format. Remember UKBB is in gr37, thus your bed file should be also in gr37. Can  
be multiple snp in a same file. Should be just in different lines.  
Second you need any of the tsv file that you have downloaded. As all the files have exactly same snps inside, it does   
not matter which of the file you use for this command.     
##### output:   
It will produce a bed file,`Present.bed` where the SNPs present in the UKBB summary file and, when the SNPs are not present, another file `closest.bed`
where the closest snp to the none-present ones are written . For now, if there is no SNP in UKBB within 1kb   
upstream or downstream, the code will ignore that snp.     
#### 10.3. Extract UKBB values for target snps   
```shell  
python 10.UKBB/Run_ExtractScore.py --bed 10.UKBB/closest.bed biomarkers-30600-both_sexes-irnt.tsv.bgz
python 10.UKBB/Run_ExtractScore.py --bed 10.UKBB/Present.bed biomarkers-30600-both_sexes-irnt.tsv.bgz 
```  
This will extract UKBB results of pval_EUR and beta_EUR for the targeted snp present in closest.bed created in previous  
step.   
As these have to be done for all the UKBB phenotypes, you can use this shell script to extract it for all the   
phenotypes automatically:  
```shell  
sh 10.UKBB/ExtractScore.sh closest.bed <UKBB_folder>  
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
#### 10.5. Extract UKBB values for all snps  
In case you want to extract UKBB values for all SNPs (which we have used for random sampling), you can use the code.   
```shell  
sh 10.UKBB/ExtractScore_AllSNP.sh <UKBB_folder>  
```  
##### input:   
You just have to give the folder path where all UKBB summary files are stored.  
##### output:  
It will create Phenotypes.pval.gz files which will be tabix indexed. Meaning you can extract any snp using:   
```shell  
tabix -h 10.UKBB/biomarkers-30600-both_sexes-irnt 1:6623020-6623020
```  
Remember all snp and all phenotype is reading and writing all of very big files. We will suggest rather than put all the  
files is same folder, break it smaller folders and then run this code independently to make it faster.    
  
#### 10.6. Enrichment for blood phenotypes
Will we use the file `pval_EUR.tsv.bgz`  generated in the precedent step to random sample 200bp (or whichever `size_windows_bp` you put in `10.UKKB/config.sk`) 
windows including at least one SNP significantly associated with at least one phenotype in the UKBB and store how many of these random 
windows include at least one SNP significantly associated with at least one phenotype of the [blood_count category](https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=100081)
of the UKBB (`blood_count_100081_list.full.name.txt`). We will compare this count distribution to the number of association with blood pheno found in our target population.
```shell  
snakemake -s UKBB/enrich.smk --profile [profile_file]
```

## 11. Archaic introgression for candidate regions for selection by hmmix (Skov et al., 2018)
You can use conda to install the dependencies necessary to run the following introgression workflow and listed in `png_introgression.yml`. 
To install conda please visit anaconda webpage.
You will also have to install hmmix following the instruction on the [hmmix github](https://github.com/LauritsSkov/Introgression-detection) 

```shell  
conda install --file png_introgression.yml
```
there are two folders, **run_hmmix** and **process_res**: <br>

#### 11.1. run_hmmix
contains the following scripts and implements the clear workflow that could be found here <br>
https://github.com/LauritsSkov/Introgression-detection <br>

1. processing_gvcf.sh # bcftools commands that were used for processing the (g)vcf files for the hmmix <br>
2. Skov_pipeline.preparing.sh # preparing specific input for the hmmix <br>
3. unite_chr_per_sample.sh # combining input for the hmmix model <br>
4. Skov_pipeline.hmm.sh # training and decoding by the hmmix model <br>
5. laucnh_Skov_pipeline.sh # launching 2-4 bullet points in sbatch mode for the regions of interest <br>

#### 11.2. **process_res** 
contains the following scripts for proccessing the hmmix output <br>

Introgression_processing.workflow.R # includes the workflow for finding the top frequent haplotype per candidate region for selection <br>
Pairwise_comparison.chr1_GBP.R # includes the counting of the pairwise distance between archaics in GBP region on chr1 <br>
