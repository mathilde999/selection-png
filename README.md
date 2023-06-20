# selection-png
Scripts developed and used in André et al., biorxiv 2022



### XP-EHH scan 
We used [selscan V2.0](https://github.com/szpiech/selscan) [[Sabeti, *et al.* Nature 2007]](https://www.nature.com/articles/nature06250)
#### create map files
`python Selscan/CreatingMapFiles.py --files_path chr/ --vcf`
#### run selscan
`selscan --xpehh --vcf High_chr{chr}.vcf.gz --vcf-ref Low_chr{chr}.vcf.gz --map High_Low_chr{chr}_GeneticMap.map --out High_chr{chr} --threads 4`


### PBS scan

#### Running PBS
`python RunPBS.py PBS.info`
#### PBS score for sliding windows
We create create sliding windows of X snp with Y snps step and give to each of these regions a PBS that is the average 
of the PBS of the X SNP that the window is made. We keep only the windows whose PBS score in the 99th percentile: output
"sliding_windows_snp_top_regions_PBS.bed". We make a file "mean_pos_PBS.bed" with the mean pos of the region (mean of 
the start and end bp) and the PBS of the region (average of the PBS of all the snp use to create the windows) 
We used windows of 20 SNPs sliding by 5 SNPs
`python sliding_windows_snp.py -b High_Low_EUR_AFR.pbs -w 20 -s 5`

We then merge top regions together when they were closer than 10kb with the closest windows  using [bedtools v2.29.2](https://bedtools.readthedocs.io/en/latest/index.html). The PBS score of the regions will be the top PBS of the regions uses to create the bigger regions
`bedtools merge -i sliding_windows_snp_top_regions_PBS.bed -d 10000 -c 4 -o max > windows_and_PBS_top_windows_final_merged.bed`

### Relate analysis
We used [Relate v1.8](https://myersgroup.github.io/relate/) [[Speidel, *et al.* Nat. Gen. 2019]](https://www.nature.com/articles/s41588-019-0484-x)

The relate snakemake pipeline is ran with the following command 
`snakemake -s Make_trees.smk -kp --profile [profile_file] `
and will create several folders in the working directory:
#### Create relate input file
The snakemake pipeline will create an input_relate folder including the *sample, *haps, *annot and *dist file used by Relate
It will also includes the file poplabel.txt that is the population files correctly ordered for the futher relate analysis
#### Make trees
Will create recombination trees for all the samples included in the vcf files you have inputed.
We used the default parameters of 1.25e-8 for the mutation rate and 3000 for the effective population size of haplotypes.
The output are *anc and *mut files in the trees_all/ directory
#### Extracting trees
Trees will be extracted for every population you included in config.sk. The extracted trees (*anc and *mut files) are 
found in the extracted_trees/ directory
#### Effective population size
Effective population size will be estimated for each of the populations you included in the config.sk file. 
We used the default recombination rate of 1.25e-8 and we specified the first time interval from 0 to 10² years ago and
then the periods between 10² and 10^7^ is split into successive bins with a step of 10^0.1^ years. We used Relate default
parameters of 28 years per generation, the number of cycles and the fractions of trees to be dropped.
Outputs (*coal) are found in the pop_size/ folder

### Clues analysis 
We used [clues v1.0](https://github.com/standard-aaron/clues/) [Stern, *et al.* Plos Gen. (2019)](https://journals.plos.org/plosgenetics/article/metrics?id=10.1371/journal.pgen.1008384) 

`snakemake -s Clues.smk --profile [profile_file]--groups create_clues_input=group1 clues=group2 --group-components group1=20 group2=10`

#### Extract local trees with Relate
The first rule of the pipeline extract the local trees for each SNPs in the list of interest for the population 
indicated in the config.sk file. The Relate rule of the snakemake pipeline will sampled the branch length for those 
trees. The branch length of the focal SNP tree is resampled 200 times (“--num_samples 200” option of the 
“SampleBranchLengths.sh” script ), and each iteration result is recorded to estimate the uncertainty in the age estimate 
of each node. To take the uncertainty in the branch length estimate into account, this can performed this several times (number_of_test in the config.sk file). 
For each SNP, we thus compute sampled branches X times. They output files *timeb can be found in the clues_input directory

#### Running clues for each focal trees
Clues is run for each of the sampled branches output are found in the clues_output directory. 

