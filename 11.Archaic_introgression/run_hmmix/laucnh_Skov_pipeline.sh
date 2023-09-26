#!/bin/bash -vx
#SBATCH -J hmmix_do
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH -t 0-10:00
#SBATCH --mem=99g
#SBATCH -p amd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=danatyermakovich@gmail.com

## upload environment
module load python/3.9.12; module load py-numpy/1.22.4; module load parallel; module load bcftools; module load vcftools
res_folder=out_gvcf_filt_phased
mkdir -p $res_folder/outG/ $res_folder/MutR/ $res_folder/observ/ $res_folder/trained/ $res_folder/res

## launch preparing steps
parallel -j 9 "bash Skov_pipeline.preparing.sh {} $res_folder" ::: {1..4} {6..10} {12..14} {17..19} 22

## concat chr per individuals by simple script (head;tail;cat commands)
bash unite_chr_per_sample.sh $res_folder # a 

## concat chr for mutatons rates
head -n1 $res_folder/MutR/chr9.mutationrate.bed > $res_folder/all_chr.mutationrate.bed
for i in {1..4} {6..10} {12..14} {17..19} 22; do tail -n+2 $res_folder/MutR/chr$i.mutationrate.bed >> $res_folder/all_chr.mutationrate.bed; done

## launch hmmix, train and decode 
parallel -j 10 -a ingroup.list "bash Skov_pipeline.hmm.sh {} $res_folder" 
# ingroup list is a list of target samples among which the analysis should be carried

