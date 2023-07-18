

i=$1 # give a chromosome
#i=12
res_folder=$2 # give a folder name

cohort=vcfs_all_positions/phasing/chr${i}_phased.YRI_PAP_gen_vcf.bcf # phased filtered gvcf
indV=individuals_yri_pap.json # json of samples

hmmix create_outgroup \
  -ind=$indV \
  -vcf=$cohort \
  -out=$res_folder/outG/chr$i.outgroup.txt \
  -ancestral=homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_$i.fa \
  -refgenome=ref_hg38/chr${i}_AC.fa
  
hmmix mutation_rate \
  -outgroup=$res_folder/outG/chr$i.outgroup.txt \
  -window_size=1000000 \
  -out $res_folder/MutR/chr$i.mutationrate.bed

hmmix create_ingroup \
  -ind=$indV \
  -vcf=$cohort \
  -out=$res_folder/observ/chr$i.obs \
  -outgroup=$res_folder/outG/chr$i.outgroup.txt \
  -ancestral=homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_$i.fa
  
#  -weights=strickmask.bed \
##  We don't use  strictmask here since the analysis was done on hg38 higv-cov gvcf genomes, 
## to which during their calling process was applied an empirical genomic mask... 
## While Archaic genomes are not used in the hmmix, and anyway were prefiltered with their own masks before the annotation process later on

