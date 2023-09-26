i=$1

## merging (g)vcf with archaics samples, to align them already and to try to use with other tool (IBDmix)

pth=combined_all-sites/full_clean

bcftools view --force-samples -S PAPs.YRI.samples $pth/chr${i}_all-sites_clean.vcf.gz -Ob -o vcfs_G/chr$i.YRI_PAP.bcf && \
  bcftools index vcfs_G/chr$i.YRI_PAP.bcf
  
bcftools annotate -x ^FORMAT/GT vcfs_G/chr$i.YRI_PAP.bcf -Ob -o vcfs_G/chr$i.YRI_PAP.no_anno.bcf && \
  bcftools index vcfs_G/chr$i.YRI_PAP.no_anno.bcf


bcftools norm -c ws --fasta-ref references/hg38/v0/Homo_sapiens_assembly38.fasta \
  archaic_genoms/four_hg38/chr$i.4_genomes.hg38.vcf.gz -Ob -o archaic_genoms/four_hg38_X/chr$i.4_genomes.hg38.fixed_ref.bcf

bcftools index archaic_genoms/four_hg38_X/chr$i.4_genomes.hg38.fixed_ref.bcf

bcftools merge archaic_genoms/four_hg38_X/chr$i.4_genomes.hg38.fixed_ref.bcf vcfs_G/chr$i.YRI_PAP.no_anno.bcf \
  -Ob -o vcfs_G/chr$i.Arch_YRI_PAP.bcf # in this file, there are REF ALT and REF . positions
  
bcftools norm -m - vcfs_G/chr$i.Arch_YRI_PAP.bcf \
  -Ob -o vcfs_G/chr$i.Arch_YRI_PAP.tmp.bcf && bcftools index vcfs_G/chr$i.Arch_YRI_PAP.tmp.bcf
  
bcftools view -v snps vcfs_G/chr$i.Arch_YRI_PAP.tmp.bcf -Ob -o vcfs_G/chr$i.Arch_YRI_PAP.spl.snps.bcf && \
  rm vcfs_G/chr$i.Arch_YRI_PAP.tmp.bcf # in this file, there are REF ALT with AC == 0 since they were subsetted from the other file
  
bcftools index vcfs_G/chr$i.Arch_YRI_PAP.spl.snps.bcf
  
bcftools filter -e 'AC == 0' vcfs_G/chr$i.Arch_YRI_PAP.spl.snps.bcf -Ob -o vcfs_G/chr$i.Arch_YRI_PAP.spl.snps.AC.bcf  #and this one will do REF ALT with some AC



## for phasing we remove archaics back

# choose them since they already was merged where somebody has alleles
bcftools view -s ^AltaiNeandertal,Vindija33.19,Chagyrskaya-Phalanx,Denisova vcfs_G/chr$i.Arch_YRI_PAP.spl.snps.AC.bcf |\
  bcftools norm -m - --fasta-ref /gpfs/space/databases/broadinstitute/references/hg38/v0/Homo_sapiens_assembly38.fasta \
   -Ob -o phasing/chr$i.YRI_PAP_Arch_pos.spl.snps.AC.realign.tmp.bcf

# -m2 -M2 -v snps = stay all;
bcftools view -m2 -M2 -v snps -i 'F_MISSING <= 0.95 & FILTER="PASS"' \
  phasing/chr$i.YRI_PAP_already_merged_Arch.spl.snps.AC.realign.tmp.bcf \
  -Ob -o phasing/chr$i.YRI_PAP_Arch_pos.spl.snps.AC.realign.FILTERED.bcf

##
shapeit4.2 \
  --input phasing/chr{i}.YRI_PAP_Arch_pos.spl.snps.AC.realign.FILTERED.bcf \
  --map maps/genetic_map_hg38_chr{i}.txt.gz \
  --output phasing/chr{i}_phased.YRI_PAP_gen_vcf.bcf \
  --region chr$i \
  --thread 10 \
  --sequencing
  
bcftools index phasing/chr{i}_phased.YRI_PAP_gen_vcf.bcf