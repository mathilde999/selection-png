#!/bin/bash
#for each bed file you need to change the name of the chr with the header of the fasta file (header in the list_headers_fa.txt)
#remplace every character in the fasta file per N before turning back to P the site that are included in the positive mask
for chr in {1..22}
do
  head -${chr} list_headers_fa.txt|tail -1 > homo_sapiens_ancestor_all_N_${chr}.fa
  zcat homo_sapiens_ancestor_${chr}.fa.gz|tail -n+2|sed -e 's/[a-z A-Z . -]/N/g' >> homo_sapiens_ancestor_all_N_${chr}.fa
  bedtools maskfasta -fi homo_sapiens_ancestor_all_N_${chr}.fa -bed full_positive_mask_chr${chr}.bed -fo full_mask_chr${chr}.fa -mc P
  gzip full_mask_chr${chr}.fa
  gzip homo_sapiens_ancestor_all_N_${chr}.fa
done