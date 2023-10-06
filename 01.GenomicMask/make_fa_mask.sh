#!/bin/bash
chr=22
#the fasta file must start with the same chr than used in the mask bed files
echo ">chr${chr}"> homo_sapiens_ancestor_all_N_${chr}.fa
#replace every character in the fasta file per N before turning back to P the site that are included in the positive mask
zcat  homo_sapiens_ancestor_${chr}.fa.gz|tail -n+2|sed -e 's/[a-z A-Z . -]/N/g' >> homo_sapiens_ancestor_all_N_${chr}.fa
bedtools maskfasta -fi homo_sapiens_ancestor_all_N_${chr}.fa -bed full_positive_mask_chr${chr}.bed -fo full_mask_chr${chr}.fa -mc P
#now we change the title of the fasta file again so that it is the same than in the original gr38 human fasta file
header=$(head -${chr} list_headers_fa.txt|tail -1)
sed -i "1s/.*/$header/" full_mask_chr${chr}.fa
gzip full_mask_chr${chr}.fa
gzip homo_sapiens_ancestor_all_N_${chr}.fa
