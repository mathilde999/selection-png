rm -fR phenotypes temp
mkdir -p temp
ls $1/*.bgz  | \
  awk -v q="'" -F "/" '{print "zcat",$0,"| head -n 1 | tr \"\\t\" \"\\n\" | \
    awk  "q"{print NR\"\\t\"$0\"\\t"$0"\\t"$NF"\"}"q"| grep \"pval_EUR\"| cut -f 1,3- "}'> temp/extract_pvalcol_index.sh
sh temp/extract_pvalcol_index.sh  > temp/path_pval_eurindex.txt
mkdir -p phenotypes
awk '{print "zcat ",$2," | cut -f "$1" >","phenotypes/"substr($3,4,length($3)-11)".pval"}' temp/path_pval_eurindex.txt \
  > temp/extract_pval.sh
sh temp/extract_pval.sh
paste \
  <(cat <(echo -e "#CHROM\tPOS\tREF\tALT") \
    <(zcat biomarkers-30600-both_sexes-irnt.tsv.bgz | cut -f 1,2,3,4 | tail -n +2| cut -f 1,2,3,4)) \
  <(cat <(ls phenotypes/*.pval |awk -F "/" '{print substr($NF,1,length($NF)-5)}'| paste -sd "\t" -) \
    <(paste phenotypes/*.pval | tail -n +2 ))  > Phenotypes.pval
bgzip Phenotypes.pval
tabix -p vcf Phenotypes.pval.gz
rm -fR phenotypes temp