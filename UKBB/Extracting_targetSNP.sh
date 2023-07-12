source ~/src/python3/MayukhPython3Source.sh 
cut -f 1,2 LL_new_Top_SNP.txt | tail -n +2 | awk '{print "chr"$1"\t"$2-1"\t"$2}' > topsnp_gr38.bed
/gpfs/space/home/mayukh/Softwares/liftover/liftOver topsnp_gr38.bed /gpfs/space/home/mayukh/References/Liftover_Chain_Files/hg38ToHg19.over.chain topsnp_hg19.bed unmmaped.bed
/gpfs/space/home/mayukh/Softwares/liftover/liftOver topsnp_hg19.bed /gpfs/space/home/mayukh/References/Liftover_Chain_Files/hg19ToGRCh37.over.chain topsnp_gr37.bed unmmaped.bed
 module load htslib/1.15
tabix -R topsnp_gr37.bed /gpfs/space/home/mayukh/xOOA2/UKBiobank/biomarkers-30600-both_sexes-irnt.tsv.bgz | cut -f 1,2 | awk '{print $1"\t"$2-1"\t"$2}' > Present.bed
cat <(cut -f 1,3 topsnp_gr37.bed ) <(cut -f 1,3 Present.bed ) | sort | uniq -c | awk '{print $1"\t"$2"\t"$3}' | grep -vw "^2" | cut -f 2,3 > notpresent.snps
cat notpresent.snps | awk -v q="'" '{print "tabix /gpfs/space/home/mayukh/xOOA2/UKBiobank/biomarkers-30600-both_sexes-irnt.tsv.bgz "$1":"$2-1000"-"$2+1000" | cut -f 2| awk "q"{print \""$1"\\t\"$1-1\"\\t\"$1\"\\t\"($1-"$2")^2}"q"| sort -k 4,4n|head -n 1 | cut -f -3"}' > closest.sh
sh closest.sh > closest.bed
cat Present.bed closest.bed| sort -k1,1n -k2,2  > All.bed 
tabix -R All.bed /gpfs/space/home/mayukh/xOOA2/UKBiobank/biomarkers-30600-both_sexes-irnt.tsv.bgz | cut -f 1,2,3,4 > AllBed.alleles
ls /gpfs/space/home/mayukh/xOOA2/UKBiobank/*.bgz | awk -F "/" '{print "python /gpfs/space/home/mayukh/src/python3/GWAS/NealeLab/Run_ExtractScore.py --bedfile All.bed "$0" --beta_percent --plot > raw/"$NF".out"}'  > Exctracting.sh 
mkdir raw
python /gpfs/space/home/mayukh/src/python3/SubmittingJobs/EBCRun.py --script Exctracting.sh  --amd --log --sh --perline --memory 10000
mkdir plot
mkdir plot_zoomed
mv *_zoomed.pdf plot_zoomed/
mv *.pdf plot/
cat raw/categorical-100940-both_sexes-100940.tsv.bgz.out | cut -f 1 | tail -n +2  > snps.list
export LC_ALL=C
cat snps.list | awk -F ":" '{print "cat <(echo -e \"phenotype\\tp_val\\tbeta\\tbetapercentile\") <(grep -w "$0" raw/*.out | cut -d \"/\" -f 2 | tr \":\" \"\\t\"  | cut -f 1,4,5,6 | sed \"s/.tsv.bgz.out//\" | sort -k2,2n ) > "$1"_"$2".tsv"}' > Out2tsv.sh
source Out2tsv.sh 
