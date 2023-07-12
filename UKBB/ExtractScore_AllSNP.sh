ls ../*.bgz  |awk -v q="'" '{print "zcat",$0,"| head -n 1 | tr \"\\t\" \"\\n\" |awk "q"{print NR\"\\t\"$0\"\\t"$0"\"}"q"| grep \"pval_EUR\"| cut -f 1,3 "}' > extract_pvalcol_index.sh
sh extract_pvalcol_index.sh  > path_pval_eurindex.txt
awk '{print "zcat ",$2," | cut -f "$1" >","phenotypes/"substr($NF,4,length($NF)-11)".pval"}' path_pval_eurindex.txt  > extract_pval.sh
python ~/src/python3/SubmittingJobs/EBCRun.py --script extract_pval.sh --perline --log --sh --amd
paste <(cat <(zcat ../biomarkers-30600-both_sexes-irnt.tsv.bgz | cut -f 1,2,3,4 | tail -n +2 | echo -e "#CHROM\tPOS\tREF\tALT"  ) <(zcat ../biomarkers-30600-both_sexes-irnt.tsv.bgz | cut -f 1,2,3,4 | tail -n +2 )) <(cat <(ls phenotypes/*.pval |awk -F "/" '{print substr($NF,1,length($NF)-5)}'| paste -sd "\t" -) <(paste phenotypes/*.pval | tail -n +2 )) | head -n 3 > test.txt
mkdir phenotypes_1
mkdir phenotypes_2
mkdir phenotypes_3
mkdir phenotypes_4
ls phenotypes/*.pval > file_list.txt
head -n 500 file_list.txt | awk '{system ("mv "$0" phenotypes_1/")}'
tail -n +501 file_list.txt| head -n500| awk '{system ("mv "$0" phenotypes_2/")}'
tail -n +1001 file_list.txt| head -n500| awk '{system ("mv "$0" phenotypes_3/")}'
tail -n +1501 file_list.txt| awk '{system ("mv "$0" phenotypes_4/")}'
cat <(ls phenotypes_1/*.pval |awk -F "/" '{print substr($NF,1,length($NF)-5)}'| paste -sd "\t" -) <(paste phenotypes_1/*.pval | tail -n +2 )  > Phenotypes_1.pval
cat <(ls phenotypes_2/*.pval |awk -F "/" '{print substr($NF,1,length($NF)-5)}'| paste -sd "\t" -) <(paste phenotypes_2/*.pval | tail -n +2 )  > Phenotypes_2.pval
cat <(ls phenotypes_3/*.pval |awk -F "/" '{print substr($NF,1,length($NF)-5)}'| paste -sd "\t" -) <(paste phenotypes_3/*.pval | tail -n +2 )  > Phenotypes_3.pval
cat <(ls phenotypes_4/*.pval |awk -F "/" '{print substr($NF,1,length($NF)-5)}'| paste -sd "\t" -) <(paste phenotypes_4/*.pval | tail -n +2 )  > Phenotypes_4.pval
module load htslib
python ~/src/python3/SubmittingJobs/EBCRun.py --command 'paste <(cat <(zcat ../biomarkers-30600-both_sexes-irnt.tsv.bgz | cut -f 1,2,3,4 | tail -n +2 | echo -e "#CHROM\tPOS\tREF\tALT"  ) <(zcat ../biomarkers-30600-both_sexes-irnt.tsv.bgz | cut -f 1,2,3,4 | tail -n +2 )) Phenotypes_1.pval Phenotypes_2.pval Phenotypes_3.pval Phenotypes_4.pval | bgzip > pval_EUR.tsv.bgz' --log --sh --amd
 python ~/src/python3/SubmittingJobs/EBCRun.py --command 'gzip phenotypes_1/*.pval' --log --sh --amd
 python ~/src/python3/SubmittingJobs/EBCRun.py --command 'gzip phenotypes_2/*.pval' --log --sh --amd
 python ~/src/python3/SubmittingJobs/EBCRun.py --command 'gzip phenotypes_3/*.pval' --log --sh --amd
 python ~/src/python3/SubmittingJobs/EBCRun.py --command 'gzip phenotypes_4/*.pval' --log --sh --amd
python ~/src/python3/SubmittingJobs/EBCRun.py --command "gzip Phenotypes_*.pval " --log --sh --amd
