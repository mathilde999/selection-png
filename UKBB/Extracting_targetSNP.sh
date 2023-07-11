mkdir score
ls *.bgz | awk -F "/" '{print "python UKBB/Run_ExtractScore.py --bed All.bed "$0" > score/"$NF".out"}'  > Extracting.sh
sh Extracting.sh
cat raw/categorical-100940-both_sexes-100940.tsv.bgz.out | cut -f 1 | tail -n +2  > snps.list
cat snps.list | awk -F ":" '{print "cat <(echo -e \"phenotype\\tp_val\\tbeta\\tbetapercentile\") <(grep -w "$0" raw/*.out | cut -d \"/\" -f 2 | tr \":\" \"\\t\"  | cut -f 1,4,5,6 | sed \"s/.tsv.bgz.out//\" | sort -k2,2n ) > "$1"_"$2".tsv"}' > Out2tsv.sh
source Out2tsv.sh 
