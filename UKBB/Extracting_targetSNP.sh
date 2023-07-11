rm -fR temp/
mkdir -p temp/score
ls $2/*.bgz | awk -F "/" -v f=$1 -v p=$3 '{print "python "p"/UKBB/Run_ExtractScore.py --bed "f" "$0" > temp/score/"$NF".out"}'  > temp/Extracting.sh
sh temp/Extracting.sh
cat $(ls temp/score/*.out| head -1) | cut -f 1 | tail -n +2  > temp/snps.list
cat temp/snps.list | awk -F ":" '{print "cat <(echo -e \"phenotype\\tp_val\\tbeta\\tbetapercentile\") <(paste <(ls temp/score/*.out | rev| cut -f 1 -d \"/\"| rev |sed \"s/.tsv.bgz.out//\")  <(grep -w "$0" temp/score/*.out | cut -d \":\" -f 2- | tr \":\" \"\\t\" )  | sort -k2,2n ) > "$1"_"$2".tsv"}' > temp/Out2tsv.sh
bash temp/Out2tsv.sh
