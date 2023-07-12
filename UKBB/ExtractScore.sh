rm -fR temp/
mkdir -p temp/score
sourcefolder=$(dirname -- "$0")
ls $2/*.bgz | \
  awk -F "/" -v f=$1 -v p=$sourcefolder \
    '{print "python "p"/Run_ExtractScore.py --bed "f" "$0" > temp/score/"$NF".out"}' >  temp/Extracting.sh
sh temp/Extracting.sh
cat $(ls temp/score/*.out| head -1) | cut -f 1 | tail -n +2  > temp/snps.list
mkdir -p scores
cat temp/snps.list | \
  awk -F ":" '{print "cat <(echo -e \"phenotype\\tpval_EUR\\tbeta_EUR\") \
    <(paste <(ls temp/score/*.out | rev| cut -f 1 -d \"/\"| rev |sed \"s/.tsv.bgz.out//\")  \
    <(grep -w "$0" temp/score/*.out | cut -d \":\" -f 2- | tr \":\" \"\\t\" |cut -f 3-)  | \
    sort -k2,2n ) > scores/"$1"_"$2".tsv"}' > \
  temp/Out2tsv.sh
bash temp/Out2tsv.sh
rm -fR temp