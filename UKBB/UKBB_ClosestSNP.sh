mkdir -p temp
tabix -R $1 $2 | cut -f 1,2 | awk '{print $1"\t"$2-1"\t"$2}' > temp/Present.bed
cat <(cut -f 1,3 $1 ) <(cut -f 1,3 temp/Present.bed ) | sort | uniq -c | awk '{print $1"\t"$2"\t"$3}' | grep -vw "^2" | cut -f 2,3 > temp/notpresent.snps
cat temp/notpresent.snps | awk -v q="'" -v f=$2 '{print "tabix "f" "$1":"$2-1000"-"$2+1000" | cut -f 2| awk "q"{print \""$1"\\t\"$1-1\"\\t\"$1\"\\t\"($1-"$2")^2}"q"| sort -k 4,4n|head -n 1 | cut -f -3"}' > temp/closest.sh
sh temp/closest.sh > temp/closest.bed
cat temp/Present.bed temp/closest.bed| sort -k1,1n -k2,2