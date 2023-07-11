tabix -R $1 $2 | cut -f 1,2 | awk '{print $1"\t"$2-1"\t"$2}' > Present.bed
cat <(cut -f 1,3 topsnp_gr37.bed ) <(cut -f 1,3 Present.bed ) | sort | uniq -c | awk '{print $1"\t"$2"\t"$3}' | grep -vw "^2" | cut -f 2,3 > notpresent.snps
cat notpresent.snps | awk -v q="'" -v q=$1 '{print "tabix "q" "$1":"$2-1000"-"$2+1000" | cut -f 2| awk "q"{print \""$1"\\t\"$1-1\"\\t\"$1\"\\t\"($1-"$2")^2}"q"| sort -k 4,4n|head -n 1 | cut -f -3"}' > closest.sh
sh closest.sh > closest.bed
cat Present.bed closest.bed| sort -k1,1n -k2,2  > All.bed