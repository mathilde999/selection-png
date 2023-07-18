
res_folder=$1

cat ingroup.list | while read name; do
head -n1 $res_folder/observ/chr9.obs.$name.txt > $res_folder/observ/obs.$name.txt
for i in {1..4} {6..10} {12..14} {17..19} 22; do 
tail -n+2 $res_folder/observ/chr$i.obs.$name.txt >> $res_folder/observ/obs.$name.txt 
done && rm $res_folder/observ/chr*.obs.$name.txt
done