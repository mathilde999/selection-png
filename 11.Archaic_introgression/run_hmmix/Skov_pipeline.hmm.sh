
name=$1

res_folder=$2

hmmix train \
  -obs=$res_folder/observ/obs.$name.txt \
  -mutrates=$res_folder/all_chr.mutationrate.bed \
  -haploid \
  -out=$res_folder/trained/trained.$name.json 
  
hmmix decode \
  -obs=$res_folder/observ/obs.$name.txt \
  -mutrates=$res_folder/all_chr.mutationrate.bed \
  -param=$res_folder/trained/trained.$name.json \
  -admixpop=archaic_genoms/four_hg38/chr*.4_genomes.hg38.fixed_ref.masked.bcf \
  -haploid \
  -out=$res_folder/res/$name.decoded
