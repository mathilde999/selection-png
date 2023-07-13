#!/bin/bash
trial=$1
RANDOM=$$
mkdir trial_${trial}

# convert ped to bed format
plink --file all_biSNP_PASS_maf_unrelated_VIFLDPruned --make-bed --out all_biSNP_PASS_maf_unrelated_VIFLDPruned

cd trial_${trial}

touch log_${trial}.out
for K in {2..10}
do
 admixture -s ${RANDOM} --cv=100 all_biSNP_PASS_maf_unrelated_VIFLDPruned.bed $K -j4|grep -w CV>> log_${trial}.out
done