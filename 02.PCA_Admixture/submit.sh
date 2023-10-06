#!/bin/bash
file=$1
for K in {2..10}
do
    for trial in {1..3}
    do
        sbatch admixture.sh ${file} ${trial} ${K} 
   done
done
