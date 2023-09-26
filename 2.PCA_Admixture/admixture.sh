#!/bin/bash
file=$1
trial=$2
K=$3

mkdir -p trial_${trial}

cd trial_${trial}

admixture -s ${trial} --cv=100 ../${file}.bed $K |grep -w CV>> log_${trial}.out

