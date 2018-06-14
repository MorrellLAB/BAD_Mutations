#!/bin/bash
#   Written by Li Lei in 20160502 in st. Paul
#   Shell script to extract the lines without NOSNP and do filtering to determine which one is deleterious mutations.

set -e
set -u
set -o pipefail
filename=/home/morrellp/llei/Deleterious_mutation_project/LTR_BAD_mutation/A_thaliana_BAD_Mutation/out_Ath/prediction/prediction_file.list
for file in $(cat $filename)
    do
       name=$(echo $file |cut -d/ -f 10 |cut -d. -f 1)
       grep -v "NOSNP" $file|grep -v "givenTree" |grep -v "global" |grep -v "Alignment order:" |grep -v "Total"|sed -e '1,4d;$d'|sed -e '$d' >/home/morrellp/llei/Deleterious_mutation_project/LTR_BAD_mutation/A_thaliana_BAD_Mutation/out_Ath/extract_prediction/${name}_extract
done
