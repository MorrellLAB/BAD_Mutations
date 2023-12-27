#!/bin/bash
#   Written by Li Lei in 20160502 in st. Paul
#   Shell script to extract the lines without NOSNP and do filtering to determine which one is deleterious mutations.

set -e
set -u
set -o pipefail
filename=/panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/out/list01/prediction_list01
for file in $(cat $filename)
    do
	       name=$(echo $file |cut -d/ -f 14 |cut -d. -f 1)
		          grep -v "NOSNP" $file|grep -v "givenTree" |grep -v "global" |grep -v "Alignment order:" |grep -v "Total"|sed -e '1,4d;$d'|sed -e '$d' >/panfs/roc/groups/9/morrellp/llei/Deleterious_mutation_project/LRT_BAD_mutation/A_thaliana_BAD_Mutation/Additional1000genes/extract_prediction/${name}_extract
				  done
