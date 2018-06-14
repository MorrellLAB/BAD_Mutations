#!/bin/bash
#   Written by Li Lei, 20160426, st.paul
#   Shell script to create the subs files

set -e
set -u
set -o pipefail
filename=at.both.snp.tsv
for geneid in $(cat 'target_gene_models.txt')
    do
    grep $geneid $filename |sort -k2,2n |cut -f 2|uniq|awk '{printf("%s\tSNP_%.0f\n", $0, 1*(NR-1))}' >./subs_files/${geneid}.subs
done
