#!/bin/bash
#   Script to remove spaces, bootstrap numbers, and other "unfriendly"
#   characters from the input phylogenetic tree and multiple sequence
#   alignment for HyPhy.
#   From Justin C. Fay

set -e
set -u
set -o pipefail

FASTA="$1"
TREE="$2"

#   '-r' (Extended regex) is a GNU option
if [ "$(uname)" == "Linux" ]
    then
        sed -i 's/[[:space:]]<unknown description>//g' ${FASTA}
        sed -i 's/[[:space:]]<unknown description>//g' ${TREE}
        sed -i -r 's/[)][0-9]\.[0-9]+/\)/g' ${TREE}
        sed -i -r "s/'//g" ${TREE}
        sed -i -r "s/\+//g" ${FASTA}
        sed -i -r "s/\+//g" ${TREE}
        sed -i -r "s/\./_/g" ${TREE}
        sed -i -r "s/\./_/g" ${FASTA}
        sed -i -e "s/:\([0-9]\)_/:\1\./g" ${TREE}
        sed -i -r '/^;+$/d' ${TREE}
#   The equivalent is -E in MacOS
elif [ "$(uname)" == "Darwin" ]
    then
        sed -i '.bak' 's/[[:space:]]<unknown description>//g' ${FASTA}
        sed -i '.bak' 's/[[:space:]]<unknown description>//g' ${TREE}
        sed -i '.bak' -E 's/[)][0-9]\.[0-9]+/\)/g' ${TREE}
        sed -i '.bak' -E "s/'//g" ${TREE}
        sed -i '.bak' -E "s/\+//g" ${FASTA}
        sed -i '.bak' -E "s/\+//g" ${TREE}
        # Remove lines with just semicolon
        sed -i '.bak' -E '/^;+$/d' ${TREE}
fi
