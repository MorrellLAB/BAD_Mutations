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

#   In place, remove the ' <unknown description>' tag from the alignment
sed -i '.bak' 's/[[:space:]]<unknown description>//g' ${FASTA}

#   In place, remove the same from the tree, as well as bootstrap values
sed -i '.bak' 's/[[:space:]]<unknown description>//g' ${TREE}
#   '-r' (Extended regex) is a GNU option
if [ "$(uname)" == "Linux" ]
    then
        sed -i -r 's/[)][0-9]\.[0-9]+/\)/g' ${TREE}
        sed -i -r "s/'//g" ${TREE}
#   The equivalent is -E in MacOS
elif [ "$(uname)" == "Darwin" ]
    then
        sed -i '.bak' -E 's/[)][0-9]\.[0-9]+/\)/g' ${TREE}
        sed -i '.bak' -E "s/'//g" ${TREE}
fi
