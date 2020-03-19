#!/bin/bash
#   Written by Thomas Kono, based off Prank_Align.sh from Paul Hoffman
#   Shell script to run Clustal-omega on an input FASTA sequence then run
#   fasttree to generate a Newick tree

set -e
set -u
set -o pipefail

CLUSTALO="${1}"
INPUT="${2}"
OUTPUT="${3}"
FASTTREE="${4}"

# Run the alignment
"${CLUSTALO}" \
    --infmt=fasta \
    -t Protein \
    --full-iter \
    --threads=1 \
    --outfmt=fasta \
    --force \
    -i "${INPUT}" \
    -o "${OUTPUT}"

# then make a tree
"${FASTTREE}" "${OUTPUT}" > "${OUTPUT/.fasta/.tre}"
