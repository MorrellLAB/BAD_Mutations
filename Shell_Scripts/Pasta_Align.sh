#!/bin/bash
#   Written by Thomas Kono, based off Prank_Align.sh from Paul Hoffman
#   Shell script to run Pasta on an input FASTA sequence

set -e
set -u
set -o pipefail

#   Path to the Pasta executable as an argument
PASTA=$1
#   Input sequence as an argument
INPUT=$2
#   Temp Dir
TEMP_DIR=$3
#   Job name
JOBNAME=$4

#   We create a variable to give the job a unique name
JOBNAME=$(date +%Y%m%d_%H%M%S%N)

$PASTA \
    -d dna\
    --no-return-final-tree-and-alignment\
    --num-cpus=1\
    --job=$JOBNAME\
    --temporaries=${TEMP_DIR}\
    -i $INPUT\
    -o ${TEMP_DIR}
