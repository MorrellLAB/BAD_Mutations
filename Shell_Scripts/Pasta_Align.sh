#!/bin/bash
#   Written by Thomas Kono, based off Prank_Align.sh from Paul Hoffman
#   Shell script to run Pasta on an input FASTA sequence

set -e
set -u
set -o pipefail

#   We have to set this environment variable to increase the java heap space,
#   else it runs out of memory someitmes and fails to finish an alignment.
export _JAVA_OPTIONS="-Xmx4g"

#   Path to the Pasta executable as an argument
PASTA=$1
#   Input sequence as an argument
INPUT=$2
#   Temp Dir
TEMP_DIR=$3
#   Job name
JOBNAME=$4

$PASTA \
    -d protein\
    --no-return-final-tree-and-alignment\
    --num-cpus=1\
    --job=$JOBNAME\
    --iter-limit=5\
    --temporaries=${TEMP_DIR}\
    -i $INPUT\
    -o ${TEMP_DIR}
