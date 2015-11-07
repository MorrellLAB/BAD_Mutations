#!/bin/bash
#   Written by Thomas Kono
#   Script to run HyPhy to predict codons
#	HYPHY script written by Justin C. Fay

set -e
set -u
set -o pipefail

HYPHY=$1
INPUT=$2
OUTPUT=$3

PREDICTION_SCRIPT="LRT.hyphy"

${HYPHY} ${PREDICTION_SCRIPT} <<< ${INPUT} > ${OUTPUT}
