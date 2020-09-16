#!/bin/bash
#   Written by Thomas Kono
#   Script to run HyPhy to predict codons
#	HYPHY script written by Justin C. Fay

set -e
set -u
set -o pipefail

HYPHY="$1"
PREDICTION_SCRIPT="$2"
INPUT="$3"
OUTPUT="$4"

${HYPHY} ${PREDICTION_SCRIPT} <<< ${INPUT} > ${OUTPUT}
