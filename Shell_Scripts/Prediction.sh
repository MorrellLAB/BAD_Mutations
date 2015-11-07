#!/bin/bash
#   Written by Thomas Kono
#   Script to run HyPhy to predict codons

set -e
set -u
set -o pipefail

HYPHY=$1
ALN=$2
SUBS=$3
