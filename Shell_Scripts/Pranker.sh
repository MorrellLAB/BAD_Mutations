#!/bin/bash
#   Written by Paul Hoffman
#   Shell script to run PRANK on an input FASTA sequence
#   To be used until a Python script for PRANK works

#   Path to the prank executable as an argument
PRANK=$1
#   Input sequence as an argument
INPUT=$2
#   Output name as an argumnet
OUTPUT=$3

$PRANK -d=$INPUT -o=$OUTPUT'_translate' -translate -F -showtree
$PRANK -d=$INPUT -o=$OUTPUT'_codon' -codon -F -showtree
