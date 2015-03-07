#!/bin/bash

#   Shell script to unzip the CDS fasta files and format them as BLAST
#   databases
#   Modified for LRT package from P. Hoffman's scipt.

#   Gain access to the NCBI executables
#   This line is specific to the Minnesota Supercomputing Institute
module load ncbi_blast+

#   Filename as argument
CDS=$1

#   Build a new filename, we replace the .gz with nothing
new_file=${CDS/.gz/}
#   Ungzip the files, and drop it into the same directory
gzip -cd $x > ${new_file}
#   Make BLAST databases out of each of the files
makeblastdb -in ${new_file} -dbtype nucl
