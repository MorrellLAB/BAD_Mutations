#!/bin/bash

#   Shell script to unzip the CDS fasta files and format them as BLAST
#   databases
#   Modified for LRT package from Paul Hoffman's scipt.

#   First argument is path to makeblastdb
MAKE_BLAST_DB=$1
#   Filename as argument
CDS=$2

#   What is today's date?
YMD=$(date +%Y%m%d)

#   Build a new filename, we replace the .gz with nothing
new_file=${CDS/.gz/}
#   Ungzip the files, and drop it into the same directory
gzip -cd $CDS > ${new_file}
#   Create a log file
LOG=${new_file}.${YMD}_makeblastdb_log
#   And an error file
ERR=${new_file}.${YMD}_makeblastdb_err
#   Make BLAST databases out of each of the files
#   We print stdout and stderr since we can handle these in python, but we
#   would also like them to be saved on disk for reference.
echo "Log file is stored in $LOG"
echo "Error file is stored in $ERR"
${MAKE_BLAST_DB} -in ${new_file} -dbtype nucl > >(tee $LOG) 2> >(tee $ERR >&2)
