#!/bin/bash

#   Shell script to unzip the CDS fasta files and format them as BLAST
#   databases
#   Modified for LRT package from P. Hoffman's scipt.

#   Gain access to the NCBI executables
#   This line is specific to the Minnesota Supercomputing Institute
#module load ncbi_blast+
#MAKE_BLAST_DB=`which makeblastdb`
#   Or, if we don't have module, specify the full path to the makeblastdb
#   executable
MAKE_BLAST_DB=/soft/ncbi_blast+/2.2.29/bin/makeblastdb

#   Filename as argument
CDS=$1

#   Build a new filename, we replace the .gz with nothing
new_file=${CDS/.gz/}
#   Ungzip the files, and drop it into the same directory
gzip -cd $CDS > ${new_file}
#   Create a log file
LOG=${new_file}.makeblastdb_log
#   And an error file
ERR=${new_file}.makeblastdb_err
#   Make BLAST databases out of each of the files
${MAKE_BLAST_DB} -in ${new_file} -dbtype nucl > $LOG 2> $ERR
