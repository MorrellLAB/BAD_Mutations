#!/bin/bash
#   Shell script to unzip the CDS fasta files and format them as BLAST
#   databases
#   Modified for LRT package from Paul Hoffman's scipt.

set -e
set -u
set -o pipefail

#   First argument is path to makeblastdb
MAKE_BLAST_DB="$1"
#   Filename as argument
CDS="$2"

#   What is today's date?
YMD=$(date +%Y%m%d)

#   Build a new filename, we replace the .gz with nothing
NEW_FILE="${CDS/.gz/}"
#   Ungzip the files, and drop it into the same directory
gzip -cd "${CDS}" > "${NEW_FILE}"
#   Create a log file
LOG="${NEW_FILE}".${YMD}_makeblastdb_log
#   And an error file
ERR="${NEW_FILE}".${YMD}_makeblastdb_err
#   Make BLAST databases out of each of the files
#   We print stdout and stderr since we can handle these in python, but we
#   would also like them to be saved on disk for reference.
echo "Log file is stored in $LOG"
echo "Error file is stored in $ERR"
#	We have to remove empty sequences from the databases, as they cause FASTA
#	parse errors. Set the record separator to > and the field separator to
#	newline. If the second field exists, print it out, otherwise, skip it.
awk 'BEGIN { RS=">"; FS="\n"; ORS="" } $2 { print ">"$0}' "${NEW_FILE}" > \
	tmp.fasta && \
	mv tmp.fasta "${NEW_FILE}"
"${MAKE_BLAST_DB}" -in "${NEW_FILE}" -dbtype nucl > >(tee "$LOG") 2> >(tee "$ERR" >&2)
