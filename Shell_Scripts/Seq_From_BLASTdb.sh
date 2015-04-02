#!/bin/bash
#   Small shell script to retrieve the FASTA sequence from 

set -e
set -u
set -o pipefail

BLASTDBCMD=$1

#   our arguments
DATABASE=$2
SEQID=$3

$BLASTDBCMD -entry $SEQID -db $DATABASE
