#! /usr/bin/env python
####20150427 by Li Lei. The aim is to extract the target fasta file cccording to the target files from the a fasta file:
from Bio import SeqIO

fasta_file = "TAIR10_cds_20101214_updated.fasta" # Input fasta file
wanted_file = "target_gene_models.txt" # Input interesting sequence IDs, one per line
result_file = "target_cds.fasta" # Output fasta file

wanted = set()
with open(wanted_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], f, "fasta")
