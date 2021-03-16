#!/usr/bin/env python3
"""Add SNP name to annotated variants. Need to read the combined reports and
substitutions files and add the SNP position to the combined report.
The inputs are taken in this order as arguments:
    1) Combined reports from running "compile" after BAD_Mutations "predict"
    2) "Long substitutions" file preferably created by running "VeP_to_Subs.py"

Code by Sam Hamman & Peter Morrell, borrowing from Tom Kono's code here:
https://github.com/MorrellLAB/BAD_Mutations/blob/dev/Supporting/VeP_to_Subs.py
"""

import sys
import os

try:
    combined = sys.argv[1]
    subs_file = sys.argv[2]
except IndexError:
    sys.stderr.write(__doc__)
    exit(1)

try:
    fp = os.path.abspath(os.path.expanduser(combined))
    h = open(combined, 'r')
    h.close()
except (OSError, PermissionError) as e:
    sys.stderr.write('Error: The combined report text file could not be read. Either the file \
        does not exist or the permissions are not open.\n')
    exit(4)

try:
    fp = os.path.abspath(os.path.expanduser(subs_file))
    h = open(subs_file, 'r')
    h.close()
except (OSError, PermissionError) as e:
    sys.stderr.write('Error: The "long substitutions" text file could not be read. Either the file \
    does not exist or the permissions are not open.\n')
    exit(4)

# Start writing the files. We will store the data in a dictionary so that we
# can make sure that the substitutions are grouped by transcript before writing
# to disk.
subs = {}
with open(subs_file, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            tmp = line.strip().split('\t')
            # We need every line of the file
            # Need to remove version information in the transcript name
            # Will also need to make transcript ID lowercase to match combined
            # report
            txid = tmp[0].replace('.v2.1', '')
            txid = txid.lower()  # transcript/gene id
            aa_pos = tmp[1]      # cds position
            aa = tmp[2]
            # Need to remove nucleotide information in the position
            # pos = tmp[3].replace([-4:], '')
            pos = tmp[3][:-4]

            key = (txid, aa_pos)
            subs.update({key: pos})

# Begin empty dict for final output
# Keys are pairs of gene id and position on said gene
# Values are the data values from the combined report
report = {}

# Grab header line from combined report
header_line = ""

with open(combined, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        elif line.startswith("VariantID"):
            header_line = line
        else:
            tmp = line.strip().split('\t')
            # We need every line of the combined report
            txid = tmp[1]        # transcript/gene id
            aa_pos = tmp[2]      # cds position

            # Grabs the substitution for the current position
            key = (txid, aa_pos)
            tmp[0] = subs[key]

            # Adds the entire adjusted row of data under
            # the unique key pair to the final dicitonary
            report.update({key: tmp})

# Write to a temporary text file, grouped first by gene id then by position
# Also note, sorted lexicographically, not numerically (e.g. 123, 32, 89, etc)
print(header_line[:-1])  # Remove newline character to avoid double-newlining
for key in sorted(report.keys()):
    out_string = '\t'.join(report[key])
    print(out_string)

sys.stderr.write('Finished converting.')
