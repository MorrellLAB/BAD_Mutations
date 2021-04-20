#!/usr/bin/env python
"""Supporting script for BAD_Mutations that will convert a report from the
Ensembl Variant Effect Predictor (VeP) to the substitutions file needed for
BAD_Mutations. Assumes that VeP report is the "tabular" report. Takes three
arguments:
    1) VeP report (gzipped)
    2) List of primary transcripts listed one per line
"""

import sys
import os
import gzip
import re

try:
    vep = sys.argv[1]
    primary_trans = sys.argv[2]
except IndexError:
    sys.stderr.write(__doc__)
    exit(1)

try:
    fp = os.path.abspath(os.path.expanduser(vep))
    h = open(vep, 'r')
    h.close()
except (OSError, PermissionError) as e:
    sys.stderr.write(
        'Error: The input VeP text file could not be read. Either the file '
        'does not exist or the permissions are not open.\n')
    exit(4)

try:
    fp = os.path.abspath(os.path.expanduser(primary_trans))
    h = open(primary_trans, 'r')
    h.close()
except (OSError, PermissionError) as e:
    sys.stderr.write(
        'Error: The input primary transcript text file could not be read. '
        'Either the file does not exist or the permissions are not open.\n')
    exit(4)

#   Read in the list of primary transcripts.
prime_trans = []
with gzip.open(primary_trans, 'rt') as f:
    for line in f:
        tmp = line.strip()
        prime_trans.append(tmp)
        uniq_trans = prime_trans

# Start writing the files. We will store the data in a dictionary so that we
# can make sure that the substitutions are grouped by transcript before writing
# to disk.

# Grab header line from combined report
comment_lines = []
header_line = ""

subs = []
with gzip.open(vep, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            tmp = line.strip()
            comment_lines.append(tmp)
        elif line.startswith('#Uploaded'):
            header_line = line.strip().split('\t')
        else:
            tmp = line.strip().split('\t')
            for lines in tmp:
                snpid = tmp[0]
                tmpid = re.split("_", snpid)
                tmpid = list([tmpid[0], tmpid[1]])
                replace_id = '_'.join(tmpid)
                tmp[0] = replace_id
                txid = tmp[4]
                impact = tmp[6]
                # if txid in uniq_trans:
                #    print(txid)
                # We want to save only the missense varaints
                if impact == 'missense_variant':
                    # if uniq_trans.count(txid) > 0:
                    if txid in uniq_trans:
                        subs.append(tmp)
                    else:
                        continue
                else:
                    continue

# Now, for each transcript:
for line in comment_lines:
    print(line)
print(*header_line, sep = "\t")

unique_subs = []
for item in subs:
    if item not in unique_subs:
        unique_subs.append(item)
for line in unique_subs:
    out_string = '\t'.join(line)
    print(out_string)

sys.stderr.write('Finished converting.\n')
