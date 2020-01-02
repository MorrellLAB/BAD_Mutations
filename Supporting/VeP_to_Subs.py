#!/usr/bin/env python
"""Supporting script for BAD_Mutations that will convert a report from the
Ensembl Variant Effect Predictor (VeP) to the substitutions file needed for
BAD_Mutations. Assumes that VeP report is the "tabular" report. Takes three
arguments:
    1) VeP report (gzipped)
    2) Output file for "long substitutions" used in the "predict" step
    3) Output directory for per-gene substitutions files
"""

import sys
import os
import gzip
import tempfile

try:
    vep = sys.argv[1]
    subs_out = sys.argv[2]
    subs_dir = sys.argv[3]
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

# Check that we can write to the output directory
try:
    subs_out_dn = os.path.abspath(
        os.path.expanduser(
            os.path.dirname(subs_out)))
    t = tempfile.TemporaryFile(dir=subs_out_dn)
    t.close()
except (OSError, PermissionError) as e:
    sys.stderr.write(
        'Error: The output directory for the long substitutions file is not '
        'writeable or does not exist.\n')
    exit(2)

try:
    subs_dir_fp = os.path.abspath(os.path.expanduser(subs_dir))
    t = tempfile.TemporaryFile(dir=subs_dir_fp)
    t.close()
except (OSError, PermissionError) as e:
    sys.stderr.write(
        'Error: The output directory for the per-gene substitutions files is '
        'not writeable or does not exist.\n')
    exit(3)

# Start writing the files. We will store the data in a dictionary so that we
# can make sure that the substitutions are grouped by transcript before writing
# to disk.
SNPID = '#Uploaded_variation'
TXID = 'Feature'
PROTPOS = 'Protein_position'
AA_STATE = 'Amino_acids'
subs = {}
with gzip.open(vep, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#Uploaded'):
            header = line.strip().split('\t')
            # A bit ugly, but slice up the columns
            save_cols = {}
            for i, h in enumerate(header):
                if h == SNPID:
                    save_cols['snp_id'] = i
                if h == TXID:
                    save_cols['tx_id'] = i
                if h == PROTPOS:
                    save_cols['prot_pos'] = i
                if h == AA_STATE:
                    save_cols['aa_state'] = i
                if h == 'Consequence':
                    save_cols['impact'] = i
        else:
            tmp = line.strip().split('\t')
            # We want to save only the missense varaints
            impact = tmp[save_cols['impact']]
            if impact == 'missense_variant':
                snpid = tmp[save_cols['snp_id']]
                txid = tmp[save_cols['tx_id']]
                cds = tmp[save_cols['prot_pos']]
                aa = tmp[save_cols['aa_state']]
                altaa = aa.split('/')[1]
                # Key them on transcript
                if txid in subs:
                    subs[txid].append((snpid, cds, altaa))
                else:
                    subs[txid] = [(snpid, cds, altaa)]
            else:
                continue

# Now, for each transcript:
lsubs_h = open(os.path.abspath(os.path.expanduser(subs_out)), 'w')
for tx in sorted(subs):
    # Open a handle to the output directory
    gsub_fp = os.path.join(subs_dir_fp, tx + '.subs')
    gsub_h = open(gsub_fp, 'w')
    for s in subs[tx]:
        # Write the subs
        gsub_h.write(s[1] + '\t' + s[0] + '\n')
        lsubs_h.write('\t'.join([tx, s[1], s[2], s[0]]) + '\n')
    gsub_h.close()
lsubs_h.close()
