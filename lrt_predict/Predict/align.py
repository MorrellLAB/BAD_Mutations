#!/usr/bin/env python

#   A script that performs the alignment and LRT prediction

#   Import standard library modules here
import tempfile
import subprocess
import os
import time
import re

#   Import Biopython modules here for sequence handling
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq

#   Import our helper scripts here
from lrt_predict.General import set_verbosity
from lrt_predict.General import check_modules


class PastaAlign(object):
    def __init__(
            self,
            pasta_path,
            clustalo_path,
            fasttree_path,
            unaligned_sequences,
            query_sequence,
            verbose):
        self.mainlog = set_verbosity.verbosity('Pasta_Align', verbose)
        #   This is file-like object
        self.input_seq = unaligned_sequences
        #   This will be populated with sequences for back-translation
        self.input_dict = {}
        self.query = query_sequence
        self.pasta_path = check_modules.check_executable(pasta_path)
        self.clustalo_path = check_modules.check_executable(clustalo_path)
        self.fasttree_path = check_modules.check_executable(fasttree_path)
        self.protein_input = None
        self.aln_out = None
        self.tree_out = None
        self.final_aln = None
        return

    def prepare_sequences(self):
        """Prepares the CDS sequences for alignment with Pasta. Checks if any
        sequences are not multiples of 3, and appends N if not. Translates
        the nucleotide sequences to amino acids for protein alignment with
        Pasta, then removes the trailing stop codon, if it is present."""
        #   First, read in the sequences and check their length
        input_seqs = list(SeqIO.parse(self.input_seq.name, 'fasta'))
        #   Start accumulating translated sequences to write into the
        #   alignment input file.
        tl_seqs = []
        for i in input_seqs:
            #   Fix the names for HyPhy. Names must be alphanumeric, with
            #   underscores
            fixed_name = re.sub('[^0-9a-zA-Z]', '_', i.id)
            #   If the length is not a multiple of 3, then we have to add Ns to
            #   make it so.
            if len(i) % 3 != 0:
                to_add = 3 - (len(i) % 3)
                #   Tack on the original sequence, with some appended Ns so we
                #   can recreate the nucleotide alignment later.
                self.input_dict[fixed_name] = str(i.seq) + to_add*'N'
                self.mainlog.debug(
                    'Length of sequence ' + i.id + ' is not a mulitple of 3. ' +
                    'Adding ' + str(to_add) + ' Ns to the end.'
                    )
                #   Create the new sequence. Pasta chokes on ambiguous
                #   amino acids that aren't X, so we replace them all with X.
                fixed_seq = Seq(str(i.seq) + to_add*'N').translate()
                sub_seq = Seq(
                    re.sub(
                        'B|Z|J|O|U',
                        'X',
                        str(fixed_seq),
                        re.I))
                new_seq = SeqRecord.SeqRecord(
                    sub_seq,
                    id=i.id,
                    description='')
            else:
                self.input_dict[fixed_name] = str(i.seq)
                sub_seq = Seq(
                    re.sub(
                        'B|Z|J|O|U',
                        'X',
                        str(i.seq.translate()),
                        re.I))
                new_seq = SeqRecord.SeqRecord(
                    sub_seq,
                    id=fixed_name,
                    description='')
            self.mainlog.debug(new_seq.id + '\t' + str(new_seq.seq))
            tl_seqs.append(new_seq)
        self.mainlog.debug('Number of species aligned: ' + str(len(tl_seqs)))
        #   Then, we have to iterate through the translated sequences and
        #   check for sequences ending in stop codons. Pasta hates these, so
        #   we will prune them. We also check for those with internal stop
        #   codons, and skip those.
        fixed_tl_seqs = []
        for i in tl_seqs:
            #   If we find an internal stop codon
            if re.match(r'.+\*[^$]', str(i.seq)):
                self.mainlog.debug(
                    'Sequence ' + i.id + ' has internal stop. Skipping.'
                    )
                continue
            else:
                #   Otherwise, just strip the stop codon off the end and save
                if i.seq.endswith('*'):
                    fixed_tl_seqs.append(
                        SeqRecord.SeqRecord(i.seq[:-1], id=i.id, description='')
                        )
                else:
                    fixed_tl_seqs.append(i)
        #   Then, we open another temporary file to hold our amino acid
        #   sequences.
        self.protein_input = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_PastaInput_',
            suffix='.fasta')
        #   And write the protein sequences into it
        SeqIO.write(fixed_tl_seqs, self.protein_input, 'fasta')
        self.protein_input.flush()
        # Return the number of sequences to align here
        return len(fixed_tl_seqs)

    def back_translate(self):
        """Back-translates from amino acid to nucleotide, using the original
        input sequences as a guide to avoid ambiguity. Assumes that a non-gap
        character in the amino acid alignment will be faithfully represented
        by a triplet in the source sequence, and will not check identity of
        translated codons."""
        #   Iterate through the aligned protein file, and start rebuilding the
        #   original nucleotide sequence.
        aln_prot = SeqIO.parse(self.aln_out, 'fasta')
        bt_seqs = []
        for rec in aln_prot:
            #   Check if we have name mismatch. This *shouldn't* happen, but
            #   this should stop some errors
            if rec.id in self.input_dict:
                in_seq = self.input_dict[rec.id]
                rebuilt_seq = ''
                pos = 0
                #   Now, iterate through the amino acid sequence and rebuild
                #   the nucleotide sequence.
                for aa in rec.seq:
                    #   If we have a gap, put in three gaps
                    if aa == '-':
                        rebuilt_seq += '---'
                    else:
                        rebuilt_seq += in_seq[pos:pos+3]
                        pos += 3
            else:
                continue
            #   Then, put in a new SeqRecord with the sequence and the name, so
            #   we can write a fasta file.
            bt_seqs.append(
                SeqRecord.SeqRecord(
                    Seq(rebuilt_seq),
                    id=rec.id,
                    description='')
                )
        #   And create a new temporary file for them
        final_seqs = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_BackTranslated_',
            suffix='.fasta',
            delete=False)
        SeqIO.write(bt_seqs, final_seqs, 'fasta')
        final_seqs.flush()
        self.final_aln = final_seqs.name
        return

    def pasta_align(self):
        """Align the amino acid sequences with Pasta."""
        #   Get the base directory of the LRT package
        lrt_path = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
        #   Then build the path to the pasta script
        pasta_script = os.path.join(
            lrt_path,
            'Shell_Scripts',
            'Pasta_Align.sh')
        #   Pasta expects a directory for output. We use the system temp dir
        pasta_out = tempfile.gettempdir()
        #   We make a job name from the time in microseconds
        #   This should be good enough...
        pasta_job = 'pastajob_' + '%.6f' % time.time()
        #   Create the command line
        cmd = [
            'bash',
            pasta_script,
            self.pasta_path,
            self.protein_input.name,
            pasta_out,
            pasta_job]
        self.mainlog.debug(' '.join(cmd))
        #   Then, we'll execute it
        p = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate()
        #   Then, build the output name
        #   The structure of the Pasta output files is
        #       Jobname.marker001.Unaligned_name.aln
        #       Jobname.tre
        aln_out = ''.join([
            pasta_out,
            os.path.sep,
            pasta_job,
            '.marker001.',
            self.protein_input.name.split(os.path.sep)[-1].replace('.fasta', ''),
            '.aln'])
        tree_out = ''.join([
            pasta_out,
            os.path.sep,
            pasta_job,
            '.tre'])
        #   And save the paths to these files as class variables
        self.aln_out = aln_out
        self.tree_out = tree_out
        return (out, err)

    def clustalo_align(self):
        """Align the amino acid sequences with Clustal-omega. We invoke this
        function in the case where a sequence only finds one homologue, which
        means we have a pairwise alignment."""
        #   Get the base directory of the LRT package
        lrt_path = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
        #   Then build the path to the pasta script
        clustalo_script = os.path.join(
            lrt_path,
            'Shell_Scripts',
            'Clustalo_Align.sh')
        # Clustal-omega needs an output filename
        clustalo_out = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_Clustalo_Out_',
            suffix='.fasta',
            delete=False)
        #   Create the command line
        cmd = [
            'bash',
            clustalo_script,
            self.clustalo_path,
            self.protein_input.name,
            clustalo_out.name,
            self.fasttree_path]
        self.mainlog.debug(' '.join(cmd))
        #   Then, we'll execute it
        p = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate()
        aln_out = clustalo_out.name
        tree_out = clustalo_out.name.replace('.fasta', '.tre')
        #   And save the paths to these files as class variables
        self.aln_out = aln_out
        self.tree_out = tree_out
        # Seek back to the beginning of the output file so it can be read again
        clustalo_out.seek(0)
        return (out, err)

    def sanitize_outputs(self):
        """Remove troublesome characters from the PASTA output files that
        cause HyPhy to crash, or the Biopython Newick parser to fail."""
        #   Get the base directory of the LRT package
        lrt_path = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
        #   Then build the path to the pasta script
        santize_script = os.path.join(
            lrt_path,
            'Shell_Scripts',
            'Prepare_HyPhy.sh')
        #   Create the command line
        cmd = [
            'bash',
            santize_script,
            self.final_aln,
            self.tree_out]
        #   Run it
        self.mainlog.debug('Sanitizing alignment in ' + self.final_aln)
        self.mainlog.debug('Sanitizing tree in ' + self.tree_out)
        p = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate()
        return(out, err)
