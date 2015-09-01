#!/usr/bin/env python

#   A script that performs the alignment and LRT prediction

#   Import standard library modules here
import tempfile
import subprocess
import os
import time

#   Import Biopython modules here for sequence handling
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq

#   Import our helper scripts here
from lrt_predict.General import parse_input
from lrt_predict.General import set_verbosity
from lrt_predict.General import check_modules


class PastaAlign:
    def __init__(self, unaligned_sequences, query_sequence, verbose):
        self.mainlog = set_verbosity.verbosity('Pasta_Align', verbose)
        #   This is file-like object
        self.input_seq = unaligned_sequences
        #   This will be populated with sequences for back-translation
        self.input_dict = {}
        self.query = query_sequence
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
            #   If the length is not a multiple of 3, then we have to add Ns to
            #   make it so.
            if len(i) % 3 != 0:
                to_add = 3 - (len(i) % 3)
                #   Tack on the original sequence, with some appended Ns so we
                #   can recreate the nucleotide alignment later.
                self.input_dict[i.id] = str(i.seq) + to_add*'N'
                new_seq = SeqRecord.SeqRecord(
                    Seq(str(i.seq) + to_add*'N').translate(),
                    id=i.id)
            else:
                self.input_dict[i.id] = str(i.seq)
                new_seq = SeqRecord.SeqRecord(
                    i.seq.translate(),
                    id=i.id)
            self.mainlog.debug(new_seq.id + '\t' + str(new_seq.seq))
            tl_seqs.append(new_seq)
        self.mainlog.debug(len(tl_seqs))
        #   Then, we have to iterate through the translated sequences and
        #   check for sequences ending in stop codons. Pasta hates these, so
        #   we will prune them.
        fixed_tl_seqs = []
        for i in tl_seqs:
            if i.seq.endswith('*'):
                fixed_tl_seqs.append(SeqRecord.SeqRecord(i.seq[:-1], id=i.id))
            else:
                fixed_tl_seqs.append(i)
        #   Then, we open another temporary file to hold our amino acid
        #   sequences.
        self.protein_input = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='LRTPredict_PastaInput_',
            suffix='.fasta')
        #   And write the protein sequences into it
        SeqIO.write(fixed_tl_seqs, self.protein_input, 'fasta')
        self.protein_input.flush()
        return

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
            #   Then, put in a new SeqRecord with the sequence and the name, so
            #   we can write a fasta file.
            bt_seqs.append(SeqRecord.SeqRecord(
                Seq(rebuilt_seq),
                id=rec.id))
        #   And create a new temporary file for them
        final_seqs = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='LRTPredict_BackTranslated_',
            suffix='.fasta')
        SeqIO.write(bt_seqs, final_seqs, 'fasta')
        final_seqs.flush()
        return final_seqs

    def pasta_align(self):
        """Align the amino acid sequences with Pasta."""
        #   Get the base directory of the LRT package
        lrt_path = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
        #   Then build the path to the pasta script
        pasta_script = os.path.join(
            lrt_path,
            'Shell_Scripts',
            'Pasta_Align.sh')
        #   Check for the presence of the pasta executable
        pasta_path = check_modules.check_executable('run_pasta.py')
        #   Pasta expects a directory for output. We use the system temp dir
        pasta_out = tempfile.gettempdir()
        #   We make a job name from the time in microseconds
        #   This should be good enough...
        pasta_job = 'pastajob_' + '%.6f' % time.time()
        #   Create the command line
        cmd = [
            'bash',
            pasta_script,
            pasta_path,
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
