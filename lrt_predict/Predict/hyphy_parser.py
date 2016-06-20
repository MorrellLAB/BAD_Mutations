#!/usr/bin/env python
"""Parse the HyPhy output and produce a report for each of the SNPs thave have
been predicted.

Originally written by Justin C. Fay, rewritten for general use in BAD_Mutations
by Thomas Kono."""

#   Import standard library modules here
import os
#   Import our helper scripts here
from lrt_predict.General import set_verbosity


class HyPhyParser(object):
    """Class to compile and evaluate the predictions for SNPs. Will read a list
    of individual HyPhy output files in a directory and pull out the relevant
    lines."""

    def __init__(self, outdir, verbose):
        self.preddir = outdir
        self.mainlog = set_verbosity.verbosity('HyPhy_Parser', verbose)
        return

    def get_prediction_files(self):
        """List the contents of the supplied directory and find all the files
        that contain predictions. These will end in '_Predictions.txt' - this
        is somewhat fragile, but the paths used throughout the package are
        consistent."""
        pred_dir_contents = os.listdir(self.preddir)
        pred_files = [
            fname
            for fname
            in pred_dir_contents
            if fname.endswith('_Predictions.txt')]
        return pred_files

    def parse_prediction(self, pred_file):
        """Parse a single prediction file. Return the lines that have prediction
        information for them, along with a gene identifier."""
        #   set a flag that says whether or not we are reading in alignment data
        #   This is somewhat messy, but since the number of species in each
        #   alignment is variable, it is the most straightforward way to do it.
        in_aln = False
        #   And an empty list to store the prediction info
        gene_preds = []
        with open(os.path.join(self.preddir, pred_file), 'r') as f:
            for line in f:
                #   The header for the actual predictions starts with 'Position'
                if line.startswith('Position'):
                    in_aln = True
                    #   Also start a counter for the CDS position we are
                    #   considering.
                    cds_pos = 0
                    continue
                #   Then, if we are in the alignment, check for predicton data
                if in_aln:
                    #   The prediction lines are ones that do NOT end in 'NOSNP'
                    if 'NOSNP' in line:
                        #   Check the third field in the 'NOSNP' line - if it is
                        #   not '-', then we increment the CDS position
                        #   counter.
                        if line.strip().split()[2] != '-':
                            cds_pos += 1
                        continue
                    else:
                        #   We first check if we've run off the edge of the
                        #   alignment data. The line right after the alignment
                        #   data starts with 'Alignment'
                        if line.startswith('Alignment'):
                            return gene_preds
                        else:
                            #   Increment the CDS position counter. There has to
                            #   be a non-gap character at the query positions
                            cds_pos += 1
                            tmp = line.strip().split()
                            #   Get the gene ID
                            #   We will assume this is the first part before the
                            #   _Predictions suffix.
                            geneid = os.path.basename(pred_file).rsplit('_', 1)[0]
                            anno = [geneid, str(cds_pos)] + tmp
                            gene_preds.append(anno)
        return gene_preds

    def compile_predictions(self, pred_data):
        """Put all the prediction data into a nice report."""
        #   Define a header for the output
        header = [
            'GeneID',
            'CDSPos',
            'AlignedPosition',
            'L0',
            'L1',
            'Constraint',
            'Chisquared',
            'P-value',
            'SeqCount',
            'Alignment',
            'ReferenceAA',
            'MaskedConstraint',
            'MaskedP-value']
        #   Define an output filename, in the predictions directory.
        self.mainlog.info(
            'Trying to use filename ' +
            self.preddir + '/Combined_Report.txt for report.'
            )
        if os.path.isfile(os.path.join(self.preddir, 'Combined_Report.txt')):
            self.mainlog.warning(
                self.preddir + '/Combined_Report.txt' +
                ' exists! Will overwrite!')
        handle = open(os.path.join(self.preddir, 'Combined_Report.txt'), 'w')
        handle.write('\t'.join(header) + '\n')
        for gene_prediction in pred_data:
            for snp_prediction in gene_prediction:
                handle.write('\t'.join(snp_prediction) + '\n')
        handle.flush()
        handle.close()
        return
