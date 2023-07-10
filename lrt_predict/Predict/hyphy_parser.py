#!/usr/bin/env python
"""Parse the HyPhy output and produce a report for each of the SNPs thave have
been predicted.

Originally written by Justin C. Fay, rewritten for general use in BAD_Mutations
by Thomas Kono."""

#   Import standard library modules here
import os
import math
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

    def logistic_p_values(self, unmasked, masked, cons, m_cons, ref, alt, aln):
        """Calculate a p-value for the test, based on a logistic regression
        reported in the accompanying manuscript. The equations are as follows:

log(p/(1-p)) = -2.407-0.2139*LRT(unmasked)-0.2056*constraint+0.07368*Rn-0.1236*An
log(p/(1-p)) = -2.453-0.1904*LRT(masked)-0.1459*constraint+0.2199*max(Rn,An)-0.2951*abs(Rn-An)

        where LRT is the log10(LRT p-value), constraint is conservation of the alignment
        column, Rn is the number of "reference" amino acids, and "An" is the
        number of alternate amino acids. We will report both the masked and
        the unmasked versions. The p-values will not be Bonferroni-adjusted.
        """
        #   Cast the variables to float
        try:
            #   Sometimes the values we pass are so small, they are 0. We set
            #   these to very small.
            if unmasked == '0':
                unmasked = 1e-16
            if masked == '0':
                masked = 1e-16
            unmasked_p = math.log10(float(unmasked))
            masked_p = math.log10(float(masked))
            con = float(cons)
            m_con = float(m_cons)
            if con > 10:
                con = 10.0
            if m_con > 10:
                m_con = 10.0
        except ValueError:
            self.mainlog.error('Non-numeric data passed to compile.')
            exit(1)
        #   Calculate the Rn and An values. These are easy, we just have to
        #   count the number of times it shows up
        rn = float(aln.count(ref))
        an = float(aln.count(alt))
        #   Then, calculate the p-values. We will calculate the right hand side
        #   of the equation first.
        m_rn_an = float(max(rn, an))
        a_rn_an = float(abs(rn - an))
        unmasked_rhs = -2.407 - (0.2139*unmasked_p) - (0.2056*con) + (0.07368*rn) - (0.1236*an)
        masked_rhs = -2.453 - (0.1904*masked_p) - (0.1459*m_con) + (0.2199*m_rn_an) - (0.2951*a_rn_an)
        #   Solve for p
        #   Protect the denominator calculation. Sometimes it overlfows. If it
        #   does, set it to Inf
        try:
            u_K = math.exp(-unmasked_rhs)
        except OverflowError:
            u_K = float('inf')
        try:
            m_K = math.exp(-masked_rhs)
        except OverflowError:
            m_K = float('inf')
        unmasked_log_p = 1 / (u_K + 1)
        masked_log_p = 1 / (m_K + 1)
        #   Return it
        return (unmasked_log_p, masked_log_p)

    def parse_prediction(self, pred_file):
        """Parse a single prediction file. Return the lines that have prediction
        information for them, along with a gene identifier."""
        #   set a flag that says whether or not we are reading in alignment data
        #   This is somewhat messy, but since the number of species in each
        #   alignment is variable, it is the most straightforward way to do it.
        in_aln = False
        #   And an empty list to store the prediction info
        gene_preds = []
        geneseq = ''
        self.mainlog.debug('Reading ' + pred_file)
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
                    self.mainlog.debug(str(cds_pos) + ' ' + line.strip())
                    #   The prediction lines are ones that do NOT end in 'NOSNP'
                    if 'NOSNP' in line:
                        geneseq += line.strip().split()[2]
                        #   Check the third field in the 'NOSNP' line - if it is
                        #   not '-', then we increment the CDS position
                        #   counter.
                        if line.strip().split()[2] == 'NA':
                            cds_pos += 1
                        continue
                    #   We first check if we've run off the edge of the
                    #   alignment data. The line right after the alignment
                    #   data starts with 'Alignment'
                    if line.startswith('Alignment'):
                        return gene_preds
                    geneseq += line.strip().split()[8]
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
        print(gene_preds)
        return gene_preds

    def add_regression(self, alt, prediction):
        """Use the information in a long format substitutions file and the
        SNP prediction information to calculate a logistic P-value. The data we
        need from the substitutions file is the alternate amino acid state."""
        #   Call the function that calculates the p-value. We send the data in
        #   the following order:
        #       Unmasked p-value
        #       Masked p-value
        #       Constraint
        #       Reference AA
        #       Alt AA
        #   First check the alt. If it's NA, we return NA as well
        if alt == 'NA':
            return prediction + ['NA', 'NA']
        u, m = self.logistic_p_values(
            prediction[7],
            prediction[12],
            prediction[5],
            prediction[11],
            prediction[10],
            alt,
            prediction[9])
        #   Just tack them on to the prediction data and return it
        u = str(u)
        m = str(m)
        return prediction + [u, m]

    def compile_predictions(self, pred_data):
        """Put all the prediction data into a nice report."""
        #   Define a header for the output
        header = [
            'VariantID',
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
            'MaskedP-value',
            'LogisticP_Unmasked',
            'LogisticP_Masked']
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
        for snp_prediction in pred_data:
            handle.write('\t'.join(snp_prediction) + '\n')
        handle.flush()
        handle.close()
        return
