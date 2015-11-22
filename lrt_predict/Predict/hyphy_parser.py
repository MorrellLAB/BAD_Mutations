#!/usr/bin/env python
"""Parse the HyPhy output and produce a report for each of the SNPs thave have
been predicted.

Originally written by Justin C. Fay, rewritten for general use in BAD_Mutations
by Thomas Kono."""

#   Import standard library modules here
import os


class HyPhyPrediction(object):
    """Holds data for a single SNP prediction."""
    pass


class HyPhyParser(object):
    """Class to compile and evaluate the predictions for SNPs. Will read a list
    of individual HyPhy output files in a directory and pull out the relevant
    lines. Will also assign a prediction of "NEUTRAL" or "DELETERIOUS"
    depending on the user-specified filters. In the case where a prediction
    cannot be made, the reason that the prediction failed is listed. A
    prediction can fail based on the number of species represented, or the
    alternate amino acid appearing in the alignment too many times."""

    def __init__(self, outdir, mseq, pval, ncodons=None):
        self.preddir = outdir
        self.minseq = mseq
        if ncodons:
            self.ntests = ncodons
        else:
            self.ntests = 0
        self.sig = pval
        return

    def get_prediction_files(self):
        """List the contents of the supplied directory and find all the files
        that contain predictions. These will end in '_Predictions.txt' - this
        is somewhat fragile, but the paths used throughout the package are
        consistent."""
        pass

    def parse_prediction(self, pred_file):
        """Parse a single prediction file. Get the constraint, the p-value, the
        number of sequences in the alignment, the runtime, and the number of
        reference and alternate amino acids in the alignment."""
        pass

    def compile_predictions(self):
        """Calculate appropriate p-value for the test based on the specified
        threshold and the number of predictions were run, and print out the
        formatted table."""
        pass
