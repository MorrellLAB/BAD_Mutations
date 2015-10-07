BAD_Mutations
============

Overview
--------
BAD_Mutations (**B**LAST-**A**lign-**D**eleterious?) performs a likelihood
ratio test (LRT) for the prediction of deleterious variants, as described in 
[Chun and Fay (2009)](http://genome.cshlp.org/content/19/9/1553.abstract). The
package is designed to identify deleterious variants in plant species.
BAD_Mutations is first used in *Kono et al. (In Prep.)*, and will have a formal
publication to follow.

Data Sources
------------
The package downloads all CDS FASTA files from
[Phytozome 10](http://phytozome.jgi.doe.gov/) and 
[Ensembl Plants](http://plants.ensembl.org). In order to run the script, you
will need a JGI Genomes account. As of October 2015, BAD_Muations sources data
from 37 genome sequences. BAD_Muatations is configured to fetch all
Angiosperm CDS sequences, but it is possible to modify the fetching scripts
to retrieve other data sets.

Method
------
The package will predict deleterious variants in a user-supplied sequence by
using [tblastx](http://blast.ncbi.nlm.nih.gov/Blast.cgi) to identify
orthologues in other plant species. The nucleotide sequences from orthologues
are then aligned using [PASTA](http://www.cs.utexas.edu/~phylo/software/pasta/)
and a phylogenetic tree estimated from the alignment. The model from 
[Chun and Fay (2009)](http://genome.cshlp.org/content/19/9/1553.abstract) is
then applied for each specified codon, and a prediction as to whether the query
variant is deleterious will be returned. Since the model predicts conservation
of codons, it can only predict coding SNPs.

Input
-----
This package accepts input as a FASTA formatted file containing a single 
sequence, and a second file specifying a list of codons with nonsynonymous
substitutions. The substitution file specifies the codon number (starting 
with codon 1), one per line. Substition positions are the residue number 
in the coding sequence (CDS) translation of the gene. Users should submit 
sequences in 5' to 3' orientation, though there is no programmatic means of
confirming this. If a gene sequence is submitted 3' to 5', the substitution
position should also be submitted in that order. The substitutions file
may optionally contain a SNP ID as a second column, separated by a tab. While
individual SNPs are reported, note that analysis is a metric of codon
conservation.

Example: The pair of files
```
Sequence_1.fasta:
>Sequence_1
ATG...

Sequence_1.subs:
4   SNP_1
10  SNP_2
25  SNP_3
100 SNP_4
```

contain four missense variants to predict in one sequence. The variants occur
at positions 4, 10, 25, and 100 in the **protein** sequence (i.e., they are
residue number 4, 10, 25, and 100 in the CDS translation), and have the IDs
``SNP_1``, ``SNP_2``, ``SNP_3``, and ``SNP_4``. Note that while
the FASTA file contains **nucleotide** sequence, the missense variant
information is given in terms of **protein** sequence, or codon number. This is
necessary to allow the software to calculate synonymous substitution rates.

Runtimes and Benchmarks
-----------------------
By far, the slowest part of BAD_Mutations is fetching CDS sequences and
converting them to BLAST databases. This may take up to several hours,
depending on your network and disk speeds. The databases and FASTA files take
up approximately 4GB, as of October 2015. As more genomes are sequnced and
annotated, this figure will increase.

For a typical barley gene (~3000bp), BAD_Mutations can generate a phylogenetic
tree and multiple sequence alignment in approximately 5-10 minutes on a dekstop
computer (Intel i7 2.8GHz). Note, however, that this figure can vary depending
on the gene you are using. Rapidly evolving genes will be much more difficult
to align and estimate phylogenies from, and will take longer. Also note that
not every gene will have an ortholog with enough sequence similarity, so not
every gene will have every species represented in the alignment and tree. This
is not a problem for BAD_Mutations.

Predictions are much slower, and are currently being benchmarked.

TODO
----
* 'setup' subcommand
    * ~~Write configuration file~~ **Done!**
    * ~~Read and parse configuration file~~ **Done!**
* 'fetch' subcommand **Done!**
    * ~~Add fetching CDS sequences from Ensembl Plants~~ **Done!**
    * ~~Integreate creation of BLAST databases~~ **Done!**
    * ~~Add MD5 checking for CDS files to avoid unnecessary downloads~~ **Done!**
    * ~~Add config file for excluding sequences from Phytozome~~ **Done!**
* 'predict' subcommand
    * Write input parsing script **Deciding on Final Input Format**
    * ~~Write BLAST searching script~~ **Done!**
    * ~~Write BLAST output parsing script~~ **Done!**
    * ~~Write PRANK sequence alignment script~~ **Done!**
    * ~~Compute phylogentic tree of species in database~~ **Not necessary, build tree from MSA**
    * ~~Add option to "prune" species from the all-species phylogenetic tree~~ **Not necessary, build tree from MSA**
    * Implement the LRT itself
    * ~~Add support for parallel execution~~ **Parallelize at the process level, not within the package**
    * Add support for substitutions in nucleotide sequence as well as protein sequence
