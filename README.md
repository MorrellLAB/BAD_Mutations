LRT Pipeline
============

Overview
--------
A package to perform a likelihood ratio test (LRT) for the prediction of
deleterious variants, as described in 
[Chun and Fay (2009)](http://genome.cshlp.org/content/19/9/1553.abstract). The
package is designed to identify deleterious mutations in plant species.

Data Sources
------------
The package downloads all CDS FASTA files from
[Phytozome 10](http://phytozome.jgi.doe.gov/). In order to run the script, you
will need a JGI Genomes account. Fetching data from [Ensembl](http://plants.ensembl.org)
is not yet implemented, but is planned.

Method
------
The package will predict deleterious variants in a user-supplied sequence by
using [tblastx](http://blast.ncbi.nlm.nih.gov/Blast.cgi) to identify likely
homologues in other plant species. The nucleotide sequences from homologues
will then be aligned using [PRANK](http://wasabiapp.org/software/prank/)
and a precomputed phylogenetic tree of the species. The model from 
[Chun and Fay (2009)](http://genome.cshlp.org/content/19/9/1553.abstract) will then
be applied for each specified codon, and a prediction as to whether the query
variant is deleterious will be returned.

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

TODO
----
* 'fetch' subcommand
    * ~~Add fetching CDS sequences from Ensembl Plants~~ **Done!**
    * ~~Integreate creation of BLAST databases~~ **Done!**
    * ~~Add MD5 checking for CDS files to avoid unnecessary downloads~~ **Done!**
    * ~~Add config file for excluding sequences from Phytozome~~ **Done!**
* 'predict' subcommand
    * Write input parsing script **Deciding on Final Input Format**
    * ~~Write BLAST searching script~~ **Done!**
    * ~~Write BLAST output parsing script~~ **Done!**
    * Write PRANK sequence alignment script
    * ~~Compute phylogentic tree of species in database~~ **Not necessary, build tree from MSA**
    * ~~Add option to "prune" species from the all-species phylogenetic tree~~ **Not necessary, build tree from MSA**
    * Implement the LRT itself
