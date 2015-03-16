LRT Pipeline
============

Overview
--------
A package to perform the a likelihood ratio test (LRT) for the prediction of
deleterious variants, as described in 
[Chun and Fay (2009)](http://genome.cshlp.org/content/19/9/1553.long). This
package is specially designed to work on plant species.

Data Sources
------------
Currently, the package downloads all CDS FASTA files from
[Phytozome 10](http://phytozome.jgi.doe.gov/). In order to run the script, you
will need a JGI Genomes account. Fetching data from 
[Ensembl](http://plants.ensembl.org) is not yet implemented, but is planned.

Method
------
The package will predict deleterious variants in a user-supplied sequence by
using [tblastx](http://blast.ncbi.nlm.nih.gov/Blast.cgi) to identify likely
homologues in each species. The homologues will then be aligned using
[PRANK](http://wasabiapp.org/software/prank/) and a precomputed phylogenetic
tree of the species. The model from 
[Chun and Fay(2009)](http://genome.cshlp.org/content/19/9/1553.long) will then
be applied over each specified codon, and a prediction as to whether the query
variant is deleterious will be returned.

Input
-----
This package accepts input into the prediction pipeline in the form of a FASTA
formatted file containing a single sequence, and another file containing a list
of affected codons, one per line. The substitutions file may optionally
contain a SNP ID as a second column, separated by a tab. 

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
at positions 4, 10, 25, and 100 in the **protein** sequence, and have the IDs
``SNP_1``, ``SNP_2``, ``SNP_3``, and ``SNP_4``, respectively. Note that while
the FASTA file contains **nucleotide** sequence, the missense variant
information is given in terms of **protein** sequence, or codon number. This is
necessary to allow the software to calculate synonymous substitution rates.

TODO
----
* 'fetch' subcommand
    * Add downloading from Ensembl Plants
    * ~~Integreate creation of BLAST databases~~ **Done!**
    * ~~Add MD5 checking for CDS files to avoid unnecessary downloads~~ **Done!**
* 'predict' subcommand
    * Write input parsing script
    * Write BLAST searching script
    * Write BLAST output parsing script
    * Write PRANK sequence alignment script
    * Compute phylogentic tree of species in database
    * Implement the LRT itself
