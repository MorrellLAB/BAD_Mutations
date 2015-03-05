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
will need a JGI Genomes account. Fetching data from Ensembl is not yet
implemented, but is planned.

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
To be written.
