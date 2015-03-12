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
This package accepts input into the prediction pipeline in the form of a FASTA-
formatted file, with affected codon information in the defline. The first piece
of qualifier information is the list of codons with nonsynonymous variants in 
them, and the second (optional) is an ID for the variants. The two qualifiers
should be space-delimited. Multi-record FASTA is supported.

Example:
```
>Sequence_1 codons=4,10,25,100 IDs=SNP_1,SNP_2,SNP_3,SNP_4
ATG...
>Sequence_2 codons=20,122 IDs=SNP_5,SNP_6
ATG...
...
...
...
```

Note that even though the defline contains codon information, the sequence in
the FASTA file should be protein-coding **nucleotide** sequence.

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
