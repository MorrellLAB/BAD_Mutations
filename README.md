BAD_Mutations
============

Overview
--------
BAD_Mutations (**B**LAST-**A**ligned-**D**eleterious?) performs a likelihood
ratio test (LRT) for the prediction of deleterious variants, as described in 
[Chun and Fay (2009)](http://genome.cshlp.org/content/19/9/1553.abstract). The
package is designed to identify deleterious variants in plant species.
BAD_Mutations is first used in *Kono et al. (In Prep.)*, with publication detailing
features of the software to follow.

For more information, see the user manual.

TODO
----
* 'setup' subcommand
    * Implement dependency downloding
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
    * ~~Implement the LRT itself`` **Done!**
    * ~~Add support for parallel execution~~ **Parallelize at the process level, not within the package**
    * Add support for substitutions in nucleotide sequence as well as protein sequence
