Objective: 
==========

Generate a flexible pipeline for predicting deleterious mutations, capable of handling 
queries from different plant species and easily updated with new plant genomes that become 
available.

<ol>
  <li>
  <p>
    Data files, sources and formats
  </p>
<p>
<b>Sequences:</b>
Multisequence FASTA file of all coding sequences (nucleotide) for each plant genome 
available. Corresponding multisequence FASTA file of translated protein sequences. This can be 
done using EMBOSS's transeq command. Output should be checked for stop codons (*) in translation 
before the end of the protein. These will be removed unless there are too many (100s) indicating 
problems with frame or A-I editing.
</p>

<p>
<b>Source:</b> We need to specify source of CDS (e.g. ensembl or phytozome, assembly/release number 
for both genome and annotation) and specific parameters used, e.g. whether to use primary transcript 
(e.g. longest) or all transcripts for a gene and any other settings used to obtain the data from BioMart.
downloads.
</p>

<p>
<b>Species:</b> We also need to specify a list of plant genomes that we intend to use. Although for Barley 
we focus on grasses and could limit database to angiosperms, I don't see harm in having other more 
distant genomes included.
</p>

<p>
Phyozome info for releases is here: http://phytozome.jgi.doe.gov/pz/portal.html#!releaseNotes
</p>

<p>
<b>Format:</b> Each protein/species dataset needs to be formated as a blast database. (Record blast 
version). As we discussed we can also try formating a single database with all species included. In this 
case, we have to ensure that each sequence name ( after > ) can be used to know which species it 
came from.
</p>

<p>
One issue we may face is speed. Computing clusters can give us 100-1000 nodes, which I suspect 
should be fine for 20k blast searches based on barley genes. But some attention to this would be 
useful. One possibility is to use BLAT (100x faster but less sensitive). I think this is possible since 
BLAT works fine for >80% ID and I don't think we want distant homologues, which is what blastp is 
designed to detect. We should also modify blast parameters to suite our project.
</p>

</li>

<li>
Alignment pipeline
<ol>
  <li>Given a query sequence, blastp to protein database(s).</li>
  <li>Filter results for only best hit for each species</li>
  <li>Extract nucleotide and protein sequences corresponding to best blast hits sequences and 
  positions (note that we may just want to extract entire sequence rather than sequences with 
  partial hits - which could be removed in filtering step).</li>
  <li>Multiple alignment based on proteins (PRANK) and then corresponding nucleotide alignment 
  based on protein alignment (EMBOSS transalign).
</ol>
</li>

<li>LRT and SIFT calling LRT requires a nucleotide alignment and will output calls for any site specified. SIFT requires a protein alignment and will output calls for any or all sites. Hopefully, we can implement 1 & 2 for all genes and then SIFT and LRT can run quickly with each user query. 1 & 2 can be updated only when new genomes/annotations are available or if someone 
wants to do predictions using another species (e.g. soy).
</li>
</ol>
