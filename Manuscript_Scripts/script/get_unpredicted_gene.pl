#!/usr/bin/perl
##by Li Lei, 20160429, in St.Paul;
#this is to extract the file fasta file path which is not predicted with BAD_Mutation!Since I did prediction for two batches of genes.
#usage:./get_unaligned_gene.pl predicted_gene_list_20160429 both_manual_fasta_file.list >both_manual_fasta_file_2.list
use strict;
use warnings;
use Data::Dumper;
my ($file1, $file2) = @ARGV;

my %gidhash;
open(GENEID,  "$file1") or die "Could not open $file1";

foreach my $row (<GENEID>){
        chomp $row;
        $gidhash{$row}=0;
        #print "$rtemp1[0]\n";
}
close (GENEID);
#print Dumper(\%gidhash);

open(IN,  "$file2") or die "Could not open $file2";

foreach my $row (<IN>){
        chomp $row;
        my @rtemp = split(/\//,$row);
        my @rtemp1 = split(/\./,$rtemp[8]);
        #print "$rtemp1[0]\n";
        if (not exists $gidhash{$rtemp1[0]}){
            print "$row\n";
        }
        #print "$rtemp1[0]\n";
}
close (IN);

