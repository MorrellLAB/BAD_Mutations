#!/usr/bin/perl
##by Li Lei, 20171001, in St.Paul;
#this is to creat the long sub file, which look like a ANNOVA table to feed BAD mutation to calculate the masked and unmasked p-value;

#usage: 
use strict;
use warnings;
#use Data::Dumper;
my $file = $ARGV[0];

#my %gidhash;
my @temp = split (/\//, $file);
my @temp1 = split (/\./,$temp[-1]);
my $gid = $temp1[0];
#print "$gid\n";
open(SUBS,  "$file") or die "Could not open $file";

#print "SNP_ID\tChromosome\tPosition\tSilent\tTranscript_ID\tCodon_Position\tRef_Base\tAlt_Base\tAA1\tAA2\tAA_Pos\tCDS_Pos\n";
foreach my $row (<SUBS>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my $SNPid = $rtemp[1];
        my $aa_pos = $rtemp[0];
        print "$SNPid\tNA\tNA\tNo\t$gid\tNA\tNA\tNA\tNA\t$aa_pos\tNA\n";     
}
close (SUBS);
#print Dumper(\%gidhash);
