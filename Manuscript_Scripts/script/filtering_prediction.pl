#!/usr/bin/perl
##by Li Lei, 20160502, in St.Paul;
#this is to do filtering to see if the variants are deleterious variants. Firstly based on the criteris:
#SeqCount>10;maskedConstraint<1;maskedP-value< 0.05/number of tested codon.
#secondly we also need to exclude variants where the derived state in barley was present in other species
#this program only takes care of the first filtering
#usage: ./filtering_prediction.pl /home/morrellp/llei/Deleterious_mutation_project/LTR_BAD_mutation/A_thaliana_BAD_Mutation/file_list/both_manual_subs_file.list  /home/morrellp/llei/Deleterious_mutation_project/LTR_BAD_mutation/A_thaliana_BAD_Mutation/file_list/prediction.extract.list
use strict;
use warnings;
use Data::Dumper;
my ($file1, $file2) = @ARGV;

my %gidhash;
open(SUBS,  "$file1") or die "Could not open $file1";

foreach my $row (<SUBS>){
        chomp $row;
        my @rtemp = split(/\//,$row);
        my @rtemp1 = split(/\./,$rtemp[8]);
        my $g_id = $rtemp1[0];
        my $count = &count_codons($row);
        $gidhash{$g_id} = $count;        
}
close (SUBS);
#print Dumper(\%gidhash);

open(IN,  "$file2") or die "Could not open $file2";

foreach my $row (<IN>){
        chomp $row;
        #print "$row\n";
        my @rtemp = split(/\//,$row);
        my @rtemp1 = split(/\./,$rtemp[9]);
        #print "$rtemp1[0]\n";
        if (exists $gidhash{$rtemp1[0]}){
            my $codons = $gidhash{$rtemp1[0]};
            my $cutoff = 0.05 / $codons;
            #print "$rtemp1[0]\t$cutoff\n";
            my $outFileName=$rtemp[9].".filtering";
            #print "$outFileName\t$rtemp1[0]\t$cutoff\n";
            &dele_filtering($row,$outFileName,$cutoff);
        }
        #print "$rtemp1[0]\n";
}
close (IN);

sub count_codons{
    my $file =shift;
    open(FILE,  $file) or die "Could not open $file";
    my $count=0;
    foreach my $row (<FILE>){
        chomp $row;
        $count++;
    }
    return $count;
    close (FILE);
}

sub dele_filtering{
    #my $file =shift;
    my ($file,$outfile,$threshold) = (@_);
    my $path = "/home/morrellp/llei/Deleterious_mutation_project/LTR_BAD_mutation/A_thaliana_BAD_Mutation/out_Ath/filter_extract_prediction";
    
    open(FILE, $file) or die "Could not open $file";
    open(OUT,  ">$path/$outfile");
    my $header=<FILE>;
       chomp $header;
    print OUT "$header\tStatus\n";
    foreach my $row (<FILE>){
        chomp $row;
        #print  OUT "$row\n";
        my @rtemp = split(/\t/,$row);
        my $adjust_P = sprintf("%.25f", $rtemp[10]);
        my $adjust_constraint = sprintf("%.25f", $rtemp[9]);
        if ($rtemp[6] > 10 && $adjust_constraint < 1 && $adjust_P < $threshold){
             print OUT "$row\tDeleterious\n";
        }
        else{
             print OUT "$row\tTolerant\n";
        }
    }
    close (FILE);
    close (OUT);
    return $outfile;
    #return print  "$row\n";

}


