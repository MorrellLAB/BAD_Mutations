#!/usr/bin/perl
#By Li Lei. For spliting the sequence and only use the gene modle id as sequence name in the individual fasta file;
my $file = shift;

open (INFILE, "< $file")or die "Can't open $file";
while (<INFILE>) {
		$line = $_;
		chomp $line;
		if ($line =~ /\>/) { #if has fasta >
			close OUTFILE;
			my $seq_name = substr($line,1);
			my @array = split(/\|/, $seq_name);
			my $new_file = $array[0];
			   $new_file =~ s/^\s+|\s+$//g; #remove all the whitespace; 
			   $new_file =~ s/\./\_/g; 
			   my $new_seq_name = $new_file;        
			      $new_file .= ".fasta";
			#print "$new_seq_name\t$new_file\n";
			open (OUTFILE, ">$new_file")or die "Can't open: $new_file $!";
			print OUTFILE ">$new_seq_name\n";
		}
		else{
		    print OUTFILE "$line\n";
	    }
}
close OUTFILE;
