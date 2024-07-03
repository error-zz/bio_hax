#!/usr/local/bin/perl -w

use strict;

my $gene_list = shift;
my $maf = shift;

# sig_genes_x179.txt
# maf_to_Rinput.top179genes.tsv


# open (FPLUS, "> ./500plus.lst");
my @gene_list;
open (GENE_IN, "$gene_list");
while (defined(my $line = <GENE_IN>)){ #get each line
	chomp $line;
	push(@gene_list, $_);
}
close GENE_IN;

# open (FPLUS, "> ./500plus.lst");
open (MAF_IN, "$maf");
while (defined(my $line = <MAF_IN>)){ #get each line
	chomp $line;
	my @xs = split(/\s+/, $line);
	print "$xs[3]\n";
}
close MAF_IN;

