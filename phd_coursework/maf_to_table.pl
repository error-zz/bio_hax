#!/usr/bin/perl -w
use strict;
no warnings;


# SUBROUTINES ###########

sub runsys {
				my $cmd = shift;
				print "++ RUN CMD : $cmd\n";
				print LOG "++ RUN CMD : $cmd\n";
				system("$cmd");
}


# MAIN ###################

my $infile =shift(@ARGV);
if (!(-s("$infile"))){ 
	print "ERROR: $infile not found\nExiting.\n";
	exit;
}

# 2 column tsv: file location - previx
my %status;
my @all_gene_ids;
open INFILE, "< $infile";
while (defined( my $line = <INFILE> )){
	chomp $line;
	my @xs = split(/\s+/, $line);
	my $patient_id = $xs[0];
	my $gene_id = $xs[3];
	my $ref_nt = $xs[8];
	my $query_a = $xs[9];
	my $query_b = $xs[10];
	
	if (($ref_nt ne $query_a) && ($ref_nt ne $query_b)){
			$status{$patient_id}{$gene_id} = 2;
	}elsif (($ref_nt ne $query_b) || ($ref_nt ne $query_a)){
		if($status{$patient_id}{$gene_id} != 2){
			$status{$patient_id}{$gene_id} = 1;
		}
	}

	push(@all_gene_ids, $gene_id);
}
my @unique_genes = do { my %seen; grep { !$seen{$_}++ } @all_gene_ids };

# PRINT OUTPUT
# HEADER - GENE IDS
#my @patients = keys(%status);
print " \t";
#foreach(sort {$a <=> $b} keys(%{$status{$patients[0]}})){ print "$_\t"; }
foreach(@unique_genes){ print "$_\t"; }
print "\n";
foreach my $patient (sort {$a <=> $b} keys(%status)){
	print "$patient\t";
	#foreach my $gene (sort {$a <=> $b} keys(%{$status{$patient}})){
	foreach my $gene (@unique_genes){
		if (!($status{$patient}{$gene})){$status{$patient}{$gene} = 0;}
		print "$status{$patient}{$gene}\t";
	}
	print "\n";
}