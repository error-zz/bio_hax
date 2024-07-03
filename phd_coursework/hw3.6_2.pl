#!/usr/bin/perl -w
# J. Craig Venter Institute
# Coded by Jamison M. McCorrison (jmccorri@jcvi.org)

use strict;
no warnings;
my $vocab = shift; # vocab.txt
my $unigram = shift; # unigram.txt
my $bigram = shift; # bigram.txt

#################################################
# DATA PREP
print "EXE : paste $vocab $unigram > 3_6_unigram_input.txt\n";
system("paste $vocab $unigram > 3_6_unigram_input.txt");
my $ttl_count = `cat $unigram | awk \'{sum += \$1} END {print sum}\'`;
chomp $ttl_count;
print "\n($ttl_count) total word occurrences\n";

open(INPUTFILE,"< 3_6_unigram_input.txt");
my $c = 1;
my %u_freqs;
my %u_probs;
my %id_key;
my %rev_id_key;
while (defined(my $line = <INPUTFILE>)){
	my @xs = split(/\s+/, $line);
	$u_freqs{$xs[0]} = $xs[1];
	$u_probs{$xs[0]} = $xs[1]/$ttl_count;
	#DEBUG# print "( $xs[0] )\t$u_freqs{$xs[0]}\t$u_probs{$xs[0]}\n";
	$id_key{$c} = $xs[0];
	$rev_id_key{$xs[0]} = $c;
	$c++;
}
close INPUTFILE;

#################################################
# 6A
print "\n\n(3A) All words starting with letter B and their probabilities:\n";
foreach(reverse sort { $u_probs{$a} <=> $u_probs{$b} } keys(%u_probs)){
	if ($_ =~ m/^B/){
		print "\t$_\t$u_probs{$_}\n";
	}
}
print "\n";

#################################################
# 6B
print "(3B) Top 10 words to follow string ONE, sorted by probability:\n";

# count instances of ONE = 17 (see debug output below)
my $ttl_count_adjacent = `cat $bigram  | awk \'{if (\$1 == 17) print}\' | awk \'{sum += \$3} END {print sum}\'`;
chomp $ttl_count_adjacent;
print "\n($ttl_count_adjacent) total occurences of words following ONE\n";
open(INPUTFILE,"< $bigram");
my %b_freqs;
my %b_probs;
my %queried_probs;
while (defined(my $line = <INPUTFILE>)){
	my @xs = split(/\s+/, $line);
	$b_freqs{ $id_key{$xs[0]} }{ $id_key{$xs[1]} } = $xs[2];
	$b_probs{ $id_key{$xs[0]} }{ $id_key{$xs[1]} } = $xs[2]/$ttl_count_adjacent;
	if ($id_key{$xs[0]} eq "ONE"){
		# print "( $xs[0] = $id_key{$xs[0]} )( $xs[1] = $id_key{$xs[1]} )\t$b_freqs{$id_key{$xs[0]}}{$id_key{$xs[1]}}\t$b_probs{$id_key{$xs[0]}}{$id_key{$xs[1]}}\n";
		$queried_probs{$id_key{$xs[1]}} = $b_probs{ $id_key{$xs[0]} }{ $id_key{$xs[1]} };
	}
}
my $cutoff = 1;
foreach(reverse sort { $queried_probs{$a} <=> $queried_probs{$b} } keys(%queried_probs)){
		print "\t$_\t$queried_probs{$_}\n";
		$cutoff++;
		if ($cutoff == 11){ last; }
}
close INPUTFILE;

#################################################
# 6C - see subroutines below
my $sentence = "The stock market fell by one hundred points last week";
print "\n$sentence\n";
my $p_u = unigrammer($sentence, %u_probs);
print "\tp_u = $p_u\n";
my $p_b = bigrammer($sentence, \%u_probs, \%b_probs);
mixturemodel($sentence, \%u_probs, \%b_probs);

$sentence = "The fourteen officials sold fire insurance";
print "\n$sentence\n";
$p_u = unigrammer($sentence, %u_probs);
print "\tp_u = $p_u\n";
my $p_b = bigrammer($sentence, \%u_probs, \%b_probs);
print "\tp_b = $p_b\n";
mixturemodel($sentence, \%u_probs, \%b_probs);


#################################################
# 6C,D - unigrammer
sub unigrammer{
	print "!--UNIGRAMMER--!\n";
	my $sentence = shift;
	my (%u_probs) = @_;
	my @ordered_words = split(/\s+/, $sentence);
	my $p = 1;
	foreach my $w (@ordered_words){
		$w = uc($w);
		$p *= $u_probs{$w};
		print "\t$w\t$u_probs{$w}\tp= $p\n";
	}
	return natural_log($p);
}
#################################################
# 6C,D - bigrammer
sub bigrammer{
	print "!--BIGRAMMER--!\n";
	my $sentence = shift;
	my ($l1, $l2) = @_;
	my (%u_probs) = %$l1;
	my (%b_probs) = %$l2;
	my @ordered_words = split(/\s+/, $sentence);
	my $p = 1;
	my $w_c = 0;
	my $prev_w = "<s>";
	foreach my $w (@ordered_words){
		$w = uc($w);
		$p *= $b_probs{$prev_w}{$w};
		print "\t$prev_w\t->\t$w\t$u_probs{$w}\tp= $p\n";
		$prev_w = $w;
	}
	if ($p == 0){ return 0; }
	else{ return natural_log($p); }
}

#################################################
# 6E - bigrammer
sub mixturemodel{
	print "!--MIXTURE MODEL--!\n";
	my $sentence = shift;
	my ($l1, $l2) = @_;
	my (%u_probs) = %$l1;
	my (%b_probs) = %$l2;
	my @ordered_words = split(/\s+/, $sentence);
	for (my $lambda = 0.01; $lambda <= 1; $lambda += 0.01){
		my $p = 1;
		my $w_c = 0;
		my $prev_w = "<s>";
		foreach my $w (@ordered_words){
			$w = uc($w);
			$p *= (1-$lambda)*$u_probs{$prev_w} + $lambda*$b_probs{$prev_w}{$w};
			#print "\t$prev_w\t->\t$w\t$u_probs{$w}\tp= $p\n";
			$prev_w = $w;
		}
		$p = natural_log($p);
		print "$lambda\t$p\n";
	}
	#if ($p == 0){ return 0; }
	#else{ return natural_log($p); }
}

#################################################
# 6C,D,E - natural_log_calculator
sub natural_log {
  my $q = shift;
  return log($q)/log(2.7182818284590452353602874713527);
}
