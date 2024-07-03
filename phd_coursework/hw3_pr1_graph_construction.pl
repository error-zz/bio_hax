#!/usr/bin/perl -w

use strict;
no warnings;

########### INPUT ##########

# hw3_pr1_graph_construction.pl <n> <weight_range_low> <weight_range_high>
my $n = shift;
	print "INPUT:\n\tN = $n\n";
my $weight_range_low = shift;
	print "\tLOW RANGE OF WEIGHT SCALE = $weight_range_low\n";
my $weight_range_high = shift;
	print "\tHIGH RANGE OF WEIGHT SCALE = $weight_range_high\n";
my $prefix = shift;
	print "\tOUTPUT PREFIX = $prefix\n";

if ( !($n) || !($weight_range_low) || !($weight_range_high) ){
	print "BAD INPUT!  SEE EXAMPLE USAGE:\n";
	print "hw3_pr1_graph_construction.pl <n> <weight_range_low> <weight_range_high>\n\n";
	exit;
}else{
	print "ALL INPUTS DETECTED, BEGINNING PROCESS.\n\n";
}

########### MAIN! ##########

# SELECT N VALUE
# UPDATE TO LOG BASE 2 AND UPPER LIMIT ON SAMPLE SIZE

my $minimum_subset_n = 2 * log2($n);
print "MINIMUM SUBSET FOR EVALUATION : $minimum_subset_n\n";
#my $sample_size = abs(int($minimum_subset_n))+1;
my $sample_size = roundup($minimum_subset_n);
print "SELECTED SAMPLE SIZE : $sample_size\n";

# INITIALIZE WITH SEQUENTIAL P-VALUES
my %edges = ();
my @vertices = ();
my @p_vals;
my @out_strings = ();

print "\nGenerating sequential edges and associated weights...\n";
push(@p_vals, $weight_range_low);
push(@p_vals, $weight_range_high);
print "BINARY WEIGHT OPTIONS: @p_vals\n";
for (my $z = 0; $z < $sample_size; $z++){
	print "( $z )\t";
	my $zz = $z+1;
	if($z != ($sample_size-1)){
		# my @ran = randarray(@p_vals);
		my $randnum = rand(1);
		# print "R: $randnum\n";
		if ($randnum > 0.5){ 
			$edges{$z}{$zz} = $p_vals[0];
		}else{ 
			$edges{$z}{$zz} = $p_vals[1];
		}
		my $tmp_str = "$z\t$zz\t$edges{$z}{$z+1}";
		print "\t$tmp_str\n";
		push(@out_strings, $tmp_str);
	}else{
		# interior vertex
		print "\tComplete.\n";
	}
}

# build all possible random weight values
my @weights = ();
for (my $w = $weight_range_low; $w <= $weight_range_high; $w++){
	push(@weights, $w);
}

print "\nGenerating non-sequential edges and associated weights...\n";
print "WEIGHT RANGE: @weights\n";

for (my $k = 0; $k < $sample_size; $k++){
	for (my $b = 0; $b < $sample_size; $b++){
		if (($k != $b) && ($k != ($b-1))){
			print "( $k ) ( $b )  ";
			my @ran = randarray(@weights);
			# print "RAND WEIGHT RANGE: @ran\n";
			$edges{$k}{$b} = $ran[0]; 
			my $tmp_str = "$k\t$b\t$edges{$k}{$b}";
			print "\t$tmp_str\n";
			push(@out_strings, $tmp_str);
		}
	}
}

my $outfile = "$prefix" . "_graph.txt";
print "\nPrinting to file $outfile ...\n\n";
				
if (-s("./tmp_out.txt")){ system("rm tmp_out.txt"); }	
open OUTOUTOUT, ">> ./tmp_out.txt";
my $tmp = $sample_size-1;
print OUTOUTOUT "0\n$tmp\n";
foreach(sort {$a <=> $b} @out_strings){
	print OUTOUTOUT "$_\n";
}
close OUTOUTOUT;
system("mv tmp_out.txt $outfile");
print "Run completed successfully.  Exiting.\n\n";
exit;


##################################################################################
sub randarray {
	my @array = @_;
	my @rand = undef;
	my $seed = $#array + 1;
	my $randnum = int(rand($seed));
	$rand[$randnum] = shift(@array);
	while (1) {
		my $randnum = int(rand($seed));
		
		if (!(defined $rand[$randnum])){
			$rand[$randnum] = shift(@array);
		}
		last if ($#array == -1);
	}
	return @rand;
}

##################################################################################
sub log2 {
    my $n = shift;
    return log($n)/log(2);
}

##################################################################################
sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}