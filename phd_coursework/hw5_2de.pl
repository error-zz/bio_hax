#!/usr/bin/perl -w

# system("source /Volumes/Mac\ HD/Applications/PDL/setup_bash");
use strict;
use PDL; # using perl data language (numpy equivalent for perl)
         # if this fails to open, install PDL and run 'setup_bash' in the install directory to
         # configure your run environment for execution!
no warnings;


my %spectX;
my $c = 0;
open (READ_IN, "spectX.txt");
while (defined(my $line = <READ_IN>)){ #get each line
	chomp $line;
	my @xs = split(/\s+/, $line);
	$spectX{$c} = \@xs;
	$c++;
}
close READ_IN;

my %spectY;
$c = 0;
open (READ_IN, "spectY.txt");
while (defined(my $line = <READ_IN>)){ #get each line
	chomp $line;
	$spectY{$c} = $line;
	$c++;
}
close READ_IN;

# initialize
my $i = 2/23;
my $p = pdl [ $i,$i,$i,$i,$i, $i,$i,$i,$i,$i, $i,$i,$i,$i,$i, $i,$i,$i,$i,$i, $i,$i,$i ];
my @p_arr = ( $i,$i,$i,$i,$i, $i,$i,$i,$i,$i, $i,$i,$i,$i,$i, $i,$i,$i,$i,$i, $i,$i,$i );

print "\n-----\nInitial p:\n";
print "p : $p\n-----\n";

my $initial_probability;
my $L = 0;
my $errors = 0;
my @keys_spectX = keys(%spectX);
my $size_spectX = @keys_spectX;
for (my $r = 0; $r < $size_spectX; $r++){

	# calculate probability (equivalent to calculation from 2.a)
	my @tmparr = @{$spectX{$r}};
	my $x_r = pdl ( [ @tmparr ] );
	my $y_r = $spectY{$r};
	my $p_r;
	foreach(@p_arr){ $p_r *= (1-$p)^$x_r }
	$initial_probability = log ( ( $y_r x (1-$p_r) ) + ( (1-$y_r) x $p_r ) );
	print "\n-----\nProbability calculation:\n";
	print "x_r : $x_r\ny_r : $y_r\np_r : $p_r\nprob : $initial_probability\n-----\n";

	# calculate log likelihood (equivalent to calcaulation from 2.b)


}
$L = $L/$size_spectX;