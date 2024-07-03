#!/usr/bin/perl -w

use strict;
use PDL; # using perl data language (numpy equivalent for perl)
no warnings;

# Parse files to hash

my %nasdaq00 = ();
my $c=0;
open (SIZE_IN, "nasdaq00.txt");
while (defined(my $line = <SIZE_IN>)){ #get each line
	chomp $line;
	$nasdaq00{$c} = $line;
	$c++;
}
my %nasdaq01 = ();
$c=0;
open (SIZE_IN, "nasdaq01.txt");
while (defined(my $line = <SIZE_IN>)){ #get each line
	chomp $line;
	$nasdaq01{$c} = $line;
	$c++;
}

# Initialize

my $aa = pdl [ [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0] ];   # 4x4 array
print "############################################\n# INITIALIZE \n############################################\n";
print "a= $aa\n";
my $bb = pdl [0,0,0,0]; # 4x1 array
print "b=\n$bb\n############################################\n";


foreach my $r (sort {$a <=> $b} keys (%nasdaq00)){
	if ($r > 3){
		print "\n\n########### $r ###########\n";
		my $x = pdl [ $nasdaq00{$r-1}, $nasdaq00{$r-2}, $nasdaq00{$r-3}, $nasdaq00{$r-4} ];
		# print "x : $x\n";
		my $xt = transpose($x);
		# print "xt : $xt\n";
		my $xtx = $xt x $x;
		# print "xtx : $xtx\n";
		$aa += $xtx;
		my $this_x = $x x $nasdaq00{$r};
		$bb += $this_x;
		# print "aa: $aa\n";
		print "bb: $bb\n";
	}
}

my $ww = ( inv($aa) x transpose($bb));
print "\nfin.\n\n############################################\n\nw= $ww\n";

# Start 5.b 
# Compare the model’s performance (in terms of mean squared error) on the data from the years 2000 and 2001
my $mse = 0;
foreach my $r (sort {$a <=> $b} keys (%nasdaq00)){
	if ($r > 3){
		my $x = pdl [ $nasdaq00{$r-1}, $nasdaq00{$r-2}, $nasdaq00{$r-3}, $nasdaq00{$r-4} ];
		print "x : $x\n";
		# multiply by precalculated weights
		my $xw = $x * $ww;
		$mse += (($xw-$x)^2);
	}
}
print "Trained on ( 2000 ) - Classifying ( 2000 ) :\n\tMSE = $mse\n;"


