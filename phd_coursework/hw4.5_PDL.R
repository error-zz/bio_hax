#!/usr/bin/perl -w

use strict;
use PDL;

open (SIZE_IN, "nasdaq00.txt");
my @nasdaq00 = ();
while (defined(my $line = <SIZE_IN>)){ #get each line
	chomp $line;
	push(@nasdaq(00), $line);
}

# How accurately can the index on one day be predicted by a linear combination of the three preceding indices?

# Initialize
my $A = pdl [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];   # 4x4 array
my $b = 

foreach my $index (@nasdaq(00)){
	
}
