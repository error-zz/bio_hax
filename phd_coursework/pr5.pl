#!/usr/bin/perl -w

use strict;

my @alpha = ("A", "C", "G", "T");
foreach (my $x = 7; $x <= 12; $x++){
    print "\>Lequals4exponent" ;
    print "$x\n";
    foreach(my $t = 0; $t < 4**$x; $t++){
        my $rando = int(rand(4));
        print "$alpha[$rando]";
    }
    print "\n";
}

# validated
# head -n 2 test.txt | tail -n 1 | wc -c
#   16385
