#!/usr/bin/perl -w

use strict;

my @alpha = ("A", "C", "G", "T");
foreach (my $x = 7; $x <= 12; $x++){
    foreach(my $r = 0; $r < 10; $r++){
        print "\>Lequals4exponent" ;
        print "$x" . "_";
        print "$r\n";
        foreach(my $t = 0; $t < 4**$x; $t++){
            my $rando = int(rand(4));
            print "$alpha[$rando]";
        }
        print "\n";
    }
}
