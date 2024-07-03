#!/usr/bin/perl -w
use strict;

my $in = shift;
print $in;

my $count = 0;

open(GC,"< $in");
my $out = $in . ".flipList";
if (-s("$out")){ system("rm -rf $out"); }
system("touch $out");
open(PRINTME,"> $out");

while (defined(my $line = <GC>)){ #get each line
    
    chomp $line;
    $count = $count + 1;
    print "Cluster Number $count : ";
    
    my @xs = split(/\s+/, $line);
    my $xcnt = @xs;
    print "$xcnt entries : ";

    foreach my $x (@xs){
        print PRINTME "$x\t$count\n";
    }
    
}
