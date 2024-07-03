#!/usr/bin/perl -w

###### VARIABLES AND WHATNOT

use strict;

my $inputfile = shift;
my $k = shift;
my $c = 0;
my $query;
open (FILE_IN, "$inputfile");
while (defined(my $line = <FILE_IN>)){ #get each line
	chomp $line;
	$query = $line;
}

###### BUILD SUFFIXES AND STORE INDICES

my %suffix_array;
for (my $i = 0; $i < length($query)-$k+1; $i++){
    my $substring = substr( $query, $i, $k );
    push( @{$suffix_array{$substring}}, $i );
    print "$i : $substring\n";
    
    
    
}

foreach my $substring ( sort(keys (%suffix_array)) ){
    print "$substring | ";
    foreach my $g ( @{$suffix_array{$substring}} ){
        print "$g ";
    }
    print "\n";
}

###### curated by hand following suffix calculation




