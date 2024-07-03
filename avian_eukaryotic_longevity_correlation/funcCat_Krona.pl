#!/usr/bin/perl -w
use strict;
open(GC,"< denovo_combined_annotations.gff");
my $count=0;
open (SEQOUT, '> out.txt');
while (defined(my $line = <GC>)){ #get each line
    chomp $line;
    if(!($line =~ m/^#/)){
        my @xs = split(/\t/, $line);
        my $ttt = $xs[8];
        my @zs = split /;/, $ttt;
        my @qs = split(/=/, $zs[1]);
        # prints the mean of the coverage histogram between the 2 coords
        my $mean = `grep $xs[0] denovo_scaffold_coverage.txt | awk \'{if ((\$2 >= $xs[3]) && (\$2 <= $xs[4])) print \$3}\' | awk '\{total += \$1} END {print total/NR}\'`;
        chomp $mean;
        #print "$xs[0]\t$xs[3]\t$xs[4]\t$qs[1]\t$mean\n";
        print SEQOUT "$xs[0]\t$mean\n";
    }
    if ($count % 500 == 0){ 
        my $perc = ($count/64434)*100;
        print "$perc\n";
    }
    $count = $count + 1
}
close SEQOUT;
