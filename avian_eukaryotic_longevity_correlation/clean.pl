#!/usr/bin/perl -w
# J. Craig Venter Institute
# Written by Jamison M. McCorrison (jmccorri@jcvi.org)
no warnings;
use strict;

my $fastq = shift;
my $prefix = shift;

my $runcmd = "grep -n -m 1 \' 2:\' $fastq | tr \':\' \' \' | awk \'{print \$1}\'";
print "1\t$runcmd\n";
my $n = `$runcmd`; chomp $n;
$n = $n - 1;

$runcmd = "head -n " . $n . " " . $fastq . " > " . $prefix . ".R1.fastq";
print "2\t$runcmd\n";
system($runcmd);

$runcmd = "awk \'NR>" . $n ."\' " . $fastq . " > " . $prefix . ".R2.fastq";
print "3\t$runcmd\n";
system($runcmd);
