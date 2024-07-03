#!/usr/bin/perl -w

# perl -MCPAN -e shell
# CPAN> install Math::Pari
use Math::Pari qw/binomial pari2num gdiv/;

my %input_hash;
$input_hash{"cases"}{"AA"} = 18;
$input_hash{"cases"}{"AG"} = 9;
$input_hash{"cases"}{"GG"} = 3;
$input_hash{"ctrls"}{"AA"} = 8;
$input_hash{"ctrls"}{"AG"} = 9;
$input_hash{"ctrls"}{"GG"} = 3;

#fishers exact

# P = 
# (a+b choose a));
my $tmp1 = $input_hash{"cases"}{"AA"} + $input_hash{"ctrls"}{"AA"};
my $r1 = binomial( $tmp1, $input_hash{"cases"}{"AA"} );
print "! $tmp1 choose $input_hash{\"cases\"}{\"AA\"} !\t";
print "( $r1 )\n";
# (c+d choose c)
$tmp1 = $input_hash{"cases"}{"AG"} + $input_hash{"ctrls"}{"AG"};
my $r2 = binomial( $tmp1, $input_hash{"cases"}{"AG"} );
print "! $tmp1 choose $input_hash{\"cases\"}{\"AG\"} !\t";
print "( $r2 )\n";
# (e+f choose e)
$tmp1 = $input_hash{"cases"}{"GG"} + $input_hash{"ctrls"}{"GG"};
my $r3 = binomial( $tmp1, $input_hash{"cases"}{"GG"} );
print "! $tmp1 choose $input_hash{\"cases\"}{\"GG\"} !\t";
print "( $r3 )\n-----------------------\n";
# (n choose a+c+e)
my $n = $input_hash{"cases"}{"AA"} + $input_hash{"cases"}{"AG"} + $input_hash{"cases"}{"GG"} + $input_hash{"ctrls"}{"AA"} + $input_hash{"ctrls"}{"AG"} + $input_hash{"ctrls"}{"GG"};
$tmp1 = $input_hash{"cases"}{"AA"} + $input_hash{"cases"}{"AG"} + $input_hash{"cases"}{"GG"};
my $d = binomial( $n, $tmp1);
print "! $n choose $tmp1 !\t";
print "( $d )\n\n=\n\n";

my $soln = ($r1 * $r2 * $r3);
print "$soln\n\n";
#my $soln2 = `expr $soln / $d`;
my $soln2 = `perl perldivide.pl $soln $d`;
chomp $soln2;
print "$soln2\n\n";




