#!/usr/bin/perl -w

use strict;
use PDL; # using perl data language (numpy equivalent for perl)
         # if this fails to open, install PDL and run 'setup_bash' in the install directory to
         # configure your run environment for execution!
no warnings;

# import from file

my $rowsize;
# import as hash of hashes
my %test3 = ();
my $x=0;
open (SIZE_IN, "newTest3.txt");
while (defined(my $line = <SIZE_IN>)){ #get each line
	chomp $line;
	my @xs = split(/\s+/, $line);
	$rowsize = @xs;
	my $y = 0;
	foreach(@xs){ 
		$test3{$x}{$y} = $_; 
		$y++;
	}
	$x++;
}
my $test3_count = $x;

# import as hash of hashes
my %test5 = ();
$x=0;
open (SIZE_IN, "newTest5.txt");
while (defined(my $line = <SIZE_IN>)){ #get each line
	chomp $line;
	my @xs = split(/\s+/, $line);
	my $y = 0;
	foreach(@xs){ 
		$test5{$x}{$y} = $_; 
		$y++;
	}
	$x++;
}
my $all_count = $test3_count + $x;

# set random weights using PDL random (perl "random" is notoriously NOT random)
my $w = random($rowsize,1);
# set error rate = 3% (as suggested at monday discussion session)
my $e = 0.03; 
# set learning rate = step size = n > 0 
my $n = 1;
# cutoff suggested in class : $n = 0.02/T
my $n_cutoff = 0.02/$all_count;

print "n = $n -> $n_cutoff\n";
print "e = $e\n";
print "RAND W :$w\n";

# the goal of maximum likelihood estimation is to update vector w <- w - n(dL/dw)
# at each stage:
# iterative update rule Phi <- Phi - n(df/dPhi)
#                              Phi - nf'
# derive L=sigmaT(yt-delta(wn*xt))xt

# evaluating test set 3
my $is3 = 1; # known truth!

while ($n > $n_cutoff){

	my $step_size = 0;
	my $L = 0;

	foreach my $x (sort {$a <=> $b} keys(%test3)){
		# for each line... pull full sequence into array, then convert to pdl matrix object
		my @tmp = ();
		foreach my $y (sort {$a <=> $b} keys $test3{$x}){ push(@tmp, $test3{$x}{$y}); }
		my $xt = pdl [@tmp];
		# print "($x)\txt:$xt\n";
		# calculate sigmoid function
		my $xtt = transpose($xt);
		# print "xtt:$xtt\n";
		my $wxtt = $w x $xtt;
		# print "wxtt:$wxtt\n";
		my $sz = 1/(1+2.7182818284590452353602874713527**(-1 - $wxtt));
		# print "($x)\tsz:$sz\n";
		# step size
		$step_size += (($is3-$sz) x $xtt);
		# log-likelihood calculation
		# with pdl operator overrides, log is natural log by default (and log10 denotes base 10 log)
		$L += ($is3 x log($sz));
		print "($x)\tstep_size:$step_size\n";
		print "($x)\tL:$L\n";
	}
	foreach my $x (sort {$a <=> $b} keys(%test3)){
		my @tmp = ();
		foreach my $y (sort {$a <=> $b} keys $test3{$x}){ push(@tmp, $test3{$x}{$y}); }
		my $xt = pdl [@tmp];
		
	}

}

