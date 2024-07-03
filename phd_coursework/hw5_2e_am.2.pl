#!/usr/bin/perl -w

system("source /Volumes/Mac\ HD/Applications/PDL/setup_bash");
use strict;
use PDL; # using perl data language (numpy equivalent for perl)
         # if this fails to open, install PDL and run 'setup_bash' in the install directory to
         # configure your run environment for execution!
no warnings;


my %spectX;
my $c = 0;
open (READ_IN, "spectX.txt");
while (defined(my $line = <READ_IN>)){ #get each line
	chomp $line;
	my @xs = split(/\s+/, $line);
	$spectX{$c} = \@xs;
	$c++;
}
close READ_IN;

my %spectY;
$c = 0;
open (READ_IN, "spectY.txt");
while (defined(my $line = <READ_IN>)){ #get each line
	chomp $line;
	$spectY{$c} = $line;
	$c++;
}
close READ_IN;

# initialize
my $i = 2/23;
my $p = pdl [ $i,$i,$i,$i,$i, $i,$i,$i,$i,$i, $i,$i,$i,$i,$i, $i,$i,$i,$i,$i, $i,$i,$i ];
my %LL;
my %errors;

print "\n-----\nInitial p:\n";
print "p : $p\n-----\n";

# perform 256 iterations
for (my $c = 0; $c <= 256; $c++){
	
	# E-STEP #

	my $L = 0;
	my $err = 0;
	my $qcount = 0;
	foreach my $r ( sort {$a <=> $b} keys %spectX){
		# DEBUG # print "\n--[ R = $r ]----------------------\n@{$spectX{$r}}\n";
		my $p_c = 1;
		foreach(@{$spectX{$r}}){ 
			$p_c *= (1-$p)**$_; 
		}
		$p_c = 1 - $p_c;
		my $p_c_1 = $p_c->slice('0,:0,');
		# DEBUG # print "\n! p_c = $p_c_1";
		my $p_y = $spectY{$r} x $p_c_1 + (1-$spectY{$r}) x (1-$p_c_1);
		# DEBUG # print "\n! p_y = $p_y";
		if ($p_y < 0.5){
			$err++;
		}
		$L += log($p_y);
		# DEBUG # print "! L running product ( $r ) = $L !\n";
		$qcount++;
	}
	my $LogLikelihood = $L / $qcount;
	if ($c == 0 || $c == 1 || $c == 2  || $c == 4  || $c == 8  || $c == 16  || $c == 32  || $c == 64  || $c == 128  || $c == 256){
		print "\n----------- (Iteration $c) -----------\n";
		print "Log Likelihood = $LogLikelihood\n";
		print "Error Count    = $err\n";
		print "P              = $p\n";
	} 

	# M-STEP # 
	my @new_p = ();

	for (my $m =0; $m < 23; $m++){
		my $init_p = 0;
		my $t = 0;
		for (my $r = 0; $r <= 256; $r++){
			my $p_c = 1;
			foreach(@{$spectX{$r}}){ 
				$p_c *= (1-$p)**$_; 
			}
			$p_c = 1 - $p_c;
			if ( ($spectY{$r} == 1) && ( $spectX{$r}[$m] == 1) ){
				my $tmpval = $p->flat->index(0);
				$init_p += $tmpval/$p_c;
			}
			if ( $spectX{$r}[$m] == 1 ){ $t++; }
		}
		if ($t > 0){
			my $pp = $init_p/$t;
			push(@new_p, $pp->flat->index($m));
		}else{
			push(@new_p, 0);
		}
	}
	# print "Updating p : $p\n";
	
	$p = pdl ([ @new_p ]);
	# undef @new_p;
	
	# print "\t-> $p\n";

	# DEBUG EXIT
	if ($c == 2){
		exit;
	}
}
