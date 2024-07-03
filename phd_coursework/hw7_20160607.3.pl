#!/usr/bin/perl -w
use strict;
use PDL; # using perl data language (numpy equivalent for perl)
         # if this fails to open, install PDL and run 'setup_bash' in the install directory to
         # configure your run environment for execution!
         # source /Volumes/Mac\ HD/Applications/PDL/setup_bash
use PDL::NiceSlice;
no warnings;

if (!(-s("./prob_a1.txt"))){ print "\tNOTE! THIS SCRIPT ASSUMES ALL INPUT FILES ARE PROVIDED IN THE CURRENT WORKING DIRECTORY\nExiting.\n"; }

print "\n\n--[ File Import. ]-------------------------------\n\n";

# A1
my %prob_a1 = (); my %prob_a1_hash = (); my $c = 0;
open (FILE_IN, "./prob_a1.txt");
while (defined(my $line = <FILE_IN>)){ 
    chomp $line; my $q = 0;
    my @xs = split(/\s+/, $line);
    $prob_a1{$c} = [ @xs ];
    foreach(@xs){ $prob_a1_hash{$xs[0]}{$xs[1]} = $xs[2]; $q++; }
    $c++;
}
my $front = defined;
my $prob_a1_pdl = undef;
foreach my $id (keys %prob_a1){
    my $row = $prob_a1{$id};
    if (defined $front){
        $front = undef;
        $prob_a1_pdl = pdl ( [@$row] );
    }
    my $p_new = pdl ( [@$row] );
    $prob_a1_pdl = $prob_a1_pdl->glue(1,$p_new);
}

# A2
my %prob_a2 = (); my %prob_a2_hash = (); my $c = 0;
open (FILE_IN, "./prob_a2.txt");
while (defined(my $line = <FILE_IN>)){ 
    chomp $line; my $q = 0;
    my @xs = split(/\s+/, $line);
    $prob_a2{$c} = [ @xs ];
    foreach(@xs){ $prob_a2_hash{$xs[0]}{$xs[1]} = $xs[2]; $q++; }
    $c++;
}
my $front = defined;
my $prob_a2_pdl = undef;
foreach my $id (keys %prob_a2){
    my $row = $prob_a2{$id};
    if (defined $front){
        $front = undef;
        $prob_a2_pdl = pdl ( [@$row] );
    }
    my $p_new = pdl ( [@$row] );
    $prob_a2_pdl = $prob_a2_pdl->glue(1,$p_new);
}

# A3
my %prob_a3 = (); my %prob_a3_hash = (); my $c = 0;
open (FILE_IN, "./prob_a3.txt");
while (defined(my $line = <FILE_IN>)){ 
    chomp $line; my $q = 0;
    my @xs = split(/\s+/, $line);
    $prob_a3{$c} = [ @xs ];
    foreach(@xs){ $prob_a3_hash{$xs[0]}{$xs[1]} = $xs[2]; $q++; }
    $c++;
}
my $front = defined;
my $prob_a3_pdl = undef;
foreach my $id (keys %prob_a3){
    my $row = $prob_a3{$id};
    if (defined $front){
        $front = undef;
        $prob_a3_pdl = pdl ( [@$row] );
    }
    my $p_new = pdl ( [@$row] );
    $prob_a3_pdl = $prob_a3_pdl->glue(1,$p_new);
}

# A4
my %prob_a4 = (); my %prob_a4_hash = (); my $c = 0;
open (FILE_IN, "./prob_a4.txt");
while (defined(my $line = <FILE_IN>)){ 
    chomp $line; my $q = 0;
    my @xs = split(/\s+/, $line);
    $prob_a4{$c} = [ @xs ];
    foreach(@xs){ $prob_a4_hash{$xs[0]}{$xs[1]} = $xs[2]; $q++; }
    $c++;
}
my $front = defined;
my $prob_a4_pdl = undef;
foreach my $id (keys %prob_a4){
    my $row = $prob_a4{$id};
    if (defined $front){
        $front = undef;
        $prob_a4_pdl = pdl ( [@$row] );
    }
    my $p_new = pdl ( [@$row] );
    $prob_a4_pdl = $prob_a4_pdl->glue(1,$p_new);
}

# reward
my %reward = (); my %reward_hash = (); my $c = 0;
open (FILE_IN, "./rewards.txt");
while (defined(my $line = <FILE_IN>)){ 
    chomp $line; my $q = 0;
    my @xs = split(/\s+/, $line);
    $reward{$c} = [ @xs ];
    foreach(@xs){ $reward_hash{$c}{$q} = $_; $q++; }
    $c++;
}
print "-\n";
my $front = defined;
my $reward_pdl = undef;
foreach my $id (keys %reward){
    my $row = $reward{$id};
    if (defined $front){
        $front = undef;
        $reward_pdl = pdl ( [@$row] );
    }
    my $p_new = pdl ( [@$row] );
    $reward_pdl = $reward_pdl->glue(1,$p_new);
}

print "\n\n--[ Beginning iterations. ]-------------------------------\n\n";

# ------------[ START MAIN. ]----------------------------------------- # 

my $v_orig = zeroes(81);
my $bool = 0;
my $counter = 0;
while($bool == 0){
	my $v_new = compute_v($v_orig);
	if(approx($v_orig,$v_new)){$bool = 1;}
	$counter++;
}
print "Finished in $counter iterations.\n\n(Successful exit.)\n";


exit;

# ------------[ END MAIN. ]------------------------------------------- # 

sub compute_v{
	my $v = shift;
	my @as = ();
	foreach(my $s=0; $s<81; $s++){ 	# iterate through all states
		my $a = compute_a($s, $v); 	# compute a
		push(@as, $a);             	# append to return vector
	}
	print "\t[ A : @as ]\n";

	# my $p_pdl = undef;
	my $p = undef;
	my @p1 = ();my @p2 = ();my @p3 = ();my @p4 = ();
	for (my $s = 0; $s < 81; $s++){
		my $a = @as[$s];
		if ($a == 1){  
			foreach(keys(%{$prob_a1_hash{$s}})){
				push(@p1,$prob_a1_hash{$s}{$_});
			}print "\n";
			print "(1) --> ptmp ($s,1) : @p1\n";
			my $p1_pdl = pdl ( [@p1] );
			print "P1 PDL : $p1_pdl";
			$p = $p->glue(1,$p1_pdl);

		}elsif ($a == 2){  
			foreach(keys(%{$prob_a2_hash{$s}})){
				push(@p2,$prob_a2_hash{$s}{$_});
			}print "\n";
			print "(2) --> ptmp ($s,2) : @p2\n";
		}elsif ($a == 3){  
			foreach(keys(%{$prob_a3_hash{$s}})){
				push(@p3,$prob_a3_hash{$s}{$_});
			}print "\n";
			print "(3) --> ptmp ($s,3) : @p3\n";
		}elsif ($a == 4){  
			foreach(keys(%{$prob_a4_hash{$s}})){
				push(@p4,$prob_a4_hash{$s}{$_});
			}print "\n";
			print "(4) --> ptmp ($s,4) : @p4\n";
		}
	}
	my $p_pdl = pdl ( [@p1], [@p2], [@p3], [@p4] );
	print "P PDL : $p_pdl";

	my $q = $p * $v;
	my $q_slice = $q->slice('0,:,');
	my $return_v = $reward_pdl + $q_slice * 0.9875;
	print "$return_v\n";
	# exit;
	return $return_v;

}

sub compute_a{
	my $s = shift;
	my $v = shift;
	my @arr= ();
	# print "S : $s\nV : $v\n";

	# 1 #
		my $total = 0;
		for (my $sprime = 0; $sprime < 81; $sprime++){
			my $str = "$sprime" . ",:,";      
            my $v_tmp = sclr($v->slice($str));
            if (!($prob_a1_hash{$s}{$sprime} =~ /^-?\d+$/)){ $prob_a1_hash{$s}{$sprime}=0; }
			$total += $prob_a1_hash{$s}{$sprime} * $v_tmp;
		}
		push(@arr, $total);

	# 2 # 
		$total = 0;
		for (my $sprime = 0; $sprime < 81; $sprime++){
			my $str = "$sprime" . ",:,";      
            my $v_tmp = sclr($v->slice($str));
            if (!($prob_a2_hash{$s}{$sprime} =~ /^-?\d+$/)){ $prob_a2_hash{$s}{$sprime}=0; }
			$total += $prob_a2_hash{$s}{$sprime} * $v_tmp;
		}
		push(@arr, $total);

	# 3 # 
		$total = 0;
		for (my $sprime = 0; $sprime < 81; $sprime++){
			my $str = "$sprime" . ",:,";      
            my $v_tmp = sclr($v->slice($str));
            if (!($prob_a3_hash{$s}{$sprime} =~ /^-?\d+$/)){ $prob_a3_hash{$s}{$sprime}=0; }
			$total += $prob_a3_hash{$s}{$sprime} * $v_tmp;
		}
		push(@arr, $total);

	# 4 # 
		$total = 0;
		for (my $sprime = 0; $sprime < 81; $sprime++){
			my $str = "$sprime" . ",:,";      
            my $v_tmp = sclr($v->slice($str));
            if (!($prob_a4_hash{$s}{$sprime} =~ /^-?\d+$/)){ $prob_a4_hash{$s}{$sprime}=0; }
			$total += $prob_a4_hash{$s}{$sprime} * $v_tmp;
		}
		push(@arr, $total);

	my $max = -999;
	foreach(@arr){ 
		if($_+1 > $max){ $max = $_+1; }
	}
	print "ARR : @arr\n";
	print "\t( s: $s )\t[ A : $max ]\n";
	return $max;
}
