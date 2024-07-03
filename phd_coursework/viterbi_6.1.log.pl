#!/usr/bin/perl -w
use strict;
use PDL; # using perl data language (numpy equivalent for perl)
         # if this fails to open, install PDL and run 'setup_bash' in the install directory to
         # configure your run environment for execution!
         # source /Volumes/Mac\ HD/Applications/PDL/setup_bash
no warnings;

#  https://courses.engr.illinois.edu/cs447/Slides/Lecture07HO.pdf

############### < INPUT PARSING > ###############

if (!(-s("./initialStateDistribution.txt"))){ print "\tNOTE! THIS SCRIPT ASSUMES ALL INPUT FILES ARE PROVIDED IN THE CURRENT WORKING DIRECTORY\nExiting.\n"; }

print "\n\nInitializing...\n";
my @alphabet = ("-","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z");

# INPUT 1 : INITIAL STATE DISTRIBUTION # 
my @initial_state_distribution = ();
open (FILE_IN, "./initialStateDistribution.txt");
while (defined(my $line = <FILE_IN>)){ chomp $line; push(@initial_state_distribution, $line); }
my $initial_state_distribution_pdl = pdl ([ @initial_state_distribution ]);
# print "\n\nINITIAL STATE DISTRIBUTION : $initial_state_distribution_pdl\n";

# INPUT 2 : EMISSION MATRIX # 
my %emissionMatrix = (); my $c = 1;
my %emissionMatrix_hash = ();
open (FILE_IN, "./emissionMatrix.txt");
while (defined(my $line = <FILE_IN>)){ 
    chomp $line; 
    my @xs = split(/\s+/, $line);
    $emissionMatrix{$c} = [ @xs ];
    $emissionMatrix_hash{$c}{0} = $xs[0];
    $emissionMatrix_hash{$c}{1} = $xs[1];
    $c++;
}
 my $front = defined;
 my $emissionMatrix_pdl = undef;
 foreach my $id (keys %emissionMatrix){
     my $row = $emissionMatrix{$id};
     if (defined $front)
     {
         $front = undef;
         $emissionMatrix_pdl = pdl ( [@$row] );
     }
     my $p_new = pdl ( [@$row] );
     $emissionMatrix_pdl = $emissionMatrix_pdl->glue(1,$p_new);
 }
print "\n\nEMISSION MATRIX : $emissionMatrix_pdl\n";

# INPUT 3 : TRANSITION MATRIX #
my %transitionMatrix = (); 
my %transitionMatrix_hash = ();
$c = 1;
open (FILE_IN, "./transitionMatrix.txt");
while (defined(my $line = <FILE_IN>)){ 
    chomp $line; 
    my @xs = split(/\s+/, $line);
    $transitionMatrix{$c} = [ @xs ];
    my $q =1;
    foreach(@xs){
        $transitionMatrix_hash{$c}{$q} = $_;
        $q++;
    }
    #$transitionMatrix_hash{$c}{0} = $xs[0];
    #$transitionMatrix_hash{$c}{1} = $xs[1];
    $c++;
}
my $front = defined;
my $transitionMatrix_pdl = undef;
foreach my $id (keys %transitionMatrix){
    my $row = $transitionMatrix{$id};
    if (defined $front)
    {
        $front = undef;
        $transitionMatrix_pdl = pdl ( [@$row] );
    }
    my $p_new = pdl ( [@$row] );
    $transitionMatrix_pdl = $transitionMatrix_pdl->glue(1,$p_new);
}
print "\n\nTRANSITION MATRIX : $transitionMatrix_pdl\n";

# INPUT 3 : OBSERVATION MATRIX #
my %observationMatrix = (); 
my @observationMatrix_arr = ();
$c = 1;
open (FILE_IN, "./observations.txt");
while (defined(my $line = <FILE_IN>)){ 
    chomp $line; 
    my @xs = split(/\s+/, $line);
    $observationMatrix{$c} = [ @xs ];
    @observationMatrix_arr = @xs;
    $c++;
}
# my $front = defined;
# my $observationMatrix_pdl = undef;
# foreach my $id (keys %observationMatrix){
#     my $row = $observationMatrix{$id};
#     if (defined $front)
#     {
#         $front = undef;
#         $observationMatrix_pdl = pdl ( [@$row] );
#     }
#     my $p_new = pdl ( [@$row] );
#     $observationMatrix_pdl = $observationMatrix_pdl->glue(1,$p_new);
# }
# print "\n\nOBSERVATION MATRIX : $observationMatrix_pdl\n";

############### < \INPUT PARSING > ##############
############### < MAIN > ########################

my %viterbi_trellis = ();
my %viterbi_backpointer = ();
print "\n\nInitialize...\n";
for (my $t = 1; $t <= 26; $t++){
    $viterbi_trellis{1}{$t} = $initial_state_distribution[$t-1] * $emissionMatrix_hash{$t}{$observationMatrix_arr[0]};
    print "\t(1)($t) $viterbi_trellis{1}{$t}\n";
}

print "\n\nRecursion...\n";
my $n = 150000;

# for (my $i = 2; $i <= $n; $i++){
#     for (my $t = 1; $t <= 26; $t++){
#         $viterbi_trellis{$i}{$t} = 1;
#     }
# }
for (my $i = 2; $i <= $n; $i++){
    if ($i % 10000 == 0){ print "\tRecusion counter : $i\n"};
    for (my $t = 1; $t <= 26; $t++){
        for (my $qt = 1; $qt <= 26; $qt++){
            my $tmp = $viterbi_trellis{$i-1}{$qt} * $transitionMatrix_hash{$qt}{$t};
            # print "\t\t$viterbi_trellis{$i-1}{$qt} * $transitionMatrix_hash{$qt}{$t} == $tmp\n";
            # print "($i)($t)->($qt) $tmp vs. $viterbi_trellis{$i}{$t}\n";
            if ($tmp > $viterbi_trellis{$i}{$t}){
                $viterbi_trellis{$i}{$t} = $tmp; 
                $viterbi_backpointer{$i}{$t} = $qt;
            #     print "* ($qt)";
            } 
            # print "\n";
        }
        # print "| t $t |  | i $i |  | obmatrix $observationMatrix_arr[$i] |\n";
        my $tmp2 = $viterbi_trellis{$i}{$t} * $emissionMatrix_hash{$t}{$observationMatrix_arr[$i]}; 
        print "! BEST !\t($i)($t)->($viterbi_backpointer{$i}{$t}) $tmp2 = $viterbi_trellis{$i}{$t} * $emissionMatrix_hash{$t}{$observationMatrix_arr[$i]} ( BP = $viterbi_backpointer{$i}{$t} = $alphabet[$viterbi_backpointer{$i}{$t}] )\n";
        $viterbi_trellis{$i}{$t} = $tmp2;
        #print "\n";
    }
    # exit;
}


foreach my $xx (sort {$a <=> $b} keys %viterbi_trellis){
     foreach my $yy (sort {$a <=> $b} keys %{$viterbi_trellis{$xx}}){
         print "$viterbi_trellis{$xx}{$yy} ";
     }
     print "\n";
 }

print "\n\nFinish...\n";
my $tmax;
my $vitmax = 0;
for (my $t = 1; $t <= 26; $t++){
    if ($viterbi_trellis{$n-1}{$t} > $vitmax){ 
        $tmax = $t;
        $vitmax = $viterbi_trellis{$n-1}{$t};
        print "\t($n-1)($t)->$viterbi_trellis{$n}{$t} (tmax = $tmax) (vitmax = $vitmax)\n";
    }
}
print "End point : tmax = $tmax , vitmax = $vitmax\n";

exit;

print "\n\nBacktrack...\n";
my $i = $n;
my %string = ();
while ($i > 0){
    $string{$i} = $tmax;
    $tmax = $viterbi_backpointer{$i}{$tmax};
    $i--;
}

print "\n\nResult:\n";
foreach (sort {$a <=> $b} keys (%string)){
    print "($_) $string{$_} -> $alphabet[$string{$_}]\n";
}
print "\n\nExit.\n";




