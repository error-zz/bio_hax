#!/usr/bin/perl -w

# j.jmccorrison - beng203 / cse283 - 4/21/2015
# example execution:
# perl ./pr1.AA_freq.pl yeast_aaseqs.txt > yeast_aaseqs.aa_frequency.txt

use strict;
no warnings;

system("cp yeast_microarraydata2.normalized_kmeans10_K_G10.kgg A.kgg");
system("cp yeast_phyloprofiles2.kmeans10_K_G10.kgg B.kgg");

# INITIAL RUN 

system("rm out.txt");
system("rm remove_class.txt");
system("rm allout.txt");
system("./pr4_classify.4.pl");
my $tmp = `cat out.txt | awk \'{print \$3}\' | sort | uniq > remove_class.txt`;
chomp $tmp;
print "$tmp\n";
system("cat out.txt >> allout.txt");


#############

my @classes;
if (-s("remove_class.txt")){
    open (FILE_IN, "remove_class.txt");
    while (defined(my $line = <FILE_IN>)){ 
        my @xs = split(/\s+/, $line);
        if ($xs[0]){ push (@classes, $xs[0]); }
    }
    close FILE_IN;
}else{
    print "NOTHING TO REMOVE - FAIL CASE!\nExit.\n";
    exit;
}

open(PST_OUT,">./newA.txt") or die "Could not open ./newA.txt!";    
open (FILE_IN, "A.kgg");
    while (defined(my $line = <FILE_IN>)){ 
        chomp $line;
        my $bool = 0;
        my @xs = split(/\s+/, $line);
        foreach(@classes){
            #print "C: $xs[0], $_\n";
            if (index($xs[0], $_) != -1){ $bool = 1; }
        }
        if ($bool == 1){
            my @ys = split(/_/, $xs[0]);
            print PST_OUT "$ys[0] $xs[1]\n";
        }else{
            print PST_OUT "$line\n";
        }
    }
close PST_OUT;
close FILE_IN;
system("mv ./newA.txt A.kgg");

open(PST_OUT,">./newB.txt") or die "Could not open ./newB.txt!";    
open (FILE_IN, "B.kgg");
    while (defined(my $line = <FILE_IN>)){ 
        chomp $line;
        my $bool = 0;
        my @xs = split(/\s+/, $line);
        foreach(@classes){
            #print "C: $xs[0], $_\n";
            if (index($xs[0], $_) != -1){ $bool = 1; }
        }
        if ($bool == 1){
            my @ys = split(/_/, $xs[0]);
            print PST_OUT "$ys[0] $xs[1]\n";
        }else{
            print PST_OUT "$line\n";
        }
    }
close PST_OUT;
close FILE_IN;
system("mv ./newB.txt B.kgg");

open (FILE_IN, "out.txt");
    my $bool = 0;
    while (defined(my $line = <FILE_IN>)){ 
        chomp $line;
        my @xs = split(/\s+/, $line);
        if ($xs[0] eq "MICRO"){ $bool = 1; }
        elsif ($xs[0] eq "MICRO"){ $bool = 2; }
        else{
            if ($bool == 1){
                #my $cmd = "cat A.kgg | awk \'{if (\$2 != \$xs[0]) print}\' > A.kgg";
                #print "CMD : $cmd\n";
                system("cat A.kgg | awk \'{if (\$2 != \$xs[0]) print}\' > A.kgg"); 
            }elsif($bool == 2){
                system("cat B.kgg | awk \'{if (\$2 != \$xs[0]) print}\' > B.kgg");
            }
        }
    }