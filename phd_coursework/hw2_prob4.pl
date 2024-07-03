#!/usr/bin/perl -w
# j.mccorrison 1/2015
use strict;
use warnings;

my $inputfile = shift;
my $african_x = shift;
my $african_y = shift;
my $euro_x = shift;
my $euro_y = shift;

my %distances_to_african_centroid;
my %distances_to_european_centroid;

open (FILE_IN, "$inputfile");
my $counter = 0;
while (defined(my $line = <FILE_IN>)){ #get each line
	chomp $line;
	if ($line){
	    my @xs = split(/,/, $line);
	    if ( ($xs[2] ne "European") && ($xs[2] ne "African") ){
	        $distances_to_african_centroid{$xs[2]} = sqrt( (abs($xs[0]-$african_x))^2 + (abs($xs[1]-$african_y))^2 );
	        print "$xs[2] : (afri) $distances_to_african_centroid{$xs[2]}\n";
	        $distances_to_european_centroid{$xs[2]} = sqrt( (abs($xs[0]-$euro_x))^2 + (abs($xs[1]-$euro_y))^2 );
	        print "$xs[2] : (euro) $distances_to_european_centroid{$xs[2]}\n";
	    }
	}
}

my @afri_cluster = ();
my @euro_cluster = ();
my @no_cluster = ();

# sort all keys


# associate with clusters
foreach my $q_a (sort { ($distances_to_african_centroid{$a} cmp $distances_to_african_centroid{$b}) || ($a cmp $b) } keys %distances_to_african_centroid) 
{
#foreach my $q_a ( sort values (%distances_to_african_centroid) ){
    my $afri_cluster_len = @afri_cluster;
    my $euro_cluster_len = @euro_cluster;
    if ( ($distances_to_african_centroid{$q_a} < $distances_to_european_centroid{$q_a}) && ($afri_cluster_len < 10) ){ 
        push(@afri_cluster, $q_a); 
    }
    elsif ( ($distances_to_european_centroid{$q_a} < $distances_to_african_centroid{$q_a}) && ($euro_cluster_len < 10) ){ 
        push(@euro_cluster, $q_a); 
    }
    elsif ( ($distances_to_european_centroid{$q_a} == $distances_to_african_centroid{$q_a}) && (($euro_cluster_len < 10) || ($afri_cluster_len < 10)) ){ 
        if (($euro_cluster_len < 10) && ($afri_cluster_len < 10)){
            # random
            my $select = rand(1);
            if ($select == 0){
                push(@afri_cluster, $q_a);
            }else{
                push(@euro_cluster, $q_a);
            }
        }elsif( $afri_cluster_len < 10 ){
            push(@afri_cluster, $q_a);
        }else{
            push(@euro_cluster, $q_a);
        }
    }else{
        push(@no_cluster, $q_a);
    }
}

print "\nAFRICAN SWALLOWS:\n(ID)\t\t(Dist. to African)\t(Dist. to European)\n";
foreach(@afri_cluster){ print "$_\t$distances_to_african_centroid{$_}\t$distances_to_european_centroid{$_}\n"; }
print "\nEUROPEAN SWALLOWS:\n(ID)\t\t(Dist. to African)\t(Dist. to European)\n";
foreach(@euro_cluster){ print "$_\t$distances_to_african_centroid{$_}\t$distances_to_european_centroid{$_}\n"; }
print "\nUNCLUSTERED SWALLOWS:\n(ID)\t\t(Dist. to African)\t(Dist. to European)\n";
foreach(@no_cluster){ print "$_\t$distances_to_african_centroid{$_}\t$distances_to_european_centroid{$_}\n"; }
print "\n";

exit;

