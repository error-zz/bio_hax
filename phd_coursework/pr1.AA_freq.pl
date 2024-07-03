#!/usr/bin/perl -w

# j.jmccorrison - beng203 / cse283 - 4/21/2015
# example execution:
# perl ./pr1.AA_freq.pl yeast_aaseqs.txt > yeast_aaseqs.aa_frequency.txt

use strict;
no warnings;

my $inputfile = shift;
my $c = 0;
open (FILE_IN, "$inputfile");

# PARSE INPUT FILE

my $cur_id;
my %FREQ_HASH;
my %LEN_HASH;

# BUILD HEADER
my @all_ids = ("A", "R", "N", "D", "B", "C", "Q", "E", "Z", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V");
print " \t";
foreach(@all_ids){ print "$_\t"; }
print "\n";

while (defined(my $line = <FILE_IN>)){ 
	chomp $line;
	
	# PARSE AND COUNT FREQUENCIES 
	
	if ($line =~m/^>/){
	
	    # get column 1 by stripping the special char from the fasta ID
	    
	    $cur_id = reverse $line;
	    chop($cur_id);
	    $cur_id = reverse $cur_id;
	    
	    # INITIALIZE HASH
	    
	    $LEN_HASH{$cur_id} = 0;
	    foreach(@all_ids){ $FREQ_HASH{$cur_id}{$_} = 0 }
	    
	}else{
	
	    # get column 2...n by capturing AA frequency
	    
	    for (my $r=0; $r < length($line); $r++){
	        my $cur_char = substr($line, $r, 1);
	        
	        if (!($FREQ_HASH{$cur_id}{$cur_char})){
	            $FREQ_HASH{$cur_id}{$cur_char} = 1;
	        }else{
	            $FREQ_HASH{$cur_id}{$cur_char}++;
	        }
	        # print "$cur_char ( $cur_id - $FREQ_HASH{$cur_id}{$cur_char} )\n";
	    }
	    $LEN_HASH{$cur_id} += length($line);
	}
	
}

# SORT ARRAY AND PRINT RATIOS TO TOTAL COUNTS

foreach my $id (keys (%FREQ_HASH)){
    print "$id\t";
    foreach my $aa (@all_ids){
        #print "($aa)";
        my $perc = $FREQ_HASH{$id}{$aa}/$LEN_HASH{$id};
        $perc = sprintf("%.3f", $perc);
        print "$perc\t";
        #print "($id ($LEN_HASH{$id}) - $aa ) : $FREQ_HASH{$id}{$aa} ($perc)\n";
    }
    print "\n";
}