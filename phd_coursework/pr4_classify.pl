#!/usr/bin/perl -w

# j.jmccorrison - beng203 / cse283 - 4/21/2015
# example execution:
# perl ./pr1.AA_freq.pl yeast_aaseqs.txt > yeast_aaseqs.aa_frequency.txt

use strict;
no warnings;

 ############################# ############################# #############################
# INITIAL DECLARATIONS ############################# #############################
my %AA_CLUSTER_FREQ;
my %MICRO_CLUSTER_FREQ;
my %PHYLO_CLUSTER_FREQ;
my %micro_ids;
my %phylo_ids;
my @classification_search_terms = ("actin", "cyto", "endo", "Golgi", "micro", "mito", "nucl", "perox", "spind", "vacu");
# my @classification_search_terms = ("actin", "endo", "Golgi", "micro", "mito", "perox", "spind", "vacu");

 ############################# ############################# #############################
# INITIALIZE ALL ############################# #############################
for(my $c = 0; $c <=9; $c++){
	foreach(@classification_search_terms){  
	    $AA_CLUSTER_FREQ{$c}{$_} = 0; 
	    $MICRO_CLUSTER_FREQ{$c}{$_} = 0; 
	    $PHYLO_CLUSTER_FREQ{$c}{$_} = 0;
	}
}

 ############################# ############################# #############################
# POPULATE ALL ############################# #############################
open (FILE_IN, "yeast_aaseqs.aa_frequency_K_G10.kgg");
while (defined(my $line = <FILE_IN>)){ 
	chomp $line;
	my @xs = split(/\s+/, $line);
	foreach(@classification_search_terms){  
	    #print "$xs[1] - $xs[0] $_ -";
	    if (index($xs[0], $_) != -1){
	        # match found
	        $AA_CLUSTER_FREQ{$xs[1]}{$_}++;
	        #print "!!!!! (  $AA_CLUSTER_FREQ{$xs[1]}{$_} )\n";
	        last;
	    }
	    #print "\n";
	}
}
open (FILE_IN, "yeast_microarraydata2.normalized_kmeans10_K_G10.kgg");
while (defined(my $line = <FILE_IN>)){ 
	chomp $line;
	my @xs = split(/\s+/, $line);
	foreach(@classification_search_terms){  
	    if (index($xs[0], $_) != -1){
	        # match found
	        $MICRO_CLUSTER_FREQ{$xs[1]}{$_}++;
	        push(@{$micro_ids{$xs[1]}}, $xs[0]);
	        last;
	    }
	}
}
close FILE_IN;
open (FILE_IN, "yeast_phyloprofiles2.kmeans10_K_G10.kgg");
while (defined(my $line = <FILE_IN>)){ 
	chomp $line;
	my @xs = split(/\s+/, $line);
	foreach(@classification_search_terms){  
	    if (index($xs[0], $_) != -1){
	        # match found
	        $PHYLO_CLUSTER_FREQ{$xs[1]}{$_}++;
	        push(@{$phylo_ids{$xs[1]}}, $xs[0]);
	        last;
	    }
	}
}
close FILE_IN;

 ############################# ############################# #############################
# DEBUG AND CALCULATE TOTALS ############################# #############################
my %micro_type;
my %phylo_type;
my %micro_local_counts;
my %phylo_local_counts;
foreach my $local (@classification_search_terms){
     $micro_local_counts{$local} = 0;
        $phylo_local_counts{$local} = 0;
}
foreach my $cluster (sort {$a <=> $b} keys (%AA_CLUSTER_FREQ)){
    print "--> CLUSTER No. $cluster\t\t|\n";
    my $sumA = 0; my $sumB = 0; my $sumC = 0;
    foreach my $local (@classification_search_terms){
            #print "$local\t| $AA_CLUSTER_FREQ{$cluster}{$local}\t| $MICRO_CLUSTER_FREQ{$cluster}{$local}\t| $PHYLO_CLUSTER_FREQ{$cluster}{$local}\t|\n";
            $sumA += $AA_CLUSTER_FREQ{$cluster}{$local};
            $sumB += $MICRO_CLUSTER_FREQ{$cluster}{$local};
            $sumC += $PHYLO_CLUSTER_FREQ{$cluster}{$local};
            $micro_local_counts{$local}+=$MICRO_CLUSTER_FREQ{$cluster}{$local};
            $phylo_local_counts{$local}+=$PHYLO_CLUSTER_FREQ{$cluster}{$local};
    }
    #print "---------------------------------\nTOTAL\t| $sumA\t| $sumB\t| $sumC\t|\n---------------------------------\n";
    
    # CONVERT TO RELATIVE FREQUENCY
    my $micro_max = 0;
    my $phylo_max = 0;
    foreach my $local (@classification_search_terms){
        if (!$micro_local_counts{$local}){ $micro_local_counts{$local} = 999; }
        if ($MICRO_CLUSTER_FREQ{$cluster}{$local}/($micro_local_counts{$local}/$sumB) > $micro_max){
            $micro_max = $MICRO_CLUSTER_FREQ{$cluster}{$local};
            $micro_type{$cluster} = $local;
        }
        if (!$phylo_local_counts{$local}){ $phylo_local_counts{$local} = 999; }
        if ($PHYLO_CLUSTER_FREQ{$cluster}{$local}/($phylo_local_counts{$local}/$sumC) > $phylo_max){ 
            $phylo_max = $PHYLO_CLUSTER_FREQ{$cluster}{$local};
            $phylo_type{$cluster} = $local;
        }
        $AA_CLUSTER_FREQ{$cluster}{$local} = sprintf("%.3f", $AA_CLUSTER_FREQ{$cluster}{$local}/$sumA);
        $MICRO_CLUSTER_FREQ{$cluster}{$local} = sprintf("%.3f", $MICRO_CLUSTER_FREQ{$cluster}{$local}/$sumB);
        $PHYLO_CLUSTER_FREQ{$cluster}{$local} = sprintf("%.3f", $PHYLO_CLUSTER_FREQ{$cluster}{$local}/$sumC);
        print "$local\t| $AA_CLUSTER_FREQ{$cluster}{$local}\t| $MICRO_CLUSTER_FREQ{$cluster}{$local}\t| $PHYLO_CLUSTER_FREQ{$cluster}{$local}\t|\n";
    }
}


 ############################# ############################# #############################
# ID CONTENT ############################# #############################

my %micro_id_match;
my %phylo_id_match;
foreach my $clusterA (sort {$a <=> $b} keys (%MICRO_CLUSTER_FREQ)){
    my $size = @{$micro_ids{$clusterA}};
    my $max = 0;
    foreach my $clusterB (sort {$a <=> $b} keys (%PHYLO_CLUSTER_FREQ)){
        my $score = 0;
        my $counter = 0;
        foreach my $a ( @{$micro_ids{$clusterA}} ){
            foreach my $b ( @{$phylo_ids{$clusterB}} ){
                if ($a eq $b){ $score++; last; }
                #$counter++;
            }
        }
        $score = $score/$size;
        $score = sprintf("%.3f", $score);
        $micro_id_match{$clusterA}{$clusterB} = $score;
        $phylo_id_match{$clusterB}{$clusterA} = $score;
        print "||| $clusterA |  $clusterB | $score |||";
        if ($score > $max){
            $max = $score;
            #$micro_id_match{$clusterA} = $clusterB;
            #$phylo_id_match{$clusterB} = $clusterA;
            print "*";
        }
        print "\n";
    }
}


 ############################# ############################# #############################
# COMPARE MICROARRAY AND PHYLO CLUSTERING ############################# #############################

my %cross_perc_AB;

my %diff_hash;
foreach my $clusterA (sort {$a <=> $b} keys (%MICRO_CLUSTER_FREQ)){
    foreach my $clusterB (sort {$a <=> $b} keys (%PHYLO_CLUSTER_FREQ)){
        $diff_hash{$clusterA}{$clusterB} = 0;
    }
}
foreach my $clusterA (sort {$a <=> $b} keys (%MICRO_CLUSTER_FREQ)){
    foreach my $clusterB (sort {$a <=> $b} keys (%PHYLO_CLUSTER_FREQ)){
        my $diff_sum = 0;
        foreach my $local (@classification_search_terms){
            my $diff = abs($MICRO_CLUSTER_FREQ{$clusterA}{$local} - $PHYLO_CLUSTER_FREQ{$clusterB}{$local});
            $diff = sprintf("%.3f", $diff);
            #my $mod_diff = (1-$MICRO_CLUSTER_FREQ{$clusterA}{$local})*$diff;
            my $mod_diff = (1-$MICRO_CLUSTER_FREQ{$clusterA}{$local})*$diff;
            $mod_diff = sprintf("%.3f", $mod_diff);
            #print "| $clusterA\t| $clusterB\t| $local\t | $diff\t| $mod_diff\n";
            $diff_hash{$clusterA}{$clusterB} += $mod_diff;
        }
    }
}
my %micro_pairs;
my %phylo_pairs;
foreach my $clusterA (sort {$a <=> $b} keys (%MICRO_CLUSTER_FREQ)){
        my $min = 9999;
    foreach my $clusterB (sort {$a <=> $b} keys (%PHYLO_CLUSTER_FREQ)){
        if ($diff_hash{$clusterA}{$clusterB} < $min){
            $micro_pairs{$clusterA} = $clusterB;
            $phylo_pairs{$clusterB} = $clusterA;
            $min = $diff_hash{$clusterA}{$clusterB};
            print "| $clusterA\t| $clusterB\t| $diff_hash{$clusterA}{$clusterB}\t| -!!! \n";
        }
    }
}

 ############################# ############################# #############################
# COMPARE AND REDUCE  ############################# #############################

my %micro_classified;
my %phylo_classified;
my @types_used;
while (@types_used < 5){

    print "\n\nMICRO -> PHYLO:\n";
    print "\t\t     | class\t| id\t|\n";
    my $min = 999;
    my $max = 0;
    my $micro_min;
    my $micro_max;
    foreach my $clusterA (sort {$a <=> $b} keys (%micro_pairs)){
        my $tmp = abs($diff_hash{$clusterA}{$micro_pairs{$clusterA}} - $micro_id_match{$clusterA}{$micro_pairs{$clusterA}});
        print "$clusterA ($micro_type{$clusterA}) -> $micro_pairs{$clusterA} ($phylo_type{$micro_pairs{$clusterA}}) | $diff_hash{$clusterA}{$micro_pairs{$clusterA}}\t| $micro_id_match{$clusterA}{$micro_pairs{$clusterA}}\t| $tmp\t|";
        if ($diff_hash{$clusterA}{$micro_pairs{$clusterA}} < $min){
            $min = $diff_hash{$clusterA}{$micro_pairs{$clusterA}};
            $micro_min = $clusterA;
            print " *";
        }
        #if ($micro_id_match{$clusterA}{$micro_pairs{$clusterA}} > $max){
        if ($tmp > $max){
            $max = $micro_id_match{$clusterA}{$micro_pairs{$clusterA}};
            $micro_max = $clusterA;
            print " !";
        }
        print "\n";
    }
    print "\n\PHYLO -> MICRO:\n";
    print "\t\t     | class\t| id\t| cmpr\t|\n";
    my $min = 999;
    my $max = 0;
    my $phylo_min;
    my $phylo_max;
    foreach my $clusterA (sort {$a <=> $b} keys (%phylo_pairs)){
        my $tmp = abs($diff_hash{$clusterA}{$phylo_pairs{$clusterA}} * $phylo_id_match{$clusterA}{$phylo_pairs{$clusterA}});
        print "$clusterA ($phylo_type{$clusterA}) -> $phylo_pairs{$clusterA} ($phylo_type{$phylo_pairs{$clusterA}}) | $diff_hash{$clusterA}{$phylo_pairs{$clusterA}}\t| $phylo_id_match{$clusterA}{$phylo_pairs{$clusterA}}\t| $tmp\t|";
        if ($diff_hash{$clusterA}{$micro_pairs{$clusterA}} < $min){
            $min = $diff_hash{$clusterA}{$phylo_pairs{$clusterA}};
            $phylo_min = $clusterA;
            print " *";
        }
        #if ($phylo_id_match{$clusterA}{$phylo_pairs{$clusterA}} > $max){
        if ($tmp > $max){
            $max = $phylo_id_match{$clusterA}{$phylo_pairs{$clusterA}};
            $phylo_max = $clusterA;
            print " !";
        }
        print "\n";
            
    }
####### UPDATE TO CONTINUE CLEANUP
    
    my $a = $micro_pairs{$micro_max};
    my $b = $phylo_pairs{$a};
    print "A : $a - B : $b\n";
    my $bool = 0;
    if ($a == $b){
        print "ideal match found (micro=$micro_max) (phylo=$b) - | $phylo_type{$micro_max} |\n"; 
        $micro_classified{$micro_max} = $phylo_type{$micro_max};
        $phylo_classified{$a} = $phylo_type{$a};
        my $zbool =0;
        foreach(@types_used){ if ($_ eq $phylo_type{$micro_max}){$zbool=1;}  }
        if ($zbool == 0){ push(@types_used, $phylo_type{$micro_max}); }
        delete $micro_pairs{$micro_max};
        delete $phylo_pairs{$a};
        $bool = 1;
    }else{
        my $a = $micro_pairs{$micro_min};
        my $b = $phylo_pairs{$a};
        print "A : $a - B : $b\n";
        if ($a == $b){
            print "secondary match found (micro=$a) (phylo=$b) - | $phylo_type{$b} |\n"; 
            $micro_classified{$micro_min} = $phylo_type{$a};
            $phylo_classified{$a} = $phylo_type{$a};
            my $zbool =0;
            foreach(@types_used){ if ($_ eq $phylo_type{$micro_max}){$zbool=1;}  }
            if ($zbool == 0){ push(@types_used, $phylo_type{$micro_max}); }
            delete $micro_pairs{$a};
            delete $phylo_pairs{$a};
            $bool = 1;
        }
    }
    if ($bool == 0){
        print "no more easy matches!\n";
        
        
        # RERUN CLASSIFICATION!
        
foreach my $local (@classification_search_terms){
    foreach(@types_used){ if ($local eq $_){next;} }
     $micro_local_counts{$local} = 0;
        $phylo_local_counts{$local} = 0;
}
foreach my $cluster (sort keys (%AA_CLUSTER_FREQ)){
    print "--> CLUSTER No. $cluster\t\t|\n";
    my $sumA = 0; my $sumB = 0; my $sumC = 0;
    foreach my $local (@classification_search_terms){
            #print "$local\t| $AA_CLUSTER_FREQ{$cluster}{$local}\t| $MICRO_CLUSTER_FREQ{$cluster}{$local}\t| $PHYLO_CLUSTER_FREQ{$cluster}{$local}\t|\n";
            $sumA += $AA_CLUSTER_FREQ{$cluster}{$local};
            $sumB += $MICRO_CLUSTER_FREQ{$cluster}{$local};
            $sumC += $PHYLO_CLUSTER_FREQ{$cluster}{$local};
            $micro_local_counts{$local}+=$MICRO_CLUSTER_FREQ{$cluster}{$local};
            $phylo_local_counts{$local}+=$PHYLO_CLUSTER_FREQ{$cluster}{$local};
    }
}
    
        exit;
    }
        
}