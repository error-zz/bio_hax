#!/usr/bin/perl -w
use strict;
no warnings;


##########################################################################
### MAIN #################################################################
##########################################################################

#############################
# INPUT HANDLING ############
#############################
my $inputfile = shift;
my $prefix = shift;

#system("cat $inputfile | tr \'-\' \' \' | tr \'>\' \' \' | tr \':\' \' \' > tmp_input.txt");

#open (FILE_IN, "tmp_input.txt");
open (FILE_IN, "$inputfile");
my %edges;
my ($source_node, $sink_node);

my $counter = 0;
while (defined(my $line = <FILE_IN>)){ #get each line
	chomp $line;
	if($counter == 0){
	    $source_node = $line;
	}elsif($counter == 1){
	    $sink_node = $line;
	 }else{
	    my @xs = split(/\s+/, $line);
	    my $starting_node = $xs[0];
	    my $receiving_node = $xs[1];
	    my $weight = $xs[2];
	    $edges{$starting_node}{$receiving_node} = $weight;
	}
	$counter++;
}

print "source : $source_node\n";
print "sink : $sink_node\n";
print "- FORWARD: -\n";
foreach my $starting_node (sort { $a <=> $b }  keys %edges){
    foreach my $receiving_node (sort { $a <=> $b }  keys($edges{$starting_node})){
        print "( $starting_node )( $receiving_node ) $edges{$starting_node}{$receiving_node}\n";
    }
}

#############################
# BEGIN #####################
#############################

my %backtracking_pointers = backpointer_edges( %edges );

my $longest_path_score = topologically_ordered_longest_path( $source_node, $sink_node, \%edges, \%backtracking_pointers );

    open OUTOUTOUT, "> ./tmp_out.txt";
    print OUTOUTOUT "$longest_path_score\n";

    my $outfile = "$prefix" . "_greedy_heaviest_simple_path.txt";
    system("mv tmp_out.txt $outfile");
    print "\n|C| = $longest_path_score\n\n";
    print "Run completed successfully.  Exiting.\n\n";
    close OUTOUTOUT;


exit;

#############################
# END #######################
#############################



##########################################################################
### topological_order_simple #############################################
##########################################################################
sub topological_order_simple{
    
    print "- TOPOLOGICAL ORDERING: -\n";
    
    my $source_node = shift;
    my %backtracking_pointers = @_;
    
    my %return;
    
    #$return{$source_node} = 0;
    #print "($source_node) $return{$source_node}\n";
    my $order_value = 1;
    
    my @all_node_ids = ();
    foreach my $rec_node_id (sort { $a <=> $b } keys(%backtracking_pointers)){ 
        foreach my $start_node_id (sort keys($backtracking_pointers{$rec_node_id})){ 
            my $checkA = 0;
            my $checkB = 0;
            foreach(@all_node_ids){
                if ($_ eq $start_node_id){ $checkA = 1; }
                elsif($_ eq $rec_node_id){ $checkB = 1; }
            }
            if ($checkA == 0){ push(@all_node_ids, $start_node_id) };
            if ($checkB == 0){ push(@all_node_ids, $rec_node_id) };
        }
    }
    
    foreach( sort { $a <=> $b } @all_node_ids ){
        $return{$_}=$order_value;
        $order_value++;
        print "($_) $return{$_}\n";
    }
    
    return %return;
    
}

##########################################################################
### topologically_ordered_longest_path ###################################
##########################################################################
sub topologically_ordered_longest_path{

    my $source_node = shift;
    my $sink_node = shift;
    print "SOURCE: $source_node\nSINK: $sink_node\n\n";
    my ( $edges_ref, $backtracking_pointers_ref ) = @_;
    my %edges = %$edges_ref;
    my %backtracking_pointers = %$backtracking_pointers_ref;
    
    my %top_order = topological_order_simple($source_node, %backtracking_pointers); 
    
    # initialize
    my %longest_path_score;
    foreach(sort keys %top_order){
        $longest_path_score{$_} = -9999;
    }
    $longest_path_score{$source_node} = 0;
    
    my $last;
    my %best_predecessors;
    foreach my $cur_node (sort { $a <=> $b } keys %top_order){
    print "-- CUR NODE $cur_node ---\n";
        if ($cur_node == $source_node){ $longest_path_score{$_} = 0; }
        else{                
            #max predecesoors of cur_node
            my $max = -9999;
            my $max_predecessor_id;
            
            if ($backtracking_pointers{$cur_node}){
                foreach(sort { $a <=> $b } keys $backtracking_pointers{$cur_node}){
                    if ( $backtracking_pointers{$cur_node}{$_} + $longest_path_score{$_} > $max ){ 
                        $max = $backtracking_pointers{$cur_node}{$_} + $longest_path_score{$_}; 
                        $max_predecessor_id = $_;
                    }
                    print "backtracking from (rec : $cur_node )<-(start: $_ ) : edge $backtracking_pointers{$cur_node}{$_} + longest path of $_ : $longest_path_score{$_} = $max (best= $max_predecessor_id)\n";
                }
            }else{
                # end of extraneous path
                #$max = $backtracking_pointers{$cur_node}{$_} - 9999;
                
                #$max = -9999;
                print "setting predecessor NULL <- $cur_node\n";
                $max_predecessor_id = "NULL";
            
            }
            $longest_path_score{$cur_node} = $max;
            if (($max_predecessor_id) && ($max_predecessor_id ne "NULL")){
                $best_predecessors{$cur_node} = $max_predecessor_id;
                print "setting predecessor $max_predecessor_id <- $cur_node\n";
            }
            
            $longest_path_score{$cur_node} = $max;
            $last = $cur_node;
         }   
    }
        
    print "\n\n$longest_path_score{$sink_node}\n";
    my $tmpstr = reverse $sink_node;
    my $str = "$tmpstr>-";
    print "$str\n";
    
    my $bool = 0;
    my $t = $sink_node;
    while ($bool == 0){
        #if ($best_predecessors{$t} == 0){ 
        if (!($best_predecessors{$t})){ 
            $tmpstr = reverse $best_predecessors{$t};
            $str .= "0"; 
            $bool = 1; 
            }
        else{ 
            $tmpstr = reverse $best_predecessors{$t};
            $str .= "$tmpstr>-"; 
        }
        $t = $best_predecessors{$t};
        print "$str\n";
    }
    #$str .= "$source_node";
    $str = reverse $str;
    print "$str\n";
    
    my $outfile = "$prefix" . "_greedy_heaviest_simple_path.txt";
    print "\nPrinting to file $outfile ...\n\n";
    if (-s("./tmp_out.txt")){ system("rm tmp_out.txt"); }   
    open OUTOUTOUT, ">> ./tmp_out.txt";
    print OUTOUTOUT "$longest_path_score{$sink_node}\n$str\n";
    close OUTOUTOUT;
    system("mv tmp_out.txt $outfile");
    #print "Run completed successfully.  Exiting.\n\n";

    #exit;
    return $longest_path_score{$last};
}


##########################################################################
### backpointer_edges #####################################
##########################################################################
sub backpointer_edges{

    my %edges = @_;
    my %backtracking_pointers;
    
    foreach my $starting_node (sort keys %edges){
        foreach my $receiving_node (sort keys($edges{$starting_node})){
            $backtracking_pointers{$receiving_node}{$starting_node} = $edges{$starting_node}{$receiving_node};
        }
    }
    
    print "- REVERSE: -\n";
    foreach my $receiving_node (sort  { $b <=> $a } keys %backtracking_pointers){
        foreach my $starting_node (sort  { $b <=> $a } keys($backtracking_pointers{$receiving_node})){
            print "( $receiving_node )( $starting_node ) $backtracking_pointers{$receiving_node}{$starting_node}\n";
        }
    }

    return %backtracking_pointers;
    
}
