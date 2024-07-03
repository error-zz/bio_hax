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
my $max_edge_weight = 0;
foreach my $starting_node (sort { $a <=> $b }  keys %edges){
    foreach my $receiving_node (sort { $a <=> $b }  keys($edges{$starting_node})){
        print "( $starting_node )( $receiving_node ) $edges{$starting_node}{$receiving_node}\n";
        if ($edges{$starting_node}{$receiving_node} > $max_edge_weight){ $max_edge_weight = $edges{$starting_node}{$receiving_node}; }
    }
}

#############################
# BEGIN MAIN  ###############
#############################

my $local_optima = simulated_annealing_heaviest_path( $source_node, $sink_node, \%edges );

my $outfile = "$prefix" . "_sim_anneal_heaviest_simple_path.txt";
system("mv tmp_out.txt $outfile");

print "\n\n|C| = $local_optima";
print "\n\nEnd of main.\n";
exit;


#############################
# END  ######################
#############################


##########################################################################
### simulated_annealing_heaviest_path ####################################
##########################################################################
sub simulated_annealing_heaviest_path{

	my $source_node = shift;
	# set sink node
    my $sink_node = shift;
    print "SOURCE: $source_node\nSINK: $sink_node\n\n";
    my ( $edges_ref ) = @_;
    my %edges = %$edges_ref;

	######	Initialize to starting node = ea
	my $cur_node = $source_node;
	
	######	Initialize heaviest path score = Smax = 0
	my $heaviest_path_score = 0;
	
	######	Initialize best path p = NULL 
	my @best_path = ();
	
	######	Back up initial weighted graph = Go
	my %initial_weighted_graph = %edges;

	######	UNTIL reach node eb OR time limit reached
	my $path_score = 0;
	my @path_contents = ();
	push(@path_contents, $cur_node);
	my $fail_counter = 10;
	until($cur_node == $sink_node){
		print "(*** CUR NODE ***) == $cur_node\n";
	
	############	FOR all vi∈ V, iterate through connected edges, C(vi) = vertices connected to vertex vi
		my @neighbors = keys( %{$edges{$cur_node}} );
		my $neighbor_count = @neighbors;
	
	##################	w(vi,vn) = weight of edge from current vertex to vertex vn … for each C(vi)
	##################	q(vi,vn) = probability of travelling to vertex vn along edge en
		my %neighbor_probabilities = {};
		my $sum = 0;
		foreach my $cur_neighbor (sort {$a <=> $b} @neighbors){
			$sum += $edges{$cur_node}{$cur_neighbor};
		}
		foreach my $cur_neighbor (sort {$a <=> $b} @neighbors){
			print "\t(*CN=$cur_node)(w= $cur_neighbor of $neighbor_count)\n";
			print "\t\t$cur_neighbor : w= $edges{$cur_node}{$cur_neighbor} ";
			# ! see biased_random_selector subroutine !
			$neighbor_probabilities{$cur_node}{$cur_neighbor} = $edges{$cur_node}{$cur_neighbor}/$sum;
			print "p= $neighbor_probabilities{$cur_node}{$cur_neighbor}\n";
		}

	############	IF edge leaving vertex vi to another neighboring vertex vn exists with q(vn) > 0, select new neighbor by selecting at pseudo-random with weighted probability q(vi) 
	##################	Modify graph (Gp) to update weight of all incoming edges to vertex vi = infinity (therefore, q(vi) = 0)
	##################	Update S = S+ w(vi,vn)
	############	ELSE IF q(vn) = 0 for all neighbors
	##################	IF length(p) = V-1
	########################	Select last neighbor = eb
	########################	Update S = S+ w(vi, eb)
	########################	IF S > Smax, set new best path
	##################  IF vn = eb (end node)
	########################	q(vi,eb) = 0 (temporarily)

		my $cur_neighbor_selection = biased_random_selector(\%neighbor_probabilities, \@path_contents);
		print "\n\tNeighbor selected : $cur_neighbor_selection\n";

	##################	ELSE (fail case)
	########################	Reset Gn to Go
	########################	Re-initialize starting node = ea
	########################	IF time limit reached, exit as fail state.  No path exists.

		if ($cur_neighbor_selection == $sink_node){
			push(@path_contents, $cur_neighbor_selection);
			my $path_size = @path_contents;
			my $hopeful_size = @neighbors + 1;
			print "~~~ PATH FOUND! ~~~ ($path_size vs. $hopeful_size) \n"; 
			foreach(@path_contents){ print "$_ -> " }
			print "fin.\n";
			if ($path_size == $hopeful_size){
				print "SUCCESS!\n"; 
				if ($path_score > $heaviest_path_score){ $heaviest_path_score = $path_score; }
				print "CURRENT PATH: $path_score\n";
				print "HEAVIEST PATH: $heaviest_path_score\n";
				$fail_counter = 10;
				last;
			}else{
				print "FAIL! TOO SHORT!\n";
				my $pop_iterator = int($fail_counter/10);
				print "POPPING BACK $pop_iterator VERTICES... (fc = $fail_counter) (pi = $pop_iterator)\n";
				if ($pop_iterator > ($path_size - 1)){
				#if ($pop_iterator == $path_size){
					$fail_counter = 10;
					$pop_iterator = 1;
				}else{
					$fail_counter++;
				}
				# print "POPPING BACK $pop_iterator VERTICES... (fc = $fail_counter)\n";
				for (my $y=0; $y < $pop_iterator; $y++){
					pop(@path_contents);
				}
				my $size = @path_contents;
				$cur_node = $path_contents[$size-1];
				#$fail_counter++;
				next;
			}
		}
	
		# DEBUG
		$cur_node = $cur_neighbor_selection;
		$fail_counter = 10;
		push(@path_contents, $cur_neighbor_selection);
	#END UNTIL
	}

	# compile score
	print "CALCULATING SCORE:\n";
	my $summary_path_score = 0;
	my $source = "NULL";
	my $sink = "NULL";
	foreach(@path_contents){
		if ($_ == $source_node){
			# source node
			$source = $_;
		#}elsif($_ == $sink_node){
			# sink node
		}else{
			# intermediate
			$sink = $_;
			$summary_path_score += $edges{$source}{$sink};
			print "$source -> $sink : $edges{$source}{$sink} += $summary_path_score\n";
		}
	}

    print "\nPrinting to file $outfile ...\n\n";
    if (-s("./tmp_out.txt")){ system("rm tmp_out.txt"); }   
    open OUTOUTOUT, ">> ./tmp_out.txt";
	print OUTOUTOUT "$summary_path_score\n";
	foreach(@path_contents){print OUTOUTOUT "$_ -> "; }
	print OUTOUTOUT "\n";

	return $summary_path_score;
} 


##################################################################################
sub biased_random_selector {
	print "\tSelecting at random...\n";
	my $pointer1 = shift;
	my %neighbor_probabilities = %$pointer1;
	my $pointer2 = shift;
	my @path_contents = @$pointer2;

	my $rand = rand(1);
	print "\t\tRAND: $rand\n";
	foreach my $cur_node (sort {$a <=> $b} keys( %neighbor_probabilities )){
		my $last;
		my $val = 0;
		print "\t\tCN : $cur_node\n";
		if (!($cur_node =~ /^\d+$/)){
			print "\t\t\t- skipping annoying hash reference!\n";
			next;
		}
		foreach my $cur_neighbor (sort {$a <=> $b} keys( %{$neighbor_probabilities{$cur_node}} )){
			
			# check for previous selection 
			print "\t\tPREV PATH: ";
			my $bool = 0;
			foreach(@path_contents){
				print "$_ ";
				if ($cur_neighbor == $_){ $bool = 1; print "<* "}
			}print "\n";
			if ($bool == 1){ next; }

			$val += $neighbor_probabilities{$cur_node}{$cur_neighbor};
			print "\t\t$cur_node -> $cur_neighbor : w= $neighbor_probabilities{$cur_node}{$cur_neighbor}\n\t\t\t$rand vs. $val\n";
			if ($neighbor_probabilities{$cur_node}{$cur_neighbor} > 0){
				if ($rand < $val){
					print "\t\t\t* selecting $cur_neighbor";
					return $cur_neighbor;
				}
			}
			$last = $cur_neighbor;

		}

		print "\t\t\t** selecting $last";
		return $last;
	}
}
