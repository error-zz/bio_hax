#!/usr/bin/perl -w
use strict;

sub outhelp(){
	print "\n***\n compare_lists.pl_rapidsort input1.txt input2.txt output_prefix\n***\n";
	exit;
}

MAIN : {
	my $i1; my $i2; my $dummy; my $out;
	`CURRENTPATH=\$PWD`;

	$i1 = shift(@ARGV);
	$i2 = shift(@ARGV);
	$out = shift(@ARGV);
	
	####### CAPTURE INPUT LINES
	
	my $line;
	my @ids_i1;
	my @ids_i2;
	
	print "\n\n++ READING INPUT 1 ( $i1 )\n  ...\n";
	open (READ_IN, "$i1");
	while (defined($line = <READ_IN>)){ #get each line
		chomp $line;
		push (@ids_i1, "$line");
	}
	@ids_i1 = sort @ids_i1;
	
	print "\n\n++ READING INPUT 2 ( $i2 )\n  ...\n";
	open (READ_IN, "$i2");
	while (defined(my $line = <READ_IN>)){ #get each line
		chomp $line;
		push (@ids_i2, "$line");
	}
	@ids_i2 = sort @ids_i2;
	
	my $outfile = $out . '.both.list';
	print "$outfile\n";
	open (SEQOUTboth, "> $outfile");

	$outfile = $out . '.1.list';
	print "$outfile\n";
	open (SEQOUT1, "> $outfile");

	$outfile = $out . '.2.list';
	print "$outfile\n";
	open (SEQOUT2, "> $outfile");
	#  'a' cmp 'b' # -1
	#  'b' cmp 'a' #  1
	#  'a' cmp 'a' #  0
	
	print "++ STRING COMPARISON STARTED ON ALL ARRAYS\n  ...\n\n";
	
	my $coord_id2 = 0;
	for my $coord_id1 (0 .. $#ids_i1) {
		
		my $bool = 0;
		if (!($ids_i2[$coord_id2])){
	    	printl("\n++ End of file 2 reached.\n\nExiting.\n");
	    	exit;
	    }
	    	
	    my $eval_a;
	    my $eval_b;	
		while ($bool==0){
	    	#compare strings
	    	$eval_a = $ids_i1[$coord_id1] cmp $ids_i2[$coord_id2];
	    	#print "\n TEST : $ids_i1[$coord_id1] vs. $ids_i2[$coord_id2] equals $eval_a ... ";
	    	
	    	if ( $eval_a == -1 ){ #id1 is lesser, iterate id1
	    		print SEQOUT1 "$ids_i1[$coord_id1]\n";
	    		#print "1a ( $ids_i1[$coord_id1] ) ";
	    		$bool = 1;
	    		
	    	}elsif( $eval_a == 1){ #id1 is lesser, iterate id2
	    		$eval_b = $ids_i1[$coord_id1] cmp $ids_i2[$coord_id2];
	    		while($eval_b == 1){
	    			print SEQOUT2 "$ids_i2[$coord_id2]\n";
	    			#print "2a ( $ids_i2[$coord_id2] ) ";
	    			$coord_id2++;
	    			if (!($ids_i2[$coord_id2])){
	    				print "\n++ End of file 2 reached.\n\nExiting.\n";
	    				# exit routine
	    					#system("mv out.both.lst $out.foundinboth.lst");
							#system("mv out.1.lst $out.onlyINin1.lst");
							#system("mv out.2.lst $out.onlyINin2.lst");
	    					exit;
	    			}
	    			$eval_b = $ids_i1[$coord_id1] cmp $ids_i2[$coord_id2];
	    		}
	    		if ($eval_b == -1){ # id1 is lesser, iterate and move forward
	    			print SEQOUT1 "$ids_i1[$coord_id1]\n";
	    			#print "1b ( $ids_i1[$coord_id1] ) ";
	    			$bool = 1;
	    		}else{ # match found, write and move forward
	    			print SEQOUTboth "$ids_i1[$coord_id1]\n";
	    			#print "Bb ( $ids_i1[$coord_id1] ) ";
	    			$bool = 1;
	    			$coord_id2++;
	    		}
	    		
	    	}else{ # match found
	    		print SEQOUTboth "$ids_i1[$coord_id1]\n";
	    		#print "Ba ( $ids_i1[$coord_id1] ) ";
	    		$bool = 1;
	    		$coord_id2++;
	    	}
		}
    	
	}
	print "\n++ End of file 1 reached. Summarizing the rest of file 2 if required\n  ...\n";
	# print "BIN 2 : sample status : ID = $coord_id2, SIZE = $#ids_i2, STRING = $ids_i2[$coord_id2]\n\n";
	if ($coord_id2 <= $#ids_i2){
		my $orig = $coord_id2;
		for my $x ($orig .. $#ids_i2) {
			print SEQOUT2 "$ids_i2[$x]\n";
			#print "2b ( $ids_i2[$x] ) ";
		}
	}
	
	#system("mv out.both.lst $out.foundinboth.lst");
	#system("mv out.1.lst $out.onlyINin1.lst");
	#system("mv out.2.lst $out.onlyINin2.lst");
	
	print "\n\nDone.\n";
}
