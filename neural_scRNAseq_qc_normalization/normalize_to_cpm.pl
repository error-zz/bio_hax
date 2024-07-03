#!/usr/bin/perl -w
# J. Craig Venter Institute
# Written by Jamison M. McCorrison (jmccorri@jcvi.org)
no warnings;
use strict;

my %Opts;
if (-e("./run.LOG.txt")){ system("rm -rf ./run.LOG.txt"); }
system("touch ./normalize.log.txt");
open LOG, ">> ./normalize.log.txt";

########################################################################################################################
# runsys : 
########################################################################################################################
sub runsys {
	my $cmd = shift;
	print "\n\t\t|| RUN CMD (16s_pipe) || $cmd\n";
	print LOG "|| RUN CMD (16s_pipe) || $cmd\n\n";
	system("$cmd");
}

########################################################################################################################
# printl : 
########################################################################################################################
sub printl {
	my $pr = shift;
	#run STDOUT print
	print $pr;
	#run LOG print (when appropriate)
	print LOG $pr;
}

########################################################################################################################
# MAIN : 
########################################################################################################################
MAIN : {	

    my $otu_table = shift;
    my $read_counts = shift;
    my $out = shift;

    print "\n|| Reading $read_counts ...\n";
    my $c = 0;
    my %lens = ();
    open(GC,"< $read_counts");
    while (defined(my $line = <GC>)){ #get each line
		chomp $line;
		my @xs = split(/\t/, $line);
		if ($c == 0){
		    # header
		}else{
		    # content
		    my $id = $xs[0];
		    $id =~ s/-/_/g;
		    $lens{$id} = $xs[1];
		    # print "LEN: $xs[0]\t$xs[1]\n";
		}
		$c++;
	}
	close GC;
    
    print "\n|| Reading $otu_table ...\n";
    open (NORMOUT, "> $out");
    $c = 0;
    my $skip_to = -1;
    my @header = ();
    open(GC,"< $otu_table");
    while (defined(my $line = <GC>)){ #get each line
		chomp $line;
		my @xs = split(/\t/, $line);
		if ($c == 0){
		    # header
		    my $d = 0;
		    foreach(@xs){
		        if ($skip_to == -1){
                    if (index($_, "ENSG") != -1){
                        $skip_to = $d;     
                        print "\nDEBUG: first count column = $skip_to ($_)\n";  
                    }  
		        }
		        $d++;
		        print NORMOUT "$_\t";
		    }
		    @header = @xs;
		    $c++;
		}else{
		    # content
		    my $d = 0;
		    foreach(@xs){
		        if ($d < $skip_to){
		            print NORMOUT "$_\t";
		        }else{
		            # CPM = [ (read.counts)/(num.input.fragments) ] * ( 10^6 )
		             my $id = $xs[0];
		            $id =~ s/-/_/g;
		            $id =~ s/\QLayer\E/layer/g;
		            my $length =  $lens{$id};
		            # print "$id : len = $length\tTC = $_\t";
		            my $CPM = ( $_ / $length ) * 1000000;
		            # print "CPM : $CPM\n";
		            print NORMOUT "$CPM\t";
		        }
		        $d++;
		    }
		}
		print NORMOUT "\n";
	}
	close GC;
	close NORMOUT;
	
    if (!(-s("$out"))){
        print "\n\n\n[ERROR] Failed to generate $out\n\n";
        exit;
    }else{
        print "\n\n\nSuccessful exit reached.\n\n";
        exit;
    }

}

