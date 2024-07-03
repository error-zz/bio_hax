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
my $f = shift;
my $prefix = shift;

print "INPUT FILE: $inputfile\n";
open (FILE_IN, "$inputfile");
my @snp_strings = ();
my %snpmatrix = ();
my $line; my $m = 0;
while (defined($line = <FILE_IN>)){ #get each line
	chomp $line;
    if ($m == 0){ $m = length($line); }
	push(@snp_strings, "$line");
}

my $n = `wc -l $inputfile | awk \'{print \$1}\'`; chomp $n;
$n++;
print "M = $m\nN = $n\n\n";
my %column_one_counter = ();
for (my $r = 0; $r < $m; $r++){
    $column_one_counter{$r} = 0;
}

### CREATE SNP MATRIX ###
my $c=0; my $d = 0;

foreach my $line (@snp_strings){
	my @chars = split(//, $line);
	$d=0;
	foreach(@chars){
        # and translate
		$snpmatrix{$d}{$c} = $_;
        if ($_ == 1){ $column_one_counter{$d}++; }
		$d++;
	}
	$c++;
}

my $c_soln = 0;
foreach(sort {$a <=> $b} keys{%column_one_counter}){
    my $ratio = $column_one_counter{$_}/$n;
    $ratio = sprintf("%.2f", $ratio);
    
    if ($ratio >= $f){
            print "YES ------ $_ : $column_one_counter{$_} ( f = $ratio )\n";
            $c_soln++;
        }else{
            print "no. ------ $_ : $column_one_counter{$_} ( f = $ratio )\n";
        }
}

print "\n|C| = $c_soln\n\nExiting.\n";

    my $outfile = "$prefix" . ".txt";
    print "\nPrinting to file $outfile ...\n\n";
    if (-s("./tmp_out.txt")){ system("rm tmp_out.txt"); }   

    open OUTOUTOUT, ">> ./tmp_out.txt";
    print OUTOUTOUT "|C| = $c_soln\n";
    close OUTOUTOUT;
    
    system("mv tmp_out.txt $outfile");
    print "Run completed successfully.  Exiting.\n\n";

# DEBUG PRINT
#printhash(%snpmatrix);

#my @columns_as_rows = ();
#foreach my $y (sort {$a <=> $b} keys(%snpmatrix)){
#    my $tmp_str = "";
#    foreach my $z (sort {$a <=> $b} keys( %{$snpmatrix{$y}} )){
#        $tmp_str .= "$snpmatrix{$y}{$z}";
#    }
#    push (@columns_as_rows, "$tmp_str");
#    print "colsasrows ($y) : $tmp_str\n";
#}



#############################################################################
sub printhash{
	# FOR DEBUG PURPOSES ONLY!

    my %hash = @_;
    
    my $q=0;
    print "   \t";
    foreach my $y (keys( %{$hash{0}} )){
    	print "$q ";
    	$q++;
    }
    print "\n";

    foreach my $y (sort {$a <=> $b} keys(%hash)){
    	print "$y :\t";
    	foreach my $z (sort {$a <=> $b} keys( %{$hash{$y}} )){
    		print "$hash{$y}{$z} ";
    	}
    	print "\n";
    }
}