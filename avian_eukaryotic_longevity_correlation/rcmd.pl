#!/usr/bin/perl -w
use strict;

#my $in = shift;
#print $in;

my $count = 0;

my $in_flip = shift;
my $in_abs = shift;
my $in_out = shift;
#open(FLIP,"< ../mkorth/out.all_vs_all.mci.I14.flipList");
#open(ABUND, "< ../mktab/22.PintailDuck.rpkm");
#open(PRINTME,"> zzzout");
open(FLIP,"< $in_flip");
open(ABUND, "< $in_abs");
open(PRINTME,"> $in_out");

my %flip;
# print "\n\nFLIP\n";
while (defined(my $line = <FLIP>)){ #get each line
    chomp $line;
    my @xs = split(/\s+/, $line);
    $flip{$xs[0]} = $xs[1];
    #print "$xs[0] : $xs[1]\n";
}
    
my %abs;
my %qid;
my %rid;
# print "\n\nABUND\n";
while (defined(my $line = <ABUND>)){ #get each line
    chomp $line;
    my @xs = split(/\s+/, $line);
    $abs{$xs[0]} = $xs[1];
    $qid{$xs[0]} = $xs[2];
    $rid{$xs[0]} = $xs[3];
    #print "$xs[0] : $xs[1] : $xs[2] : $xs[3]\n";
}
#exit;

# open(TtoP,"< ../mkorth/Anas_platyrhynchos.transcript2prot");
my %ttop;
my @ref_array = ("Anas_platyrhynchos","Calypte_anna","Charadrius_vociferus","Columba_livia","Corvus_brachyrhynchos","Gallus_gallus","Meleagris_gallopavo","Nipponia_nippon","Picoides_pubescens","Struthio_camelus","Taeniopygia_guttata","Falcons_peregrine");
foreach my $r (@ref_array){
    my $inlocation = '../mkorth/' . $r . '.transcript2prot';
    open(TtoP,"< $inlocation");
    # print "\n\nTtoP : $r\n";
    while (defined(my $line = <TtoP>)){ #get each line
        chomp $line;
        my @xs = split(/\s+/, $line);
        $ttop{$r}{$xs[1]} = $xs[0];
        #print "$r : $xs[1] : $xs[0]\n";
    }
    close TtoP;
}
# exit;

my %report_out;
my @col_ids = ("Transcript_ID","Protein_ID","Abundance","Cluster_ID","Reference_ID","Query_ID",@ref_array); # ,"Num_Uniq_QuerySpecies_Represented");
foreach my $p ( keys(%flip) ){
    foreach my $ri ( @ref_array ){
        $report_out{$p}{$ri} = 0;
    }
    $report_out{$p}{"Abundance"} = 0;
}



# foreach protein id in fliplist
# foreach my $p ( keys(%flip) ){
foreach my $p ( keys(%flip) ){  # ENSMGAP00000013893.1	15363 #
    #$p = 'Apl_R001668';
    # print("\n\n$p :\n");
    # check if protein id is in each reference's transcript array
    # if( grep {$_ eq $p} keys(%abs) ){
        # speed up check by using strings we know
        my $s = substr( $p, 0, 3 );
        #print ($s);
        # print "$s\n";
        if($s eq "ENS"){
            # ncbi ref OR falcon

            print "\np _ $p\n";
            
            foreach my $ref ("Taeniopygia_guttata","Meleagris_gallopavo","Gallus_gallus","Falcons_peregrine"){
                if( exists( $ttop{$ref}{$p} ) ){
                    $report_out{$p}{$ref}++;
                    $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                    $report_out{$p}{"Protein_ID"} = $p;
                    $report_out{$p}{"Cluster_ID"} = $flip{$p};
                    if(exists($abs{$ttop{$ref}{$p}})){
                        $report_out{$p}{"Abundance"} += $abs{$ttop{$ref}{$p}};
                        $report_out{$p}{"Reference_ID"} = $rid{$ttop{$ref}{$p}};
                        $report_out{$p}{"Query_ID"} = $qid{$ttop{$ref}{$p}};
                        #exit;
                    }else{
                        $report_out{$p}{"Reference_ID"} = $ref;
                        $report_out{$p}{"Query_ID"} = "Undef";
                    }
                    $report_out{$p}{$ref}++;
                    print"(A) $ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";                    
                    last;
                #}elsif( ($input_pass=~/^$/) ){
                }else{
                    # fix big 
                     my $queryResponse = `grep $p $in_flip | head -n 1`; chomp $queryResponse;
                     # search for isoform specific version of query
                     if(length($queryResponse)>0){
                        print "(B) $ref : $queryResponse\n"; # $queryResponse\n";
                        my @xs = split(/\s+/, $queryResponse);
                        $p = $xs[0];
                        if( exists( $ttop{$ref}{$p} ) ){
                            print "$queryResponse : $p\n";
                            $report_out{$p}{$ref}++;
                            $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                            $report_out{$p}{"Protein_ID"} = $p;
                            $report_out{$p}{"Cluster_ID"} = $flip{$p};
                            print "\tLooking for $abs{$ttop{$ref}{$p}}\n";
                            if(exists($abs{$ttop{$ref}{$p}})){
                                $report_out{$p}{"Abundance"} += $abs{$ttop{$ref}{$p}};
                                $report_out{$p}{"Reference_ID"} = $rid{$ttop{$ref}{$p}};
                                $report_out{$p}{"Query_ID"} = $qid{$ttop{$ref}{$p}};
                                #exit;
                            }else{
                                $report_out{$p}{"Reference_ID"} = $ref;
                                $report_out{$p}{"Query_ID"} = "Undef";
                            }
                            $report_out{$p}{$ref}++;
                            print"(B) $ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";               
                            last;
                        }else{
                            $report_out{$p}{"Reference_ID"} = $ref;
                            $report_out{$p}{"Query_ID"} = "Undef";
                        }
                     }
                }
                #}else{
                if(!(exists($report_out{$p}{"Query_ID"} ))) {
                    print "Queried $ref and found none. (ENS)\n";
                }
                #}
                #print "\n\n\n\nDEBUG EXIT!\n\n";
                #exit;
            }
            
        }elsif($s eq "Apl"){
            my $ref = "Anas_platyrhynchos";
            if( exists( $ttop{$ref}{$p} ) ){
                $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                $report_out{$p}{"Protein_ID"} = $p;
                $report_out{$p}{"Cluster_ID"} = $flip{$p};
                if(exists($abs{$p})){
                    $report_out{$p}{"Abundance"} += $abs{$p};
                    $report_out{$p}{"Reference_ID"} = $rid{$p};
                    $report_out{$p}{"Query_ID"} = $qid{$p};
                }else{
                    $report_out{$p}{"Reference_ID"} = $ref;
                    $report_out{$p}{"Query_ID"} = "Undef";
                    # $report_out{$p}{"Query_ID"} = $qid{$p};
                }
                $report_out{$p}{$ref}++;
                #print"$ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";
                #exit;
            }else{
                # print "Queried $ref and found none. (ENS)\n";
            }
            #exit;
        }elsif($s eq "Aan"){
            my $ref = "Calypte_anna";
            if( exists( $ttop{$ref}{$p} ) ){
                $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                $report_out{$p}{"Protein_ID"} = $p;
                $report_out{$p}{"Cluster_ID"} = $flip{$p};
                if(exists($abs{$p})){
                    $report_out{$p}{"Abundance"} += $abs{$p};
                    $report_out{$p}{"Reference_ID"} = $rid{$p};
                    $report_out{$p}{"Query_ID"} = $qid{$p};
                }else{
                    $report_out{$p}{"Reference_ID"} = $ref;
                    $report_out{$p}{"Query_ID"} = "Undef";
                }
                $report_out{$p}{$ref}++;
                #print"$ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";
            }else{
                # print "Queried $ref and found none. (ENS)\n";
            }
        }elsif($s eq "Cvo"){
            my $ref = "Charadrius_vociferus";
            if( exists( $ttop{$ref}{$p} ) ){
                $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                $report_out{$p}{"Protein_ID"} = $p;
                $report_out{$p}{"Cluster_ID"} = $flip{$p};
                if(exists($abs{$p})){
                    $report_out{$p}{"Abundance"} += $abs{$p};
                    $report_out{$p}{"Reference_ID"} = $rid{$p};
                    $report_out{$p}{"Query_ID"} = $qid{$p};
                }else{
                    $report_out{$p}{"Reference_ID"} = $ref;
                    $report_out{$p}{"Query_ID"} = "Undef";
                }
                $report_out{$p}{$ref}++;
                #print"$ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";
            }else{
                # print "Queried $ref and found none. (ENS)\n";
            }
        }elsif($s eq "Cli"){
            my $ref = "Columba_livia";
            if( exists( $ttop{$ref}{$p} ) ){
                $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                $report_out{$p}{"Protein_ID"} = $p;
                $report_out{$p}{"Cluster_ID"} = $flip{$p};
                if(exists($abs{$p})){
                    $report_out{$p}{"Abundance"} += $abs{$p};
                    $report_out{$p}{"Reference_ID"} = $rid{$p};
                    $report_out{$p}{"Query_ID"} = $qid{$p};
                }else{
                    $report_out{$p}{"Reference_ID"} = $ref;
                    $report_out{$p}{"Query_ID"} = "Undef";
                }
                $report_out{$p}{$ref}++;
                #print"$ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";
            }else{
                # print "Queried $ref and found none. (ENS)\n";
            }
        }elsif($s eq "Cbr"){
            my $ref = "Corvus_brachyrhynchos";
            if( exists( $ttop{$ref}{$p} ) ){
                $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                $report_out{$p}{"Protein_ID"} = $p;
                $report_out{$p}{"Cluster_ID"} = $flip{$p};
                if(exists($abs{$p})){
                    $report_out{$p}{"Abundance"} += $abs{$p};
                    $report_out{$p}{"Reference_ID"} = $rid{$p};
                    $report_out{$p}{"Query_ID"} = $qid{$p};
                }else{
                    $report_out{$p}{"Reference_ID"} = $ref;
                    $report_out{$p}{"Query_ID"} = "Undef";
                }
                $report_out{$p}{$ref}++;
                #print"$ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";
            }else{
                # print "Queried $ref and found none. (ENS)\n";
            }
        }elsif($s eq "Nni"){
            my $ref = "Nipponia_nippon";
            if( exists( $ttop{$ref}{$p} ) ){
                $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                $report_out{$p}{"Protein_ID"} = $p;
                $report_out{$p}{"Cluster_ID"} = $flip{$p};
                if(exists($abs{$p})){
                    $report_out{$p}{"Abundance"} += $abs{$p};
                    $report_out{$p}{"Reference_ID"} = $rid{$p};
                    $report_out{$p}{"Query_ID"} = $qid{$p};
                }else{
                    $report_out{$p}{"Reference_ID"} = $ref;
                    $report_out{$p}{"Query_ID"} = "Undef";
                }
                $report_out{$p}{$ref}++;
                #print"$ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";
            }else{
                # print "Queried $ref and found none. (ENS)\n";
            }
        }elsif($s eq "Ppu"){
            my $ref = "Picoides_pubescens";
            if( exists( $ttop{$ref}{$p} ) ){
                $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                $report_out{$p}{"Protein_ID"} = $p;
                $report_out{$p}{"Cluster_ID"} = $flip{$p};
                if(exists($abs{$p})){
                    $report_out{$p}{"Abundance"} += $abs{$p};
                    $report_out{$p}{"Reference_ID"} = $rid{$p};
                    $report_out{$p}{"Query_ID"} = $qid{$p};
                }else{
                    $report_out{$p}{"Reference_ID"} = $ref;
                    $report_out{$p}{"Query_ID"} = "Undef";
                }
                $report_out{$p}{$ref}++;
                #print"$ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";
            }else{
                # print "Queried $ref and found none. (ENS)\n";
            }
        }elsif($s eq "Sca"){
           my $ref = "Struthio_camelus";
            if( exists( $ttop{$ref}{$p} ) ){
                $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                $report_out{$p}{"Protein_ID"} = $p;
                $report_out{$p}{"Cluster_ID"} = $flip{$p};
                if(exists($abs{$p})){
                    $report_out{$p}{"Abundance"} += $abs{$p};
                    $report_out{$p}{"Reference_ID"} = $rid{$p};
                    $report_out{$p}{"Query_ID"} = $qid{$p};
                }else{
                    $report_out{$p}{"Reference_ID"} = $ref;
                    $report_out{$p}{"Query_ID"} = "Undef";
                }
                $report_out{$p}{$ref}++;
                #print"$ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";
            }else{
                # print "Queried $ref and found none. (ENS)\n";
            }
        }elsif($s eq "per"){
            my $ref = "Falcons_peregrine";
            if( exists( $ttop{$ref}{$p} ) ){
                $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                $report_out{$p}{"Protein_ID"} = $p;
                $report_out{$p}{"Cluster_ID"} = $flip{$p};
                if(exists($abs{$p})){
                    $report_out{$p}{"Abundance"} += $abs{$p};
                    $report_out{$p}{"Reference_ID"} = $rid{$p};
                    $report_out{$p}{"Query_ID"} = $qid{$p};
                }else{
                    $report_out{$p}{"Reference_ID"} = $ref;
                    $report_out{$p}{"Query_ID"} = "Undef";
                }
                $report_out{$p}{$ref}++;
                #print"$ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";
            }else{
                # print "Queried $ref and found none. (ENS)\n";
            }
        }else{
            foreach my $ref (@ref_array){
                if( exists( $ttop{$ref}{$p} ) ){
                    $report_out{$p}{"Transcript_ID"} = $ttop{$ref}{$p};
                    $report_out{$p}{"Protein_ID"} = $p;
                    $report_out{$p}{"Cluster_ID"} = $flip{$p};
                    if(exists($abs{$p})){
                        $report_out{$p}{"Abundance"} += $abs{$p};
                        $report_out{$p}{"Reference_ID"} = $rid{$p};
                        $report_out{$p}{"Query_ID"} = $qid{$p};
                    }else{
                        $report_out{$p}{"Reference_ID"} = $ref;
                        $report_out{$p}{"Query_ID"} = "Undef";
                    }
                    $report_out{$p}{$ref}++;
                    #print"$ref\t: $report_out{$p}{\"Transcript_ID\"} $report_out{$p}{\"Protein_ID\"} $report_out{$p}{\"Abundance\"} $report_out{$p}{\"Cluster_ID\"} $report_out{$p}{\"Reference_ID\"} $report_out{$p}{\"Query_ID\"} $report_out{$p}{$ref}\n";
                    last;
                }else{
                    # print "Queried $ref and found none. (ALL)\n";
                }
            }
        }
    # }
}



foreach my $x ( keys(%report_out) ){
    if(exists($report_out{$x}{"Protein_ID"})){
    if(exists($report_out{$x}{"Query_ID"})){
    if($report_out{$x}{"Query_ID"} ne "Undef"){
        for my $y ( @col_ids ){
            print PRINTME "$report_out{$x}{$y}\t";
        }
        print PRINTME "\n";
    }
    }
    }
}


# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N ortholog_report -cwd sh run.sh
# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N ortholog_report -cwd perl rcmd.pl ../mkorth/out.all_vs_all.mci.I14.flipList ../mktab/22.PintailDuck.rpkm 22.PintailDuck.I14

# perl rcmd.pl ../mkorth/out.all_vs_all.mci.I14.flipList ../mktab/15.BarnSwallow.rpkm 15.BarnSwallow.I14

# perl rcmd.pl ../mkorth/out.all_vs_all.mci.I14.flipList ../mktab/61.Mute_Swan.rpkm 61.Mute_Swan.I14
 
# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N ortholog_report -cwd perl rcmd.pl ../mkorth/out.all_vs_all.mci.I14.flipList all.rpkm all.I14
 
# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N ortholog_report -cwd perl rcmd.pl ../mkorth/out.all_vs_all.mci.I17.flipList all.rpkm all.I17

# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N ortholog_report -cwd perl rcmd.pl ../mkorth/out.all_vs_all.mci.I20.flipList all.rpkm all.I20

# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N ortholog_report -cwd perl rcmd.pl ../mkorth/out.all_vs_all.mci.I25.flipList all.rpkm all.I25

# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N ortholog_report -cwd perl rcmd.pl ../mkorth/out.all_vs_all.mci.I30.flipList all.rpkm all.I30

# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N ortholog_report -cwd perl rcmd.pl ../mkorth/out.all_vs_all.mci.I40.flipList all.rpkm all.I40
 
