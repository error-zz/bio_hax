#!/usr/bin/perl -w
# J. Craig Venter Institute
# Written by Jamison M. McCorrison (jmccorri@jcvi.org)
no warnings;
use strict;

my $in = shift;

my $run = 'cat Sample_' . $in . '.interleaved_pairs.fastq | tr \' \' \'_\' > ' . $in . '.all.fastq';
print("$run\n");
system("$run");
my $expected = $in . '.all.fastq';
if(!-s("$expected")){ print("$expected not found"); exit; }

my $run = '/usr/local/devel/BCIS/assembly/tools/NeatFreq/lib/fastq_ID_fixer.pl -i ' . $in . '.all.fastq > ' . $in . '.idfix.fastq';
print("$run\n");
system("$run");
my $expected = $in . '.idfix.fastq';
if(!-s("$expected")){ print("$expected not found"); exit; }

my $run = 'grep \'@HISEQ\' ' . $in . '.all.fastq | grep \'_1\' > tmp.tmp.' . $in . '.R1.ids.txt';
print("$run\n");
system("$run");
my $expected = $in . '.R1.ids.txt';
if(!-e("$expected")){ print("$expected not found");}

my $run = 'grep \'@HISEQ\' ' . $in . '.all.fastq | grep \'_2\' > tmp.tmp.' . $in . '.R2.ids.txt';
print("$run\n");
system("$run");
my $expected = $in . '.R2.ids.txt';
if(!-e("$expected")){ 
    print("$expected not found");
}

my $run = 'cat tmp.tmp.' . $in . '.R1.ids.txt | tr \'_\' \'\\t\' | awk \'{print $1}\' | sort | uniq > tmp.tmp.' . $in . '.R1.uniq_ids.txt';
print("$run\n");
system("$run");
my $expected = $in . '.R1.uniq_ids.txt';
if(!-e("$expected")){ print("$expected not found\n");}

my $run = 'cat tmp.tmp.' . $in . '.R2.ids.txt | tr \'_\' \'\\t\' | awk \'{print $1}\' | sort | uniq > tmp.tmp.' . $in . '.R2.uniq_ids.txt';
print("$run\n");
system("$run");
my $expected = 'tmp.tmp.' . $in . '.R2.uniq_ids.txt';
if(!-e("$expected")){ print("$expected not found\n"); }

my $run = 'perl compare_lists_rapidsort.pl tmp.tmp.' . $in . '.R1.uniq_ids.txt tmp.tmp.' . $in . '.R2.uniq_ids.txt tmp' . $in;
 print("$run\n");
system("$run");
my $expected = 'tmp' . $in . '.both.list';
if(!-e("$expected")){ print("$expected not found\n");}

my $run = 'wc -l tmp' . $in . '.both.list | awk \'{print $1}\'';
my $tmp0 = `$run`; chomp $tmp0;
my $run = 'wc -l tmp' . $in . '.1.list | awk \'{print $1}\'';
my $tmp1 = `$run`; chomp $tmp1;
my $run = 'wc -l tmp' . $in . '.2.list | awk \'{print $1}\'';
my $tmp2 = `$run`; chomp $tmp2;
my $run = 'wc -l tmp.tmp.' . $in . '.R1.uniq_ids.txt | awk \'{print $1}\'';
my $tmp3 = `$run`; chomp $tmp3;
my $run = 'wc -l tmp.tmp.' . $in . '.R2.uniq_ids.txt | awk \'{print $1}\'';
my $tmp4 = `$run`; chomp $tmp4;

print "$tmp0 $tmp1 $tmp2 $tmp3 $tmp4\n";
if($tmp0 == 0 && $tmp1 == 0 && $tmp2 == 0){

    print "the script failed again you idiot\n";
    exit;
    
}elsif($tmp0 == $tmp3 && $tmp3 == $tmp4){

    print("SAME LENGTH\n");
    my $run = '/usr/local/devel/BCIS/Illumina_tools/splitSequences_fastq.pl Sample_' . $in . '.interleaved_pairs.fastq Q.' . $in . '.Final.R1.fastq Q.' . $in . '.Final.R2.fastq';
   print("$run\n");
    system("$run");
    my $expected = 'Q.' . $in . '.Final.R2.fastq';
    if(!-e("$expected")){ print("$expected not found\n"); exit; }
    
}else{

    print("DIFF LENGTH\n");
    my $run = 'cat tmp' . $in . '.both.list | awk \'{print \$1,\"_1\"}\' | tr -d \' \' | grep -f - ' . $in . '.all.fastq > ' . $in . '.r1.ids.txt';
 print("$run\n");
system("$run");
    my $run = 'cat tmp' . $in . '.both.list | awk \'{print \$1,\"_2\"}\' | tr -d \' \' | grep -f - ' . $in . '.all.fastq > ' . $in . '.r2.ids.txt';
 print("$run\n");
system("$run");
    my $run = 'cat tmp' . $in . '.1.list tmp' . $in . '.2.list | grep -f - ' . $in . '.all.fastq > ' . $in . '.frag.ids.txt';
 print("$run\n");
system("$run");
    my $run = '/usr/local/devel/BCIS/assembly/tools/NeatFreq/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile ' . $in . '.idfix.fastq -name ' . $in . '.r1.ids.txt -outfile ' . $in . '.Final.R1.fastq';
 print("$run\n");
system("$run");
    my $run = '/usr/local/devel/BCIS/assembly/tools/NeatFreq/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile ' . $in . '.idfix.fastq -name ' . $in . '.r2.ids.txt -outfile ' . $in . '.Final.R2.fastq';
 print("$run\n");
system("$run");
    my $run = '/usr/local/devel/BCIS/assembly/tools/NeatFreq/lib/mira_3.1.15_dev_linux-gnu_x86_64_static/scripts/fastqselect.tcl -infile ' . $in . '.idfix.fastq -name ' . $in . '.frag.ids.txt -outfile ' . $in . '.Final.fragments.fastq';
 print("$run\n");
system("$run");
    
}

# cd /usr/local/scratch/CORE/jmccorri/2018/query_raw
# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N x.68.L003 -cwd perl readqc.pl 68.L003
# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N x.68.L004 -cwd perl readqc.pl 68.L004
# qsub -P 9043 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N x.68.L005 -cwd perl readqc.pl 68.L005

