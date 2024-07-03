#!/usr/bin/perl -w
# J. Craig Venter Institute
# Written by Jamison M. McCorrison (jmccorri@jcvi.org)
no warnings;
use strict;

my $queryprefix = shift;
my $refprefix = shift; 

# TRIM
my $runcmd = "/usr/local/bin/java -jar /usr/local/devel/BCIS/pratap/bin/trimmomatic_install/Trimmomatic-0.35/trimmomatic-0.35.jar PE -threads 4 -phred33 /usr/local/scratch/CORE/jmccorri/2018/query_raw/" . $queryprefix . ".R1.fastq /usr/local/scratch/CORE/jmccorri/2018/query_raw/" . $queryprefix . ".R2.fastq /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.R1.fastq /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.S1.fastq /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.R2.fastq /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.S2.fastq ILLUMINACLIP:/usr/local/scratch/CORE/pratap/hlorenzi/eag-164/miRNA/human/sequences/TruSeq_Small_RNA_Sample_Prep_Kits.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:12 MINLEN:25";
print "\n\n$queryprefix     |\n\t$runcmd";
system("$runcmd");

my $runcmd = "/usr/local/bin/java -jar /usr/local/devel/BCIS/pratap/bin/trimmomatic_install/Trimmomatic-0.35/trimmomatic-0.35.jar SE -threads 4 -phred33 /usr/local/scratch/CORE/jmccorri/2018/query_raw/" . $queryprefix . ".Fragments.fastq /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.Frag.fastq ILLUMINACLIP:/usr/local/scratch/CORE/pratap/hlorenzi/eag-164/miRNA/human/sequences/TruSeq_Small_RNA_Sample_Prep_Kits.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:12 MINLEN:25";
print "\n\n$queryprefix     |\n\t$runcmd";
system("$runcmd");
  
my $runcmd = "cat /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.S*q /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.Frag.fastq > /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.fragments.fastq";
print "\n\n$queryprefix     |\n\t$runcmd";
system("$runcmd");
 
my $runcmd = "/usr/local/devel/BCIS/pratap/bin/hisat2 --threads 4 --time --dta --dta-cufflinks -x /usr/local/scratch/CORE/jmccorri/2018/ref_cds/" . $refprefix . " -1 /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.R1.fastq -2 /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.R2.fastq -U /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/trimmed/" . $queryprefix . "_trimmed.fragments.fastq -S " . $queryprefix . ".sam";
print "\n\n$queryprefix     |\n\t$runcmd";
system("$runcmd");
 
my $runcmd = "/usr/local/devel/BCIS/pratap/bin/samtools view -Sbo " . $queryprefix . ".bam " . $queryprefix . ".sam";
print "$\n\nqueryprefix     |\n\t$runcmd";
system("$runcmd");
 
my $runcmd = "/usr/local/devel/BCIS/pratap/bin/samtools sort " . $queryprefix . ".bam -o " . $queryprefix . "_simplSort.bam";
print "\n\n$queryprefix     |\n\t$runcmd";
system("$runcmd");
 # samtools sort 01_Common_Grackle.bam -o 01_Common_Grackle_simplSort.bam
 
#my $runcmd = "/usr/local/devel/BCIS/pratap/bin/samtools sort " . $queryprefix . "_trimmed.sam -n -T " . $queryprefix . ".sortedByName -o " . $queryprefix . ".sortedByName.bam";
#print "$queryprefix     |\n\t$runcmd";
#system("$runcmd");
 
my $runcmd = "/usr/local/devel/BCIS/pratap/bin/stringtie " . $queryprefix . "_simplSort.bam -o ./" . $queryprefix . ".gtf -p 4 -A " . $queryprefix . ".gtf.gene_abundance.txt";
print "\n\n$queryprefix     |\n\t$runcmd";
system("$runcmd");
 
my $runcmd = "samtools view " . $queryprefix . ".sortedByName.bam | /usr/local/devel/BCIS/pratap/bin/htseq-count -s no - ./" . $queryprefix . ".gtf >|" . $queryprefix . ".bam.counts.txt";
print "\n\n$queryprefix     |\n\t$runcmd";
system("$runcmd");
 
 # cd /usr/local/scratch/CORE/jmccorri/2018/AlignAndCounts.StringTie/61_MuteSwan
 # qsub -P 700020 -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N hs61 -cwd perl run.pl Mute_Swan Anas_platyrhynchos
 
 
 
 

