#!/usr/bin/perl -w
use strict;
no warnings;

# EXAMPLE EXE:
# perl x.pl PREFIX /usr/local/scratch/CORE/jmccorri/BCIS481_16sOTU/dataset1_mothurexample/F3D0_S188_L001_R1_001.fastq /usr/local/scratch/CORE/jmccorri/BCIS481_16sOTU/dataset1_mothurexample/F3D0_S188_L001_R2_001.fastq NULL NULL /usr/local/scratch/CORE/jmccorri/BCIS481_16sOTU/TESTEXE_BASIC/z01_fuzznuc/utaxref/rdp_16s_trainset15/fasta/refdb.fa /usr/local/scratch/CORE/jmccorri/BCIS481_16sOTU/TESTEXE_BASIC/z01_fuzznuc/utaxref/rdp_16s_trainset15/taxconfs/full_length.tc /usr/local/scratch/CORE/jmccorri/BCIS481_16sOTU/TESTEXE_BASIC/z01_fuzznuc/silva.bacteria/silva.bacteria.fasta /usr/local/scratch/CORE/jmccorri/BCIS481_16sOTU/TESTEXE_BASIC/z01_fuzznuc/silva.bacteria/silva.bacteria.silva.tax

my $prefix = shift;
my $in1 = shift;
my $in2 = shift;
# Allow user to specify percent of samples (OTU_per_sample_min_cov) 
my $OTU_per_sample_min_cov = shift;
# with minimum coverage (OTU_min_cov) required for pipeline to continue. 
my $OTU_min_cov = shift;
# If percent_of_sample and OTU_min_coverage are NA then option will not be applied.
my $sixteenSreferenceDB = shift;
my $taxconfs = shift;
my $bacterial_reference = shift;
my $bacterial_tax = shift; 

print "\n\n-----------------------------------------------------\n-----------------------------------------------------\n\n";
print "\n\n██████╗ ███████╗   ██████╗     ██████╗ ████████╗██╗   ██╗\n";
print "██╔══██╗██╔════╝   ╚════██╗   ██╔═══██╗╚══██╔══╝██║   ██║\n";
print "██████╔╝█████╗      █████╔╝   ██║   ██║   ██║   ██║   ██║\n";
print "██╔═══╝ ██╔══╝     ██╔═══╝    ██║   ██║   ██║   ██║   ██║\n";
print "██║     ███████╗██╗███████╗██╗╚██████╔╝   ██║   ╚██████╔╝\n";
print "╚═╝     ╚══════╝╚═╝╚══════╝╚═╝ ╚═════╝    ╚═╝    ╚═════╝\n ";
print "\n\n-----------------------------------------------------\n-----------------------------------------------------\n\n";
print " ██╗\n███║\n╚██║\n ██║\n ██║\n ╚═╝\n"; # ONE # 
print "Step 1: Entering Uparse pipeline....\n";
print "Assigning reads to operational taxonomic units.\n";

print "\n\n-----------------------------------------------------\n\n";
print "1. (1/5) Merging overlapping paired end sequences...\n\n";
runsys("/usr/local/bin/usearch -fastq_mergepairs $in1 -reverse $in2 -fastqout merged.fq -fastq_maxdiffs 6 -fastq_truncqual 2 -fastq_minlen 36 -fastq_minmergelen 50 -fastq_allowmergestagger");
if (!(-s("merged.fq"))){ print "\n\n!!! FAIL CASE !!!\nMerge failed to generate merged.fq.  See command above."; exit; }

print "\n\n-----------------------------------------------------\n\n";
print "1. (2/5) Filtering...\n\n";
runsys("/usr/local/bin/usearch -fastq_filter merged.fq -fastq_maxee 1.0 -relabel FILT_ -fastaout filtered.fa");
if (!(-s("filtered.fa"))){ print "\n\n!!! FAIL CASE !!!\nFilter failed to generate filtered.fa.  See command above."; exit; }

print "\n\n-----------------------------------------------------\n\n";
print "1. (3/5) Determining unique sequences (OTUs)...\n\n";
runsys("/usr/local/bin/usearch -derep_fulllength filtered.fa -relabel UNIQ_ -sizeout -fastaout uniques.fa");
if (!(-s("uniques.fa"))){ print "\n\n!!! FAIL CASE !!!\nFilter failed to generate uniques.fa.  See command above."; exit; }

# print "\n\n-----------------------------------------------------\n\n";
# print "1. (4/6) Sorting unique OTUs...\n\n";
# runsys("/usr/local/bin/usearch -sortbysize otus.fa -fastaout sorted.fa -minsize 2");
# if (!(-s("sorted.fa"))){ print "\n\n!!! FAIL CASE !!!\nFilter failed to generate sorted.fa.  See command above."; exit; }

print "\n\n-----------------------------------------------------\n\n";
print "1. (4/5) Clustring OTUs...\n\n";
runsys("/usr/local/bin/usearch -cluster_otus uniques.fa -otus otus.fa -relabel OTU_ -sizeout");
if (!(-s("otus.fa"))){ print "\n\n!!! FAIL CASE !!!\nFilter failed to generate otus.fa.  See command above."; exit; }

print "\n\n-----------------------------------------------------\n\n";
print "1. (5/5) Printing high quality matches between non-unique reads (merged.fa) and otus (otus.fa) to generate OTU table...\n\n";
#runsys("/usr/local/bin/usearch -usearch_global merged.fq -db otus.fa -strand plus -id 0.97 -otutabout otutab.txt -biomout otutab.json");
runsys("/usr/local/bin/usearch -usearch_global merged.fq -db otus.fa -strand plus -id 0.97 -alnout otutab.txt -blast6out hits.b6 -maxaccepts 8 -maxrejects 256");
if (!(-s("otutab.txt"))){ print "\n\n!!! FAIL CASE !!!\nFilter failed to generate otutab.txt.  See command above."; exit; }
# ^ NOTE : Possible update to assign taxos early
# If the OTU sequences have tax=xxx; annotations, these will be included as an extra column in the tabbed file or as taxonomy metadata in the BIOM file. These annotations can be generated using the -fastaout option of the utax command, e.g.:

# ^ NOTE : 
# [ http://www.drive5.com/usearch/manual/mapreadstootus.html  ]
# [ http://www.drive5.com/usearch/manual/uparse_cmds.html     ]
# [ http://www.drive5.com/usearch/manual/uparse_pipeline.html ]

print "\nStep 1 completed successfully (otutab.txt and hits.b6 generated successfully.)";

print "\n\n-----------------------------------------------------\n-----------------------------------------------------\n\n";

print "██████╗ \n╚════██╗\n █████╔╝\n██╔═══╝ \n███████╗\n╚══════╝\n"; # TWO # 
#print "Step 2: Entering mothur pipeline....\n";
print "Step 2. Assigning taxonomy using mothur\n\n";

# ^ NOTE : "Most taxonomy prediction algorithms don't provide a confidence estimate, including 
# GAST, 
# the default QIIME method (assign_taxonomy.py ‑m uclust) 
# and the mothur Classify_seqs command with method=knn."
# 					[ http://www.drive5.com/usearch/manual/tax_conf.html ]
# mothur > classify.seqs(fasta=abrecovery.fasta, template=nogap.bacteria.fasta, taxonomy=silva.bacteria.silva.tax, method=knn)

print "\n\n-----------------------------------------------------\n\n";
my $template = "/usr/local/scratch/CORE/jmccorri/BCIS481_16sOTU/TESTEXE_BASIC/z01_fuzznuc/silva.bacteria/silva.bacteria.ng.fasta";
print "2. (1/2) De-gapping Bacterial Reference\n\n";
if (!(-s("/usr/local/scratch/CORE/jmccorri/BCIS481_16sOTU/TESTEXE_BASIC/z01_fuzznuc/silva.bacteria/silva.bacteria.ng.fasta"))){
	runsys("/usr/local/bin/mothur <<< \"degap.seqs(fasta=$bacterial_reference)\"");
	if (!(-s("/usr/local/scratch/CORE/jmccorri/BCIS481_16sOTU/TESTEXE_BASIC/z01_fuzznuc/silva.bacteria/silva.bacteria.ng.fasta"))){  print "\n\n!!! FAIL CASE !!!\nFailed to generate non-gapped bacterial reference (silva.bacteria.ng.fasta)\n"; exit; }
	else{ $template = "/usr/local/scratch/CORE/jmccorri/BCIS481_16sOTU/TESTEXE_BASIC/z01_fuzznuc/silva.bacteria/silva.bacteria.ng.fasta"; }
}else{
	print "\tSkipping step.  (File exists.)\n\n"
}

print "\n\n-----------------------------------------------------\n\n";
print "2. (2/2) Classifying OTU Sequences with Bacterial Reference\n";
runsys("/usr/local/bin/mothur <<< \"classify.seqs(fasta=otus.fa, template=$template, taxonomy=$bacterial_tax, method=knn))\"");
# otus.silva.knn.taxonomy
# otus.silva.knn.tax.summary
# otus.silva.knn.flip.accnos
if (!(-s("otus.silva.knn.taxonomy"))){ print "\n\n!!! FAIL CASE !!!\nFailed to classify OTU sequences (and generate otus.silva.knn.taxonomy)\n\n"; exit; }
# ^ NOTE : [  wget http://www.mothur.org/w/images/9/98/Silva.bacteria.zip ]
# alternate run methods?
# mothur > classify.seqs(fasta=stability.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.files.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, reference=trainset14_032015.pds.fasta, taxonomy=trainset14_032015.pds.tax, cutoff=80)

print "\nStep 2 completed successfully (otus.silva.knn.taxonomy, otus.silva.knn.tax.summary, and otus.silva.knn.flip.accnos generated successfully.)";

print "\n\n-----------------------------------------------------\n-----------------------------------------------------\n\n";


# ALTERNATE METHOD FOR 2. 
# ASSIGN TAXONOMY WITH UTAX (option makeudb_utax does not work)

# 1. Download the 16S reference files and extract them from the archive using tar -zxvf. The archive contains a FASTA database with full-length sequences plus pre-trained taxconfs files for shorter segmenets (500nt, 250nt and 120nt).
# print "Using 16s reference database : $sixteenSreferenceDB\n";
# ^ # [ http://drive5.com/utax/data/utax_rdp_16s_tainset15.tar.gz ]
# 2. Choose the taxconfs file for length closest to your sequences.
# print "Using taxconfs config file : $taxconfs\n\t";
# 3. Make the UDB file using the makeudb_utax command (it is ok to use full-length reference sequences with parameters for shorter sequences):
# runsys("usearch -makeudb_utax $sixteenSreferenceDB -output refdb.udb -taxconfsin $taxconfs");
# 4. Run the utax command or cluster_otus_utax command:
# print "2. (1/2) Printing high quality matches between non-unique reads (merged.fa) and otus (otus.fa) to generate OTU table...\n\n";
# runsys("usearch -utax reads.fastq -db tax.udb -utaxout utax.txt -strand both");

# print "\n\n-----------------------------------------------------\n\n";
# print "\n4. (1/2) Taxonomy assignment to OTUs.\n\n";
# runsys("/usr/local/bin/usearch -utax otus.fa -db rdp_16s -strand both -taxconfs rdp_16s_short.tc -utaxout tax.txt");
# ^ NOTE : It is not possible to generate your own .tc files...

# print "\n\n-----------------------------------------------------\n\n";
# print "\n4. (2/2) Generating OTU table.\n\n";
# runsys("python ~/py/uc2otutab.py map.uc > otu_table.txt");

# ^ NOTE : 
# [ http://www.drive5.com/usearch/manual/utax_16s.html       ]


# print "\n\n-----------------------------------------------------\n-----------------------------------------------------\n\n";

# print "██████╗ \n╚════██╗\n █████╔╝\n ╚═══██╗\n██████╔╝\n╚═════╝\n"; # THREE # 
# print "Step 3: Filtering by match count\n";

# print "\n\n\n\t ! TBA !\n\n\n\n";

# print "\n\n-----------------------------------------------------\n-----------------------------------------------------\n\n";
 
# print "██╗  ██╗\n██║  ██║\n███████║\n╚════██║\n     ██║\n     ╚═╝\n"; # FOUR # 

# print "\n\n-----------------------------------------------------\n-----------------------------------------------------\n\n";

print "███████╗\n██╔════╝\n███████╗\n╚════██║\n███████║\n╚══════╝\n"; # FIVE
print "Step 5. Visualization\n";
print "\n\n!! TBA !!\n\nExiting.";

exit;

sub runsys {
				my $cmd = shift;
				print "!|| SYS_EXE ||\t$cmd\n";
				system("$cmd");
}
