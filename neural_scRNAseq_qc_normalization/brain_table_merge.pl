#!/usr/bin/perl -w
use strict;
no warnings;

#
# cd /usr/local/scratch/CORE/jmccorri/BRAINS_merge_tables
#
# ln -s /usr/local/projdata/0672/projects/Illumina_Data/RUN_1/analysis_data/reports_and_counts/AIBS_RSEM-countMatrices_AND_RSEM-TC-QC-reports_03312016/AIBS_run1_FC1_tpm_matrix.txt /usr/local/projdata/0672/projects/Illumina_Data/RUN_1/analysis_data/reports_and_counts/AIBS_RSEM-countMatrices_AND_RSEM-TC-QC-reports_03312016/AIBS_run1_FC2_tpm_matrix.txt /usr/local/projdata/0672/projects/Illumina_Data/RUN_2/analysis_data/reports_and_counts/AIBS_RSEM-countMatrices_AND_RSEM-TC-QC-reports_05262016/run2_FC1_tpm_matrix.txt /usr/local/projdata/0672/projects/Illumina_Data/RUN_2/analysis_data/reports_and_counts/AIBS_RSEM-countMatrices_AND_RSEM-TC-QC-reports_05262016/run2_FC2_tpm_matrix.txt /usr/local/projdata/0672/projects/Illumina_Data/RUN_3/analysis_data/reports_and_counts/AIBS_RSEM-countMatrices_AND_RSEM-TC-QC-reports_06262016/AIBS_run3_FC1_tpm_matrix.txt /usr/local/projdata/0672/projects/Illumina_Data/RUN_3/analysis_data/reports_and_counts/AIBS_RSEM-countMatrices_AND_RSEM-TC-QC-reports_06262016/AIBS_run3_FC2_tpm_matrix.txt .
# 
# perl brain_table_merge covariate_table.txt 
# 
# perl brain_table_merge.pl covariate_table.txt AIBS_run1_FC1_tpm_matrix.txt,run2_FC1_tpm_matrix.txt,AIBS_run3_FC1_tpm_matrix.txt,AIBS_run1_FC2_tpm_matrix.txt,run2_FC2_tpm_matrix.txt,AIBS_run3_FC2_tpm_matrix.txt  
# 

my $covariate_table = shift;
my $counts_tables = shift;

print "\nReading in covariate table...\n";
open INFILE, "< $covariate_table";
my $c = 0;
my %covariate_matrix;
my @covariate_header = ();
while (defined( my $line = <INFILE> )){
	my @xs = split(/\s+/, $line);
	if ($c == 0){
		# header
		@covariate_header = @xs;
	}else{
		$covariate_matrix{ $xs[0] }{ $covariate_header[$c] } = $_;
	}
	$c++;

}
close INFLIE;

print "\nReading in counts tables...\n";
my %counts_matrix;
my @ct_tbls = split(/,/, $counts_tables);
my $c = 0;
my @counts_header = ();
my @all_genes_list = ();
foreach my $ct_tbl (@ct_tbls){
	open INFILE, "< $ct_tbl";
	my $i = 0;
	while (defined( my $line = <INFILE> )){
		my @xs = split(/\s+/, $line);
		if ($i == 0){
			# header
			@counts_header = @xs;
		}else{
			push(@all_genes_list, $xs[0]);
			my $r = 0;
			foreach my $count (@xs){
				if ($r > 0){
					$counts_matrix{ $counts_header[$i] }{ $xs[0] } = $count;
				}
				$r++;
			}
		}
		$i++;
	}
	close INFILE;
}


print "\nReducing to unique gene array...\n";
# reduce gene list to unique array
my @nonuniq = @all_genes_list;
@all_genes_list = ();
my $prev = "qqq";
foreach(sort {$a <=> $b} @nonuniq){
	if ($_ ne $prev){ push (@all_genes_list, $_); print "| \t$_\t|\n"; }
	$prev = $_;
}


# initialize
print "\nInitializing report matrix...\n";
my %report_matrix = ();
foreach my $sample_id (sort {$a <=> $b} keys(%covariate_matrix)){
	print "\t| $sample_id\t|\n";
	foreach my $col_id (@covariate_header){
		$report_matrix{$sample_id}{$col_id} = "null";
	}
	foreach my $col_id (@all_genes_list){
		$report_matrix{$sample_id}{$col_id} = 0;
	}
}

# populate
print "\nPopulate report matrix...\n";
foreach my $sample_id (sort {$a <=> $b} keys( %covariate_matrix) ){
	print "\t| $sample_id\t|\n";
	foreach my $col_id (keys( %{$covariate_matrix{$sample_id}} )){
		$report_matrix{$sample_id}{$col_id} = $covariate_matrix{$sample_id}{$col_id};
	}
	foreach my $col_id (keys( %{$counts_matrix{$sample_id}} )){
		$report_matrix{$sample_id}{$col_id} = $counts_matrix{$sample_id}{$col_id};
	}
}

# print
print "\nPrinting report matrix...\n";
if (-e("brain_table_merge_output.txt")){ system("rm -rf brain_table_merge_output.txt"); }
system("touch brain_table_merge_output.txt");
open PLOTME, ">> brain_table_merge_output.txt";
print PLOTME "SampleID\t";
foreach my $col_id (@covariate_header){ print PLOTME "$col_id\t"; }
foreach my $col_id (@all_genes_list){ print PLOTME "$col_id\t"; }
print PLOTME "\n";

foreach my $sample_id (sort {$a <=> $b} keys(%report_matrix)){
	foreach my $col_id (@covariate_header){
		print PLOTME "$report_matrix{$sample_id}{$col_id}\t";
	}
	foreach my $col_id (@all_genes_list){
		print PLOTME "$report_matrix{$sample_id}{$col_id}\t";
	}
	print PLOTME "\n";
}
close PLOTME;

