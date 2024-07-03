setwd ("~/Desktop/exe.DrugLongevity")
#install.packages('jsonlite', dependencies=TRUE, repos='http://cran.rstudio.com/')
library(read_json)

# # # # # # # # # # # # # # # 

# >> We will work with the following excel file, focusing on the significant variants. 
#  "Druggable_Genome_GWAS_2017.xlsx":
#  - 4479 gene symbols
# - Includes chromosome and coords but don't use!  use the ensembl gene id based coordinates instead.

Druggable_Genome_GWAS_2017    <- read.table("Druggable_Genome_GWAS_2017.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)


us <- unique(  Druggable_Genome_GWAS_2017$ensembl_gene_id )
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
huref_genes <- getBM(attributes=
                       c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), 
                     filters =
                       'ensembl_gene_id', 
                     values =us, mart = ensembl
)
tIDs <- huref_genes[,2]
names(tIDs) <- huref_genes[,1]
huref_genes <- unique(huref_genes[,-2])




snps.df.colnames <-c( "ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position",
                      #allcolnames)
                      "id","consequence_type","strand","seq_region_name","assembly_name","feature_type","source","end","start",
                      "alleles",
                      "alleles1","alleles2","alleles3","alleles4","alleles5","alleles6","alleles7","alleles8","alleles9","alleles10",
                      "alleles11","alleles12","alleles13","alleles14","alleles15","alleles16","alleles17","alleles18","alleles19","alleles20",
                      "alleles21","alleles22","alleles23","alleles24","alleles25","alleles26","alleles27","alleles28","alleles29","alleles30",
                      "alleles31","alleles32","alleles33","alleles34","alleles35","alleles36","alleles37","alleles38","alleles39","alleles40",
                      "alleles41","alleles42","alleles43","alleles44","alleles45","alleles46","alleles47","alleles48","alleles49","alleles50",
                      "clinical_significance",
                      "clinical_significance1","clinical_significance2","clinical_significance3","clinical_significance4","clinical_significance5","clinical_significance6","clinical_significance7","clinical_significance8","clinical_significance9","clinical_significance10",
                      "clinical_significance11","clinical_significance12","clinical_significance13","clinical_significance14","clinical_significance15","clinical_significance16","clinical_significance17","clinical_significance18","clinical_significance19","clinical_significance20",
                      "clinical_significance21","clinical_significance22","clinical_significance23","clinical_significance24","clinical_significance25","clinical_significance26","clinical_significance27","clinical_significance28","clinical_significance29","clinical_significance30"
) 


#save.image("xxx3396_incomplete.Rdata")
load("xxx3396_incomplete.Rdata")


#all.snps.df <- as.data.frame(matrix( ncol=length(snps.df.colnames), nrow=1 ))
#colnames(all.snps.df) <- snps.df.colnames
#a <- 1
#for( h in 1:dim(huref_genes)[1]){
  #for( h in 393:dim(huref_genes)[1]){
#for( h in 2997:dim(huref_genes)[1]){
#for( h in 3396:dim(huref_genes)[1]){
for( h in 3982:dim(huref_genes)[1]){
  print( paste(a, (h/dim(huref_genes)[1])*100 ) )
  x.start <- huref_genes[h,"start_position"]-2500
  x.end   <- huref_genes[h,"start_position"]+2500
  query <- paste0(  "https://rest.ensembl.org/overlap/region/human/",
                    huref_genes[h,"chromosome_name"], #7
                    ":",
                    x.start, #140424943
                    "-",
                    x.end, #140624564
                    "?content-type=application/json;feature=variation")
  xmlout <- paste0(  
    huref_genes[h,"chromosome_name"], #7
    ".",
    huref_genes[h,"start_position"], #140424943
    ".",
    huref_genes[h,"end_position"], #140624564
    ".xml")
  zzz <- 1
  while(!file.exists(xmlout)){
    if(zzz > 1){ print("restart!") }
    snps <- download.file(query, xmlout)
    zzz <- zzz + 1
  }
  xmlin <- read_json(xmlout)
  
  snps.df <- matrix( ncol=length(snps.df.colnames), nrow=length(xmlin) )
  colnames(snps.df) <- snps.df.colnames
  for(x in 1:length(xmlin)){
    if(x %% 5000 == 0){
      print( paste0( huref_genes[h,"hgnc_symbol"] , " ",  length(xmlin), " ---[read]---> ", (x/length(xmlin))*100, "%"  ) )
    }
    snps.df[ x , 1:5 ] <-  unname(unlist(huref_genes[h,]))
    for(n in names(unlist(xmlin[x]))){
      snps.df[ x , n ] <- unlist(xmlin[x])[n]
    }
  }
  all.snps.df <- rbind(all.snps.df, snps.df)
  a <- a + 1
}
 
#save.image("xxx3396_incomplete.Rdata")
#save.image("xxx_complete.Rdata")
 
setwd ("~/Desktop/exe.DrugLongevity")
load("xxx_complete.Rdata")
dim(unique(data.frame(A=all.snps.df$id, B=all.snps.df$ensembl_gene_id)))
write.table( all.snps.df, "all.snps.dataframe.20190516.tsv", sep="\t" )
 
# all.snps_uniq.df <- unique(all.snps.df)
 