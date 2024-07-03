setwd ("~/Desktop/exe.DrugLongevity")
#library(UniProt.ws)
#library(biomaRt)

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


#dgidb_export_2019-05-11.tsv
dgidb_export_2019    <- read.table("dgidb_export_2019-05-11.tsv", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)

tmp <- read.table("dgidb_export_2019-05-11.tsv", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
dim(unique(data.frame(A=tmp$gene, B=tmp$drug)))


head(dgidb_export_2019)
 
dgidb_export_2019$ensembl_gene_id <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$hgnc_symbol     <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$chromosome_name <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$start_position  <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$end_position    <- rep(NA, dim(dgidb_export_2019)[1])

#dgidb_export_2019$

for(d in 1:dim(dgidb_export_2019)[1]){
  if( sum( huref_genes$hgnc_symbol == dgidb_export_2019[d,]$gene ) > 0 ){
    dgidb_export_2019[ d , colnames(huref_genes) ] <- huref_genes[ huref_genes$hgnc_symbol == dgidb_export_2019[d,]$gene , ]
  }
}

save.image("q20190514.Rdata")
 
#
dgidb_export_2019$druggability_tier <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$hgnc_names     <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$chr_b37 <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$start_b37  <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$end_b37    <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$strand    <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$description    <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$no_of_gwas_regions    <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$small_mol_druggable    <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$bio_druggable    <- rep(NA, dim(dgidb_export_2019)[1])
dgidb_export_2019$adme_gene    <- rep(NA, dim(dgidb_export_2019)[1])


for(d in 1:dim(dgidb_export_2019)[1]){
  if( sum(Druggable_Genome_GWAS_2017$ensembl_gene_id == dgidb_export_2019$ensembl_gene_id[1]) > 0){
    dgidb_export_2019[ d , colnames(Druggable_Genome_GWAS_2017)[2:12] ] <- Druggable_Genome_GWAS_2017[ Druggable_Genome_GWAS_2017$ensembl_gene_id == Druggable_Genome_GWAS_2017[d,]$ensembl_gene_id , 2:12 ]
  }
}


#save.image("q20190515.Rdata")


setwd ("~/Desktop/exe.DrugLongevity")
#biocLite("BiocUpgrade")
 #remove.packages("BiocInstaller")
#install.packages("BiocInstaller", repos="http://bioconductor.org/packages/2.13/bioc")
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
load("q20190515.Rdata")
 
sum_all_GTex_eQTL    <- read.table("GTEx_Analysis_v7_eQTL/sum_all_GTex_eQTL.tsv", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
colnames(sum_all_GTex_eQTL)[1]
ensembl_split <- strsplit(sum_all_GTex_eQTL$gene_id, '\\.')
ensembl_gene_simple  <- list()
ensembl_gene_isoform <- list()
ei <- 1
for(e in ensembl_split){
  ensembl_gene_simple[ei] <- ensembl_split[[ei]][1]
  ensembl_gene_isoform[ei] <- ensembl_split[[ei]][2]
  ei <- ei + 1
}
sum_all_GTex_eQTL$ensembl_gene_simple  <- unlist(ensembl_gene_simple)
sum_all_GTex_eQTL$ensembl_gene_isoform <- unlist(ensembl_gene_isoform)


dim(unique(data.frame(A=sum_all_GTex_eQTL$gtex_tissue, B=sum_all_GTex_eQTL$rs_id_dbSNP147_GRCh37p13)))


# now compare on gene simple

# FIX DUPE
colnames(dgidb_export_2019)[ colnames(dgidb_export_2019) == "strand" ] <- "dgidb.strand"
colnames(sum_all_GTex_eQTL)[1] <- "gtex_tissue"

# #Adipose_Subcutaneous
# dgidb_export_2019$gtex_tissue <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$gene_id <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$gene_name <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$gene_chr <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$gene_start <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$gene_end <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$strand <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$num_var <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$beta_shape1 <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$beta_shape2 <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$true_df <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$pval_true_df <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$variant_id <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$tss_distance <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$chr <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$pos <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$ref <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$alt <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$num_alt_per_site <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$rs_id_dbSNP147_GRCh37p13 <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$minor_allele_samples <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$minor_allele_count <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$maf <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$ref_factor <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$pval_nominal <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$slope <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$slope_se <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$pval_perm <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$pval_beta <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$qval <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$pval_nominal_threshold <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$log2_aFC <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$log2_aFC_lower <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$log2_aFC_upper <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$ensembl_gene_simple <- rep(NA, dim(dgidb_export_2019)[1])
# dgidb_export_2019$ensembl_gene_isoform <- rep(NA, dim(dgidb_export_2019)[1])
# sum_all_GTex_eQTL$search_term <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$match_term <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$match_type <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$gene <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$drug <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$interaction_types <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$sources <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$pmids <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$ensembl_gene_id <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$hgnc_symbol <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$chromosome_name <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$start_position <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$end_position <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$druggability_tier <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$hgnc_names <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$chr_b37 <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$start_b37 <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$end_b37 <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$dgidb.strand <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$description <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$no_of_gwas_regions <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$small_mol_druggable <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$bio_druggable <- rep(NA, dim(sum_all_GTex_eQTL)[1])
# sum_all_GTex_eQTL$adme_gene <- rep(NA, dim(sum_all_GTex_eQTL)[1])

p <- unique(dgidb_export_2019[ dgidb_export_2019$ensembl_gene_id == dgidb_export_2019$ensembl_gene_id[1] , ])
q <- sum_all_GTex_eQTL[ sum_all_GTex_eQTL$ensembl_gene_simple == dgidb_export_2019$ensembl_gene_id[1] ,  ]
summary <- as.data.frame(matrix( ncol=60, nrow=1 ))
colnames(summary) <- c(colnames(q),colnames(p))

#for(d in 1:dim(dgidb_export_2019)[1]){
for(d in 3225:dim(dgidb_export_2019)[1]){
  
  if(d %%25 == 0){ 
    print(paste((d/dim(dgidb_export_2019)[1])*100,"%",dim(summary)[1])) 
    write.table( summary , file = paste0("summary.x",d,".tsv"), sep = "\t" )
    summary <- as.data.frame(matrix( ncol=60, nrow=1 ))
    colnames(summary) <- c(colnames(q),colnames(p))
  }
  if( dgidb_export_2019$ensembl_gene_id[d] %in% sum_all_GTex_eQTL$ensembl_gene_simple ){
    
    #dgidb_export_2019$ensembl_gene_id[d]
    
    p <- unique(dgidb_export_2019[ dgidb_export_2019$ensembl_gene_id == dgidb_export_2019$ensembl_gene_id[d] , ])
    p <- p[ rowSums(!is.na(p)) != 0 , ]
    
    q <- sum_all_GTex_eQTL[ sum_all_GTex_eQTL$ensembl_gene_simple == dgidb_export_2019$ensembl_gene_id[d] ,  ]
    
    nu <- matrix(nrow = dim(p)[1]*dim(q)[1] , ncol = dim(p)[2]+dim(q)[2] )
    colnames(nu) <- c(colnames(q),colnames(p))
    c <- 1
    for( x in 1:dim(q)[1]){
      for( r in 1:dim(p)[1]){
        nu[c, 1:dim(q)[2] ]                         <- unname(unlist( q[x,1:dim(q)[2]] ))
        nu[c, (dim(q)[2]+1):(dim(q)[2]+dim(p)[2]) ] <- unname(unlist( p[r,1:dim(p)[2]] ))
        c <- c + 1
      }
    }
    
    summary <- rbind(summary, nu)
    
    #dgidb_export_2019[d,"gtex_tissue"] <- sum_all_GTex_eQTL[ sum_all_GTex_eQTL$ensembl_gene_simple == dgidb_export_2019$ensembl_gene_id[d] , ]
  }
}
 
# save.image("www565_incomplete.Rdata")
save.image("www_complete.Rdata")
# load("www565_incomplete.Rdata")

#write.table( summary , file = "summary.3225.tsv", sep = "\t" )
