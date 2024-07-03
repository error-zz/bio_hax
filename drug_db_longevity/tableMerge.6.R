install.packages("UniProt.ws")
 
setwd ("~/Desktop/exe.DrugLongevity")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

TTD_with_ChEMBL    <- read.table("TTD_with_ChEMBL_Targets_TD_032519.txt", sep="\t", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)
t.parse <- TTD_with_ChEMBL[-1,] 
colnames(t.parse) <- t.parse[1,]
t.parse <- t.parse[-1,1:25] 

Girke_drugage      <- read.table("Girke_drugage_id_mapping.txt",          sep="\t", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)
g.parse <- Girke_drugage[-1,] 
colnames(g.parse) <- g.parse[1,]
g.parse <- g.parse[-1,] 
g.parse <- g.parse[ , colSums(!is.na(g.parse)) != 0 ]
 
m.parse            <- read.table("Merged_Drugs_022519.MergedNClean.txt",  sep="\t", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)
colnames(m.parse) <- m.parse[1,]
m.parse <- m.parse[-1,] 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# source("https://bioconductor.org/biocLite.R")
# biocLite("UniProt.ws") #, type = "source", dependencies=TRUE)
# library(UniProt.ws)
# up <- UniProt.ws(taxId=9606)
# #Current Taxonomy ID:
# #  9606
# #Current Species name:
# #  Homo sapiens
# 
# select(up, "P43166", "SYMBOL", "UNIPROT")
#   #org.Rn.eg.db, head(t.parse$UniProt_ID), "SYMBOL", "UNIPROT")
#    # head(t.parse$UniProt_ID)

# F THIS GARBAGE LIBRARY IM JUST USING THE WEB PORTAL WHAT A MASSIVE WASTE OF TIME
# https://www.uniprot.org/uploadlists/

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# I. add symbols
# https://www.uniprot.org/uploadlists/

t.symbols <- read.table("t.parse.uniprotIDs_toSymbols.txt",  sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
t.parse$Symbols <- rep(NA, dim(t.parse)[1])
for( q in 1:dim(t.parse)[1] ){
  
  syms <- list()
  si <- 1
  for(s in strsplit( t.parse$UniProt_ID[ q ] , "," )){
    if(length(s) > 0){
      if(!is.na(s)){
        if(length(t.symbols[ t.symbols$From == s , "To" ] )>0){
          syms[si] <- t.symbols[ t.symbols$From == s , "To" ] 
          si <- si + 1
        }
      }
    }
  }
  if(length(paste(unique(unlist(syms)),sep=","))>0){
    t.parse[q,"Symbols"] <- paste(unique(unlist(syms)),sep=",")
  }
  
}

# add ensembl IDs too

t.symbols <- read.table("toENSG.txt",  sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
t.parse$EnsemblID <- rep(NA, dim(t.parse)[1])
for( q in 1:dim(t.parse)[1] ){
  
  syms <- list()
  si <- 1
  for(s in strsplit( t.parse$UniProt_ID[ q ] , "," )){
    if(length(s) > 0){
      if(!is.na(s)){
        if(length(t.symbols[ t.symbols$From == s , "To" ] )>0){
          syms[si] <- t.symbols[ t.symbols$From == s , "To" ] 
          si <- si + 1
        }
      }
    }
  }
  if(length(paste(unique(unlist(syms)),sep=","))>0){
    t.parse[q,"EnsemblID"] <- paste(unique(unlist(syms)),sep=",")
  }
  
}




dim(unique(data.frame(A=t.parse$pref_name, B=t.parse$EnsemblID)))




# add the coordinates and chromosome for each ensembl id 

us <- unique(  t.parse[,"EnsemblID"] )
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
huref_genes <- getBM(attributes=
                        c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), 
                    filters =
                         'ensembl_gene_id', 
                    values =us, mart = ensembl
                   )
# ditch the transcript IDs but store them separately
tIDs <- huref_genes[,2]
names(tIDs) <- huref_genes[,1]
huref_genes <- unique(huref_genes[,-2])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# IIa. https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0%20-O%20variants.tsv.bgz

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
all.snps.df <- as.data.frame(matrix( ncol=length(snps.df.colnames), nrow=1 ))
colnames(all.snps.df) <- snps.df.colnames
a <- 1
for( h in 1:dim(huref_genes)[1]){
#for( h in 393:dim(huref_genes)[1]){
  print( paste(a, (h/dim(huref_genes)[1])*100 ) )
  query <- paste0(  "https://rest.ensembl.org/overlap/region/human/",
                    huref_genes[h,"chromosome_name"], #7
                    ":",
                    huref_genes[h,"start_position"], #140424943
                    "-",
                    huref_genes[h,"end_position"], #140624564
                    "?content-type=application/json;feature=variation")
  xmlout <- paste0(  
    huref_genes[h,"chromosome_name"], #7
    ".",
    huref_genes[h,"start_position"], #140424943
    ".",
    huref_genes[h,"end_position"], #140624564
    ".xml")
  if(!file.exists(xmlout)){
    snps <- download.file(query, xmlout)
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

#save.image( "20190505.Rdata" )
load("20190505.Rdata")
all.snps.df <- all.snps.df[ , colnames(all.snps.df)[ colSums(!is.na(all.snps.df)) != 0 ] ]

# ALL OF THE METHODS BELOW DID NOT WORK
# I HAD TO WRITE AN API QUERY ENGINE - FML

#biocLite("xml2")
#library(biomaRt)
#library(xml2)
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# variation = useEnsembl(biomart="snp", dataset="hsapiens_snp")
#uk_bioank_snps_overlapping_with_ttd_target_gene_ensembl_ids <- as.data.frame(matrix(ncol=7,nrow=1))
#c <- 1
#ui <- 1
#' for( u in us ) ){
#'   print(paste(ui,u))
#'   try( report <- getBM(attributes = c(
#'     "ensembl_gene_stable_id",
#'     'refsnp_id'), #,
#'     #'chr_name',
#'     #'chrom_start',
#'     #'chrom_end',
#'     #'minor_allele',
#'     #'minor_allele_freq'),
#'     filters = 'ensembl_gene',
#'     values = u,
#'     mart = variation
#'   ) , TRUE)
#'   print(.Last.value)
  # #result <- try(myFunc(a), silent=TRUE)
  #[1] "3 ENSG00000133019"
  #Error in curl::curl_fetch_memory(url, handle = handle) : 
  #  Timeout was reached: Operation timed out after 300004 milliseconds with 839501 bytes received
  
  # https://www.biostars.org/p/363530/
  # Since this query is too large for BioMart, I would suggest using the Ensembl REST API Overlap endpoint, restricting the query using the 'feature=variation' parameter: https://rest.ensembl.org/documentation/info/overlap_region
  
  # in other words, if it fails i have to parse the chr:coordStart-coordEnd and feed it to a web search 
  # UGHHHHH!
  # example : 
  # https://rest.ensembl.org/overlap/region/human/7:140424943-140624564?content-type=application/json;feature=variation
  # 
  # print(paste(ui,u,dim(report)[1]))
  # if(dim(report)[1] > 0){
  #   for( r in 1:dim(report)[1] ){
  #     uk_bioank_snps_overlapping_with_ttd_target_gene_ensembl_ids[ c  , ] <- report[ r , ]
  #     c <- c + 1
  #   }
  #   print(dim(uk_bioank_snps_overlapping_with_ttd_target_gene_ensembl_ids))
  # }
#  
#  ui <- ui + 1
#}
 
#look up a single gene and get SNP data


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# IIb. now merge on 

genage    <- read.table("genage.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)

t.sym       <- unique(toupper(t.parse$Symbols))
g.sym       <- toupper(genage$symbol)
g.sym_alias <- toupper(unlist(strsplit(genage$aliases, ' ')))

#length( unique( c( t.sym[ t.sym %in% g.sym ] , t.sym[ t.sym %in% g.sym_alias ] ) ) )

t.sym_gen <- t.parse
t.sym_gen[ , paste0("genage.",colnames(genage)) ] <- rep(NA, dim(t.sym_gen)[1])
for( u in unique( unlist(c( t.sym[ t.sym %in% g.sym ] , t.sym[ t.sym %in% g.sym_alias ] ))) ){
  
  if( u %in% g.sym){
    tt <- toupper(t.parse$Symbols) == u
    tt[is.na(tt)] <- FALSE
    t.sym_gen[ tt
               , 
               paste0("genage.",colnames(genage)) 
             ] <- unname(unlist(genage[ toupper(genage$symbol) == u , ]))
    print(paste("A", u))
  }else if( u %in% g.sym_alias){
    print("C")
    match <- 0
    gi <- 1
    for(g in genage$aliases){
      if( u %in% toupper(unlist(strsplit(g, ' '))) ){
        match <- gi
      }
      gi <- gi + 1
    }
    t.sym_gen[match
              , 
              paste0("genage.",colnames(genage)) 
              ] <- unname(unlist(genage[ match , ]))
    print(paste("B", u))
  }
  
}
   

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# III. Merge on pref_name
# 
# tg.parse <- t.parse
# for(c in colnames(g.parse)){
#   c <- paste0("GIRKE.",c)
#   tg.parse[,c] <- rep(NA, dim(t.parse)[1])
# }
# 
# for( q in 1:dim(t.parse)[1] ){
#   
# 
#   
# }
# 
# 
# 
# 
# 
















































# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# ***** 





L100_mrna_profiling_assay.csv    <- read.table("L100_mrna_profiling_assay.csv", sep=",", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)

disease_signatures.v    <- read_json("disease_signatures-v1.0.json")
# disease_signatures.p    <- read_json("disease_signatures-p1.0.json")
# Error in parse_con(txt, bigint_as_char) : 
#   lexical error: invalid char in json text.
# n_genes": [["5330432E05Rik", NaN], ["1700080E11Rik", NaN], [
#                     (right here) ------^

single_drug_perturbations.v    <- read_json("single_drug_perturbations-v1.0.json")
single_drug_perturbations.dm    <- read_json("single_drug_perturbations-DM.json")
single_drug_perturbations.p    <- read_json("single_drug_perturbations-p1.0.json")

single_gene_perturbations.v    <- read_json("single_gene_perturbations-v1.0.json")
single_gene_perturbations.p    <- read_json("single_gene_perturbations-p1.0.json")
 
DSigDB_All.gmt

DSigDB_All_detailed.txt






























# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# SKIP TO 









