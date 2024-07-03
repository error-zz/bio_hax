setwd ("~/Desktop/exe.DrugLongevity")

####################################################################################################################################################################################################

####################################################################################################################################################################################################

####################################################################################################################################################################################################

####################################################################################################################################################################################################

##########################
# load all ttd genes

TTD_with_ChEMBL    <- read.table("TTD_with_ChEMBL_Targets_TD_032519.txt", sep="\t", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)
t.parse <- TTD_with_ChEMBL[-1,] 
colnames(t.parse) <- t.parse[1,]
t.parse <- t.parse[-1,1:25] 

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

##########################
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
# 4862    2


write.table(t.parse, file="t.parse.tsv", sep="\t")
 
##########################
# add Thomas' conversions

# Thomas took the TTD gene list and brought in some columns from other tables.  These are necessary given other conversions.
# n1 <- read.table( "Nik_Gene_Protein_Mappings_12May19.Nik_Gene_Protein_Mappings_12May19.tsv", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
# n2 <- read.table( "Nik_GeneIDs_DrugTargetBioact_12May19.Nik_GeneIDs_DrugTargetBioact_12May19.tsv", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
# n3 <- read.table( "Nik_GeneIDs_DrugTargetAnnotation_12May19.Nik_GeneIDs_DrugTargetAnnotation_12May19.tsv", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)






####################################################################################################################################################################################################

####################################################################################################################################################################################################

####################################################################################################################################################################################################

####################################################################################################################################################################################################



Girke_drugage      <- read.table("Girke_drugage_id_mapping.txt",          sep="\t", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)
g.parse <- Girke_drugage[-1,] 
colnames(g.parse) <- g.parse[1,]
g.parse <- g.parse[-1,] 
g.parse <- g.parse[ , colSums(!is.na(g.parse)) != 0 ]


####################################################################################################################################################################################################

####################################################################################################################################################################################################

####################################################################################################################################################################################################

####################################################################################################################################################################################################



m.parse            <- read.table("Merged_Drugs_022519.MergedNClean.txt",  sep="\t", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)
colnames(m.parse) <- m.parse[1,]
m.parse <- m.parse[-1,] 


####################################################################################################################################################################################################

####################################################################################################################################################################################################

####################################################################################################################################################################################################

####################################################################################################################################################################################################




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
