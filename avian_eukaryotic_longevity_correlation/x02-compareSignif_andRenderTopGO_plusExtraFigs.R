setwd("C:\\Users\\jamis\\Desktop\\AvianPostPHD")
library(ape)
library(colorRamps)
library(caper)
library(gplots)
library(DescTools)
library(vegan)
library(bcv)
library(adephylo)
library(MDMR)
library(topGO)
library(biomaRt)
library(ALL)
library(Rgraphviz)

# pull in environment from x01_modelExe_20200612.R
load("deveL-complete-202200613.Rdata")

library(ggplot2)
library(vegan)
library(Hmisc)
library(gplots)

### i) McCorrison et. al. (Avian Transcriptomics, fibroblasts)
AvianAndMa.all <- read.table("C:\\Users\\jamis\\Desktop\\mbp\\exe.Avian\\rankedLists\\OGsignif.orderedOn.align_rank.reslogMLS.tsv", sep="\t", stringsAsFactors = FALSE)
colnames(AvianAndMa.all) <- AvianAndMa.all[1,]
AvianAndMa.all <- AvianAndMa.all[-1,]
rownames(AvianAndMa.all) <- AvianAndMa.all[,1]
AvianAndMa.all <- AvianAndMa.all[,-1]
colSums(AvianAndMa.all[ , c("MLS.p_value.Align","MLS.p_value.DeNovo","resMLS.p_value.Align",
                            "resMLS.p_value.DeNovo","reslogMLS.p_value.Align","reslogMLS.p_value.DeNovo","reslogMLS_binary.p_value.Align","reslogMLS_binary.p_value.DeNovo",
                            "reslogMLS_sd2.p_value.Align","reslogMLS_sd2.p_value.DeNovo","reslogMLS_sd1.p_value.Align","reslogMLS_sd1.p_value.DeNovo") ] < 0.005, na.rm = TRUE)
Ma.all <- AvianAndMa.all[ rowSums(is.na(AvianAndMa.all[,c("MLS.p_value.Ma","resMLS.p_value.Ma")])) < 2 , c("huref_gene_symbol","MLS.p_value.Ma","resMLS.p_value.Ma") ]
Ma.lt0.05 <- Ma.all # [ rowSums( Ma.all[,c("MLS.p_value.Ma","resMLS.p_value.Ma")] > 0.05 ) < 2 , c("huref_gene_symbol","MLS.p_value.Ma","resMLS.p_value.Ma") ]
  
### iii) Orwoll et. al (???)
# Which manuscript? https://ohsu.pure.elsevier.com/en/persons/eric-orwoll/publications/
# iiia) 1115 highlighted as associated with bone loss, longevity, mortality, or obesity.
# iiib) 25 "Longevity-associated"
# Boneloss	Longevity	Mortality	Obesity
Orwoll.all <- read.table("C:\\Users\\jamis\\Desktop\\mbp\\exe.Avian\\rankedLists\\Orwoll.eric_pep_Ncmps.txt", sep="\t", stringsAsFactors = FALSE)
colnames(Orwoll.all) <- Orwoll.all[1,]
Orwoll.all <- Orwoll.all[-1,]
rownames(Orwoll.all) <- Orwoll.all$Symbol
Orwoll.boneloss  <-  Orwoll.all[,"Boneloss"][ Orwoll.all[,"Boneloss"] != "" ]
names(Orwoll.boneloss) <- rownames(Orwoll.all)[ Orwoll.all[,"Boneloss"] != "" ]
Orwoll.longevity  <-  Orwoll.all[,"Longevity"][ Orwoll.all[,"Longevity"] != "" ]
names(Orwoll.longevity) <- rownames(Orwoll.all)[ Orwoll.all[,"Longevity"] != "" ]
Orwoll.mortality  <-  Orwoll.all[,"Mortality"][ Orwoll.all[,"Mortality"] != "" ]
names(Orwoll.mortality) <- rownames(Orwoll.all)[ Orwoll.all[,"Mortality"] != "" ]
Orwoll.obesity  <-  Orwoll.all[,"Obesity"][ Orwoll.all[,"Obesity"] != "" ]
names(Orwoll.obesity) <- rownames(Orwoll.all)[ Orwoll.all[,"Obesity"] != "" ]

###  iv) Paola et. al (SNP association, Centenarian study)
# 621 highlighted as "top hits" and most immediately highlight huref gene symbols using https://biit.cs.ut.ee/gsnpense (more than 900 avalable total)
# From these 2 studies:
#   2012: https://www.ncbi.nlm.nih.gov/pubmed/22279548
# 2013: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3808698/
Paola.symbols <- read.table("C:\\Users\\jamis\\Desktop\\mbp\\exe.Avian\\rankedLists\\Paola.topHits.GeneSymbols.csv", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
colnames(Paola.symbols) <- 	c("rank","rs-code","ChrPos","GeneNames","EnsemblGeneIDs","missense variant","splice region variant","synonymous variant","3 prime UTR variant","non coding transcript exon variant","intron variant","NMD transcript variant","non coding transcript variant","upstream gene variant","downstream gene variant")
Paola.topHits <- read.table("C:\\Users\\jamis\\Desktop\\mbp\\exe.Avian\\rankedLists\\Paola.top.hits.meta.analysis.pm5.txt", sep="\t", stringsAsFactors = FALSE)
rownames(Paola.topHits) <- Paola.topHits[,1]
Paola.pop <- Paola.symbols[ Paola.symbols$EnsemblGeneIDs != "" , ]
Paola.pop$pValue <- Paola.topHits[ Paola.pop$`rs-code` , 3 ]
Paola.pop$effect <- Paola.topHits[ Paola.pop$`rs-code` , 2 ]

### v) Peters et. al (Human transcriptomics, blood)
# "identify 1,497 genes that are differentially expressed with chronological age"
# 2015: https://www.ncbi.nlm.nih.gov/pubmed/26490707
Peters.brief <- read.table("C:\\Users\\jamis\\Desktop\\mbp\\exe.Avian\\rankedLists\\Peters.brief.txt", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
colnames(Peters.brief) <- Peters.brief[1,]
Peters.brief <- Peters.brief[-1,]
Peters.brief[ Peters.brief$`NEW-Gene-ID` == "DENND1B\xca" , ]$`NEW-Gene-ID` <- "DENND1B"
#Peters.all <- read.table("C:\\Users\\jamis\\Desktop\\mbp\\exe.Avian\\rankedLists\\Peters.johnson-2015-SD1.txt", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
#colnames(Peters.all) <- Peters.all[3,]
#Peters.all <- Peters.all[-1,]
#Peters.all <- Peters.all[-1,]
#Peters.all <- Peters.all[-1,]

### vi) Parrot longevity genes
# 31 genes highlighted with positive sleection analysis and human reference orthologs (see tab E)
# 2018: https://www.sciencedirect.com/science/article/pii/S0960982218314179
Parrot.all <- read.table("C:\\Users\\jamis\\Desktop\\mbp\\exe.Avian\\rankedLists\\ParrotPositiveLRTSelect.txt", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
colnames(Parrot.all) <- Parrot.all[1,]
Parrot.all <- Parrot.all[-1,]
Parrot.all <- Parrot.all[-1,]

### vii) Sood et al 
Sood.all <- read.table("C:\\Users\\jamis\\Desktop\\mbp\\exe.Avian\\rankedLists\\Sood.EtAl2015.tsv", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
colnames(Sood.all) <- Sood.all[1,]
Sood.all <- Sood.all[-1,]


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

# gather all unique names

allSymbols <- c(
  # ii) Ma et. al. (Rodent transcriptomics)
  Ma.lt0.005$huref_gene_symbol,
  # iii) Orwoll et. al (???)
  Orwoll.all$Symbol,
  #  iv) Paola et. al (SNP association, Centenarian study)
  unlist(strsplit(Paola.pop$EnsemblGeneIDs,",")),
  # v) Peters et. al (Human transcriptomics, blood)
  Peters.brief$`NEW-Gene-ID`,
  #Peters.all$`NEW-Gene-ID`,
  # vi) Parrot longevity genes
  Parrot.all$Gene,
  # vii) Sood et al
  Sood.all$`Gene Symbol`
)
length(allSymbols)
# 15989
length(unique(allSymbols))
# 12033
allSymbols <- allSymbols[!is.na(allSymbols)]
allSymbols <- allSymbols[allSymbols != ""]

###############################################################################################################################




qqqqqqqqqqqqqqqqqq <- as.data.frame(matrix(ncol=12,nrow=length(cc[37:96])))


signif_cmpr <- as.data.frame(matrix())


colnam <- c("Ma.MLSandMLSW",
            
            "Orwall.boneloss.longevity.mortality.obesity",
            
            "Paola.SNP",
            
            "Peters.HurefBlood.DISC.REPL.META",
            
            "Parrot.LRT",
            
            "Sood.P")

colnames(qqqqqqqqqqqqqqqqqq) <- colnam
 
for(c in cc[37:96] ){
  
      ty <- grepl( 'FDR' , c )
      if( ty == TRUE){
        prefix <- substr(c,1,nchar(c)-10)
        yylab <- '-log( FDR-adjusted(p) )'
      }else{
        prefix <- substr(c,1,nchar(c)-6)
        yylab <- '-log( p )'
      }
      pValue   <- paste0( prefix , '.p_val' )
      fdrValue <- paste0( pValue , ".FDR" )
      xmax <- max(model_performance[ , c][ !is.na(model_performance[ , c])  ])
      
      estrhocols <- colnames(model_performance)[ grepl('est.rho',colnames(model_performance)) ]
      iter <- strsplit(c,'\\.')[[1]][1]
      
      model_performance.pos <- rownames(model_performance)[ model_performance[,c] > 0 ]
      model_performance.neg <- rownames(model_performance)[ model_performance[,c] < 0 ]
      
      cutoffs <- c(0.005,0.0001,0.0005,0.00001,0.00005, 0.01,0.05,0.001)
      print(paste(nrow(model_performance)," total"))
      
      for(cutoff in cutoffs){
        
        model_performance.cut  <- model_performance[ model_performance[,pValue] < cutoff ,  ]
        model_performance.cut_fdr  <- model_performance[ model_performance[,fdrValue] < cutoff ,  ]
        
        model_performance.cut.pos <- model_performance.cut[ model_performance.cut[ , c ] > 0 ,  ]
        model_performance.cut.neg <- model_performance.cut[ model_performance.cut[ , c ] < 0 ,  ]
        model_performance.cut_fdr.pos <- model_performance.cut_fdr[ model_performance.cut_fdr[ , c ] > 0 ,  ]
        model_performance.cut_fdr.neg <- model_performance.cut_fdr[ model_performance.cut_fdr[ , c ] < 0 ,  ]
        
        model_performance.cut.pos.sym <- unique(strsplit(paste(unlist(unname(cmpq[ rownames(model_performance.cut.pos) ]))[ unlist(unname(cmpq[ rownames(model_performance.cut.pos) ])) != "" ], collapse=","), split = ","))
        model_performance.cut.neg.sym <- unique(strsplit(paste(unlist(unname(cmpq[ rownames(model_performance.cut.neg) ]))[ unlist(unname(cmpq[ rownames(model_performance.cut.neg) ])) != "" ], collapse=","), split = ","))
        model_performance.cut_fdr.pos.sym <- unique(strsplit(paste(unlist(unname(cmpq[ rownames(model_performance.cut_fdr.pos) ]))[ unlist(unname(cmpq[ rownames(model_performance.cut_fdr.pos) ])) != "" ], collapse=","), split = ","))
        model_performance.cut_fdr.neg.sym <- unique(strsplit(paste(unlist(unname(cmpq[ rownames(model_performance.cut_fdr.neg) ]))[ unlist(unname(cmpq[ rownames(model_performance.cut_fdr.neg) ])) != "" ], collapse=","), split = ","))
      
        sum_mat <- as.data.frame(matrix(ncol=6,nrow=4))
        sum_mat.str <- as.data.frame(matrix(ncol=6,nrow=4))
        colnames(sum_mat) <- colnames(sum_mat.str) <- colnam
        rownames(sum_mat) <- rownames(sum_mat.str) <- c("pos","neg","fdr.pos","fdr.neg")
        
        # ii) Ma et. al. (Rodent transcriptomics)model_performance.cut_fdr.pos
        sum_mat[1,1] <- length( Ma.lt0.05$huref_gene_symbol[ Ma.lt0.05$huref_gene_symbol                   %in% model_performance.cut.pos.sym[[1]] ] )
        # iii) Orwoll et. al (???)
        sum_mat[1,2] <- length( Orwoll.all$Symbol[            Orwoll.all$Symbol                              %in% model_performance.cut.pos.sym[[1]] ] )
        #  iv) Paola et. al (SNP association, Centenarian study)
        sum_mat[1,3] <- length( unlist(strsplit(Paola.pop$EnsemblGeneIDs,","))[ unlist(strsplit(Paola.pop$EnsemblGeneIDs,",")) %in% model_performance.cut.pos.sym[[1]] ] )
        # v) Peters et. al (Human transcriptomics, blood)
        sum_mat[1,4] <- length( Peters.brief$`NEW-Gene-ID`[ Peters.brief$`NEW-Gene-ID`                       %in% model_performance.cut.pos.sym[[1]] ] )
        # vi) Parrot longevity genes
        sum_mat[1,5] <- length( Parrot.all$Gene[ Parrot.all$Gene                                             %in% model_performance.cut.pos.sym[[1]] ] )
        # vii) Sood et al
        sum_mat[1,6] <- length( Sood.all$`Gene Symbol`[ Sood.all$`Gene Symbol`                               %in% model_performance.cut.pos.sym[[1]] ] )
        
        # ii) Ma et. al. (Rodent transcriptomics)
        sum_mat[2,1] <- length( Ma.lt0.05$huref_gene_symbol[ Ma.lt0.05$huref_gene_symbol                   %in% model_performance.cut.neg.sym[[1]] ] )
        # iii) Orwoll et. al (???)
        sum_mat[2,2] <- length( Orwoll.all$Symbol[            Orwoll.all$Symbol                              %in% model_performance.cut.neg.sym[[1]] ] )
        #  iv) Paola et. al (SNP association, Centenarian study)
        sum_mat[2,3] <- length( unlist(strsplit(Paola.pop$EnsemblGeneIDs,","))[ unlist(strsplit(Paola.pop$EnsemblGeneIDs,",")) %in% model_performance.cut.neg.sym[[1]] ] )
        # v) Peters et. al (Human transcriptomics, blood)
        sum_mat[2,4] <- length( Peters.brief$`NEW-Gene-ID`[ Peters.brief$`NEW-Gene-ID`                       %in% model_performance.cut.neg.sym[[1]] ] )
        # vi) Parrot longevity genes
        sum_mat[2,5] <- length( Parrot.all$Gene[ Parrot.all$Gene                                             %in% model_performance.cut.neg.sym[[1]] ] )
        # vii) Sood et al
        sum_mat[2,6] <- length( Sood.all$`Gene Symbol`[ Sood.all$`Gene Symbol`                               %in% model_performance.cut.neg.sym[[1]] ] )
        
        # ii) Ma et. al. (Rodent transcriptomics)
        sum_mat[3,1] <- length( Ma.lt0.05$huref_gene_symbol[ Ma.lt0.05$huref_gene_symbol                   %in% model_performance.cut_fdr.pos.sym[[1]] ] )
        # iii) Orwoll et. al (???)
        sum_mat[3,2] <- length( Orwoll.all$Symbol[            Orwoll.all$Symbol                              %in% model_performance.cut_fdr.pos.sym[[1]] ] )
        #  iv) Paola et. al (SNP association, Centenarian study)
        sum_mat[3,3] <- length( unlist(strsplit(Paola.pop$EnsemblGeneIDs,","))[ unlist(strsplit(Paola.pop$EnsemblGeneIDs,",")) %in% model_performance.cut_fdr.pos.sym[[1]] ] )
        # v) Peters et. al (Human transcriptomics, blood)
        sum_mat[3,4] <- length( Peters.brief$`NEW-Gene-ID`[ Peters.brief$`NEW-Gene-ID`                       %in% model_performance.cut_fdr.pos.sym[[1]] ] )
        # vi) Parrot longevity genes
        sum_mat[3,5] <- length( Parrot.all$Gene[ Parrot.all$Gene                                             %in% model_performance.cut_fdr.pos.sym[[1]] ] )
        # vii) Sood et al
        sum_mat[3,6] <- length( Sood.all$`Gene Symbol`[ Sood.all$`Gene Symbol`                               %in% model_performance.cut_fdr.pos.sym ] )
        
        # ii) Ma et. al. (Rodent transcriptomics)
        sum_mat[4,1] <- length( Ma.lt0.05$huref_gene_symbol[ Ma.lt0.05$huref_gene_symbol                   %in% model_performance.cut_fdr.neg.sym[[1]] ] )
        # iii) Orwoll et. al (???)
        sum_mat[4,2] <- length( Orwoll.all$Symbol[            Orwoll.all$Symbol                              %in% model_performance.cut_fdr.neg.sym[[1]] ] )
        #  iv) Paola et. al (SNP association, Centenarian study)
        sum_mat[4,3] <- length( unlist(strsplit(Paola.pop$EnsemblGeneIDs,","))[ unlist(strsplit(Paola.pop$EnsemblGeneIDs,",")) %in% model_performance.cut_fdr.neg.sym[[1]] ] )
        # v) Peters et. al (Human transcriptomics, blood)
        sum_mat[4,4] <- length( Peters.brief$`NEW-Gene-ID`[ Peters.brief$`NEW-Gene-ID`                       %in% model_performance.cut_fdr.neg.sym[[1]] ] )
        # vi) Parrot longevity genes
        sum_mat[4,5] <- length( Parrot.all$Gene[ Parrot.all$Gene                                             %in% model_performance.cut_fdr.neg.sym[[1]] ] )
        # vii) Sood et al
        sum_mat[4,6] <- length( Sood.all$`Gene Symbol`[ Sood.all$`Gene Symbol`                               %in% model_performance.cut_fdr.neg.sym[[1]] ] )
        
        # print(paste( prefix , cutoff ) )
        # print(sum_mat)
        
        
        # ii) Ma et. al. (Rodent transcriptomics)model_performance.cut_fdr.pos
        sum_mat.str[1,1] <- paste( Ma.lt0.05$huref_gene_symbol[ Ma.lt0.05$huref_gene_symbol                   %in% model_performance.cut.pos.sym[[1]] ] , collapse = "," )
        # iii) Orwoll et. al (???)
        sum_mat.str[1,2] <- paste( Orwoll.all$Symbol[            Orwoll.all$Symbol                              %in% model_performance.cut.pos.sym[[1]] ] , collapse = "," )
        #  iv) Paola et. al (SNP association, Centenarian study)
        sum_mat.str[1,3] <- paste( unlist(strsplit(Paola.pop$EnsemblGeneIDs,","))[ unlist(strsplit(Paola.pop$EnsemblGeneIDs,",")) %in% model_performance.cut.pos.sym[[1]] ] , collapse = "," )
        # v) Peters et. al (Human transcriptomics, blood)
        sum_mat.str[1,4] <- paste( Peters.brief$`NEW-Gene-ID`[ Peters.brief$`NEW-Gene-ID`                       %in% model_performance.cut.pos.sym[[1]] ] , collapse = "," )
        # vi) Parrot longevity genes
        sum_mat.str[1,5] <- paste( Parrot.all$Gene[ Parrot.all$Gene                                             %in% model_performance.cut.pos.sym[[1]] ] , collapse = "," )
        # vii) Sood et al
        sum_mat.str[1,6] <- paste( Sood.all$`Gene Symbol`[ Sood.all$`Gene Symbol`                               %in% model_performance.cut.pos.sym[[1]] ] , collapse = ",")
        
        # ii) Ma et. al. (Rodent transcriptomics)
        sum_mat.str[2,1] <- paste( Ma.lt0.05$huref_gene_symbol[ Ma.lt0.05$huref_gene_symbol                   %in% model_performance.cut.neg.sym[[1]] ] , collapse = ",")
        # iii) Orwoll et. al (???)
        sum_mat.str[2,2] <- paste( Orwoll.all$Symbol[            Orwoll.all$Symbol                              %in% model_performance.cut.neg.sym[[1]] ] , collapse = ",")
        #  iv) Paola et. al (SNP association, Centenarian study)
        sum_mat.str[2,3] <- paste( unlist(strsplit(Paola.pop$EnsemblGeneIDs,","))[ unlist(strsplit(Paola.pop$EnsemblGeneIDs,",")) %in% model_performance.cut.neg.sym[[1]] ] , collapse = ",")
        # v) Peters et. al (Human transcriptomics, blood)
        sum_mat.str[2,4] <- paste( Peters.brief$`NEW-Gene-ID`[ Peters.brief$`NEW-Gene-ID`                       %in% model_performance.cut.neg.sym[[1]] ] , collapse = "," )
        # vi) Parrot longevity genes
        sum_mat.str[2,5] <- paste( Parrot.all$Gene[ Parrot.all$Gene                                             %in% model_performance.cut.neg.sym[[1]] ] , collapse = ",")
        # vii) Sood et al
        sum_mat.str[2,6] <- paste( Sood.all$`Gene Symbol`[ Sood.all$`Gene Symbol`                               %in% model_performance.cut.neg.sym[[1]] ] , collapse = ",")
        
        # ii) Ma et. al. (Rodent transcriptomics)
        sum_mat.str[3,1] <- paste( Ma.lt0.05$huref_gene_symbol[ Ma.lt0.05$huref_gene_symbol                   %in% model_performance.cut_fdr.pos.sym[[1]] ] , collapse = ",")
        # iii) Orwoll et. al (???)
        sum_mat.str[3,2] <- paste( Orwoll.all$Symbol[            Orwoll.all$Symbol                              %in% model_performance.cut_fdr.pos.sym[[1]] ] , collapse = ",")
        #  iv) Paola et. al (SNP association, Centenarian study)
        sum_mat.str[3,3] <- paste( unlist(strsplit(Paola.pop$EnsemblGeneIDs,","))[ unlist(strsplit(Paola.pop$EnsemblGeneIDs,",")) %in% model_performance.cut_fdr.pos.sym[[1]] ], collapse = "," )
        # v) Peters et. al (Human transcriptomics, blood)
        sum_mat.str[3,4] <- paste( Peters.brief$`NEW-Gene-ID`[ Peters.brief$`NEW-Gene-ID`                       %in% model_performance.cut_fdr.pos.sym[[1]] ] , collapse = ",")
        # vi) Parrot longevity genes
        sum_mat.str[3,5] <- paste( Parrot.all$Gene[ Parrot.all$Gene                                             %in% model_performance.cut_fdr.pos.sym[[1]] ] , collapse = ",")
        # vii) Sood et al
        sum_mat.str[3,6] <- paste( Sood.all$`Gene Symbol`[ Sood.all$`Gene Symbol`                               %in% model_performance.cut_fdr.pos.sym ] , collapse = ",")
        
        # ii) Ma et. al. (Rodent transcriptomics)
        sum_mat.str[4,1] <- paste( Ma.lt0.05$huref_gene_symbol[ Ma.lt0.05$huref_gene_symbol                   %in% model_performance.cut_fdr.neg.sym[[1]] ] , collapse = ",")
        # iii) Orwoll et. al (???)
        sum_mat.str[4,2] <- paste( Orwoll.all$Symbol[            Orwoll.all$Symbol                              %in% model_performance.cut_fdr.neg.sym[[1]] ] , collapse = ",")
        #  iv) Paola et. al (SNP association, Centenarian study)
        sum_mat.str[4,3] <- paste( unlist(strsplit(Paola.pop$EnsemblGeneIDs,","))[ unlist(strsplit(Paola.pop$EnsemblGeneIDs,",")) %in% model_performance.cut_fdr.neg.sym[[1]] ] , collapse = ",")
        # v) Peters et. al (Human transcriptomics, blood)
        sum_mat.str[4,4] <- paste( Peters.brief$`NEW-Gene-ID`[ Peters.brief$`NEW-Gene-ID`                       %in% model_performance.cut_fdr.neg.sym[[1]] ] , collapse = ",")
        # vi) Parrot longevity genes
        sum_mat.str[4,5] <- paste( Parrot.all$Gene[ Parrot.all$Gene                                             %in% model_performance.cut_fdr.neg.sym[[1]] ] , collapse = ",")
        # vii) Sood et al
        sum_mat.str[4,6] <- paste( Sood.all$`Gene Symbol`[ Sood.all$`Gene Symbol`                               %in% model_performance.cut_fdr.neg.sym[[1]] ] , collapse = ",")
        
        # print( sum_mat.str )
        
        write.table(sum_mat, paste0('overlaps/geneHurefSym_overlapFrequency.',prefix,'.',cutoff,'.tsv'), sep='\t')
        write.table(sum_mat.str, paste0('overlaps/geneHurefSym_overlapGeneNames.',prefix,'.',cutoff,'.tsv'), sep='\t')
        
      }
    
}
 

