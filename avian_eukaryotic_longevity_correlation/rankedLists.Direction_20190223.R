#setwd ("/Users/jmccorri/Desktop/exe.Avian")

library(ggplot2)
library(vegan)
library(ggplot2)
library(Hmisc)
library(gplots)
setwd ("~/Desktop/exe.Avian/rankedLists")

### i) McCorrison et. al. (Avian Transcriptomics, fibroblasts)
AvianAndMa.all <- read.table("OGsignif.orderedOn.align_rank.reslogMLS.tsv", sep="\t", stringsAsFactors = FALSE)
colnames(AvianAndMa.all) <- AvianAndMa.all[1,]
AvianAndMa.all <- AvianAndMa.all[-1,]
rownames(AvianAndMa.all) <- AvianAndMa.all[,1]
AvianAndMa.all <- AvianAndMa.all[,-1]
colSums(AvianAndMa.all[ , c("MLS.p_value.Align","MLS.p_value.DeNovo","resMLS.p_value.Align",
                            "resMLS.p_value.DeNovo","reslogMLS.p_value.Align","reslogMLS.p_value.DeNovo","reslogMLS_binary.p_value.Align","reslogMLS_binary.p_value.DeNovo",
                            "reslogMLS_sd2.p_value.Align","reslogMLS_sd2.p_value.DeNovo","reslogMLS_sd1.p_value.Align","reslogMLS_sd1.p_value.DeNovo") ] < 0.005, na.rm = TRUE)
keepers <- list()
ki <- 1
for (c in c("MLS.p_value.Align","MLS.p_value.DeNovo","resMLS.p_value.Align",
                                  "resMLS.p_value.DeNovo","reslogMLS.p_value.Align","reslogMLS.p_value.DeNovo","reslogMLS_binary.p_value.Align","reslogMLS_binary.p_value.DeNovo",
                                  "reslogMLS_sd2.p_value.Align","reslogMLS_sd2.p_value.DeNovo","reslogMLS_sd1.p_value.Align","reslogMLS_sd1.p_value.DeNovo")){
    subA <- AvianAndMa.all[ !is.na(AvianAndMa.all[,c]) , ]
    for( t in rownames(subA[ as.double(subA[,c]) < 0.005 , ]) ){
      keepers[ki] <- t
      ki <- ki + 1
    }
}
AvianAndMa.AboveCutoff <- AvianAndMa.all[ unique(unlist(keepers)) , ]

#apply cutoffs
# ia) MLS
# ib) MLS.denovo
# id) MLSW.denovo
# ic) MLSW 
# ie) MLSLW
# if) MLSLW.deNovo
# ig) MLSLW.dropout
# ih) MLSLW.std.dev.1.extreme.longevity
# ii) MLSLW.std.dev.1POS.extreme.longevity
# ij) MLSLW.std.dev.2.extreme.longevity
# ik) MLSLW.std.dev.2POS.extreme.longevity
# ii) Ma et. al. (Rodent transcriptomics)
# Significant cutoff to use?
#  https://elifesciences.org/articles/19130


### iii) Orwoll et. al (???)
# Which manuscript? https://ohsu.pure.elsevier.com/en/persons/eric-orwoll/publications/
# iiia) 1115 highlighted as associated with bone loss, longevity, mortality, or obesity.
# iiib) 25 "Longevity-associated"
# Boneloss	Longevity	Mortality	Obesity
Orwoll.all <- read.table("Orwoll.eric_pep_Ncmps.txt", sep="\t", stringsAsFactors = FALSE)
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
Paola.symbols <- read.table("Paola.topHits.GeneSymbols.csv", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
colnames(Paola.symbols) <- 	c("rank","rs-code","ChrPos","GeneNames","EnsemblGeneIDs","missense variant","splice region variant","synonymous variant","3 prime UTR variant","non coding transcript exon variant","intron variant","NMD transcript variant","non coding transcript variant","upstream gene variant","downstream gene variant")
Paola.topHits <- read.table("Paola.top.hits.meta.analysis.pm5.txt", sep="\t", stringsAsFactors = FALSE)
rownames(Paola.topHits) <- Paola.topHits[,1]
Paola.pop <- Paola.symbols[ Paola.symbols$EnsemblGeneIDs != "" , ]
Paola.pop$pValue <- Paola.topHits[ Paola.pop$`rs-code` , 3 ]
Paola.pop$effect <- Paola.topHits[ Paola.pop$`rs-code` , 2 ]

### v) Peters et. al (Human transcriptomics, blood)
# "identify 1,497 genes that are differentially expressed with chronological age"
# 2015: https://www.ncbi.nlm.nih.gov/pubmed/26490707
Peters.brief <- read.table("Peters.brief.txt", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
colnames(Peters.brief) <- Peters.brief[1,]
Peters.brief <- Peters.brief[-1,]
Peters.brief[ Peters.brief$`NEW-Gene-ID` == "DENND1B\xca" , ]$`NEW-Gene-ID` <- "DENND1B"
#Peters.all <- read.table("Peters.johnson-2015-SD1.txt", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
#colnames(Peters.all) <- Peters.all[3,]
#Peters.all <- Peters.all[-1,]
#Peters.all <- Peters.all[-1,]
#Peters.all <- Peters.all[-1,]
 
### vi) Parrot longevity genes
# 31 genes highlighted with positive sleection analysis and human reference orthologs (see tab E)
# 2018: https://www.sciencedirect.com/science/article/pii/S0960982218314179
Parrot.all <- read.table("ParrotPositiveLRTSelect.txt", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
colnames(Parrot.all) <- Parrot.all[1,]
Parrot.all <- Parrot.all[-1,]
Parrot.all <- Parrot.all[-1,]
   
### vii) Sood et al 
Sood.all <- read.table("Sood.EtAl2015.tsv", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
colnames(Sood.all) <- Sood.all[1,]
Sood.all <- Sood.all[-1,]


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

# gather all unique names

allSymbols <- c(
# i) McCorrison et. al. (Avian Transcriptomics, fibroblasts)
# ii) Ma et. al. (Rodent transcriptomics)
AvianAndMa.AboveCutoff$huref_gene_symbol,

# ma in
toupper(MaMLS$V10)[ toupper(MaMLS$V10) != "Symbol" ][!is.na( toupper(MaMLS$V10)[ toupper(MaMLS$V10) != "Symbol" ] )],
toupper(MaMLSW$V10)[ toupper(MaMLSW$V10) != "Symbol" ][!is.na( toupper(MaMLSW$V10)[ toupper(MaMLSW$V10) != "Symbol" ] )],

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


MaMLS$V10 <- toupper(MaMLS$V10)
MaMLSW$V10 <- toupper(MaMLSW$V10)
 
colnam <- c("Avian.MLS","Avian.MLS.denovo",
            "Avian.MLSW","Avian.MLSW.denovo",
            "Avian.MLSLW","Avian.MLSLW.deNovo",
            "Avian.MLSLW.dropout","Avian.MLSLW.dropout.denovo",
            "Avian.MLSLW.std.dev.1.extreme.longevity","Avian.MLSLW.std.dev.1POS.extreme.longevity",
            "Avian.MLSLW.std.dev.2.extreme.longevity","Avian.MLSLW.std.dev.2POS.extreme.longevity",
            
            "Ma.MLS","Ma.MLSW",
            
            "Orwall.boneloss", "Orwall.longevity", "Orwall.mortality", "Orwall.obesity",
            
            "Paola.SNP",
            
            "Peters.HurefBlood.DISC","Peters.HurefBlood.REPL","Peters.HurefBlood.META",
            
            "Parrot.LRT",
            
            "Sood.P"
)

pvalTab <- matrix(nrow = length(unique(allSymbols)),
                 ncol = length(colnam))
colnames(pvalTab) <- colnam
rownames(pvalTab) <- unique(allSymbols)

for(x in rownames(pvalTab)){
if(!is.na(x)){
  # i,ii)
  if(x %in% AvianAndMa.AboveCutoff$huref_gene_symbol){
    pvalTab[x,c(1:12)] <- as.double(unname(unlist(
                              apply(
                                AvianAndMa.AboveCutoff[ AvianAndMa.AboveCutoff$huref_gene_symbol == x , 
                                                     c("MLS.p_value.Align","MLS.p_value.DeNovo","resMLS.p_value.Align","resMLS.p_value.DeNovo","reslogMLS.p_value.Align","reslogMLS.p_value.DeNovo",
                                                       "reslogMLS_binary.p_value.Align","reslogMLS_binary.p_value.DeNovo",
                                                       "reslogMLS_sd1.p_value.Align","reslogMLS_sd1.p_value.DeNovo","reslogMLS_sd2.p_value.Align","reslogMLS_sd2.p_value.DeNovo"   
                                                      ) ]
                              , 2, FUN=min)
                          )))
    # pvalTab[x,c(13:14)] <- as.double(unname(unlist(AvianAndMa.AboveCutoff[ AvianAndMa.AboveCutoff$huref_gene_symbol == x , c("MLS.p_value.Ma","resMLS.p_value.Ma") ][1,])))
  }
  
  if(x %in% MaMLS$V10){
    pvalTab[x,13] <- MaMLS[ MaMLS$V10 == x , 6 ]
  }
  if(x %in% MaMLSW$V10){
    pvalTab[x,14] <- MaMLSW[ MaMLSW$V10 == x , 6 ]
  }
  
  # iii) orwall
  if(x %in% names(Orwoll.boneloss)){
    pvalTab[x,15] <- 1
  }
  if(x %in% names(Orwoll.longevity)){
    pvalTab[x,16] <- 1
  }
  if(x %in% names(Orwoll.mortality)){
    pvalTab[x,17] <- 1
  }
  if(x %in% names(Orwoll.obesity)){
    pvalTab[x,18] <- 1
  }
  # iv)  Paola et. al (SNP association, Centenarian study)
  if(sum(grepl(x, Paola.pop$EnsemblGeneIDs)) > 0){
    pvalTab[x,19] <-min(as.double(Paola.pop[ grepl(x, Paola.pop$EnsemblGeneIDs) , ]$pValue))
  }
  
  # v)   Peters et. al (Human transcriptomics, blood)
  if(x %in% Peters.all$`NEW-Gene-ID`){
    pvalTab[x,c(20:22)] <- as.double(unname(unlist(Peters.all[Peters.all$`NEW-Gene-ID` == x,c(11,16,20)][1,])))
  }
  
  # vi)  Parrot longevity genes
  if(x %in% Parrot.all$Gene){
    pvalTab[x,23] <-as.double(Parrot.all[ Parrot.all$Gene == x , 2 ])
  }
  
  # vii) Sood et al
  if(x %in% Sood.all$`Gene Symbol`){
    pvalTab[x,24] <- 1
  }
  
}
}

write.table(pvalTab, file = "pvalTab.20190208.tsv", sep = "\t")

# save.image(file=paste0("pvalTab.20190206.Rdata"))
# load(file=paste0("pvalTab.20190206.Rdata"))


rowSums(!is.na(pvalTab))
 
# this subroutine is why R is the worst
pq <- matrix(nrow=dim(pvalTab)[1],ncol=dim(pvalTab)[2])
rownames(pq) <- rownames(pvalTab)
colnames(pq) <- colnames(pvalTab)
for(x in rownames(pvalTab)){
  for(y in colnames(pvalTab)){
    pq[x,y] <- as.double(pvalTab[x,y])
  }
}
pq[is.na(pq)] <- -1

png("20190208.top1000overlaps.png", width = 1000, height=1000)
heatmap.2( 
  
            pq[ names(sort(rowSums(!is.na(pq))))[ 1:1000 ] , ],
            trace='none',
            col=c(rep("black",50),redblue(50)),
            margins=c(22,10),
            keysize=0.7, key.par = list(cex=0.5),
            rowV = FALSE
           
        )
dev.off()

png("20190208.top100overlaps.png", width = 1000, height=1000)
heatmap.2( 
  
  pq[ names(sort(rowSums(!is.na(pq))))[ 1:100 ] , ],
  trace='none',
  col=c(rep("black",50),redblue(50)),
  margins=c(22,10),
  keysize=0.7, key.par = list(cex=0.5),

)
dev.off()


png("20190208.top500overlaps.less.png", width = 1000, height=1000)
heatmap.2( 
  
  pq[ names(rev(sort(rowSums(pq[  ,
                                    c("Avian.MLSLW","Ma.MLS","Ma.MLSW",
                                      "Orwall.boneloss","Orwall.longevity","Orwall.mortality","Orwall.obesity",
                                      "Paola.SNP","Peters.HurefBlood.META","Parrot.LRT","Sood.P")
                                   ] > -1))))[ 1:500 ] , c("Avian.MLSLW","Avian.MLSLW.deNovo","Ma.MLS","Ma.MLSW",
                                                      "Orwall.boneloss","Orwall.longevity","Orwall.mortality","Orwall.obesity",
                                                      "Paola.SNP","Peters.HurefBlood.META","Parrot.LRT","Sood.P") ],
  trace='none',
  col=c(rep("black",50),redblue(50)),
  margins=c(22,10),
  keysize=0.7, key.par = list(cex=0.5),
  
)
dev.off()



# names(!is.na(pvalTab[ , "Ma.MLS" ]) == TRUE) %in% MaMLS$V10
# FIXED ! 
# rownames(pvalTab[ !is.na(pvalTab[ , "Ma.MLS" ] ) , ]) %in% MaMLS$V10






#avian in

alignment.CAPER.MLS_v_binaryAbs  <- read.table("../RERUN.alignment.CAPER.MLS_v_binaryAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

alignment.CAPER.MLS_v_relAbs  <- read.table("../RERUN.alignment.CAPER.MLS_v_relAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

alignment.CAPER.logwMLSresid_v_binaryAbs  <- read.table("../RERUN.alignment.CAPER.logwMLSresid_v_binaryAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

alignment.CAPER.logwMLSresid_v_relAbs  <- read.table("../RERUN.alignment.CAPER.logwMLSresid_v_relAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

alignment.CAPER.wMLSresid_v_binaryAbs  <- read.table("../RERUN.alignment.CAPER.wMLSresid_v_binaryAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

alignment.CAPER.wMLSresid_v_relAbs  <- read.table("../RERUN.alignment.CAPER.wMLSresid_v_relAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

denovo_highConfidence.CAPER.MLS_v_binaryAbs  <- read.table("../RERUN.denovo_highConfidence.CAPER.MLS_v_binaryAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

denovo_highConfidence.CAPER.MLS_v_relAbs  <- read.table("../RERUN.denovo_highConfidence.CAPER.MLS_v_relAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

denovo_highConfidence.CAPER.logwMLSresid_v_binaryAbs  <- read.table("../RERUN.denovo_highConfidence.CAPER.logwMLSresid_v_binaryAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

denovo_highConfidence.CAPER.logwMLSresid_v_relAbs  <- read.table("../RERUN.denovo_highConfidence.CAPER.logwMLSresid_v_relAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

denovo_highConfidence.CAPER.wMLSresid_v_binaryAbs  <- read.table("../RERUN.denovo_highConfidence.CAPER.wMLSresid_v_binaryAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

denovo_highConfidence.CAPER.wMLSresid_v_relAbs  <- read.table("../RERUN.denovo_highConfidence.CAPER.wMLSresid_v_relAbs.report.tsv", sep="\t", stringsAsFactors = FALSE)

denovo_highConfidence.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport  <- read.table("../denovo_highConfidence.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport.tsv", sep="\t", stringsAsFactors = FALSE)

denovo_highConfidence.20.extremeLongevity_sd2_hi.logwMLSresid_v_relAbsreport  <- read.table("../denovo_highConfidence.20.extremeLongevity_sd2_hi.logwMLSresid_v_relAbsreport.tsv", sep="\t", stringsAsFactors = FALSE)

denovo_highConfidence.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport  <- read.table("../denovo_highConfidence.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport.tsv", sep="\t", stringsAsFactors = FALSE)

denovo_highConfidence.26.extremeLongevity_sd1andsd2_hi.logwMLSresid_v_relAbsreport  <- read.table("../denovo_highConfidence.26.extremeLongevity_sd1andsd2_hi.logwMLSresid_v_relAbsreport.tsv", sep="\t", stringsAsFactors = FALSE)

alignment.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport  <- read.table("../alignment.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport.tsv", sep="\t", stringsAsFactors = FALSE)

alignment.20.extremeLongevity_sd2_hi.logwMLSresid_v_relAbsreport  <- read.table("../alignment.20.extremeLongevity_sd2_hi.logwMLSresid_v_relAbsreport.tsv", sep="\t", stringsAsFactors = FALSE)

alignment.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport  <- read.table("../alignment.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport.tsv", sep="\t", stringsAsFactors = FALSE)

alignment.26.extremeLongevity_sd1andsd2_hi.logwMLSresid_v_relAbsreport  <- read.table("../alignment.26.extremeLongevity_sd1andsd2_hi.logwMLSresid_v_relAbsreport.tsv", sep="\t", stringsAsFactors = FALSE)





dirTab <- matrix(nrow = length(unique(allSymbols)),
                  ncol = length(colnam))
colnames(dirTab) <- colnam
rownames(dirTab) <- unique(allSymbols)

for(x in rownames(dirTab)){
  if(!is.na(x)){
    # i,ii)
    if(x %in% AvianAndMa.AboveCutoff$huref_gene_symbol){
      #pvalTab[x,c(1:12)] <- as.double(unname(unlist(
      #  apply(
      #    AvianAndMa.AboveCutoff[ AvianAndMa.AboveCutoff$huref_gene_symbol == x , 
      #                            c("MLS.p_value.Align","MLS.p_value.DeNovo","resMLS.p_value.Align","resMLS.p_value.DeNovo","reslogMLS.p_value.Align","reslogMLS.p_value.DeNovo",
      #                              "reslogMLS_binary.p_value.Align","reslogMLS_binary.p_value.DeNovo",
      #                              "reslogMLS_sd1.p_value.Align","reslogMLS_sd1.p_value.DeNovo","reslogMLS_sd2.p_value.Align","reslogMLS_sd2.p_value.DeNovo"   
      #                            ) ]
      #    , 2, FUN=min)
      #)))
      #pvalTab[x,c(13:14)] <- as.double(unname(unlist(AvianAndMa.AboveCutoff[ AvianAndMa.AboveCutoff$huref_gene_symbol == x , c("MLS.p_value.Ma","resMLS.p_value.Ma") ][1,])))
      
      
      #x <- "KLHL26"
      clustx <- rownames(AvianAndMa.AboveCutoff[ AvianAndMa.AboveCutoff$huref_gene_symbol == x , ])
      clustx <- paste0("orthoMCL_",strsplit(clustx,'r')[[1]][2])
      if(clustx %in% rownames(alignment.CAPER.MLS_v_relAbs)){
        dirTab[x,1] <- alignment.CAPER.MLS_v_relAbs[ clustx , 'slopes' ]
      }
      if(clustx %in% rownames(alignment.CAPER.wMLSresid_v_relAbs)){
        dirTab[x,3] <- alignment.CAPER.wMLSresid_v_relAbs[ clustx , 'slopes' ]
      }
      if(clustx %in% rownames(alignment.CAPER.logwMLSresid_v_relAbs)){
        dirTab[x,5] <- alignment.CAPER.logwMLSresid_v_relAbs[ clustx , 'slopes' ]
      }
      if(clustx %in% rownames(alignment.CAPER.logwMLSresid_v_binaryAbs)){
        dirTab[x,7] <- alignment.CAPER.logwMLSresid_v_binaryAbs[ clustx , 'slopes' ]
      }
      if(clustx %in% rownames(alignment.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport)){
        dirTab[x,9] <- alignment.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport[ clustx , 'slopes' ]
      }
      if(clustx %in% rownames(alignment.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport)){
        dirTab[x,11] <- alignment.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport[ clustx , 'slopes' ]
      }
      
      
      
      clustx <- rownames(AvianAndMa.AboveCutoff[ AvianAndMa.AboveCutoff$huref_gene_symbol == x , ])
      clustnum <- strsplit(clustx,'r')[[1]][2]
      pad <- 5 - nchar(strsplit(clustx,'r')[[1]][2])
      clustx <- paste0("MCL_",paste0(rep(0,pad),collapse =""),clustnum)
      clustx
      #clustx <- paste0("orthoMCL_",strsplit(clustx,'r')[[1]][2])
      if(clustx %in% rownames(denovo_highConfidence.CAPER.MLS_v_relAbs)){
        dirTab[x,2] <- denovo_highConfidence.CAPER.MLS_v_relAbs[ clustx , 'slopes' ]
      }
      if(clustx %in% rownames(denovo_highConfidence.CAPER.wMLSresid_v_relAbs)){
        dirTab[x,4] <- denovo_highConfidence.CAPER.wMLSresid_v_relAbs[ clustx , 'slopes' ]
      }
      if(clustx %in% rownames(denovo_highConfidence.CAPER.logwMLSresid_v_relAbs)){
        dirTab[x,6] <- denovo_highConfidence.CAPER.logwMLSresid_v_relAbs[ clustx , 'slopes' ]
      }
      if(clustx %in% rownames(denovo_highConfidence.CAPER.logwMLSresid_v_binaryAbs)){
        dirTab[x,8] <- denovo_highConfidence.CAPER.logwMLSresid_v_binaryAbs[ clustx , 'slopes' ]
      }
      if(clustx %in% rownames(denovo_highConfidence.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport)){
        dirTab[x,10] <- denovo_highConfidence.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport[ clustx , 'slopes' ]
      }
      if(clustx %in% rownames(denovo_highConfidence.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport)){
        dirTab[x,12] <- denovo_highConfidence.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport[ clustx , 'slopes' ]
      }
      
      
    }
    
    
    #x <- "POLR3H"
    if(x %in% MaMLS$V10){
      dirTab[x,13] <- MaMLS[ MaMLS$V10 == x , 3 ][1]
    }
    if(x %in% MaMLSW$V10){
      dirTab[x,13] <- MaMLSW[ MaMLSW$V10 == x , 3 ][1]
    }
    
    # iii) orwall
    if(x %in% names(Orwoll.boneloss)){
      dirTab[x,15] <- Orwoll.boneloss[x] # <- SAME
    }
    if(x %in% names(Orwoll.longevity)){
      dirTab[x,16] <- Orwoll.longevity[x] # <- SAME
    }
    if(x %in% names(Orwoll.mortality)){
      dirTab[x,17] <- Orwoll.mortality[x] # <- SAME
    }
    if(x %in% names(Orwoll.obesity)){
      dirTab[x,18] <- Orwoll.obesity[x] # <- SAME
    }
    # iv)  Paola et. al (SNP association, Centenarian study)
    if(sum(grepl(x, Paola.pop$EnsemblGeneIDs)) > 0){
      matchp <- Paola.pop[ grepl(x, Paola.pop$EnsemblGeneIDs) , ]
      if(dim(matchp)[1] > 1){
        dirTab[x,19] <-  as.double(matchp[ as.double(matchp$pValue) == min(as.double(matchp$pValue)) , ]$effect)[1]
      }else{
        dirTab[x,19] <-  as.double(matchp$effect)
      }
    }
    
    # v)   Peters et. al (Human transcriptomics, blood)
    if(x %in% Peters.all$`NEW-Gene-ID`){
      sigs <- unlist(Peters.all[Peters.all$`NEW-Gene-ID` == x,c(10,15,19)][1,])
      sigs[sigs == "-"] <- -1
      sigs[sigs == "+"] <- 1
      dirTab[x,c(20:22)] <- sigs
    }
    
    # vi)  Parrot longevity genes
    if(x %in% Parrot.all$Gene){
      dirTab[x,23] <-as.double(Parrot.all[ Parrot.all$Gene == x , 2 ])
    }
    
    # vii) Sood et al
    if(x %in% Sood.all$`Gene Symbol`){
      if( Sood.all[ Sood.all$`Gene Symbol` == x , ]$`Ratio of Y:0 muscle` == "up"){
        dirTab[x,24] <- 1
      }else if( Sood.all[ Sood.all$`Gene Symbol` == x , ]$`Ratio of Y:0 muscle` == "down"){
        dirTab[x,24] <- -1
      }
    }
    
  }
}


dim(pvalTab)
dim(dirTab)
colnames(dirTab) <- paste0("DIRECTION.",colnames(dirTab))
 
allTab <- t(rbind(t(pvalTab),t(dirTab)))


colSums(is.na(allTab))

save.image(file="20190209.Rdata")
  
n <- is.na(allTab[,"DIRECTION.Avian.MLSLW.deNovo"])
nn <- n[n==TRUE]
n <- is.na(allTab[,"Avian.MLSLW.deNovo"])
nnn <- n[n==TRUE]



write.table(allTab, file="pvalue_and_direction.20180209.tsv", sep="\t")
 





##########################



projlist <- c(
  "Ma.MLS"       ,                              "Ma.MLSW"        ,                            "Orwall.boneloss"        ,                    "Orwall.longevity"     ,                     
  "Orwall.mortality"                ,           "Orwall.obesity"       ,                      "Paola.SNP"                ,                  "Peters.HurefBlood.DISC"   ,                 
  "Peters.HurefBlood.REPL"        ,             "Peters.HurefBlood.META"    ,                 "Parrot.LRT"                       ,          "Sood.P"  ,
  "Avian.MLS"                            ,      "Avian.MLS.denovo"        ,                   "Avian.MLSW"                 ,                "Avian.MLSW.denovo" ,                        
  "Avian.MLSLW"                         ,       "Avian.MLSLW.deNovo"             ,            "Avian.MLSLW.dropout"         ,               "Avian.MLSLW.dropout.denovo"             ,   
  "Avian.MLSLW.std.dev.1.extreme.longevity"  ,  "Avian.MLSLW.std.dev.1POS.extreme.longevity" ,"Avian.MLSLW.std.dev.2.extreme.longevity"   , "Avian.MLSLW.std.dev.2POS.extreme.longevity" )


outmat <- matrix(nrow=dim(allTab)[1],ncol=7)
rownames(outmat) <- rownames(allTab)
colnames(outmat) <- c("signif_avian_models",
              "signif_models_overall",
              "signif_studies",
              "signifUP_avian_models",
              "signifUP_models_overall",
              "signifDOWN_avian_models",
              "signifDOWN_models_overall")

signif_avian_models       <- list()
signif_models_overall     <- list()
signif_studies <- list()
signifUP_avian_models       <- list()
signifUP_models_overall     <- list()
signifDOWN_avian_models       <- list()
signifDOWN_models_overall     <- list()

ri <- 1
for( r in rownames(allTab) ){
  
  #r <- rownames(allTab) [1]
  sum(allTab[ r , projlist[ grepl("Avian", projlist) ] ] < 0.005, na.rm=TRUE)
  n <- allTab[ r , projlist[ grepl("Avian", projlist) ] ] < 0.005
  na <-names( n[ is.na(n) ])
  n <- n[ !is.na(n)]
  nt <- n[ n == TRUE]
  nf <- n[ n == FALSE]
  na
  names(nt)
  names(nf)
  
  outmat[ r , "signif_avian_models" ] <- length(nt)
  
  if(length(nt) > 0){
    dirs <- table(as.numeric( allTab[ r , paste0("DIRECTION.", names(nt)) ] > 0 ))
    outmat[ r , "signifUP_avian_models"]    <- table(as.numeric( allTab[ r , paste0("DIRECTION.", names(nt)) ] > 0 ))[1]
    outmat[ r , "signifDOWN_avian_models"]  <- table(as.numeric( allTab[ r , paste0("DIRECTION.", names(nt)) ] > 0 ))[2]
  }
  
  sum(allTab[ r , projlist[ !grepl("Avian", projlist) ] ] < 0.005, na.rm=TRUE)
  n <- allTab[ r , projlist[ !grepl("Avian", projlist) ] ] #[ !is.na(allTab[ r , projlist[ !grepl("Avian", projlist) ] ] ) ] 
  #allTab[ r , projlist[ !grepl("Avian", projlist) ] ] < 0.005
  na <- names( n[ is.na(n) ])
  n <- n[ !is.na(n)]
  #nt <- n[ n == TRUE]
  #nf <- n[ n == FALSE]
  na
  names(n)
  #names(nt)
  #names(nf)
  
  outmat[ r , "signif_models_overall"] <- length(n)
  
  if(length(n) > 0){
    dirs <- as.numeric( allTab[ r , paste0("DIRECTION.", names(n)) ] > 0 )
    outmat[ r , "signifUP_models_overall"]    <- table(as.numeric( allTab[ r , paste0("DIRECTION.", names(n)) ] > 0 ))[1]
    outmat[ r , "signifDOWN_models_overall"]  <- table(as.numeric( allTab[ r , paste0("DIRECTION.", names(n)) ] > 0 ))[2]
  }
  
  sum <- 0
  for(c in c("Ma.","Orwall.","Paola.","Peters.","Parrot.","Sood.","Avian.")){
    if(c == "Avian"){
      x <- sum(allTab[ r , projlist[ grepl(c, projlist) ] ] < 0.005, na.rm=TRUE)
    }else{
      x <- sum(!is.na(allTab[ r , projlist[ grepl(c, projlist) ] ]))
    }
    if(x > 0){ sum <- sum + 1 }
  }
  outmat[ r , "signif_studies"] <- sum 
 
  ri <- ri + 1 
}
   
outmat[ is.na(outmat) ] <- 0

dim(allTab)
dim(outmat)
allTab2 <- t(rbind(t(allTab),t(outmat)))
dim(allTab2)

write.table(allTab2, file="pvalue_direction_and_projectSum.20190209.tsv", sep="\t")
 
#colSums(is.na(allTab))

#save.image(file="20190208.QQQ.Rdata")
# load(file="20190208.QQQ.Rdata")
 
#n <- is.na(allTab[,"DIRECTION.Avian.MLSLW.deNovo"])
#nn <- n[n==TRUE]
#n <- is.na(allTab[,"Avian.MLSLW.deNovo"])
#nnn <- n[n==TRUE]





rr <- allTab2[ rev(order(allTab2[,"signif_studies"])), ]

pq2 <- pq
for( w in colnames(pq[ , grepl("Avian",colnames(pq)) ]) ){
  pq2[ pq2[,w] > 0.005 , w ] <- -1
}

pq3 <- pq2
for( w in colnames(pq[ , grepl("Peters",colnames(pq)) ]) ){
  pq3[ pq3[,w] > 0.005 , w ] <- -1
}

pq3[ pq3[,"Sood.P"] != -1 , "Sood.P"] <- 0.0001
pq3[ pq3[,"Orwall.boneloss"] != -1 , "Orwall.boneloss"] <- 0.0001
pq3[ pq3[,"Orwall.obesity"] != -1 , "Orwall.obesity"] <- 0.0001
pq3[ pq3[,"Orwall.mortality"] != -1 , "Orwall.mortality"] <- 0.0001
pq3[ pq3[,"Orwall.longevity"] != -1 , "Orwall.longevity"] <- 0.0001
pq3[pq3 == -1] <- -0.00995

png("signif_studies_500.png", width = 1000, height=2000)
heatmap.2( 
  
  pq2[ rownames(rr)[ 1:100 ] , ],
  
  trace='none',
  col=c(rep("black",50),redblue(50)),
  margins=c(22,10),
  keysize=0.7, key.par = list(cex=0.5),
  rowV = FALSE
  
)
dev.off()



# 1. Fill out the table with the studies and the information with them, importantly the number of individuals 
   # (or species) used and the criteria for significance

# 2. Make one heat map slide showing overall consistency with the genes across the studies (provide a legend for each slide)

# 3. Make one heat map slide showing overall consistency with genes *positively* associated across the studies, for studies 
   # where this is doable (include Sebastiani GWAS even though you can’t get direction of effect)

# 4. Make one heat map slide showing overall consistency with genes *negatively* associated across the studies

# 5. Take the genes with 1 study showing significance and run a pathway analysis

# 6. Take the genes with 2 studies showing significance and run a pathway analysis

# 7. Take the genes with 3 studies showing significance and run a pathway analysis


projectSum_heatmap <- matrix()



for(c in c("Ma.","Orwall.","Paola.","Peters.","Parrot.","Sood.","Avian.")){
  
  if(c == "NULL"){
    
    #for(c in colnames(pq3)[grepl(c,colnames(pq3))]){
    #  m <- rownames(pq3[ sum(pq3[ , grepl(c,colnames(pq3)) ] < 0.005 & pq3[ , grepl(c,colnames(pq3)) ] > 0) > 0 , ])
    #}
    
  }else{
    
    proj_syms <- list()
    for(d in colnames(pq3)[grepl(c,colnames(pq3))]){
      
      if(c == "Avian."){
        m <- rownames( pq3[  pq3[ , d ] > 0 & pq3[ , d ] < 0.005, ] )
      }else{
        m <- rownames( pq3[ pq3[ , d ] > 0, ] )
      }
      pqm <- pq3[ m , ]
      sorted_symbols <- names(sort(pqm[,d]))
      proj_syms <- append(proj_syms, sorted_symbols)
    }
    
    proj_syms <- unique(unlist(proj_syms))
    
    direction <- dirTab[ unlist(proj_syms) , grepl(c,colnames(dirTab)) ]
    #dirTab[ rownames(pq3) , grepl(c,colnames(dirTab)) ]
    if(!c %in% c("Paola.","Parrot.","Sood.")){
      direction <- direction[ unlist(proj_syms) , ]
    }#}else{
    #  direction <- direction[ unlist(proj_syms) , ]
    #}
    direction[ is.na(direction) ] <- 0
    if(c == "Peters."){
      direction[,1] <- as.numeric(direction[,1])
      direction[is.na(direction[,1]),1] <- 0
      direction[,2] <- as.numeric(direction[,2])
      direction[is.na(direction[,2]),2] <- 0
      direction[,3] <- as.numeric(direction[,3])
      direction[is.na(direction[,3]),3] <- 0
    }
    
    if(c == "Parrot."){
      direction <- rep(1,length(direction))
    }
    
    #direction[ direction=="" ] <- 0
    if(!c %in% c("Paola.","Parrot.","Sood.")){
      dirlist <- rep("None",dim(direction)[1])
    }else{
      dirlist <- rep("None",length(direction))
    }
    if(!c %in% c("Paola.","Parrot.","Sood.")){
      for(d in 1:dim(direction)[1]){
        negd <- sum(as.numeric(direction[d,]) < 0)
        posd <- sum(as.numeric(direction[d,]) > 0)
        if(negd > 0 & posd > 0){
          if(negd > posd){ dirlist[d] <- "Negative_Mixed" }
          if(negd < posd){ dirlist[d] <- "Positive_Mixed"  }
          if(negd == posd){ dirlist[d] <- "Mixed"  }
        }else{
          if(negd > 0){ dirlist[d] <- "Negative" }
          if(posd > 0){ dirlist[d] <- "Positive" }
        }
      }
    }else{
      for(d in 1:length(direction)){
        negd <- sum(as.numeric(direction[d]) < 0)
        posd <- sum(as.numeric(direction[d]) > 0)
        if(negd > 0 & posd > 0){
          if(negd > posd){ dirlist[d] <- "Negative_Mixed" }
          if(negd < posd){ dirlist[d] <- "Positive_Mixed"  }
          if(negd == posd){ dirlist[d] <- "Mixed"  }
        }else{
          if(negd > 0){ dirlist[d] <- "Negative" }
          if(posd > 0){ dirlist[d] <- "Positive" }
        }
      }
    }
    dircols <- dirlist
    dircols[dirlist == "Negative"] <- "brown4"
    dircols[dirlist == "Negative_Mixed"] <- "rosybrown3"
    dircols[dirlist == "Positive"] <- "darkolivegreen2"
    dircols[dirlist == "Positive_Mixed"] <- "green3"
    dircols[dirlist == "Mixed"] <- "gray"
    names(dircols)
    
    png(paste0(c,"all.overlap.png"), width = 1000, height=1500)
    heatmap.2( 
      pq3[ unlist(proj_syms) , c(1,3,5,2,4,6,7:24) ],
      trace='none',
      col=c(rep("black",50),redblue(50)),
      margins=c(22,10),
      keysize=0.7, key.par = list(cex=0.5),
      #Rowv = FALSE,
      Colv = FALSE,
      ColSideColors = c("slateblue4","slateblue4","slateblue4","slateblue3","slateblue3","slateblue3","slateblue2","slateblue2","slateblue1","slateblue1","slateblue1","slateblue1",
                        "darkgoldenrod3","darkgoldenrod3",
                        "springgreen4","springgreen4","springgreen4","springgreen4",
                        "navy",
                        "darksalmon","darksalmon","darksalmon",
                        "dodgerblue4",
                        "darkred"),
      RowSideColors = dircols,#[grepl("Negative",dirlist)],
      main=paste("\n\n",c,"\nAll Significant\n",length(unlist(proj_syms)),"significant symbols\nmax pval = ",max(pq3[ unlist(proj_syms) , c(13,14) ]))
    )
    legend("topright",
           pch=16,
           legend = c(paste0("(",  dim(pq3[ unlist(proj_syms) , c(1,3,5) ][ rowSums(pq3[ unlist(proj_syms) , c(1,3,5) ] > 0) > 0 , ])[1]  ,") Avian : Alignment"),
                      paste0("(",  dim(pq3[ unlist(proj_syms) , c(2,4,6) ][ rowSums(pq3[ unlist(proj_syms) , c(2,4,6) ] > 0) > 0 , ])[1]  ,") Avian : De novo"),
                      paste0("(",  dim(pq3[ unlist(proj_syms) , c(7,8) ][ rowSums(pq3[ unlist(proj_syms) , c(7,8) ] > 0) > 0 , ])[1]  ,") Avian : Drop-out"),
                      paste0("(",  dim(pq3[ unlist(proj_syms) , c(9,10,11,12) ][ rowSums(pq3[ unlist(proj_syms) , c(9,10,11,12) ] > 0) > 0 , ])[1]  ,") Avian : Extreme Longevity"),
                      paste0("(",  dim(pq3[ unlist(proj_syms) , c(13,14) ][ rowSums(pq3[ unlist(proj_syms) , c(13,14) ] > 0) > 0 , ])[1]  ,") Rodent"),
                      paste0("(",  dim(pq3[ unlist(proj_syms) , c(15,16,17,18) ][ rowSums(pq3[ unlist(proj_syms) , c(15,16,17,18) ] > 0) > 0 , ])[1]  ,") Human Muscle"),
                      paste0("(",  sum(pq3[ unlist(proj_syms) , c(19) ] > 0)  ,") Human SNP"),
                      paste0("(",  dim(pq3[ unlist(proj_syms) , c(20,21,22) ][ rowSums(pq3[ unlist(proj_syms) , c(20,21,22) ] > 0) > 0 , ])[1]  ,") Human Blood"),
                      paste0("(",  sum(pq3[ unlist(proj_syms) , c(23) ] > 0)  ,") Parrot Extreme Longevity"),
                      paste0("(",  sum(pq3[ unlist(proj_syms) , c(24) ] > 0)  ,") Human CHiP-seq")
           ),
           col = c("slateblue4","slateblue3","slateblue2","slateblue1","darkgoldenrod3","springgreen4","navy","darksalmon","dodgerblue4","darkred") ,
           cex = 0.8,
           pt.cex=2,
           #inset=c(0,-0.03)
    )
    
    legend("bottomleft",
           legend=c("Negative","Negative_Mixed","Mixed","Positive_Mixed","Positive"),
           col=c("brown4","rosybrown3","gray","green3","darkolivegreen2"),
           pch=16,
           cex = 0.8,
           pt.cex=2,
           title = "Count of Overlapping Significant Genes")
    dev.off()
    
    if(c != "Parrot."){
      png(paste0(c,"pos.overlap.png"), width = 1000, height=1500)
      heatmap.2( 
        pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(1,3,5,2,4,6,7:24) ],
        trace='none',
        col=c(rep("black",50),redblue(50)),
        margins=c(22,10),
        keysize=0.7, key.par = list(cex=0.5),
        #Rowv = FALSE,
        Colv = FALSE,
        ColSideColors = c("slateblue4","slateblue4","slateblue4","slateblue3","slateblue3","slateblue3","slateblue2","slateblue2","slateblue1","slateblue1","slateblue1","slateblue1",
                          "darkgoldenrod3","darkgoldenrod3",
                          "springgreen4","springgreen4","springgreen4","springgreen4",
                          "navy",
                          "darksalmon","darksalmon","darksalmon",
                          "dodgerblue4",
                          "darkred"),
        RowSideColors = dircols[grepl("Positive",dirlist)],
        #main=paste(c,"\nPositive"),#[grepl("Negative",dirlist)],
        main=paste("\n\n",c,"\nPositive Significant\n",length(unlist(proj_syms)[grepl("Positive",dirlist)] ),"significant symbols\nmax pval = ",max(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(13,14) ]))
      )
      legend("topright",
             pch=16,
             legend = c(paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(1,3,5) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(1,3,5) ] > 0) > 0 , ])[1]  ,") Avian : Alignment"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(2,4,6) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(2,4,6) ] > 0) > 0 , ])[1]  ,") Avian : De novo"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(7,8) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(7,8) ] > 0) > 0 , ])[1]  ,") Avian : Drop-out"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(9,10,11,12) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(9,10,11,12) ] > 0) > 0 , ])[1]  ,") Avian : Extreme Longevity"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(13,14) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(13,14) ] > 0) > 0 , ])[1]  ,") Rodent"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(15,16,17,18) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(15,16,17,18) ] > 0) > 0 , ])[1]  ,") Human Muscle"),
                        paste0("(",  sum(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(19) ] > 0)  ,") Human SNP"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(20,21,22) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(20,21,22) ] > 0) > 0 , ])[1]  ,") Human Blood"),
                        paste0("(",  sum(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(23) ] > 0)  ,") Parrot Extreme Longevity"),
                        paste0("(",  sum(pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(24) ] > 0)  ,") Human CHiP-seq")
             ),
             col = c("slateblue4","slateblue3","slateblue2","slateblue1","darkgoldenrod3","springgreen4","navy","darksalmon","dodgerblue4","darkred") ,
             cex = 0.8,
             pt.cex=2
             #inset=c(0,-0.03)
      )
      
      legend("bottomleft",
             legend=c("Negative","Negative_Mixed","Mixed","Positive_Mixed","Positive"),
             col=c("brown4","rosybrown3","gray","green3","darkolivegreen2"),
             pch=16,
             cex = 0.8,
             pt.cex=2)
      dev.off()
      
      
      png(paste0(c,"neg.overlap.png"), width = 1000, height= 1500)
      heatmap.2( 
        pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(1,3,5,2,4,6,7:24) ],
        trace='none',
        col=c(rep("black",50),redblue(50)),
        margins=c(22,10),
        keysize=0.7, key.par = list(cex=0.5),
        #Rowv = FALSE,
        Colv = FALSE,
        ColSideColors = c("slateblue4","slateblue4","slateblue4","slateblue3","slateblue3","slateblue3","slateblue2","slateblue2","slateblue1","slateblue1","slateblue1","slateblue1",
                          "darkgoldenrod3","darkgoldenrod3",
                          "springgreen4","springgreen4","springgreen4","springgreen4",
                          "navy",
                          "darksalmon","darksalmon","darksalmon",
                          "dodgerblue4",
                          "darkred"),
        RowSideColors = dircols[grepl("Negative",dirlist)],
        #main=paste(c,"\nPositive"),#[grepl("Negative",dirlist)],
        main=paste("\n\n",c,"\nNegative Significant\n",length(unlist(proj_syms)[grepl("Negative",dirlist)] ),"significant symbols\nmax pval = ",max(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(13,14) ]))
      )
      legend("topright",
             pch=16,
             legend = c(paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(1,3,5) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(1,3,5) ] > 0) > 0 , ])[1]  ,") Avian : Alignment"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(2,4,6) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(2,4,6) ] > 0) > 0 , ])[1]  ,") Avian : De novo"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(7,8) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(7,8) ] > 0) > 0 , ])[1]  ,") Avian : Drop-out"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(9,10,11,12) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(9,10,11,12) ] > 0) > 0 , ])[1]  ,") Avian : Extreme Longevity"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(13,14) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(13,14) ] > 0) > 0 , ])[1]  ,") Rodent"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(15,16,17,18) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(15,16,17,18) ] > 0) > 0 , ])[1]  ,") Human Muscle"),
                        paste0("(",  sum(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(19) ] > 0)  ,") Human SNP"),
                        paste0("(",  dim(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(20,21,22) ][ rowSums(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(20,21,22) ] > 0) > 0 , ])[1]  ,") Human Blood"),
                        paste0("(",  sum(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(23) ] > 0)  ,") Parrot Extreme Longevity"),
                        paste0("(",  sum(pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(24) ] > 0)  ,") Human CHiP-seq")
             ),
             col = c("slateblue4","slateblue3","slateblue2","slateblue1","darkgoldenrod3","springgreen4","navy","darksalmon","dodgerblue4","darkred") ,
             cex = 0.8,
             pt.cex=2
             #inset=c(0,-0.03)
      )
      
      legend("bottomleft",
             legend=c("Negative","Negative_Mixed","Mixed","Positive_Mixed","Positive"),
             col=c("brown4","rosybrown3","gray","green3","darkolivegreen2"),
             pch=16,
             cex = 0.8,
             pt.cex=2)
      dev.off()
    }
  }
  
}

# save.image(file=paste0("pvalTab.20190209.6pm.Rdata"))


# run for Ma. 




for(c in c("Ma.","Orwall.","Paola.","Peters.","Parrot.","Sood.","Avian.")){
  
  proj_syms <- list()
  for(d in colnames(pq3)[grepl(c,colnames(pq3))]){
    
    if(c == "Avian."){
      m <- rownames( pq3[  pq3[ , d ] > 0 & pq3[ , d ] < 0.005, ] )
    }else{
      m <- rownames( pq3[ pq3[ , d ] > 0, ] )
    }
    pqm <- pq3[ m , ]
    sorted_symbols <- names(sort(pqm[,d]))
    proj_syms <- append(proj_syms, sorted_symbols)
  }
  
  proj_syms <- unique(unlist(proj_syms))
  
  direction <- dirTab[ unlist(proj_syms) , grepl(c,colnames(dirTab)) ]
  #dirTab[ rownames(pq3) , grepl(c,colnames(dirTab)) ]
  if(!c %in% c("Paola.","Parrot.","Sood.")){
    direction <- direction[ unlist(proj_syms) , ]
  }#}else{
  #  direction <- direction[ unlist(proj_syms) , ]
  #}
  direction[ is.na(direction) ] <- 0
  if(c == "Peters."){
    direction[,1] <- as.numeric(direction[,1])
    direction[is.na(direction[,1]),1] <- 0
    direction[,2] <- as.numeric(direction[,2])
    direction[is.na(direction[,2]),2] <- 0
    direction[,3] <- as.numeric(direction[,3])
    direction[is.na(direction[,3]),3] <- 0
  }
  
  if(c == "Parrot."){
    direction <- rep(1,length(direction))
  }
  
  #direction[ direction=="" ] <- 0
  if(!c %in% c("Paola.","Parrot.","Sood.")){
    dirlist <- rep("None",dim(direction)[1])
  }else{
    dirlist <- rep("None",length(direction))
  }
  if(!c %in% c("Paola.","Parrot.","Sood.")){
    for(d in 1:dim(direction)[1]){
      negd <- sum(as.numeric(direction[d,]) < 0)
      posd <- sum(as.numeric(direction[d,]) > 0)
      if(negd > 0 & posd > 0){
        if(negd > posd){ dirlist[d] <- "Negative_Mixed" }
        if(negd < posd){ dirlist[d] <- "Positive_Mixed"  }
        if(negd == posd){ dirlist[d] <- "Mixed"  }
      }else{
        if(negd > 0){ dirlist[d] <- "Negative" }
        if(posd > 0){ dirlist[d] <- "Positive" }
      }
    }
  }else{
    for(d in 1:length(direction)){
      negd <- sum(as.numeric(direction[d]) < 0)
      posd <- sum(as.numeric(direction[d]) > 0)
      if(negd > 0 & posd > 0){
        if(negd > posd){ dirlist[d] <- "Negative_Mixed" }
        if(negd < posd){ dirlist[d] <- "Positive_Mixed"  }
        if(negd == posd){ dirlist[d] <- "Mixed"  }
      }else{
        if(negd > 0){ dirlist[d] <- "Negative" }
        if(posd > 0){ dirlist[d] <- "Positive" }
      }
    }
  }
  dircols <- dirlist
  dircols[dirlist == "Negative"] <- "brown4"
  dircols[dirlist == "Negative_Mixed"] <- "rosybrown3"
  dircols[dirlist == "Positive"] <- "darkolivegreen2"
  dircols[dirlist == "Positive_Mixed"] <- "green3"
  dircols[dirlist == "Mixed"] <- "gray"
  names(dircols)
  
  q1_all <- rownames( pq3[ unlist(proj_syms) , c(1,3,5,2,4,6,7:24) ])
  
  z <- pq3[ unlist(proj_syms) , c(1,3,5,2,4,6,7:24) ]
  z <- z > 0
  z <- z[ , !grepl(c,colnames(z)) ]
  z <- z[rowSums(z) > 0,]
  q1_shared <- rownames(z)
  
  q2_all <- rownames( pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(1,3,5,2,4,6,7:24) ] )
  z <- pq3[ unlist(proj_syms)[grepl("Positive",dirlist)] , c(1,3,5,2,4,6,7:24) ]
  z <- z > 0
  z <- z[ , !grepl(c,colnames(z)) ]
  z <- z[rowSums(z) > 0,]
  q2_shared <- rownames(z)
  
  if(c != "Parrot."){
    q3_all <- rownames( pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(1,3,5,2,4,6,7:24) ] )
    z <- pq3[ unlist(proj_syms)[grepl("Negative",dirlist)] , c(1,3,5,2,4,6,7:24) ]
    z <- z > 0
    z <- z[ , !grepl(c,colnames(z)) ]
    z <- z[rowSums(z) > 0,]
    q3_shared <- rownames(z)
  }
  
  outz <- matrix(nrow=max(c(length(q1_all),length(q1_shared),length(q2_all),length(q2_shared),length(q3_all),length(q3_shared))) , ncol=6)
  outz[1:length(q1_all),   1] <- q1_all
  outz[1:length(q1_shared),2] <- q1_shared
  outz[1:length(q2_all),   3] <- q2_all
  outz[1:length(q2_shared),4] <- q2_shared
  if(c != "Parrot."){
    outz[1:length(q3_all),   5] <- q3_all
    outz[1:length(q3_shared),6] <- q3_shared
  }
  
  write.table(outz,file=paste0(c,".sharedlist.tsv"),sep = "\t")
  
  if(c == "Ma."){
    outz.ma <- outz
  }else if(c =="Orwall."){
    outz.orwall <- outz
  }else if(c == "Paola."){
    outz.paola <- outz
  }else if(c == "Peters."){
    outz.peters <- outz
  }else if(c == "Parrot."){
    outz.parrot <- outz
  }else if(c == "Sood."){
    outz.sood <- outz
  }else if(c == "Avian."){
    outz.avian <- outz
  }
} 
 

distMat <- matrix(nrow=7,ncol=7)
colnames(distMat) <- rownames(distMat) <- c("Ma.","Orwall.","Paola.","Peters.","Parrot.","Sood.","Avian.")
stringMat <-  matrix(nrow=7,ncol=7)
colnames(stringMat) <- rownames(stringMat) <- c("Ma.","Orwall.","Paola.","Peters.","Parrot.","Sood.","Avian.")

m <- unique(unlist(as.list(unlist(outz.ma))))
m <- m[!is.na(m)]
m <- m[m != "N/A"]
distMat["Ma.","Ma."] <- length( m )
distMat["Ma.","Orwall."] <- sum( m %in% outz.orwall )
distMat["Ma.","Paola."] <- sum( m %in% outz.paola )
distMat["Ma.","Peters."] <- sum( m %in% outz.peters )
distMat["Ma.","Parrot."] <- sum( m %in% outz.parrot )
distMat["Ma.","Sood."] <- sum( m %in% outz.sood )
distMat["Ma.","Avian."] <- sum( m %in% outz.avian )

stringMat["Ma.","Ma."]      <- paste(m,collapse=",")
stringMat["Ma.","Orwall."]  <- paste(m[m %in% outz.orwall],collapse=",")
stringMat["Ma.","Paola."]   <- paste(m[m %in% outz.paola],collapse=",")
stringMat["Ma.","Peters."]  <- paste(m[m %in% outz.peters],collapse=",")
stringMat["Ma.","Parrot."]  <- paste(m[m %in% outz.parrot],collapse=",")
stringMat["Ma.","Sood."]    <- paste(m[m %in% outz.sood],collapse=",")
stringMat["Ma.","Avian."]   <- paste(m[m %in% outz.avian],collapse=",")
 
allm <- unique(c( m[m %in% outz.orwall] , m[m %in% outz.paola] , m[m %in% outz.peters] , m[m %in% outz.parrot] , m[m %in% outz.sood] , m[m %in% outz.avian]) )
w <- matrix(nrow=length(allm),ncol=6)
rownames(w) <- allm
colnames(w) <- c("Orwall.","Paola.","Peters.","Parrot.","Sood.","Avian.") #"Ma.",
for(x in allm){
  w[x,1:6] <- c( sum( x %in% outz.orwall ) , sum( x %in% outz.paola ) , sum( x %in% outz.peters ) , sum( x %in% outz.parrot ) , sum( x %in% outz.sood ) , sum( x %in% outz.avian ) )
}
png(filename = "ma.shared_signif.png",height=2800,width=800)
heatmap.2(w, trace='none', col=c("black","red"), margins=c(10,10), keysize=0.7, key.par = list(cex=0.5))
dev.off()
w <- as.data.frame(w)
w$sum <- rowSums(w)
write.table(w[ rev(order(w$sum)) , ], file="ma.shared_signif.tsv",sep="\t") 

m <- unique(unlist(as.list(unlist(outz.orwall))))
m <- m[!is.na(m)]
m <- m[m != "N/A"]
distMat["Orwall.","Orwall."] <- length( m )
distMat["Orwall.","Paola."] <- sum( m %in% outz.paola )
distMat["Orwall.","Peters."] <- sum( m %in% outz.peters )
distMat["Orwall.","Parrot."] <- sum( m %in% outz.parrot )
distMat["Orwall.","Sood."] <- sum( m %in% outz.sood )
distMat["Orwall.","Avian."] <- sum( m %in% outz.avian )

stringMat["Orwall.","Orwall."]  <- paste(m,collapse=",")
stringMat["Orwall.","Paola."]   <- paste(m[m %in% outz.paola],collapse=",")
stringMat["Orwall.","Peters."]  <- paste(m[m %in% outz.peters],collapse=",")
stringMat["Orwall.","Parrot."]  <- paste(m[m %in% outz.parrot],collapse=",")
stringMat["Orwall.","Sood."]    <- paste(m[m %in% outz.sood],collapse=",")
stringMat["Orwall.","Avian."]   <- paste(m[m %in% outz.avian],collapse=",")

allm <- unique(c( m[m %in% outz.ma] , m[m %in% outz.paola] , m[m %in% outz.peters] , m[m %in% outz.parrot] , m[m %in% outz.sood] , m[m %in% outz.avian]) )
w <- matrix(nrow=length(allm),ncol=6)
rownames(w) <- allm
colnames(w) <- c("Ma.","Paola.","Peters.","Parrot.","Sood.","Avian.") #"Ma.",
for(x in allm){
  w[x,1:6] <- c( sum( x %in% outz.ma ) , sum( x %in% outz.paola ) , sum( x %in% outz.peters ) , sum( x %in% outz.parrot ) , sum( x %in% outz.sood ) , sum( x %in% outz.avian ) )
}
png(filename = "orwall.shared_signif.png",height=800,width=800)
heatmap.2(w, trace='none', col=c("black","red"), margins=c(10,10), keysize=0.7, key.par = list(cex=0.5))
dev.off()
w <- as.data.frame(w)
w$sum <- rowSums(w)
write.table(w[ rev(order(w$sum)) , ], file="orwall.shared_signif.tsv",sep="\t") 

m <- unique(unlist(as.list(unlist(outz.paola))))
m <- m[!is.na(m)]
m <- m[m != "N/A"]
distMat["Paola.","Paola."] <- length( unique(unlist(as.list(unlist(outz.paola)))) )
distMat["Paola.","Peters."] <- sum( unique(unlist(as.list(unlist(outz.paola)))) %in% outz.peters )
distMat["Paola.","Parrot."] <- sum( unique(unlist(as.list(unlist(outz.paola)))) %in% outz.parrot )
distMat["Paola.","Sood."] <- sum( unique(unlist(as.list(unlist(outz.paola)))) %in% outz.sood )
distMat["Paola.","Avian."] <- sum( unique(unlist(as.list(unlist(outz.paola)))) %in% outz.avian )

stringMat["Paola.","Paola."]   <- paste(m,collapse=",")
stringMat["Paola.","Peters."]  <- paste(m[m %in% outz.peters],collapse=",")
stringMat["Paola.","Parrot."]  <- paste(m[m %in% outz.parrot],collapse=",")
stringMat["Paola.","Sood."]    <- paste(m[m %in% outz.sood],collapse=",")
stringMat["Paola.","Avian."]   <- paste(m[m %in% outz.avian],collapse=",")

allm <- unique(c( m[m %in% outz.ma] , m[m %in% outz.orwall] , m[m %in% outz.peters] , m[m %in% outz.parrot] , m[m %in% outz.sood] , m[m %in% outz.avian]) )
w <- matrix(nrow=length(allm),ncol=6)
rownames(w) <- allm
colnames(w) <- c("Ma.","Orwall.","Peters.","Parrot.","Sood.","Avian.") #"Ma.",
for(x in allm){
  w[x,1:6] <- c( sum( x %in% outz.ma ) , sum( x %in% outz.orwall ) , sum( x %in% outz.peters ) , sum( x %in% outz.parrot ) , sum( x %in% outz.sood ) , sum( x %in% outz.avian ) )
}
png(filename = "paola.shared_signif.png",height=800,width=800)
heatmap.2(w, trace='none', col=c("black","red"), margins=c(10,10), keysize=0.7, key.par = list(cex=0.5))
dev.off()
w <- as.data.frame(w)
w$sum <- rowSums(w)
write.table(w[ rev(order(w$sum)) , ], file="paola.shared_signif.tsv",sep="\t") 

m <- unique(unlist(as.list(unlist(outz.peters))))
m <- m[!is.na(m)]
m <- m[m != "N/A"]
distMat["Peters.","Peters."] <- length( unique(unlist(as.list(unlist(outz.peters)))) )
distMat["Peters.","Parrot."] <- sum( unique(unlist(as.list(unlist(outz.peters)))) %in% outz.parrot )
distMat["Peters.","Sood."] <- sum( unique(unlist(as.list(unlist(outz.peters)))) %in% outz.sood )
distMat["Peters.","Avian."] <- sum( unique(unlist(as.list(unlist(outz.peters)))) %in% outz.avian )

stringMat["Peters.","Peters."]  <- paste(m,collapse=",")
stringMat["Peters.","Parrot."]  <- paste(m[m %in% outz.parrot],collapse=",")
stringMat["Peters.","Sood."]    <- paste(m[m %in% outz.sood],collapse=",")
stringMat["Peters.","Avian."]   <- paste(m[m %in% outz.avian],collapse=",")

allm <- unique(c( m[m %in% outz.ma] , m[m %in% outz.orwall] , m[m %in% outz.paola] , m[m %in% outz.parrot] , m[m %in% outz.sood] , m[m %in% outz.avian]) )
w <- matrix(nrow=length(allm),ncol=6)
rownames(w) <- allm
colnames(w) <- c("Ma.","Orwall.","Paola.","Parrot.","Sood.","Avian.") #"Ma.",
for(x in allm){
  w[x,1:6] <- c( sum( x %in% outz.ma ) , sum( x %in% outz.orwall ) , sum( x %in% outz.paola ) , sum( x %in% outz.parrot ) , sum( x %in% outz.sood ) , sum( x %in% outz.avian ) )
}
png(filename = "peters.shared_signif.png",height=2800,width=800)
heatmap.2(w, trace='none', col=c("black","red"), margins=c(10,10), keysize=0.7, key.par = list(cex=0.5))
dev.off()
w <- as.data.frame(w)
w$sum <- rowSums(w)
write.table(w[ rev(order(w$sum)) , ], file="peters.shared_signif.tsv",sep="\t") 

m <- unique(unlist(as.list(unlist(outz.parrot))))
m <- m[!is.na(m)]
m <- m[m != "N/A"]
distMat["Parrot.","Parrot."] <- length( unique(unlist(as.list(unlist(outz.parrot)))) )
distMat["Parrot.","Sood."] <- sum( unique(unlist(as.list(unlist(outz.parrot)))) %in% outz.sood )
distMat["Parrot.","Avian."] <- sum( unique(unlist(as.list(unlist(outz.parrot)))) %in% outz.avian )

stringMat["Parrot.","Parrot."]  <- paste(m,collapse=",")
stringMat["Parrot.","Sood."]    <- paste(m[m %in% outz.sood],collapse=",")
stringMat["Parrot.","Avian."]   <- paste(m[m %in% outz.avian],collapse=",")

allm <- unique(c( m[m %in% outz.ma] , m[m %in% outz.orwall] , m[m %in% outz.paola] , m[m %in% outz.peters] , m[m %in% outz.sood] , m[m %in% outz.avian]) )
w <- matrix(nrow=length(allm),ncol=6)
rownames(w) <- allm
colnames(w) <- c("Ma.","Orwall.","Paola.","Peters.","Sood.","Avian.") #"Ma.",
for(x in allm){
  w[x,1:6] <- c( sum( x %in% outz.ma ) , sum( x %in% outz.orwall ) , sum( x %in% outz.paola ) , sum( x %in% outz.peters ) , sum( x %in% outz.sood ) , sum( x %in% outz.avian ) )
}
png(filename = "parrot.shared_signif.png",height=800,width=800)
heatmap.2(w, trace='none', col=c("black","red"), margins=c(10,10), keysize=0.7, key.par = list(cex=0.5))
dev.off()
w <- as.data.frame(w)
w$sum <- rowSums(w)
write.table(w[ rev(order(w$sum)) , ], file="parrot.shared_signif.tsv",sep="\t") 
 
m <- unique(unlist(as.list(unlist(outz.sood))))
m <- m[!is.na(m)]
m <- m[m != "N/A"]
distMat["Sood.","Sood."] <- length( unique(unlist(as.list(unlist(outz.sood)))) )
distMat["Sood.","Avian."] <- sum( unique(unlist(as.list(unlist(outz.sood)))) %in% outz.avian )

stringMat["Sood.","Sood."]    <- paste(m,collapse=",")
stringMat["Sood.","Avian."]   <- paste(m[m %in% outz.avian],collapse=",")

allm <- unique(c( m[m %in% outz.ma] , m[m %in% outz.orwall] , m[m %in% outz.paola] , m[m %in% outz.peters] , m[m %in% outz.parrot] , m[m %in% outz.avian]) )
w <- matrix(nrow=length(allm),ncol=6)
rownames(w) <- allm
colnames(w) <- c("Ma.","Orwall.","Paola.","Peters.","Parrot.","Avian.") #"Ma.",
for(x in allm){
  w[x,1:6] <- c( sum( x %in% outz.ma ) , sum( x %in% outz.orwall ) , sum( x %in% outz.paola ) , sum( x %in% outz.peters ) , sum( x %in% outz.parrot ) , sum( x %in% outz.avian ) )
}
png(filename = "sood.shared_signif.png",height=800,width=800)
heatmap.2(w, trace='none', col=c("black","red"), margins=c(10,10), keysize=0.7, key.par = list(cex=0.5))
dev.off()
w <- as.data.frame(w)
w$sum <- rowSums(w)
write.table(w[ rev(order(w$sum)) , ], file="sood.shared_signif.tsv",sep="\t") 

distMat["Avian.","Avian."]    <- length(unique(unlist(as.list(unlist(outz.avian)))))
stringMat["Avian.","Avian."]   <- paste(m[m %in% outz.avian],collapse=",")

m <- unique(unlist(as.list(unlist(outz.avian))))
m <- m[!is.na(m)]
m <- m[m != "N/A"]
allm <- unique(c( m[m %in% outz.ma] , m[m %in% outz.orwall] , m[m %in% outz.paola] , m[m %in% outz.peters] , m[m %in% outz.parrot] , m[m %in% outz.sood]) )
w <- matrix(nrow=length(allm),ncol=6)
rownames(w) <- allm
colnames(w) <- c("Ma.","Orwall.","Paola.","Peters.","Parrot.","Sood.") #"Ma.",
for(x in allm){
  w[x,1:6] <- c( sum( x %in% outz.ma ) , sum( x %in% outz.orwall ) , sum( x %in% outz.paola ) , sum( x %in% outz.peters ) , sum( x %in% outz.parrot ) , sum( x %in% outz.sood ) )
}
png(filename = "avian.shared_signif.png",height=2800,width=800)
heatmap.2(w, trace='none', col=c("black","red"), margins=c(10,10), keysize=0.7, key.par = list(cex=0.5))
dev.off()
w <- as.data.frame(w)
w$sum <- rowSums(w)
write.table(w[ rev(order(w$sum)) , ], file="avian.shared_signif.tsv",sep="\t") 





png(filename = "distMat.png", height=800, width=800)
distMat[is.na(distMat)] <- -1
distMat.l <- log(distMat)
distMat.l[ 2, 5:6 ] <- 0
#distMat[d]
heatmap.2(distMat.l,
          trace='none',
          col=c(redgreen(50)),
          margins=c(10,10),
          keysize=0.7, key.par = list(cex=0.5),
          #scale = "column",
          Rowv = FALSE, Colv = FALSE,
          cellnote=distMat,
          notecol="white",
          notecex = 1.5,
          main="Overlapping Signficant Gene Count"
) 
dev.off()
  
write.table( stringMat , file = "stringMat.tsv", sep = "\t" )

          #rowV = FALSE)


m <- unlist(unique(unlist(as.list(unlist(outz.ma)))))
m <- m[ m %in% outz.orwall ]
m <- m[ !is.na(m) ]
 
distMat["Ma.","Paola."] <- sum( unique(unlist(as.list(unlist(outz.ma)))) %in% outz.paola )
distMat["Ma.","Peters."] <- sum( unique(unlist(as.list(unlist(outz.ma)))) %in% outz.peters )
distMat["Ma.","Parrot."] <- sum( unique(unlist(as.list(unlist(outz.ma)))) %in% outz.parrot )
distMat["Ma.","Sood."] <- sum( unique(unlist(as.list(unlist(outz.ma)))) %in% outz.sood )
distMat["Ma.","Avian."] <- sum( unique(unlist(as.list(unlist(outz.ma)))) %in% outz.avian )

















if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO", version = "3.8")

