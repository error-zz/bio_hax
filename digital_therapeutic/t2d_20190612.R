#!/usr/bin/env Rscript
#install.packages("beanplot",repos='http://cran.us.r-project.org')
library(gplots)
library(colorRamps)
library(ape)
library(vegan)
library(caper)
library(plot3D)
#library(rgl)
library(heatmap.plus)
library(beanplot)
library(edgeR)

setwd ("~/Desktop/exe.Marcelo/qq3")

#load previous environment?
#load("20181012star.Rdata")

############################################################################################################################################################
############################################################################################################################################################
# IMPORT #
############################################################################################################################################################
############################################################################################################################################################

RvE1 <- read.table("RvE1.out.annotated.txt",                 sep="\t", stringsAsFactors = FALSE, header=TRUE)
RvE1 <- RvE1[ , c( "Ensembl_GeneID","Gene_name","Biotype", colnames(RvE1)[grep('norm', colnames(RvE1))] ) ]
rownames(RvE1) <- RvE1[,1]
 
HvsT2D <- read.table("HvsT2D.all.out.annotated.txt", sep="\t", header=TRUE, fill=TRUE)
HvsT2D <- HvsT2D[ , c( "Ensembl_GeneID","Gene_name","Biotype", colnames(HvsT2D)[grep('norm', colnames(HvsT2D))] ) ]
rownames(HvsT2D) <- HvsT2D[,1]
 
#metadata_corrected <- read.table("metadata_corrected_020719_SEK.tab1.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
metadata_corrected <- read.table("metadata_corr_020719_sort.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
#metadata_corrected <- metadata_corrected[-1,]

T2DvsHealthy.M.W.corrected          <- read.table("T2DvsHealthy.M.W.corrected.txt", sep="\t", header=TRUE, fill=TRUE)
rownames(T2DvsHealthy.M.W.corrected) <- T2DvsHealthy.M.W.corrected[,1]
sortnames.pval.T2DvsHealthy.M.W.corrected <- rownames(T2DvsHealthy.M.W.corrected[ order(T2DvsHealthy.M.W.corrected$PValue) , ])
sortnames.fdr.T2DvsHealthy.M.W.corrected  <- rownames(T2DvsHealthy.M.W.corrected[ order(T2DvsHealthy.M.W.corrected$FDR) , ])
   
#T2DvsHealthy.M.W.corrected.selected <- read.table("T2DvsHealthy.M.W.corrected.selected.txt", sep="\t", header=TRUE, fill=TRUE)
#rownames(T2DvsHealthy.M.W.corrected.selected) <- T2DvsHealthy.M.W.corrected.selected[,1]
#sortnames.pval.T2DvsHealthy.M.W.corrected.selected <- rownames(T2DvsHealthy.M.W.corrected.selected[ order(T2DvsHealthy.M.W.corrected.selected$PValue) , ])
#sortnames.fdr.T2DvsHealthy.M.W.corrected.selected  <- rownames(T2DvsHealthy.M.W.corrected.selected[ order(T2DvsHealthy.M.W.corrected.selected$FDR) , ])

perturbation                        <- read.table("perturbation.txt", sep="\t", header=TRUE, fill=TRUE)
 



############################################################################################################################################################
############################################################################################################################################################
# 1) TEST YANMEI METHODS

# dim(HvsT2D)
# 57992    25
count_matrix <- HvsT2D[ , colnames(HvsT2D) != "MF.5_Healthy_F_35_S_W_norm"]
# count_matrix <- HvsT2D[ , colnames(HvsT2D) != "MF.7_Healthy_F_35_S_W_norm"]
count_matrix <- count_matrix[ , colnames(count_matrix) != "MF.7_Healthy_F_35_S_W_norm"]

count_matrix <- count_matrix[ , colnames(count_matrix) != "MF.1_Healthy_M_60_NS_W_norm"]
count_matrix <- count_matrix[ , colnames(count_matrix) != "MF.9_Healthy_F_28_N_B_norm"]
count_matrix <- count_matrix[,4:dim(count_matrix)[2]]
#beanplot(count_matrix[ , colnames(count_matrix) == "MF.1_Healthy_M_60_NS_W_norm"], HvsT2D[ , colnames(count_matrix) == "MF.48_Healthy_M_59_NS_W_norm"], 
#         names = c("MF.1_Healthy_M_60_NS_W_norm","MF.48_Healthy_M_59_NS_W_norm"))

rel_test <- t(decostand(t(count_matrix), method="total"))
 

count_matrix <- count_matrix[ rowSums(count_matrix[,1:18] >= 4) >= 4 , ]
# dim(count_matrix)
# 17792    24

#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#install.packages('edgeR', dependencies=TRUE, repos='http://cran.rstudio.com/')

#count_cmpr <- count_matrix[,4:24]
count_cmpr <- count_matrix
count_grps <- c("Healthy","T2D")[ as.factor(grepl("T2D",colnames(count_cmpr))) ]
d <- DGEList(counts=count_cmpr,group=factor(count_grps))

d <- calcNormFactors(d)

png(filename = "20190311.bcv.dist.png", width = 1000, height=1000)
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", rev(as.character(unique(d$samples$group))), col=1:3, pch=20)
dev.off()

# fishers exact test

d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)
names(d1)

png("20190311.BCVplot.png", width = 1000, height=1000)
plotBCV(d1)
title("A decreasing estimate for the coefficient of variation vs. CPM is a good sign!")
dev.off()

png("20190311.FishersLMFit.png", width = 1000, height=1000)
et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
topTags(et12, n=10)
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)
de1tags12 <- rownames(d1)[as.logical(de1)] 
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
dev.off()

write.table(x = et12$table, file = "20190311.fishers_pval.tsv", sep = "\t")
 
# switch to glm

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat, method="deviance",robust=TRUE)
d2 <- estimateGLMTrendedDisp(d2,design.mat)
d2 <- estimateGLMTagwiseDisp(d2,design.mat)

png("20190311.BCVplot2.png", width = 1000, height=1000)
plotBCV(d2)
title("A decreasing estimate for the coefficient of variation vs. CPM is a good sign!")
dev.off()

fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast=c(1,-1))
tab <- lrt12$table
tab$fdr <- p.adjust(lrt12$table$PValue,"fdr")
tab$GeneSymbol <- HvsT2D[ rownames(tab) , 2 ]
write.table(x = tab, file = "20190311.glm_pval.fdr.tsv", sep = "\t")
 

png("20190311.GLMFit.png", width = 1000, height=1000)
de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05)
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")
dev.off()

###############################  redo with rel distance
###############################  redo with rel distance
###############################  redo with rel distance
###############################  redo with rel distance
###############################  redo with rel distance


# dim(count_matrix)
# 17792    24

#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#library(edgeR)

#count_cmpr <- count_matrix[,4:24]
rel_matrix <- t(decostand(t(count_cmpr), method="total"))
count_grps <- c("Healthy","T2D")[ as.factor(grepl("T2D",colnames(count_cmpr))) ]
d <- DGEList(counts=count_cmpr,group=factor(count_grps))

d <- calcNormFactors(d)

png(filename = "rel.bcv.dist.png", width = 1000, height=1000)
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", rev(as.character(unique(d$samples$group))), col=1:3, pch=20)
dev.off()

# fishers exact test

d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)
names(d1)

png("rel.BCVplot.png", width = 1000, height=1000)
plotBCV(d1)
title("A decreasing estimate for the coefficient of variation vs. CPM is a good sign!")
dev.off()

png("rel.FishersLMFit.png", width = 1000, height=1000)
et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
topTags(et12, n=10)
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)
de1tags12 <- rownames(d1)[as.logical(de1)] 
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
dev.off()

write.table(x = et12$table, file = "rel.fishers_pval.tsv", sep = "\t")

# switch to glm

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat, method="deviance",robust=TRUE)
d2 <- estimateGLMTrendedDisp(d2,design.mat)
d2 <- estimateGLMTagwiseDisp(d2,design.mat)

png("rel.BCVplot2.png", width = 1000, height=1000)
plotBCV(d2)
title("A decreasing estimate for the coefficient of variation vs. CPM is a good sign!")
dev.off()

fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast=c(1,-1))
tab <- lrt12$table
tab$fdr <- p.adjust(lrt12$table$PValue,"fdr")
tab$GeneSymbol <- HvsT2D[ rownames(tab) , 2 ]
write.table(x = tab, file = "rel.glm_pval.fdr.tsv", sep = "\t")


png("rel.GLMFit.png", width = 1000, height=1000)
de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05)
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")



save.image("evaloQ.Rdata")


#  FIX DIRECTIONALITY ON THE OBSERVED GLM ASSOCIATIONS WITH T2D

library(beanplot)

png("NOduplicates.genesOfinterest.x13.beanplot.png", width=1350,height=1700)
plot.new()
par(mfrow=c(5,3))
#for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
for(g in c("NECTIN2","FAM3B","PLPP3","KCNH4","VWF","C3","SLC9A4","RBPMS2","ENPP7P8","MX1","C1orf61","GNPDA1","IRAK2","LILRB5","AKR1C1")){
  geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
  print(paste(g,geneID))
  beanplot( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))],# [rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] != 0],
            rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))],
            col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
            names = c("Healthy","T2D"),
            bw="nrd0",
            main=g
  )
}
dev.off()


par(mfrow=c(1,1))
#for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
for(g in c("NECTIN2","FAM3B","PLPP3","KCNH4","VWF","C3","SLC9A4","RBPMS2","ENPP7P8","MX1","C1orf61","GNPDA1","IRAK2","LILRB5","AKR1C1")){
  
  png(paste0("NOduplicates.genesOfinterest.",g,".beanplot.png"), width=1600,height=800)
  #par(mar=c(3,8,3,8))
  par(mar=c(4,4,4,4))
  plot.new()
  
  geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
  print(paste(g,geneID))
  beanplot( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))],# [rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] != 0],
            rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))],
            col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
            names = c("Healthy","T2D"),
            bw="nrd0",
            main=g,
            xlim = c(0.2,2.8),
            cex.axis=2,
            cex.main=2,
            cex.names=2
  )
  text( x=0.3, y=unname(rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] ), names(rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] ), cex=0.6)
  text( x=2.7, y=unname(rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] ),     names(rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] ),     cex=0.6)
  dev.off()
  
  
  png(paste0("NOduplicates.genesOfinterest.",g,".beanplotNOtext.png"), width=800,height=800)
  par(mar=c(4,4,4,4))
  beanplot( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))],# [rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] != 0],
            rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))],
            col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
            names = c("Healthy","T2D"),
            bw="nrd0",
            main=g,
            cex.axis=2,
            cex.main=2,
            cex.names=2)
  dev.off()
  
}

#dev.off()



# 
# png("NOduplicates.side.genesOfinterest.x13.beanplot.png", width=800,height=1700)
# plot.new()
# par(mfrow=c(5,3))
# #for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
# for(g in c("NECTIN2","FAM3B","PLPP3","KCNH4","VWF","C3","SLC9A4","RBPMS2","ENPP7P8","MX1","C1orf61","GNPDA1","IRAK2","LILRB5","AKR1C1")){
#   geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
#   print(paste(g,geneID))
#   beanplot( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))],# [rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] != 0],
#             rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))],
#             col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
#             names = c("Healthy","T2D"),
#             bw="nrd0",
#             side="both",
#             main=g
#   )
# }
# dev.off()

# png("NOduplicates.LOG.genesOfinterest.x13.beanplot.png", width=800,height=1700)
# plot.new()
# par(mfrow=c(5,3))
# #for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
# for(g in c("NECTIN2","FAM3B","PLPP3","KCNH4","VWF","C3","SLC9A4","RBPMS2","ENPP7P8","MX1","C1orf61","GNPDA1","IRAK2","LILRB5","AKR1C1")){
#   geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
#   print(paste(g,geneID))
#   beanplot( log(rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))]+0.0000001),# [rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] != 0],
#             log(rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))]+0.0000001),
#             col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
#             names = c("Healthy","T2D"),
#             bw="nrd0",
#             main=g
#   )
# }
# dev.off()
# # 


par(mfrow=c(1,1))
#for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
for(g in c("NECTIN2","FAM3B","PLPP3","KCNH4","VWF","C3","SLC9A4","RBPMS2","ENPP7P8","MX1","C1orf61","GNPDA1","IRAK2","LILRB5","AKR1C1")){
  
  png(paste0("NOduplicates.genesOfinterest.",g,".LOGbeanplot.png"), width=1600,height=800)
  #par(mar=c(3,8,3,8))
  par(mar=c(4,4,4,4))
  plot.new()
  
  geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
  print(paste(g,geneID))
  beanplot( log( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] + 0.0000001 ),# [rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] != 0],
            log( rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] + 0.0000001 ),
            col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
            names = c("Healthy","T2D"),
            bw="nrd0",
            main=g,
            xlim = c(0.2,2.8),
            cex.axis=2,
            cex.main=2,
            cex.names=2
  )
  text( x=0.3, y= log( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] + 0.0000001 ), names(rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] ), cex=0.6)
  text( x=2.7, y=log( rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] + 0.0000001 )     ,     names(rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] ),     cex=0.6)
  dev.off()
  
  
  png(paste0("NOduplicates.genesOfinterest.",g,".LOGbeanplotNOtext.png"), width=800,height=800)
  par(mar=c(4,4,4,4))
  beanplot( log( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] + 0.0000001 ),
            log( rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] + 0.0000001 ),
            col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
            names = c("Healthy","T2D"),
            bw="nrd0",
            main=g,
            cex.axis=2,
            cex.main=2,
            cex.names=2)
  dev.off()
  
}

# png("NOduplicates.sideLOG.genesOfinterest.x13.beanplot.png", width=800,height=1700)
# plot.new()
# par(mfrow=c(5,3))
# #for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
# for(g in c("NECTIN2","FAM3B","PLPP3","KCNH4","VWF","C3","SLC9A4","RBPMS2","ENPP7P8","MX1","C1orf61","GNPDA1","IRAK2","LILRB5","AKR1C1")){
#   geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
#   print(paste(g,geneID))
#   beanplot( log(rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))]+0.0000001),# [rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] != 0],
#             log(rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))]+0.0000001),
#             col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
#             names = c("Healthy","T2D"),
#             bw="nrd0",
#             side="both",
#             main=g
#   )
# }
# dev.off()


save.image("eval0.Rdata")
load("eval0.Rata")
 
 
# exampleIDs <- c("BAX","BCL2","CASP3")

# heatmap.2( -log(viz.relabund)[ exampleIDs, ] )

# fdr[ fdr$Gene_name %in% exampleIDs , 1]

#   
#   nuNam , 
#                                nuOrd 
#                                ]),            
#            labCol = paste(metadata_corrected[metadata_corrected$Sample != 5,]$Treatment,metadata_corrected[metadata_corrected$Sample != 5,]$ID,metadata_corrected[metadata_corrected$Sample != 5,]$Sample,sep="\t"),
#            col=c(redblue(150)),           trace="none",           margins = c(45,15),           keysize=0.4, key.par = list(cex=0.5),           
#            labRow = HvsT2D[nuNam,]$Gene_name,
#            ColSideColors = c("darkgreen","maroon")[ as.factor(metadata_corrected[metadata_corrected$Sample != 5,]$Disease) ], 
#            colsep = c(c(1:3),6,10,11,c(15:19),c(21:25),29,33,c(37:43)), sepwidth = 0.2,
#            cexRow = 1.4, cexCol = 1.8,
#            RowSideColors = c("red","green")[as.numeric(as.factor(rowMeans( viz.relabund[ nuNam , colnames(viz.relabund)[ grepl("T2D",colnames(viz.relabund)) ]  ] ) - rowMeans( viz.relabund[ nuNam , colnames(viz.relabund)[ !grepl("T2D",colnames(viz.relabund)) ]  ] ) > 0))],           
#            Colv = FALSE
# )
# 
# 
# 
# 
# 
# 
# 




































































#old code
 
#et13 <- exactTest(d2, pair=c(1,3)) # compare groups 1 and 3
#et23 <- exactTest(d2, pair=c(2,3))

############################################################################################################################################################
############################################################################################################################################################
# 1) MOST SIGNIF HEATMAP #
# 
# 1.
# For the Healthy vs Diabetes analysis (T2DvsHealthy.M.W.corrected.selected), there are 236 genes with an FDR-adjusted p<0.05 (they are also in T2DvsHealthy.M.W.corrected). Is it possible to generate a heatmap on all of those 236. I know you’ve already got heatmaps for the top 44, 100, and 1000 genes, but I didn’t see a heatmap for the 236 (if I’ve missed it, sorry, just point me in that direction!).
# 
############################################################################################################################################################
############################################################################################################################################################

# request: 

fdrEnsemblIds <- rownames(tab[ tab$fdr < 0.05, ])
fdGeneIds <- tab[ tab$fdr < 0.05, 6]






par(mfrow=c(1,1))
#for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
ci <- 1
for(g in fdGeneIds){
  
  png(paste0("fdr.",ci,".genesOfinterest.",g,".LOGbeanplot.png"), width=1600,height=800)
  #par(mar=c(3,8,3,8))
  par(mar=c(4,4,4,4))
  plot.new()
  
  geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
  print(paste(g,geneID))
  beanplot( log( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] + 0.0000001 ),# [rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] != 0],
            log( rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] + 0.0000001 ),
            col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
            names = c("Healthy","T2D"),
            bw="nrd0",
            main=g,
            xlim = c(0.2,2.8),
            cex.axis=2,
            cex.main=2,
            cex.names=2
  )
  text( x=0.3, y= log( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] + 0.0000001 ), names(rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] ), cex=0.6)
  text( x=2.7, y=log( rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] + 0.0000001 )     ,     names(rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] ),     cex=0.6)
  dev.off()
  
  
  png(paste0("NOduplicates.genesOfinterest.",g,".LOGbeanplotNOtext.png"), width=800,height=800)
  par(mar=c(4,4,4,4))
  beanplot( log( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] + 0.0000001 ),
            log( rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] + 0.0000001 ),
            col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
            names = c("Healthy","T2D"),
            bw="nrd0",
            main=g,
            cex.axis=2,
            cex.main=2,
            cex.names=2)
  dev.off()
  
  ci <- ci + 1
  
}


png(paste0("fdr.ALL.genesOfinterest.LOGbeanplot.png"), width=3000,height=4500)
par(mfrow=c(9,6))
#for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
ci <- 1
for(g in fdGeneIds){
  
  geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
  print(paste(g,geneID))
  par(mar=c(4,4,4,4))
  beanplot( log( rel_matrix[geneID,grepl("Healthy",colnames(rel_matrix))] + 0.0000001 ),
            log( rel_matrix[geneID,grepl("T2D",colnames(rel_matrix))] + 0.0000001 ),
            col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
            names = c("Healthy","T2D"),
            bw="nrd0",
            main=g,
            cex.axis=2,
            cex.main=2,
            cex.names=2)
  
}
dev.off()



#fdr <- T2DvsHealthy.M.W.corrected[ sortnames.fdr.T2DvsHealthy.M.W.corrected , ]
#fdr <- fdr[fdr$FDR<=0.05,]
 
#abund    <- HvsT2D[ ,  colnames(HvsT2D)[grep('norm', colnames(HvsT2D))] ]
abund <- count_matrix[  , colnames(count_matrix)[grep('norm', colnames(count_matrix))] ]
viz.abund <- abund
viz.abund[viz.abund == 0] <- 0.001

relabund <- t(decostand( t(abund), method = 'total'))
viz.relabund <- relabund
viz.relabund[viz.relabund == 0] <- 0.0000000001

#fdr.abund <- HvsT2D[ fdr$Ensembl_GeneID ,  colnames(HvsT2D)[grep('norm', colnames(HvsT2D))] ]

ord <-  c( colnames(viz.relabund)[ grepl("T2D",colnames(viz.relabund)) ] , colnames(viz.relabund)[ !grepl("T2D",colnames(viz.relabund)) ] )

png("20190612.fdrSignif_0.05_geneRelativeAbundance_heatmap_sidePval.png", height=2000, width=2000)
heatmap.2( -log(viz.relabund[ fdrEnsemblIds , ord ]), 
          col=c(redblue(150)),          
          trace="none",          
          margins = c(45,15),          
          keysize=0.4, key.par = list(cex=0.5),          
          labRow = fdGeneIds,
          cexRow = 1.4, cexCol = 1.4,
          ColSideColors = c("darkgreen","maroon")[ as.factor(grepl("T2D",ord)) ]#,
          #RowSideColors = c(redgreen(dim(fdr)[1]))[ order(fdr$FDR) ]
)
legend("bottomright",legend = c("Healthy","T2D"),pch=16,col=c("darkgreen","maroon"), cex = 1.4, inset = c(0.2,0), title = "Patient Status")
#legend("bottomleft",col = c("green","black","red"),pch=16,legend=c("more","...","less"), cex = 1.4, title = "Rank of Significance (FDR)")
dev.off()
 

png("20190612noSort.fdrSignif_0.05_geneRelativeAbundance_heatmap_sidePval.png", height=2000, width=2000)
heatmap.2( -log(viz.relabund[ fdrEnsemblIds , ord ]), 
           col=c(redblue(150)),          
           trace="none",          
           margins = c(45,15),          
           keysize=0.4, key.par = list(cex=0.5),          
           labRow = as.character(fdGeneIds),
           cexRow = 1.4, cexCol = 1.4,
           ColSideColors = c("darkgreen","maroon")[ as.factor(grepl("T2D",ord)) ],
           Colv = FALSE#,
           #RowSideColors = c(redgreen(dim(fdr)[1]))[ order(fdr$FDR) ]
)
legend("bottomright",legend = c("Healthy","T2D"),pch=16,col=c("darkgreen","maroon"), cex = 1.4, inset = c(0.2,0), title = "Patient Status")
#legend("bottomleft",col = c("green","black","red"),pch=16,legend=c("more","...","less"), cex = 1.4, title = "Rank of Significance (FDR)")
dev.off()


ipathway_brief_sum <- read.table("../ipathway_brief_sum.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE)
rownames(ipathway_brief_sum) <- ipathway_brief_sum[,1]
ipathway_brief_sum <- ipathway_brief_sum[,-1]
ipathway_brief_sum[ is.na(ipathway_brief_sum) ] <- 0
 
colnames(ipathway_brief_sum) %in% fdGeneIds
row_label_matrix <- matrix(nrow=length(fdGeneIds),ncol=dim(ipathway_brief_sum)[1])
rownames(row_label_matrix) <- fdGeneIds
colnames(row_label_matrix) <- rownames( ipathway_brief_sum )
row_label_matrix[ is.na(row_label_matrix) ] <- "white"
for(f in fdGeneIds){
  if(f %in% colnames(ipathway_brief_sum)){
    row_label_matrix[ f,  ] <- ipathway_brief_sum[ , f ]
  }
} 
row_label_matrix[ row_label_matrix == 0 ] <- "white"
row_label_matrix[ row_label_matrix == 1 ] <- "purple"

png("20190612noSort.labels.png", height=2000, width=2000)
heatmap.plus( 
              -log(viz.relabund[ fdrEnsemblIds , ord ]),
              col=c(redblue(150)),          
              margins = c(30,15),        
              labRow = fdGeneIds,
              RowSideColors = row_label_matrix,
              cexRow = 1.4, cexCol = 1.4,
              main="50 Significant @ FDR-adjusted P-value < 0.05"
            )
dev.off()
 

ipathwayEnsembl <- list()
ii <- 1
query <- colnames(ipathway_brief_sum)
query[ query == "ABCG1.1" ] <- "ABCG1"
for(i in query){
  if( i %in% HvsT2D$Gene_name ){
   ipathwayEnsembl[ii] <- HvsT2D[ HvsT2D$Gene_name == i,  ]
   ii <- ii + 1
  }else{
    print(i)
  }
}
ipathwayEnsembl <- as.character(unlist(ipathwayEnsembl))
    
ipathway_brief_sum_row_label_matrix <- ipathway_brief_sum
ipathway_brief_sum_row_label_matrix[ ipathway_brief_sum_row_label_matrix == 0] <- "white"
ipathway_brief_sum_row_label_matrix[ ipathway_brief_sum_row_label_matrix == 1] <- "purple"

png("ipathway.labels.png", height=2000, width=2000)

  heatmap.plus( 
    -log(viz.relabund[ ipathwayEnsembl , ord ]),
    col=c(redblue(150)),          
    margins = c(30,15),        
    labRow = fdGeneIds,
    RowSideColors = t(as.matrix(ipathway_brief_sum_row_label_matrix)),
    cexRow = 1.4, cexCol = 1.4,
    main="Diff. Exp Ingenuity Pathway Genes"
  )
  dev.off()

  
  keggList_sum <- read.table("../keggList_sum.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE)
  rownames(keggList_sum) <- keggList_sum[,1]
  keggList_sum <- keggList_sum[,-1]
  keggList_sum[ is.na(keggList_sum) ] <- 0
  keggList_sum <- keggList_sum[ , colnames(keggList_sum) != "osmE" ]
  
  keggList_Ensembl <- list()
  ii <- 1
  query <- colnames(keggList_sum)
  #query[ query == "osmE" ] <- "OSME"
  for(i in query){
    if( i %in% HvsT2D$Gene_name ){
      keggList_Ensembl[ii] <- HvsT2D[ HvsT2D$Gene_name == i,  ]
      ii <- ii + 1
    }else{
      print(i)
    }
  }
  keggList_Ensembl <- as.character(unlist(keggList_Ensembl))
  
  keggList_label_matrix <- keggList_sum
  keggList_label_matrix[ keggList_sum == 0] <- "white"
  keggList_label_matrix[ keggList_sum == 1] <- "purple"
  
  kegg.pop <- keggList_Ensembl[ keggList_Ensembl %in% rownames(viz.relabund) ]
  
  png("kegglist.labels.png", height=2000, width=2000)
  heatmap.plus( 
    -log(viz.relabund[ kegg.pop , ord ]),
    col=c(redblue(150)),          
    margins = c(30,15),        
    labRow = fdGeneIds,
    RowSideColors = t(as.matrix(keggList_label_matrix[ ,  keggList_Ensembl %in% rownames(viz.relabund)  ])),
    cexRow = 1.4, cexCol = 1.4,
    main="All Genes From Highlighted Kegg Pathways"
  ) 
   dev.off()
    
  
   
   
   
   
   
   
   
   metadat <- read.table("Metadata4REAL.txt", sep="\t", header=FALSE, fill=TRUE)
   metadat <- metadat[2:19,1:32]
   #metadat <- metadat[ !metadat$Sample %in% c(5,7) , ]
   
   pca_box_hash <- prcomp( t(rel_matrix[ fdrEnsemblIds , ord ]) , center = TRUE) 
   
   plot(pca_box_hash$x[,1:2],
        xlim=c(-0.001,0.007),
        pch=19, cex=4,
        #xlim=c(-0.001,0.003),
        col=c("maroon","darkgreen")[as.factor(grepl("Healthy",rownames(pca_box_hash$x)))],
        main="PCA Over 50 Genes at FDR < 0.05"
   )
   text(pca_box_hash$x[,1:2], 
        rownames(pca_box_hash$x), cex=0.5,
        xlim=c(-0.001,0.003)
   )
   
   
   
   biplot(pca_box_hash, main="PCA + Eigenvectors (Top 50 FDR-adjusted, FDR < 0.05)", main.cex = 0.5, 
          xlim=c(-0.3,1.2),
          cex=0.5,
          col=c("gray","red")
   )
   
   
   
   
   
   
   
   
   
  
  r2 <- RvE1[,4:25]
  r.abund    <- r2[ ,  colnames(r2)[grep('norm', colnames(r2))] ]
  viz.r.abund <- r.abund
  viz.r.abund[viz.r.abund == 0] <- 0.001
  
  r.relabund <- t(decostand( t(r.abund), method = 'total'))
  
  
  png("20190612.genesOfinterest.x5.perturb_beanplot.png", width=500,height=1200)
  plot.new()
  par(mfrow=c(5,1))
  for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
    geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
    beanplot( r.relabund[geneID, colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))])  ]] ,
              r.relabund[geneID, colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))])  ]] ,
              r.relabund[geneID, colnames(r.relabund)[grepl("_10nM",colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_10nM", colnames(r.relabund))])  ]] ,
              r.relabund[geneID, colnames(r.relabund)[grepl("_100nM",colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_100nM", colnames(r.relabund))])  ]] ,
              r.relabund[geneID, colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))])  ]] ,
              r.relabund[geneID, colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))])  ]] ,
              r.relabund[geneID, colnames(r.relabund)[grepl("_10nM",colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_10nM", colnames(r.relabund))])  ]] ,
              r.relabund[geneID, colnames(r.relabund)[grepl("_100nM",colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_100nM", colnames(r.relabund))])  ]] ,
              col = list(c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),
                         c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3")
              ), 
              names = c("Healthy (0 nM)","Healthy (1 nM)","Healthy (10 nM)","Healthy (100 nM)","T2D (0 nM)","T2D (1 nM)","T2D (10 nM)","T2D (100 nM)"),
              las=2,
              main=paste0("\n\n",g)
    )
  }
  dev.off()
  
  
  png("20190612.LOG.genesOfinterest.x5.perturb_beanplot.png", width=500,height=1200)
  plot.new()
  par(mfrow=c(5,1))
  for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
    geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
    beanplot( log( r.relabund[geneID, colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))])  ]] ),
              log( r.relabund[geneID, colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))])  ]] ),
              log( r.relabund[geneID, colnames(r.relabund)[grepl("_10nM",colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_10nM", colnames(r.relabund))])  ]] ),
              log( r.relabund[geneID, colnames(r.relabund)[grepl("_100nM",colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_100nM", colnames(r.relabund))])  ]] ),
              log( r.relabund[geneID, colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))])  ]] ),
              log( r.relabund[geneID, colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))])  ]] ),
              log( r.relabund[geneID, colnames(r.relabund)[grepl("_10nM",colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_10nM", colnames(r.relabund))])  ]] ),
              log( r.relabund[geneID, colnames(r.relabund)[grepl("_100nM",colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_100nM", colnames(r.relabund))])  ]] ),
              col = list(c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),
                         c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3")
              ), 
              names = c("Healthy (0 nM)","Healthy (1 nM)","Healthy (10 nM)","Healthy (100 nM)","T2D (0 nM)","T2D (1 nM)","T2D (10 nM)","T2D (100 nM)"),
              las=2,
              main=paste0("\n\n",g)
    )
  }
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  






png("20190612.fdrSignif_0.05_geneRelativeAbundance_heatmap_sideDirection.png", height=4000, width=2000)
heatmap.2( -log(viz.relabund[ fdr$Ensembl_GeneID , ord ]),            col=c(redblue(150)),           trace="none",           margins = c(45,15),           keysize=0.4, key.par = list(cex=0.5),           labRow = fdr$Gene_name, cexRow = 1.4, cexCol = 2.4,
           ColSideColors = c("darkgreen","maroon")[ as.factor(grepl("T2D",ord)) ], #c("darkgreen","maroon")[ as.factor(grepl("T2D",colnames(viz.relabund))) ],
           RowSideColors = c("red","green")[as.numeric(as.factor(rowMeans( viz.relabund[ fdr$Ensembl_GeneID , colnames(viz.relabund)[ grepl("T2D",colnames(viz.relabund)) ]  ] ) - rowMeans( viz.relabund[ fdr$Ensembl_GeneID , colnames(viz.relabund)[ !grepl("T2D",colnames(viz.relabund)) ]  ] ) > 0))],           
)
legend("topright",legend = c("Healthy","T2D"),pch=16,col=c("darkgreen","maroon"), cex = 2, inset = c(0.2,0), title = "Patient Status")
legend("topright",col = c("green","red"),pch=16,legend=c("Positive","Negative"), cex = 2, title = "Association with T2D Status")
dev.off()

png("20190612.fdrSignif_0.05_geneRelativeAbundance_heatmap_sidePval.sort.png", height=4000, width=2000)
heatmap.2( -log(viz.relabund[ fdr$Ensembl_GeneID , ord ]), 
           col=c(redblue(150)),          trace="none",          margins = c(45,15),          keysize=0.4, key.par = list(cex=0.5),          labRow = fdr$Gene_name, cexRow = 1.4, cexCol = 2.4,
           ColSideColors = c("darkgreen","maroon")[ as.factor(grepl("T2D",ord)) ],
           RowSideColors = c(redgreen(dim(fdr)[1]))[ order(fdr$FDR) ],
           Colv = FALSE
)
legend("topright",legend = c("Healthy","T2D"),pch=16,col=c("darkgreen","maroon"), cex = 2, inset = c(0.2,0), title = "Patient Status")
legend("topright",col = c("green","black","red"),pch=16,legend=c("more","...","less"), cex = 2, title = "Rank of Significance (FDR)")
dev.off()

png("20190612.fdrSignif_0.05_geneRelativeAbundance_heatmap_sideDirection.sort.png", height=4000, width=2000)
heatmap.2( -log(viz.relabund[ fdr$Ensembl_GeneID , ord ]),            col=c(redblue(150)),           trace="none",           margins = c(45,15),           keysize=0.4, key.par = list(cex=0.5),           labRow = fdr$Gene_name, cexRow = 1.4, cexCol = 2.4,
           ColSideColors = c("darkgreen","maroon")[ as.factor(grepl("T2D",ord)) ], #c("darkgreen","maroon")[ as.factor(grepl("T2D",colnames(viz.relabund))) ],
           RowSideColors = c("red","green")[as.numeric(as.factor(rowMeans( viz.relabund[ fdr$Ensembl_GeneID , colnames(viz.relabund)[ grepl("T2D",colnames(viz.relabund)) ]  ] ) - rowMeans( viz.relabund[ fdr$Ensembl_GeneID , colnames(viz.relabund)[ !grepl("T2D",colnames(viz.relabund)) ]  ] ) > 0))],           
           Colv = FALSE
)
legend("topright",legend = c("Healthy","T2D"),pch=16,col=c("darkgreen","maroon"), cex = 2, inset = c(0.2,0), title = "Patient Status")
legend("topright",col = c("green","red"),pch=16,legend=c("Positive","Negative"), cex = 2, title = "Association with T2D Status")
dev.off()


# include perturbaton counts

nuNam <- c( as.character(fdr$Ensembl_GeneID) , as.character(unique(perturbation[ perturbation$FDR < 0.05 , "Ensembl_GeneID" ])) )

#r2 <- RvE1[sortnames.fdr.T2DvsHealthy.M.W.corrected,4:25]
r2 <- RvE1[,4:25]
r.abund    <- r2[ ,  colnames(r2)[grep('norm', colnames(r2))] ]
viz.r.abund <- r.abund
viz.r.abund[viz.r.abund == 0] <- 0.001

r.relabund <- t(decostand( t(r.abund), method = 'total'))
viz.r.relabund <- r.relabund
viz.r.relabund[viz.r.relabund == 0] <- 0.0000000001

numID <- list()
n <- 1
r.ord <- sort(c(ord,colnames(viz.r.relabund)))
for(r in r.ord){
  qqq <- strsplit(r,"\\.")[[1]][2]
  qq2 <- strsplit(qqq[1],"_")[[1]][1]
  numID[n] <- as.numeric(qq2)
  n <- n + 1
}
names(numID) <- r.ord
IDnum <- names(numID)
names(IDnum) <- numID
nuOrd <- list()
for(w in 1:dim(metadata_corrected)[1]){
  if( metadata_corrected$Sample[w] != 5){
    nuOrd[w] <- IDnum[names(IDnum) == metadata_corrected$Sample[w]]
  }
}
nuOrd <- unlist(nuOrd)

both.relabund <- matrix(nrow=dim(viz.r.relabund)[1],ncol=length(nuOrd))
rownames(both.relabund) <- rownames(viz.r.relabund)
colnames(both.relabund) <- nuOrd
for(r in rownames(viz.r.abund)){
  both.relabund[r,ord                     ] <- viz.relabund[r,]
  both.relabund[r,colnames(viz.r.relabund)] <- viz.r.relabund[r,]
}

#metadata_col_matrix <- 
png("noDUPE.allAbundance.allSignif.png", height=4000, width=2000)
heatmap.2( -log(both.relabund[ nuNam , 
                               nuOrd 
            ]),            
           labCol = paste(metadata_corrected[metadata_corrected$Sample != 5,]$Treatment,metadata_corrected[metadata_corrected$Sample != 5,]$ID,metadata_corrected[metadata_corrected$Sample != 5,]$Sample,sep="\t"),
           col=c(redblue(150)),           trace="none",           margins = c(45,15),           keysize=0.4, key.par = list(cex=0.5),           
           labRow = HvsT2D[nuNam,]$Gene_name,
           ColSideColors = c("darkgreen","maroon")[ as.factor(metadata_corrected[metadata_corrected$Sample != 5,]$Disease) ], 
           colsep = c(c(1:3),6,10,11,c(15:19),c(21:25),29,33,c(37:43)), sepwidth = 0.2,
           cexRow = 1.4, cexCol = 1.8,
           RowSideColors = c("red","green")[as.numeric(as.factor(rowMeans( viz.relabund[ nuNam , colnames(viz.relabund)[ grepl("T2D",colnames(viz.relabund)) ]  ] ) - rowMeans( viz.relabund[ nuNam , colnames(viz.relabund)[ !grepl("T2D",colnames(viz.relabund)) ]  ] ) > 0))],           
           Colv = FALSE
)
legend("topright",legend = c("Healthy","T2D"),pch=16,col=c("darkgreen","maroon"), cex = 2, inset = c(0.3,0), title = "Patient Status")
legend("topright",col = c("green","red"),pch=16,legend=c("Positive","Negative"), cex = 2, title = "Baseline Association with T2D Status")
dev.off()
 

png("allAbundance.allSignif.hierarchical.png", height=4000, width=2000)
heatmap.2( -log(both.relabund[ nuNam , 
                               nuOrd 
                               ]),            
           labCol = paste(metadata_corrected$Sample,metadata_corrected$Treatment,metadata_corrected$ID,sep="\t"),
           col=c(redblue(150)),           trace="none",           margins = c(45,15),           keysize=0.4, key.par = list(cex=0.5),           
           labRow = HvsT2D[nuNam,]$Gene_name,
           ColSideColors = c("darkgreen","maroon")[ as.factor(metadata_corrected$Disease) ], 
           colsep = c(1:10,14:15,19:31,35,38,39,40,41), sepwidth = 0.2,
           cexRow = 1.4, cexCol = 1.8,
           RowSideColors = c("red","green")[as.numeric(as.factor(rowMeans( viz.relabund[ nuNam , colnames(viz.relabund)[ grepl("T2D",colnames(viz.relabund)) ]  ] ) - rowMeans( viz.relabund[ nuNam , colnames(viz.relabund)[ !grepl("T2D",colnames(viz.relabund)) ]  ] ) > 0))]
)
legend("topright",legend = c("Healthy","T2D"),pch=16,col=c("darkgreen","maroon"), cex = 2, inset = c(0.3,0), title = "Patient Status")
legend("topright",col = c("green","red"),pch=16,legend=c("Positive","Negative"), cex = 2, title = "Baseline Association with T2D Status")
dev.off()



############################################################################################################################################################
############################################################################################################################################################
# 2.
# I know you’ve been waiting on top genes to generate violin plots. Can we generate them for:
#   H vs. Diabetes: NECTIN2, LILRB4
# Perturbation: Within diabetic cohort (all T2D comparisons): RN7SL2
# Perturbation: Diabetic vs Healthy (all T2D vs healthy comparisons): LILRB5, AKR1C1
############################################################################################################################################################
############################################################################################################################################################

library(beanplot)

#beanplot( as.numeric(t( l[ q , grepl("Healthy",colnames(l)) ] ))[ !is.na(as.numeric(t( l[ q , grepl("Healthy",colnames(l)) ] ))) ]  ,
#          as.numeric(t( l[ q , !grepl("Healthy",colnames(l)) ]))[ !is.na(as.numeric(t( l[ q , !grepl("Healthy",colnames(l)) ]))) ],
#          bw="nrd0" ,
#          col = list(c("darkgreen","green","black","gray3"),c("red","black","black","gray3")), names = c("Healthy","T2D") #,
#)

png("NOduplicates.genesOfinterest.x13.beanplot.png", width=900,height=1000)
plot.new()
par(mfrow=c(5,3))
#for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
for(g in c("NECTIN2","FAM3B","PLPP3","KCNH4","VWF","C3","SLC9A4","RBPMS2","ENPP7P8","MX1","C1orf61","GNPDA1","IRAK2")){
  geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
  beanplot( relabund[geneID,grepl("Healthy",colnames(relabund))],
            relabund[geneID,grepl("T2D",colnames(relabund))],
            col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
            names = c("Healthy","T2D"),
            main=g
  )
}
dev.off()

png("genesOfinterest.x238.beanplot.png", width=2000,height=6200)
plot.new()
par(mfrow=c(30,8))
nuID <- paste0("\n\n",HvsT2D[ nuNam , ]$Gene_name)
ni <- 1
for(geneID in nuNam){
  beanplot( relabund[geneID,grepl("Healthy",colnames(relabund))],
            relabund[geneID,grepl("T2D",colnames(relabund))],
            col = list(c("darkgreen","green","black","gray3"),c("maroon","pink","black","gray3")), 
            names = c("Healthy","T2D"),
            main=nuID[ni]
  )
  ni <- ni + 1
}
dev.off()

png("genesOfinterest.x5.perturb_beanplot.png", width=500,height=1200)
plot.new()
par(mfrow=c(5,1))
for(g in c("NECTIN2","LILRB4","RN7SL2","LILRB5","AKR1C1")){
  geneID <- as.character(HvsT2D[ HvsT2D$Gene_name==g , "Ensembl_GeneID" ])
  beanplot( both.relabund[geneID, colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))])  ]] ,
            both.relabund[geneID, colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))])  ]] ,
            both.relabund[geneID, colnames(r.relabund)[grepl("_10nM",colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_10nM", colnames(r.relabund))])  ]] ,
            both.relabund[geneID, colnames(r.relabund)[grepl("_100nM",colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_100nM", colnames(r.relabund))])  ]] ,
            both.relabund[geneID, colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))])  ]] ,
            both.relabund[geneID, colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))])  ]] ,
            both.relabund[geneID, colnames(r.relabund)[grepl("_10nM",colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_10nM", colnames(r.relabund))])  ]] ,
            both.relabund[geneID, colnames(r.relabund)[grepl("_100nM",colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_100nM", colnames(r.relabund))])  ]] ,
            col = list(c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),
                       c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3")
                       ), 
            names = c("Healthy (0 nM)","Healthy (1 nM)","Healthy (10 nM)","Healthy (100 nM)","T2D (0 nM)","T2D (1 nM)","T2D (10 nM)","T2D (100 nM)"),
            las=2,
            main=paste0("\n\n",g)
  )
}
dev.off()

png("genesOfinterest.x238.perturb_beanplot.png", width=5000,height=9000)
plot.new()
par(mfrow=c(30,8))
nuID <- paste0("\n\n",HvsT2D[ nuNam , ]$Gene_name)
ni <- 1
for(geneID in nuNam){
  if(!geneID %in% c("ENSG00000197705","ENSG00000185479")){
    beanplot( both.relabund[geneID, colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))])  ]] ,
              both.relabund[geneID, colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))])  ]] ,
              both.relabund[geneID, colnames(r.relabund)[grepl("_10nM",colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_10nM", colnames(r.relabund))])  ]] ,
              both.relabund[geneID, colnames(r.relabund)[grepl("_100nM",colnames(r.relabund))][  grepl("Healthy", colnames(r.relabund)[grepl("_100nM", colnames(r.relabund))])  ]] ,
              both.relabund[geneID, colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_0nM", colnames(r.relabund))])  ]] ,
              both.relabund[geneID, colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_1nM", colnames(r.relabund))])  ]] ,
              both.relabund[geneID, colnames(r.relabund)[grepl("_10nM",colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_10nM", colnames(r.relabund))])  ]] ,
              both.relabund[geneID, colnames(r.relabund)[grepl("_100nM",colnames(r.relabund))][  grepl("T2D", colnames(r.relabund)[grepl("_100nM", colnames(r.relabund))])  ]] ,
              col = list(c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),c("darkgreen","green","black","gray3"),
                         c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3"),c("maroon","pink","black","gray3")
              ), 
              names = c("Healthy (0 nM)","Healthy (1 nM)","Healthy (10 nM)","Healthy (100 nM)","T2D (0 nM)","T2D (1 nM)","T2D (10 nM)","T2D (100 nM)"),
              las=2,
              main=nuID[ni]
    )
    ni <- ni + 1
  }
}
dev.off()

############################################################################################################################################################
############################################################################################################################################################
# 3.
# PC plots: can we redo the ordinal variables by bins?
#   Age by decades (</=30, 31-40, 41-50, 51-60, 61-70)
# BMI by 5 (</=20, >20-25, >25-30, >30-35, >35-40, >40)
# HbA1c levels (by clinical relevance): <6.4 (pre-diabetes), 6.5-8 (diabetes), >8 (uncontrolled diabetes)
# Cholesterol (by clinical relevance): <200 (good), 200-239 (borderline high), >240 (high) (I know most of our data is <200, so we might have to break that down further, but I figure starting with the clinical measures is best)
# Glucose (by clinical relevance): <140 (normal), >/=140 (high) (*I need to check if this is fasting or not or this may change, but I would suspect it’s not fasting glucose based on the metadata range)
# Neutrophils (I’m just picking based on the metadata ranges): <50, 50-99,100-149,150-200, >200 (1 individual) 
# Monocytes (also based on metadata ranges): <50, 50-99,100-149,150-200 (1 individual)
# Additional feedback:
#   Update 1 (above): When adding 236 genes, experiment with text readability and/or which genes to include
# Heatmaps of T2D vs healthy gene abundance:
#   Change to highlight significance of effect between case and control by re-ordering and comparing
# Volcano plots: update for both standard and perturbation analyses to show gene IDs above newly suggested cutoff.
# Also generate new beanplots for all signfiicant genes in standard and perturbation analyses
# For FDR-adjusted perturbation tests, convert presentation to a bubble plot a la attached michelet2018.pdf using bubblePlot
############################################################################################################################################################
############################################################################################################################################################





