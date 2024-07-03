library(colorRamps)
library(vegan)

#library(vegan) 
library(ape)
library(gtools)
library(colorRamps)
library(fpc)
library(cluster)
library(dplyr)
library(ade4)
library(tidyr)
library(ggplot2)
library(broom)
library(MASS)

setwd ("/Users/apple/Desktop/DEV.gcid")
# NEBnext removed
samp_assign <- c("x29UC",                    "x68UC",               "x33UC",                     "x43UC",                     "x68UC",                     "x59UC",                                            "x68UC",             "x73UC",             "x29UC",              "x29UC",                "x33UC",               "x43UC");
prep_assign <- c("BloodAgar",               "META",               "BloodAgar",                "BloodAgar",                "BloodAgar",                "ORI",                                     "ORI",              "ORI",              "ORI",                "META",               "META",               "META") 
days_after  <- c("0",                        "0",                   "0",                          "21",                         "0",                          "84",                "0",                                  "0",               "91",                    "0",                  "0",                    "0",                  "21");                                   

#######################################################################################################################################################################################
# AMR ###########################################################################################################################################################################################
#################################################################################################################################################################################################

amr_tab <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/amr_count_report.tsv", sep="\t", header=TRUE)
tmp_names <- amr_tab$orf_id
t <- length(colnames(amr_tab))
amr_tab <- amr_tab[,2:t]
amr_tab <- mapply(amr_tab, FUN=as.numeric)
amr_tab <- amr_tab[complete.cases(amr_tab),]
rownames(amr_tab) <- tmp_names


# remove nebnext
amr_tab <- amr_tab[,colnames(amr_tab)!="X67.68UC.ORI.NEBNEXT.PE.P51.1_S13_001"]

amr_tab <- amr_tab[2:length(rownames(amr_tab)),]
amr_tab <- decostand(amr_tab,"total",2,na.rm=TRUE)
amr_tab <- -log(amr_tab+0.0000001)

library(gplots)
png(filename=paste0("NU_08a_Only_AMR_Transcripts_Heatmap_SQUARE.png"), width = 1000, height = 1000)
plot.new()
test <- heatmap.2(amr_tab)
heatmap.2(amr_tab,margins=c(22,5),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ as.factor(samp_assign) ],cexRow = 0.1)
par(new=TRUE)
legend("bottomleft",c(unique(samp_assign)),border="white", bty = "n",pch=22,pt.bg=c(colorRamps::matlab.like2(length(unique(samp_assign)))))
dev.off()

write.table(amr_tab, file="NU_08a_Only_AMR_Transcripts_Heatmap_SQUARE.tsv", quote=FALSE, sep='\t')

# debug
# assign color
h <- heatmap.2(amr_tab,margins=c(22,5),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ as.factor(samp_assign) ],cexRow = 0.1)
colz <- c("red","blue")[as.factor(cutree( as.hclust(h$rowDendrogram),k = 2 ) == 1)]
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue")[cutree( as.hclust(h$rowDendrogram), 2 ) ] )

clust.1.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 2 ) == 1)[cutree( as.hclust(h$rowDendrogram),k = 2 ) == 1]
clust.2.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 2 ) == 2)[cutree( as.hclust(h$rowDendrogram),k = 2 ) == 2]

#################################
#################################
#################################

best_hit_species <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/annotation_20170710.filtered.summary_best_hit_species.tab", sep="\t", header=TRUE)
inlist <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/inlist.txt", sep="\t", header=TRUE)
inlist.v <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/inlist_viruses.txt", sep="\t", header=TRUE)

inlist <- unname(unlist(inlist)[1:136])
inlist.v <- unname(unlist(inlist.v))

species_counts <- best_hit_species[4:16]
colnames(species_counts) <- colnames(best_hit_species)[4:16]
rownames(species_counts) <- best_hit_species$best_hit_species

# v4 - remove NebNext
species_counts <- species_counts[,colnames(species_counts)!="X67.68UC.ORI.NEBNEXT.PE.P51.1_S13_001"]

species_counts.bact <- species_counts[unique(grep( paste(inlist,collapse='|'), rownames(species_counts), value=TRUE )), ]
species_counts.virus <- species_counts[unique(grep( paste(inlist.v,collapse='|'), rownames(species_counts), value=TRUE )), ]










########################################
# all

r_c <- 1
genus_id <- list()
for(r in rownames(best_hit_species)){
  genus_id[r_c] <- strsplit(as.character(best_hit_species$best_hit_species),' ')[[r_c]][1]
  r_c <- r_c + 1
}

unique_genus_id <- unique(unlist(genus_id))

genus_abundance_table <- data.frame(matrix(nrow=length(unique_genus_id),ncol=length(colnames(species_counts))))
rownames(genus_abundance_table) <- unique_genus_id
colnames(genus_abundance_table) <- colnames(species_counts)
for(g in unique_genus_id){
  genus_abundance_table[g,] <- colSums(species_counts[grepl(g,unlist(rownames(species_counts))),])
}

genus_abundance_table.rel <- t(decostand(t(genus_abundance_table), method="total"))
genus_abundance_table.rel.sorted <- genus_abundance_table.rel[rev(order(rowSums(genus_abundance_table.rel))),]
# .4 updates
genus_abundance_table.rel.sorted <- genus_abundance_table.rel.sorted[!grepl(']',rownames(genus_abundance_table.rel.sorted)),]

png(filename=paste0("NU_Genus_Heatmap.top50.samp.png"), width = 1000, height = 1000)
heatmap.2(-log(as.matrix(genus_abundance_table.rel.sorted[1:50,])+0.000001),margins=c(18,10),
               ColSideColors = c(matlab.like(length(unique(samp_assign))))[as.factor(samp_assign)],
               cexRow = 1.2)
dev.off()

write.table(genus_abundance_table.rel.sorted, file="NU_Genus_Heatmap.all.tsv", quote=FALSE, sep='\t')


png(filename=paste0("NU_Genus_Heatmap.top50.prep.png"), width = 1000, height = 1000)
heatmap.2(-log(as.matrix(genus_abundance_table.rel.sorted[1:50,])+0.000001),margins=c(18,10),
          ColSideColors = c(primary.colors(length(unique(prep_assign))))[as.factor(prep_assign)],
          cexRow = 1.2)
dev.off()

NU_Genus_Heatmap.top50.samp.png

#h <- heatmap.2(amr_tab,margins=c(22,5),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ as.factor(samp_assign) ],cexRow = 0.1)
colz <- c("red","blue")[as.factor(cutree( as.hclust(h$rowDendrogram),k = 2 ) == 1)]

plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue")[cutree( as.hclust(h$rowDendrogram), 2 ) ] )
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue","green")[cutree( as.hclust(h$rowDendrogram), 3 ) ] )
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c(rainbow(8))[cutree( as.hclust(h$rowDendrogram), 8 ) ] )


clust.1.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 3 ) == 1)[cutree( as.hclust(h$rowDendrogram),k = 3 ) == 1]
clust.2.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 3 ) == 2)[cutree( as.hclust(h$rowDendrogram),k = 3 ) == 2]
clust.3.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 3 ) == 3)[cutree( as.hclust(h$rowDendrogram),k = 3 ) == 3]








########################################
# bacterial

r_c <- 1
genus_id <- list()
for(r in rownames(best_hit_species)){
  genus_id[r_c] <- strsplit(as.character(best_hit_species$best_hit_species),' ')[[r_c]][1]
  r_c <- r_c + 1
}

unique_genus_id <- unique(unlist(genus_id))

genus_abundance_table <- data.frame(matrix(nrow=length(unique_genus_id),ncol=length(colnames(species_counts.bact))))
rownames(genus_abundance_table) <- unique_genus_id
colnames(genus_abundance_table) <- colnames(species_counts.bact)
for(g in unique_genus_id){
  genus_abundance_table[g,] <- colSums(species_counts.bact[grepl(g,unlist(rownames(species_counts.bact))),])
}

genus_abundance_table.rel <- t(decostand(t(genus_abundance_table), method="total"))
genus_abundance_table.rel.sorted <- genus_abundance_table.rel[rev(order(rowSums(genus_abundance_table.rel))),]
genus_abundance_table.rel.sorted <- genus_abundance_table.rel.sorted[!grepl(']',rownames(genus_abundance_table.rel.sorted)),]

png(filename=paste0("BACT_Genus_Heatmap.top50.samp.png"), width = 1000, height = 1000)
heatmap.2(-log(as.matrix(genus_abundance_table.rel.sorted[1:50,])+0.000001),margins=c(18,14),
          ColSideColors = c(matlab.like(length(unique(samp_assign))))[as.factor(samp_assign)],
          cexRow = 1.2)
dev.off()


write.table(genus_abundance_table.rel.sorted, file="NU_Genus_Heatmap.bact.tsv", quote=FALSE, sep='\t')


png(filename=paste0("BACT_Genus_Heatmap.top50.prep.png"), width = 1000, height = 1000)
heatmap.2(-log(as.matrix(genus_abundance_table.rel.sorted[1:50,])+0.000001),margins=c(18,14),
          ColSideColors = c(primary.colors(length(unique(prep_assign))))[as.factor(prep_assign)],
          cexRow = 1.2)
dev.off()


#h <- heatmap.2(amr_tab,margins=c(22,5),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ as.factor(samp_assign) ],cexRow = 0.1)
colz <- c("red","blue")[as.factor(cutree( as.hclust(h$rowDendrogram),k = 2 ) == 1)]

h <- heatmap.2(-log(as.matrix(genus_abundance_table.rel.sorted[1:50,])+0.000001))
              
png(filename=paste0("BACT_phylo.n2.png"), width = 500, height = 1000)
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue")[cutree( as.hclust(h$rowDendrogram), 2 ) ] )
dev.off()

png(filename=paste0("BACT_phylo.n16.png"), width = 500, height = 1000)
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c(primary.colors(14))[cutree( as.hclust(h$rowDendrogram), 16 ) ] )
dev.off()
   
clust.1.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 3 ) == 1)[cutree( as.hclust(h$rowDendrogram),k = 3 ) == 1]
clust.2.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 3 ) == 2)[cutree( as.hclust(h$rowDendrogram),k = 3 ) == 2]
clust.3.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 3 ) == 3)[cutree( as.hclust(h$rowDendrogram),k = 3 ) == 3]





plotit <- -log(as.matrix(genus_abundance_table.rel.sorted[1:26,])+0.000001)
plotit <- plotit[rownames(plotit)!="[unknown]",]
png(filename=paste0("BACT_Genus_Heatmap.top26_withoutUnknown.samp.png"), width = 1000, height = 1000)
heatmap.2(as.matrix(plotit),margins=c(18,14),
          ColSideColors = c(matlab.like(length(unique(samp_assign))))[as.factor(samp_assign)],
          labCol = days_after,
          cexRow = 1.2)
dev.off()
png(filename=paste0("BACT_Genus_Heatmap.top26_withoutUnknown.prep.png"), width = 1000, height = 1000)
heatmap.2(as.matrix(plotit),margins=c(18,14),
          ColSideColors = c(primary.colors(length(unique(prep_assign))))[as.factor(prep_assign)],
          labCol = days_after,
          cexRow = 1.2)
dev.off()


h <- heatmap.2(as.matrix(plotit))

png(filename=paste0("BACT_Genus_Heatmap.top26_withoutUnknown.n2.png"), width = 500, height = 1000)
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue")[cutree( as.hclust(h$rowDendrogram), 2 ) ] )
dev.off()

png(filename=paste0("BACT_Genus_Heatmap.top26_withoutUnknown.n5.png"), width = 500, height = 1000)
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue","darkgreen","purple","black")[cutree( as.hclust(h$rowDendrogram), 5 ) ] )
dev.off()
 







########################################
# viral

r_c <- 1
genus_id <- list()
for(r in rownames(best_hit_species)){
  genus_id[r_c] <- strsplit(as.character(best_hit_species$best_hit_species),' ')[[r_c]][1]
  r_c <- r_c + 1
}

unique_genus_id <- unique(unlist(genus_id))

genus_abundance_table <- data.frame(matrix(nrow=length(unique_genus_id),ncol=length(colnames(species_counts.virus))))
rownames(genus_abundance_table) <- unique_genus_id
colnames(genus_abundance_table) <- colnames(species_counts.virus)
for(g in unique_genus_id){
  genus_abundance_table[g,] <- colSums(species_counts.virus[grepl(g,unlist(rownames(species_counts.virus))),])
}

genus_abundance_table.rel <- t(decostand(t(genus_abundance_table), method="total"))
genus_abundance_table.rel.sorted <- genus_abundance_table.rel[rev(order(rowSums(genus_abundance_table.rel))),]
genus_abundance_table.rel.sorted <- genus_abundance_table.rel.sorted[!grepl(']',rownames(genus_abundance_table.rel.sorted)),]

png(filename=paste0("VIRAL_Genus_Heatmap.top50.samp.png"), width = 1000, height = 1000)
heatmap.2(-log(as.matrix(genus_abundance_table.rel.sorted[1:50,])+0.000001),margins=c(18,10),
          ColSideColors = c(matlab.like(length(unique(samp_assign))))[as.factor(samp_assign)],
          cexRow = 1.2)
dev.off()

write.table(genus_abundance_table.rel.sorted, file="NU_Genus_Heatmap.viral.tsv", quote=FALSE, sep='\t')

png(filename=paste0("VIRAL_Genus_Heatmap.top50.prep.png"), width = 1000, height = 1000)
heatmap.2(-log(as.matrix(genus_abundance_table.rel.sorted[1:50,])+0.000001),margins=c(18,10),
          ColSideColors = c(primary.colors(length(unique(prep_assign))))[as.factor(prep_assign)],
          cexRow = 1.2)
dev.off()

#h <- heatmap.2(amr_tab,margins=c(22,5),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ as.factor(samp_assign) ],cexRow = 0.1)
colz <- c("red","blue")[as.factor(cutree( as.hclust(h$rowDendrogram),k = 2 ) == 1)]

plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue")[cutree( as.hclust(h$rowDendrogram), 2 ) ] )
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue","green")[cutree( as.hclust(h$rowDendrogram), 3 ) ] )
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c(rainbow(8))[cutree( as.hclust(h$rowDendrogram), 8 ) ] )


clust.1.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 3 ) == 1)[cutree( as.hclust(h$rowDendrogram),k = 3 ) == 1]
clust.2.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 3 ) == 2)[cutree( as.hclust(h$rowDendrogram),k = 3 ) == 2]
clust.3.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 3 ) == 3)[cutree( as.hclust(h$rowDendrogram),k = 3 ) == 3]





png(filename=paste0("VIRAL_Genus_Heatmap.top10.samp.png"), width = 1000, height = 1000)
heatmap.2(-log(as.matrix(genus_abundance_table.rel.sorted[1:10,])+0.000001),margins=c(18,15),
          ColSideColors = c(matlab.like(length(unique(samp_assign))))[as.factor(samp_assign)],
          cexRow = 1.2)
dev.off()

png(filename=paste0("VIRAL_Genus_Heatmap.top10.prep.png"), width = 1000, height = 1000)
heatmap.2(-log(as.matrix(genus_abundance_table.rel.sorted[1:10,])+0.000001),margins=c(18,15),
          ColSideColors = c(primary.colors(length(unique(prep_assign))))[as.factor(prep_assign)],
          cexRow = 1.2)
dev.off()



h <- heatmap.2(-log(as.matrix(genus_abundance_table.rel.sorted[1:13,])+0.000001))
png(filename=paste0("VIRAL_phylo.n2.png"), width = 1000, height = 1000)
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue")[cutree( as.hclust(h$rowDendrogram), 2 ) ] )
dev.off()
 






















allcounts.rgi.rel <- as.matrix(decostand(allcounts.rgi,"total",2,na.rm=TRUE))












h <- heatmap.2(amr_tab,margins=c(22,5),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ as.factor(samp_assign) ],cexRow = 0.1)
colz <- c("red","blue")[as.factor(cutree( as.hclust(h$rowDendrogram),k = 2 ) == 1)]
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue")[cutree( as.hclust(h$rowDendrogram), 2 ) ] )



clust.1.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 2 ) == 1)[cutree( as.hclust(h$rowDendrogram),k = 2 ) == 1]

clust.2.orfids <- names(cutree( as.hclust(h$rowDendrogram),k = 2 ) == 2)[cutree( as.hclust(h$rowDendrogram),k = 2 ) == 2]
































##############################################
# 0. curate output table #############( from amrTXT.4.20170927_byCarbapenemResType.R )###########
##############################################

allcounts <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/annot3.tab", sep="\t", header=TRUE, fill=TRUE)
rownames(allcounts) <- allcounts$orf_id
allcounts.rgi <- allcounts[,34:46]
allcounts.rgi.rel <- as.matrix(decostand(allcounts.rgi,"total",2,na.rm=TRUE))
orfhits.t <- read.table("/Users/apple/Desktop/DEV.gcid/format2.assembly.orf.strict.txt", sep="\t", header=TRUE)
delivery_csv.t <- orfhits.t[,c("Best_Hit_evalue","Best_Hit_ARO","ARO","ARO_name","SNP","Best_Hit_ARO_category","ARO_category")]
rownames(delivery_csv.t) <-orfhits.t$ORF_ID
delivery.df <- as.data.frame(delivery_csv.t)
delivery_csv.df.sort.Best_Hit_ARO <- delivery.df[order(delivery.df$Best_Hit_ARO),]
defined.rgi.rel <- delivery_csv.df.sort.Best_Hit_ARO[rownames(delivery_csv.df.sort.Best_Hit_ARO) %in% rownames(allcounts.rgi.rel),]

abundance.1 <- allcounts[names(cutree( as.hclust(h$rowDendrogram),k = 2 ))[cutree( as.hclust(h$rowDendrogram),k = 2 ) == 1] %in% rownames(allcounts),]
allcounts.1.rgi <- abundance.1[,34:46]
allcounts.1.rgi.rel <- as.matrix(decostand(allcounts.1.rgi,"total",2,na.rm=TRUE))

#png(filename=paste0("NU_occurrence_histo.png"), width = 1000, height = 1000)

symbol_abundance_table.1 <- as.data.frame(table(as.factor((unlist(strsplit(as.character((abundance.1$best_hit_annotation)),'; '))))))
symbol_abundance_table.1 <- as.data.frame(table(as.factor((unlist(strsplit(as.character((abundance.1$best_hit_annotation)),', '))))))
symbol_abundance_table.1$Freq <- as.numeric(symbol_abundance_table.1$Freq)
symbol_abundance_table.1 <- symbol_abundance_table.1[order(-symbol_abundance_table.1$Freq),]
par(mar = c(4, 25, 4, 2))

#dev.off()

abundance.2 <- allcounts[names(cutree( as.hclust(h$rowDendrogram),k = 2 ))[cutree( as.hclust(h$rowDendrogram),k = 2 ) == 2] %in% rownames(allcounts),]
allcounts.2.rgi <- abundance.2[,34:46]
allcounts.2.rgi.rel <- as.matrix(decostand(allcounts.2.rgi,"total",2,na.rm=TRUE))

symbol_abundance_table.2 <- as.data.frame(table(as.factor((unlist(strsplit(as.character((abundance.2$best_hit_annotation)),'; '))))))
symbol_abundance_table.2 <- as.data.frame(table(as.factor((unlist(strsplit(as.character((abundance.2$best_hit_annotation)),', '))))))
symbol_abundance_table.2$Freq <- as.numeric(symbol_abundance_table.2$Freq)
symbol_abundance_table.2 <- symbol_abundance_table.2[order(-symbol_abundance_table.2$Freq),]


write.table(symbol_abundance_table.1, file="symbol_abundance_table.1.tsv", quote=FALSE, sep='\t')
write.table(symbol_abundance_table.2, file="symbol_abundance_table.2.tsv", quote=FALSE, sep='\t')


png(filename=paste0("hclust_n2_top25annotations.png"), width = 850, height = 1100)

par(mfrow = c(2,2))
par(mar = c(4, 25, 4, 2))

barplot(rev(head(symbol_abundance_table.2$Freq, n = 25)), names.arg = gsub('PF','\nPF',rev(head(symbol_abundance_table.2$Var1, n=25))),horiz = TRUE,las=2,main="Cluster 2\nwith\nHypothetical Protein",col.main="navyblue",cex.names=0.8,col=c(replicate(24,"gray"),"red"))
barplot(rev(head(symbol_abundance_table.2$Freq, n = 26))[1:25], names.arg = gsub('PF','\nPF',rev(head(symbol_abundance_table.2$Var1, n=26)))[1:25],horiz = TRUE,las=2,main="Cluster 2\nwithout\nHypothetical Protein",col.main="navyblue",cex.names=0.8)

barplot(rev(head(symbol_abundance_table.1$Freq, n = 25)), names.arg = gsub('PF','\nPF',rev(head(symbol_abundance_table.1$Var1, n=25))),horiz = TRUE,las=2,main="Cluster 1\nwith\nHypothetical Protein",col.main="tomato4",cex.names=0.8,col=c(replicate(24,"gray"),"red"))
barplot(rev(head(symbol_abundance_table.1$Freq, n = 26))[1:25], names.arg = gsub('PF','\nPF',rev(head(symbol_abundance_table.1$Var1, n=26)))[1:25],horiz = TRUE,las=2,main="Cluster 1\nwithout\nHypothetical Protein",col.main="tomato4",cex.names=0.8)

dev.off()

# Return to ORF matrix, remove all ORFs with ARO_Category = "efflux pump*" 



symbol_abundance_table.1limited <- defined.rgi.rel$ARO_name[grepl("fflux", defined.rgi.rel$ARO_category)]
symbol_abundance_table.1limited.sat <- d



symbol_abundance_table.1limited <- as.data.frame(table(as.factor((unlist(strsplit(as.character((symbol_abundance_table.1limited$best_hit_annotation)),'; '))))))
symbol_abundance_table.1limited <- as.data.frame(table(as.factor((unlist(strsplit(as.character((symbol_abundance_table.1limited$best_hit_annotation)),', '))))))
symbol_abundance_table.1limited$Freq <- as.numeric(symbol_abundance_table.1limited$Freq)


















png(filename=paste0("hclust_n2_top25dendro.png"), width = 900, height = 1100)
par(mfrow = c(1,1))
par(mar = c(2, 2, 2, 2))
plot( as.phylo(as.hclust(h$rowDendrogram)), tip.colo=c("tomato4","navyblue")[cutree( as.hclust(h$rowDendrogram), 2 ) ] )
dev.off()
 
alluniq <- unique(sort(c(as.character(symbol_abundance_table.1$Var1),as.character(symbol_abundance_table.2$Var1))))
dff <- as.data.frame(matrix(nrow=length(alluniq),ncol=2))
rownames(dff) <- alluniq
ai <- 0
for(a in alluniq){
  ai<-ai+1
  if (ai %% 100 == 0){print(ai)}
  if (a %in% as.character(symbol_abundance_table.1$Var1)){
    dff[a,1] <- as.numeric(symbol_abundance_table.1[as.character(symbol_abundance_table.1$Var1) == a,]$Freq)
  }else{
    dff[a,1] <- 0
  }
  if (a %in% as.character(symbol_abundance_table.2$Var1)){
    dff[a,2] <- as.numeric(symbol_abundance_table.2[as.character(symbol_abundance_table.2$Var1) == a,]$Freq)
  }else{
    dff[a,2] <- 0
  }
}

# ABC
length(intersect(symbol_abundance_table.1$Var1, symbol_abundance_table.2$Var1))
58574

# symbol_abundance_table.1)
68596     
# 2
66318 
install.packages("Vennerable", repos="http://R-Forge.R-project.org")

