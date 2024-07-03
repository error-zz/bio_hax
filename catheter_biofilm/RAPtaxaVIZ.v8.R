library(colorRamps)
library(vegan)
library(ape)
library(gtools)
library(fpc)
library(cluster)
library(dplyr)
library(ade4)
library(tidyr)
library(ggplot2)
library(MASS)

setwd ("/Users/apple/Desktop/DEV.gcid/genus")

best_hit_species <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/annotation_20170710.filtered.summary_best_hit_species.tab", sep="\t", header=TRUE)
rownames(best_hit_species) <- best_hit_species$best_hit_species

count_ABS_table <- best_hit_species[,4:16]
rownames(count_ABS_table) <- best_hit_species$best_hit_species

count_REL_table <- count_ABS_table
colsums <- colSums(best_hit_species[,4:16])
r_i <- 1
for(r in best_hit_species[,4:16]){
  print(r)
  rel <- r/colsums[r_i]
  count_REL_table[,r_i] <- rel
  r_i <- r_i + 1
}
rownames(count_REL_table) <- best_hit_species$best_hit_species
colSums(count_REL_table)

#samp_assign <- c("x29UC",                    "x68UC",               "x33UC",                     "x43UC",                     "x68UC",                     "x59UC",             "x68UC",                               "x68UC",             "x73UC",             "x29UC",              "x29UC",                "x33UC",               "x43UC");
samp_assign <- c("x29UC",                    "x68UC",               "x33UC",                     "x43UC",                     "x68UC",                     "x59UC",                                                     "x68UC",             "x73UC",             "x29UC",              "x29UC",                "x33UC",               "x43UC");
#prep_assign <- c("BloodAgar",                "META",                "BloodAgar",                 "BloodAgar",                 "BloodAgar",                 "ORI",               "ORI",                                 "ORI",               "ORI",               "ORI",                "META",                "META",                 "META") 
prep_assign <- c("BloodAgar",                "META",                "BloodAgar",                 "BloodAgar",                 "BloodAgar",                 "ORI",                                                      "ORI",               "ORI",               "ORI",                "META",                "META",                 "META") 

best_hit_species <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/annotation_20170710.filtered.summary_best_hit_species.tab", sep="\t", header=TRUE)

# v4 - remove NebNext
best_hit_species <- best_hit_species[,colnames(best_hit_species)!="X67.68UC.ORI.NEBNEXT.PE.P51.1_S13_001"]
inlist <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/inlist.txt", sep="\t", header=TRUE)
inlist.v <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/inlist_viruses.txt", sep="\t", header=TRUE)

inlist <- unname(unlist(inlist)[1:136])
inlist.v <- unname(unlist(inlist.v))

species_counts <- best_hit_species[4:16]
colnames(species_counts) <- colnames(best_hit_species)[4:16]
rownames(species_counts) <- best_hit_species$best_hit_species

species_counts <- species_counts[,colnames(species_counts)!="X67.68UC.ORI.NEBNEXT.PE.P51.1_S13_001"]
species_counts.bact <- species_counts[unique(grep( paste(inlist,collapse='|'), rownames(species_counts), value=TRUE )), ]
r_c <- 1
genus_id <- list()
for(r in rownames(species_counts.bact)){
  genus_id[r_c] <- strsplit(as.character(rownames(species_counts.bact)),' ')[[r_c]][1]
  r_c <- r_c + 1
}
#unique_genus_id <- sort(unique(unlist(genus_id)))
unique_genus_id <- unique(unlist(genus_id))
all_cols <- c("yellowgreen","red","turquoise4","goldenrod1","slateblue4","lawngreen","royalblue1","orangered4","purple1","navyblue","mediumblue","lemonchiffon3","indianred1","lightblue","sienna1","steelblue4","lightcyan",'lightgreen','pink','yellow',"khaki","violetred2","honeydew","orangered2","navy","snow3")
 
species_counts.bact.ORIGINAL <- species_counts.bact

#####################################################

species_counts.bact.29 <- species_counts.bact[,c("X69.29UC.ORI._S12_001","X1.29UCBLOODAGAR.META_S1_001","X7.29UCLB.META_S5_001")]
colnames(species_counts.bact.29) <- c("ORI","BLOOD_AGAR","LB_AGAR")
species_counts.bact.68 <- species_counts.bact[,c("X67.68UC.ORI_S10_001","X4.68UCBLOODAGAR.META_S4_001","X10.68UCLB.META_S8_001")]
colnames(species_counts.bact.68) <- c("ORI","BLOOD_AGAR","LB_AGAR")
species_counts.bact.33 <- species_counts.bact[,c("X2.33UCBLOODAGAR.META_S2_001","X8.33UCLB.META_S6_001")]
colnames(species_counts.bact.33) <- c("BLOOD_AGAR","LB_AGAR")
species_counts.bact.43 <- species_counts.bact[,c("X3.43UCBLOODAGAR.META_S3_001","X9.43UCLB.META_S7_001")]
colnames(species_counts.bact.43) <- c("BLOOD_AGAR","LB_AGAR")
species_counts.bact.59 <- as.data.frame(as.matrix(species_counts.bact[,"X66.59UC.ORI_S9_001"]))
rownames(species_counts.bact.59) <- rownames(species_counts.bact)
colnames(species_counts.bact.59) <- "ORI"
species_counts.bact.73 <- as.data.frame(as.matrix(species_counts.bact[,"X68.73UC.ORI_S11_001"]))
rownames(species_counts.bact.73) <- rownames(species_counts.bact)
colnames(species_counts.bact.73) <- "ORI"

r_c <- 1
genus_id <- list()
for(r in rownames(species_counts.bact)){
  genus_id[r_c] <- strsplit(as.character(rownames(species_counts.bact)),' ')[[r_c]][1]
  r_c <- r_c + 1
}
unique_genus_id.all <- unique(unlist(genus_id))

genus.df.all <- as.data.frame(matrix(ncol=12,nrow=length(unique_genus_id.all)))
colnames(genus.df.all) <- c("29UC_ORI","29UC_BA","29UC_LB","68UC_ORI","68UC_BA","68UC_LB","33UC_BA","33UC_LB","43UC_BA","43UC_LB","59UC_ORI","73UC_ORI")
#colnames(genus.df.all) <- 
rownames(genus.df.all) <- unique_genus_id.all
for(u in unique_genus_id.all){
  
  genus.df.all[u,"29UC_ORI"] <- sum(species_counts.bact.29[grepl(u,rownames(species_counts.bact.29)),]$ORI)
  genus.df.all[u,"29UC_BA"]  <- sum(species_counts.bact.29[grepl(u,rownames(species_counts.bact.29)),]$BLOOD_AGAR)
  genus.df.all[u,"29UC_LB"]  <- sum(species_counts.bact.29[grepl(u,rownames(species_counts.bact.29)),]$LB_AGAR)
  
  genus.df.all[u,"68UC_ORI"] <- sum(species_counts.bact.68[grepl(u,rownames(species_counts.bact.68)),]$ORI)
  genus.df.all[u,"68UC_BA"]  <- sum(species_counts.bact.68[grepl(u,rownames(species_counts.bact.68)),]$BLOOD_AGAR)
  genus.df.all[u,"68UC_LB"]  <- sum(species_counts.bact.68[grepl(u,rownames(species_counts.bact.68)),]$LB_AGAR)
  
  genus.df.all[u,"33UC_BA"]  <- sum(species_counts.bact.33[grepl(u,rownames(species_counts.bact.33)),]$BLOOD_AGAR)
  genus.df.all[u,"33UC_LB"]  <- sum(species_counts.bact.33[grepl(u,rownames(species_counts.bact.33)),]$LB_AGAR)
  
  genus.df.all[u,"43UC_BA"]  <- sum(species_counts.bact.43[grepl(u,rownames(species_counts.bact.43)),]$BLOOD_AGAR)
  genus.df.all[u,"43UC_LB"]  <- sum(species_counts.bact.43[grepl(u,rownames(species_counts.bact.43)),]$LB_AGAR)
  
  genus.df.all[u,"59UC_ORI"]  <- sum(as.matrix(species_counts.bact.59)[grepl(u,rownames(species_counts.bact.59))])
  
  genus.df.all[u,"73UC_ORI"]  <- sum(as.matrix(species_counts.bact.73)[grepl(u,rownames(species_counts.bact.73))])
  
}

par(mfrow = c(1,1))
par(mar = c(10, 4, 8, 2))
barplot(as.matrix(genus.df.all),space = c(0.2,0.2,0.2,0.8,0.2,0.2,0.8,0.2,0.8,0.2,0.8,0.8), col=all_cols[as.factor(rownames(genus.df.all))], horiz = FALSE, las=3)

par(mfrow = c(1,1))
par(mar = c(20, 8, 8, 11))
genus.df.all.rel <- t(decostand(t(as.matrix(genus.df.all)),method="total"))
barplot(as.matrix(genus.df.all.rel),space = c(0.2,0.2,0.2,0.8,0.2,0.2,0.8,0.2,0.8,0.2,0.8,0.8), col=all_cols[as.factor(rownames(genus.df.all))], horiz = FALSE, las=3)


png(filename=paste0("abundance.A_all.png"), width = 1054, height = 1769)
genus.df.all.rel.DerrickSort <- genus.df.all.rel[,rev(c("73UC_ORI","68UC_ORI","59UC_ORI","29UC_ORI","68UC_BA","43UC_BA","33UC_BA","29UC_BA","68UC_LB","43UC_LB","33UC_LB","29UC_LB"))]
order.samp <- c("29UC","33UC","43UC",       "68UC",
                "29UC","33UC","43UC",       "68UC",
                "29UC",              "59UC","68UC","73UC")
order.subj <- c("6",        "2","8","8",
                "6","2","2",    "8",
                "6","2","2",    "8")
order.date <- c("0",        "63","0","91",
                "0","0","21",    "0",
                "0","0","21",    "0")
par(mfrow = c(1,1))
par(mar = c(20, 8, 8, 11))
barplot(as.matrix(genus.df.all.rel.DerrickSort)*100,space = rev(c(0.2,0.2,0.2,0.8,0.2,0.2,0.2,0.8,0.2,0.2,0.2,0.2)), col=all_cols[as.factor(rownames(genus.df.all))],las=2,
        main="Genus Level Relative Abundance\n(D) Without Unknown and Dominant Taxa",names.arg=paste(order.date,"      ",order.subj,"      ",order.samp," "),cex.names = 2)
par(new=TRUE)
par(mar = c(2, 2, 2, 2))
legendtxt <- rownames(genus.df.all.rel)
legend("topright", inset=c(-0.65,0),legend=legendtxt,pch=22,pt.bg=all_cols[as.factor(rownames(genus.df.all))],cex=1.2,border="white", bty = "n", fill=NULL, bg="white")
legend("bottomleft", inset=c(-0.13,-0.065), legend=c("\n\n\n\n\n\nSample ID\n\n\n\n\nSubject ID\n\n\n\nDate from\nFirst Sample"),pch=22,pt.bg="white",cex=1.4,border="white", bty = "n", fill=NULL, bg="white")
legend("bottomleft", inset=c(0,-0.11), legend=c("                 ---------                  ---------                            ---------"),pch=22,pt.bg="white",cex=2,col="white",border="white", bty = "n", fill=NULL, bg="white")
legend("bottomleft", inset=c(0,-0.02), legend=c("-------- ORI --------        -- BLOOD AGAR --       ------ LB AGAR ------"),pch=22,pt.bg="white",cex=2,col="white",,border="white", bty = "n", fill=NULL, bg="white")
dev.off()

write.table(genus.df.all.rel.DerrickSort, file="abundance.A_all.tsv", quote=FALSE, sep='\t')


plot.new()
legend("center", inset=c(-0.65,0),legend=legendtxt,pch=22,pt.bg=all_cols[as.factor(rownames(genus.df.all))],cex=1.2,border="white", bty = "n", fill=NULL, bg="white")
 

png(filename=paste0("abundance.B_reduced.png"), width = 1054, height = 1769)
par(mfrow = c(1,1))
par(mar = c(20, 8, 8, 11))
plotit <- t(decostand(t(as.matrix(genus.df.all.rel.DerrickSort[!rownames(genus.df.all.rel.DerrickSort) %in% c("Klebsiella","Enterococcus","Proteus","Staphylococcus","Veillonella"),])),method="total"))
barplot(plotit,space = rev(c(0.2,0.2,0.2,0.8,0.2,0.2,0.2,0.8,0.2,0.2,0.2,0.2)), col=all_cols[as.factor(rownames(genus.df.all))][6:26],las=2,
  main="Genus Level Relative Abundance\n(D) Without Unknown and Dominant Taxa",names.arg=paste(order.date,"      ",order.subj,"      ",order.samp," "),cex.names = 2)
par(new=TRUE)
par(mar = c(2, 2, 2, 2))
legendtxt <- rownames(genus.df.all.rel)
legend("topright", inset=c(-0.65,0),legend=legendtxt,pch=22,pt.bg=all_cols[as.factor(rownames(genus.df.all))][6:26],cex=1.2,border="white", bty = "n", fill=NULL, bg="white")
legend("bottomleft", inset=c(-0.13,-0.065), legend=c("\n\n\n\n\n\nSample ID\n\n\n\n\nSubject ID\n\n\n\nDate from\nFirst Sample"),pch=22,pt.bg="white",cex=1.4,border="white", bty = "n", fill=NULL, bg="white")
legend("bottomleft", inset=c(0,-0.11), legend=c("                 ---------                  ---------                            ---------"),pch=22,pt.bg="white",cex=2,col="white",border="white", bty = "n", fill=NULL, bg="white")
legend("bottomleft", inset=c(0,-0.02), legend=c("-------- ORI --------        -- BLOOD AGAR --       ------ LB AGAR ------"),pch=22,pt.bg="white",cex=2,col="white",,border="white", bty = "n", fill=NULL, bg="white")
dev.off()
write.table(plotit, file="abundance.B_all.tsv", quote=FALSE, sep='\t')


plot.new()
legend("center", inset=c(-0.65,0),legend=legendtxt[6:26],pch=22,pt.bg=all_cols[as.factor(rownames(genus.df.all))][6:26],cex=1.2,border="white", bty = "n", fill=NULL, bg="white")

par(mar = c(20, 8, 8, 11))

barplot(as.matrix(genus.df.all[!rownames(genus.df.all.rel.DerrickSort) %in% c("Klebsiella","Enterococcus","Proteus","Staphylococcus","Veillonella"),]),space = c(0.2,0.2,0.2,0.8,0.2,0.2,0.8,0.2,0.8,0.2,0.8,0.8), col=all_cols[as.factor(rownames(genus.df.all))][6:26], horiz = FALSE, las=3)
