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

best_hit_species <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/annotation_20170710.filtered.summary_best_hit_species.tab", sep="\t", header=TRUE)

png(filename=paste0("NU_00_AbsAbs_All.png"), width = 1000, height = 1000)
plot.new()
par(fig=c(0.0, 1.0, 0.15, 1.0), new = TRUE)
boxplot(best_hit_species[,4:end],xaxt='n',col=c(colorRamps::matlab.like2(length(unique(samp_assign))))[as.factor(samp_assign)], varwidth=TRUE,notch=TRUE,outline=FALSE, margins = c(0, 0), range=2, ylab = "Absolute Abundance Within Sample", main="SPECIES ABUNDANCE (OVERALL)")
# col=unlist(full_species_color_assigned)
par(new=TRUE)
label_list <- colnames(log_best_hit_species)[4:end]
axis(1, at=1:length(label_list), labels=label_list, las=3, cex.axis=0.8)
legend("topleft",c(unique(samp_assign)),border="white", bty = "n",pch=22,pt.bg=c(colorRamps::matlab.like2(length(unique(samp_assign)))))
dev.off()

####################################################################################################################################################################################
################### SUMMARIZE BY SPECIES ###########################################################################################################################################
####################################################################################################################################################################################

# conver to rel abundance
end <- dim(best_hit_species)[2] - 2
best_hit_species[,4:end] <- decostand(best_hit_species[,4:end],"total",2,na.rm=TRUE)

# and log
log_best_hit_species <- best_hit_species
log_best_hit_species[4:end] <- log(best_hit_species[4:end])

# and assign type
#            c(1-29UCBLOODAGAR-META_S1_001	10-68UCLB-META_S8_001	2-33UCBLOODAGAR-META_S2_001	3-43UCBLOODAGAR-META_S3_001	4-68UCBLOODAGAR-META_S4_001	66-59UC-ORI_S9_001	67-68UC-ORI-NEBNEXT-PE-P51-1_S13_001	67-68UC-ORI_S10_001	68-73UC-ORI_S11_001	69-29UC-ORI-_S12_001	7-29UCLB-META_S5_001	8-33UCLB-META_S6_001	9-43UCLB-META_S7_001)

# full prep assignment
#prep_assign <- c("Blood_agar",               "LB_agar",            "Blood_agar",              "Blood_agar",               "Blood_agar",               "ori_biofilm",      "ori_biofilm_nebnext",                "ori_biofilm",      "ori_biofilm",      "ori_biofilm",       "LB_agar",             "LB_agar",            "LB_agar");
# sample_id
samp_assign <- c("x29UC",                    "x68UC",               "x33UC",                     "x43UC",                     "x68UC",                     "x59UC",             "x68UC",                               "x68UC",             "x73UC",             "x29UC",              "x29UC",                "x33UC",               "x43UC");

# JAMISON'S PREP ASSIGNMENT
#samp_assign <- c("BloodAgar",               "META",               "BloodAgar",                "BloodAgar",                "BloodAgar",                "ORI",              "ORI-NEBNEXT",                        "ORI",              "ORI",              "ORI",                "META",               "META",               "META") 
#prep_assign <- c("BloodAgar",               "META",               "BloodAgar",                "BloodAgar",                "BloodAgar",                "ORI",              "ORI-NEBNEXT",                        "ORI",              "ORI",              "ORI",                "META",               "META",               "META") 
prep_assign <- c("BloodAgar",               "META",               "BloodAgar",                "BloodAgar",                "BloodAgar",                "ORI",              "ORI",                        "ORI",              "ORI",              "ORI",                "META",               "META",               "META") 

id_assign   <- c(1,                         10,                   2,                          3,                          4,                          66,                 67,                                   67,                 68,                  69,                   7,                   8,                    9)
is.na(log_best_hit_species) <- sapply(log_best_hit_species, is.infinite)

# LOG

plot.new()
par(fig=c(0.0, 1.0, 0.15, 1.0), new = TRUE)
boxplot(-log_best_hit_species[,4:end],xaxt='n',col=c(colorRamps::matlab.like2(length(unique(samp_assign))))[as.factor(samp_assign)],varwidth=TRUE,notch=TRUE,outline=FALSE, margins = c(0, 0), range=2, ylab = "-log(relativeAbundanceWithinSample)", main="LOG(SPECIES ABUNDANCE)")
# col=unlist(full_species_color_assigned)
par(new=TRUE)
label_list <- colnames(log_best_hit_species)[4:end]
axis(1, at=1:length(label_list), labels=label_list, las=3, cex.axis=0.8)
legend("topleft",c(unique(samp_assign)),border="white", bty = "n",pch=22,pt.bg=c(colorRamps::matlab.like2(length(unique(samp_assign)))))
dev.off()


####################################################################################################################################################################################
################### SUMMARIZE BY TAXA ##############################################################################################################################################
####################################################################################################################################################################################


plot.new()
par(fig=c(0.0, 1.0, 0.15, 1.0), new = TRUE)

log_best_hit_taxa <- t(log_best_hit_species)[4:end,]
species_id_vector <- log_best_hit_species$best_hit_species
log_best_hit_taxa <- matrix(mapply(log_best_hit_taxa, FUN=as.numeric), nrow=dim(log_best_hit_taxa)[1], ncol=dim(log_best_hit_taxa)[2])
log_best_hit_taxa[is.na(log_best_hit_taxa)] <- 0
rownames(log_best_hit_taxa) <- colnames(best_hit_species)[4:end]
colnames(log_best_hit_taxa) <- best_hit_species$best_hit_species

# boxplot(log_best_hit_taxa,xaxt='n',varwidth=TRUE,notch=TRUE,outline=FALSE, margins = c(0, 0), range=2, ylab = "log(relativeAbundanceAcrossAllSamples)", main="LOG(UNIQUE ORF ABUNDANCE WITHIN ALL SAMPLES)")
#### col=unlist(full_species_color_assigned)
# par(new=TRUE)
label_list <- colnames(log_best_hit_taxa)
# axis(1, at=1:length(label_list), labels=label_list, las=3, cex.axis=0.8)

# Reduce to most abundant (all assigned)
genus_count <- 0
species_count <- 0
subpecies_count <- 0
unique_genus_list <- list()
unique_species_list <- list()
unique_subpecies_list <- list()
for (taxa_id in colnames(log_best_hit_taxa)){
  cur_taxa <- as.list(unlist(strsplit(taxa_id, ' ')))
  print(cur_taxa[1])
  if(!(cur_taxa[1] %in% unique_genus_list)){
    unique_genus_list[genus_count] <- cur_taxa[1]
    genus_count <- genus_count + 1
  }
  if(!(paste(cur_taxa[1]," ",cur_taxa[2]) %in% unique_species_list)){
    unique_species_list[species_count] <- paste(cur_taxa[1]," ",cur_taxa[2])
    species_count <- species_count + 1
  }
  if(!(taxa_id %in% unique_subpecies_list)){
    unique_subpecies_list[subpecies_count] <- taxa_id
    subpecies_count <- subpecies_count + 1
  }
}

unique_genus_list <- unique(unique_genus_list)
unique_species_list <- unique(unique_species_list)
unique_subpecies_list <- unique(unique_subpecies_list)

genus_count_table <- matrix(nrow=dim(log_best_hit_taxa)[1],ncol=length(unique_genus_list))
genus_sum <- list()
species_count_table <- matrix(nrow=dim(log_best_hit_taxa)[1],ncol=length(unique_species_list))
species_sum <- list()
subspecies_count_table <- matrix(nrow=dim(log_best_hit_taxa)[1],ncol=length(unique_subpecies_list))
subspecies_sum <- list()

##########################################################################################################################################################################################
# MOST PREVALENT GENUS ###################################################################################################################################################################
##########################################################################################################################################################################################

gi <- 1
for(g in unique_genus_list){
  genus_sum[gi] <- sum(log_best_hit_taxa[,grepl(g,colnames(log_best_hit_taxa))])
  gi <- gi + 1
}
names(genus_sum) <- unique_genus_list

genus_box_hash <- matrix(nrow=length(unique_genus_list),ncol=length(colnames(best_hit_species))-5)
rownames(genus_box_hash) <- unique_genus_list
colnames(genus_box_hash) <- colnames(best_hit_species)[4:end]

mean_of_means <- list()
g_i <- 1
for (genus in unique_genus_list){
  
  genus_sub_table <- log_best_hit_taxa[,grepl(genus,colnames(log_best_hit_taxa))]
  means <- list()
  m_i <- 1
  if(length(rownames(genus_sub_table))>0){
    for(r in rownames(genus_sub_table)){
      if( sum(genus_sub_table[r,][ as.numeric(genus_sub_table[r,]) != 0 ] > 0)  ){
        means[m_i] <- 0
      }else{
        means[m_i] <- mean(genus_sub_table[r,][ as.numeric(genus_sub_table[r,]) != 0 ])
      }
      m_i <- m_i + 1
    }
    mean_of_means[g_i] <- mean(as.numeric(means),na.rm=TRUE)
  }else{
    mean_of_means[g_i] <- 0
    means <- rep(0, length(names(genus_sub_table)))
  }
  
  #for(subspecies in colnames(genus_sub_table)){
  genus_box_hash[genus,] <- unlist(as.numeric(means))
  #}
  
  g_i <- g_i + 1
}

genus_box_hash <- replace(genus_box_hash, is.na(genus_box_hash), 0)

# remove 0s
nz_genus_box_hash <- genus_box_hash[rowSums(genus_box_hash)!=0,]
nz_genus_box_hash <- nz_genus_box_hash[order(rowMeans(nz_genus_box_hash)),]

png(filename=paste0("NU_03_Genus_RelAbs_All.png"), width = 3200, height = 800)
par(fig=c(0.0, 1.0, 0.15, 1.0), new = TRUE)
boxplot(-t(nz_genus_box_hash),xaxt='n',varwidth=TRUE,notch=FALSE,outline=FALSE, margins = c(25, 0), col="yellow", range=2, ylab = "-log(relativeAbundance)\nAs Defined By Mean of all Subspecies\nSorted by Mean", main="LOG(GENUS subspecies ABUNDANCE over all SAMPLE)")
par(new=TRUE)
label_list <- rownames(nz_genus_box_hash)
axis(1, at=1:length(label_list), labels=label_list, las=3, cex.axis=0.8)
dev.off()

boxplot(-t(nz_genus_box_hash),xaxt='n',varwidth=TRUE,notch=FALSE,outline=FALSE, margins = c(25, 0), col="yellow", range=2, ylab = "-log(relativeAbundance)\nAs Defined By Mean of all Subspecies\nSorted by Mean", main="LOG(GENUS subspecies ABUNDANCE over all SAMPLE)")


png(filename=paste0("NU_04_Genus_x_Sample_RelAbs_All.png"), width = 1000, height = 1000)
boxplot(-nz_genus_box_hash,xaxt='n',col=c(colorRamps::matlab.like2(length(unique(samp_assign))))[as.factor(samp_assign)],varwidth=TRUE,notch=TRUE,outline=FALSE, margins = c(0, 0), range=2, ylab = "-log(relativeAbundanceWithinSample)", main="LOG(GENUS subspecies ABUNDANCE over all SAMPLE)")
par(new=TRUE)
label_list <- colnames(genus_box_hash)
axis(1, at=1:length(label_list), labels=label_list, las=3, cex.axis=0.8)
legend("topleft",c(unique(samp_assign)),border="white", bty = "n",pch=22,pt.bg=c(colorRamps::matlab.like2(length(unique(samp_assign)))))
dev.off()

write.table(-t(nz_genus_box_hash), file=paste("GENUS.negLogRelAbs.tsv", sep = ""), quote=FALSE, sep='\t')

top25 <- sapply(nz_genus_box_hash, as.numeric)

library(gplots)
png(filename=paste0("NU_05_GENUS.heatmap.png"), width = 1000, height = 1000)
plot.new()
test <- heatmap.2(nz_genus_box_hash)
B <- unique(samp_assign)
A <- samp_assign
nu_order <- list()
a_i <- 1
for (a in A){
  b_i <- 1
  for (b in B){
    if (a == b){  
      nu_order[a_i] <- b_i
    }
    b_i <- b_i + 1
  }
  a_i <- a_i + 1
}
heatmap.2(nz_genus_box_hash,margins=c(22,5),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ unlist(nu_order) ])
par(new=TRUE)
legend("bottomleft",c(unique(samp_assign)),border="white", bty = "n",pch=22,pt.bg=c(colorRamps::matlab.like2(length(unique(samp_assign)))))
dev.off()

####################################################################################################################################################################################
# MOST SIGNIF 25 ###################################################################################################################################################################
####################################################################################################################################################################################

pca.nz_genus_box_hash <- prcomp(t(nz_genus_box_hash), center = TRUE) 

covar.pvals <- data.frame(matrix(ncol=3))
#covar.model_pvals <- data.frame(matrix(ncol=1))
for(i in rownames(nz_genus_box_hash)){
  reg <- lm( as.numeric(nz_genus_box_hash[i,]) ~ as.numeric(pca.nz_genus_box_hash$x[,1]) + as.numeric(pca.nz_genus_box_hash$x[,2]), na.action=na.omit )
  
  if (t(coef(summary(reg)))[2,1] > 0 && is.nan(t(coef(summary(reg)))[2,1]) == 0){
    cpus.lm2 <- stepAIC(reg, trace = FALSE)
    
    fstat <- summary(reg)$fstatistic 
    #covar.model_pvals[i] <- unname(p f(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
    covar.pvals[i,] <- unname(t(as.matrix(coef(summary(reg))[,4])))
    
    covar.deviance <- data.frame() #data.frame(matrix(ncol=dim(t(cpus.lm2$anova))[2],nrow=length(rownames(qc))))
    covar.resid_deviance <- data.frame() #data.frame(matrix(ncol=dim(t(cpus.lm2$anova))[2],nrow=length(rownames(qc))))
    
    L <- dim(t(cpus.lm2$anova))[2]
    if (L > 1){
      for (r in 2:L){
        rr <- t(cpus.lm2$anova)[,r][1]
        covar.deviance[i,rr] <- as.numeric(t(cpus.lm2$anova)[,r][3])
        covar.resid_deviance[i,rr] <- as.numeric(t(cpus.lm2$anova)[,r][5])
      }
    }else{
      covar.deviance[i,rr] <- 0
      covar.resid_deviance[i,rr] <- 0
    }
    
  }else{
    
    #covar.model_pvals[i] <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
    covar.pvals[i,] <- matrix(0,1,3)
    
  }
}
covar.pvals.sorted <- covar.pvals[with(covar.pvals,order(X1)),]
covar.pvals.sorted <- covar.pvals.sorted[-1,]
top25_signifIDs <- rownames(covar.pvals.sorted)[1:25]
#top25_signifIDs[1] <- "Eubacterium"
#top25_signifIDs[2] <- "Clostridium"
#top25_signifIDs[3] <- "Unknown"
#covar.pvals
pca.nz_genus_box_hash <- prcomp(t(nz_genus_box_hash[top25_signifIDs,]), center = TRUE) 
biplot(pca.nz_genus_box_hash, main="Genus Frequencies (25 Most Significant Predictor of PC Orientation)", cex=0.8, pch=20, expand=7, xlim=c(-1.1,2.5), ylim=c(-1.2,1.7))





#############################################################################################################################################################################
# SUM TOP 25 ABUNDANCE BY PREP ##############################################################################################################################################
#############################################################################################################################################################################

#top25_signifIDs <- unique(top25_signifIDs)

species.counts <- best_hit_species$total_reads
genus_ids <- unique(unlist(genus_id))
countz <- data.frame(matrix(nrow=length(genus_ids),ncol=1))
rownames(countz) <- genus_ids
s_i <- 1
for (s in species.counts){
  countz[unlist(genus_id[s_i]),1] <- 0
  s_i <- s_i + 1
}
s_i <- 1
for (s in species.counts){
  countz[unlist(genus_id[s_i]),1] <- countz[unlist(genus_id[s_i]),1] + s
  s_i <- s_i + 1
}

#countz[order(countz),]

r_ids <- list()
#for(r in top25_signifIDs){
for(r in rev(rownames(countz)[order(countz)])){
  print(r)
  print (r_ids)
  r_ids <- c(r_ids, as.character(best_hit_species[grepl(r,unlist(best_hit_species$best_hit_species)),]$best_hit_species))
}

top25_signifIDs <- rev(rownames(countz)[order(countz)])

r_ids <- unique(r_ids)
rownames(best_hit_species) <- best_hit_species$best_hit_species
signif_abundance_table <- best_hit_species[unlist(unique(r_ids)),]
rownames(signif_abundance_table) <- signif_abundance_table$best_hit_species

r_c <- 1
genus_id <- list()
for(r in rownames(signif_abundance_table)){
  genus_id[r_c] <- strsplit(as.character(signif_abundance_table$best_hit_species),' ')[[r_c]][1]
  r_c <- r_c + 1
}

plot.new()
par(mfrow = c(5,1))
signif_abundance_table$genus <- genus_id
genus_abundance_table <- data.frame(matrix(nrow=3*length(top25_signifIDs),ncol=4))
listof <- list()
li <- 1
#for (g in unique(genus_id)){
row_ids <- list()

#top25_signifIDs[21] <- "Clostridium"

for (g in top25_signifIDs){
  
  #nu_entry <- list()
  #nu_entry[1] <- g
  tmp <- signif_abundance_table[genus_id %in% g,]
  
  tmp.blood <- colSums(tmp[4:16][,prep_assign=="BloodAgar"])/colSums(best_hit_species[,4:16][,prep_assign=="BloodAgar"])
  #nu_entry <- c(nu_entry, unname(tmp.blood))
  genus_abundance_table[(li*3)-2,] <- unname(tmp.blood)
  row_ids[(li*3)-2] <- paste0(g,".BloodAgar")

  tmp.meta <- colSums(tmp[4:16][,prep_assign=="META"])/colSums(best_hit_species[,4:16][,prep_assign=="META"])
  #nu_entry <- c(nu_entry, unname(tmp.meta))
  genus_abundance_table[(li*3)-1,] <- unname(tmp.meta)
  row_ids[(li*3)-1] <- paste0(g,".META")
  
  tmp.ori <- colSums(tmp[4:16][,prep_assign=="ORI"])/colSums(best_hit_species[,4:16][,prep_assign=="ORI"])
  genus_abundance_table[(li*3),] <- unname(tmp.meta)
  row_ids[(li*3)] <- paste0(g,".ORI")
  
  boxplot(tmp.blood,tmp.meta,tmp.ori,main=paste0(li," : ",g),col=c("red","green","blue"))
  legend("topright",legend=c("Blood Agar","Meta","Ori"),pch=22,pt.bg=c("red","green","blue"))
  li <- li + 1
  
}

png(filename=paste0("NU_Fig3_GenusAbundance_SideBySide.png"), width = 1400, height = 900)
plot.new()
par(mfrow = c(1,1))
label_ids <- row_ids 
label_ids[!grepl('BloodAgar',label_ids)] <- ""
par(mar = c(12, 2, 4, 2))
boxplot(-log(t(genus_abundance_table)),col=c("red","green","lightblue"),names=c(unlist(label_ids)),las=2,main="Top 25 (Max Total Abundance)\n-log(RelAbs)")
legend("topright",legend=c("Blood Agar","Meta","Ori"),pch=22,pt.bg=c("red","green","lightblue"))
dev.off()

###############
# sum by prep #
###############

genus_abundance_table.colSums <- rowSums(genus_abundance_table)
prep_dfx <- data.frame(matrix(nrow=length(top25_signifIDs),ncol=3))
c <- 1
for(r in 0:(length(top25_signifIDs)-1)){
  rr <- r + 1
  prep_dfx[rr,1] <- genus_abundance_table.colSums[c]
  c <- c + 1
  prep_dfx[rr,2] <- genus_abundance_table.colSums[c]
  c <- c + 1
  prep_dfx[rr,3] <- genus_abundance_table.colSums[c]
  c <- c + 1
}
#prep_dfx <- prep_dfx[-35,]
colnames(prep_dfx) <- c("BloodAgar","Meta","Ori")
rownames(prep_dfx) <- top25_signifIDs

png(filename=paste0("NU_F3_MaxGenusAbundance_neglogRelAbundance.png"), width = 1400, height = 900)
par(mfrow = c(1,1))
par(mar = c(15, 2, 4, 2))
barplot(as.matrix(t(-log(prep_dfx[-21,]))),las=2,main="Genus Level Match\n-Log(Relative Abundance)",col=c("red","green","lightblue"))
legend("topright",legend=c("Blood Agar","Meta","Ori"),pch=22,pt.bg=c("red","green","lightblue"))
dev.off()

#################################################################################################################################################################################################
# PATHWAY  ######################################################################################################################################################################################
#################################################################################################################################################################################################


#pathway_tab <- read.table("annotation_20170710.filtered.summary_KO_pathway.tab", sep="\t", header=TRUE)
pathway_tab <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/annotation_20170710.filtered.summary_KO_pathway_reformat.txt", sep="\t", header=TRUE)
tmp_names <- pathway_tab$KO_pathway
t <- length(colnames(pathway_tab))-2
pathway_tab <- pathway_tab[,2:t]
pathway_tab <- mapply(pathway_tab, FUN=as.numeric)
pathway_tab <- pathway_tab[complete.cases(pathway_tab),]
rownames(pathway_tab) <- tmp_names
  # unlist(tmp_names[1:168])
 
pathway_tab <- decostand(pathway_tab,"total",2,na.rm=TRUE)

png(filename=paste0("NU_06_KO_Pathway_Heatmap.png"), width = 1000, height = 3300)
plot.new()
test <- heatmap.2(pathway_tab)
B <- unique(samp_assign)
A <- samp_assign
nu_order <- list()
a_i <- 1
for (a in A){
  b_i <- 1
  for (b in B){
    if (a == b){  
      nu_order[a_i] <- b_i
    }
    b_i <- b_i + 1
  }
  a_i <- a_i + 1
}
heatmap.2(pathway_tab,margins=c(22,28),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ unlist(nu_order) ],cexRow = 1.1)
par(new=TRUE)
legend("bottomleft",c(unique(samp_assign)),border="white", bty = "n",pch=22,pt.bg=c(colorRamps::matlab.like2(length(unique(samp_assign)))))
dev.off()




#######################################################################################################################################################
# CLASS  ##############################################################################################################################################
#######################################################################################################################################################

class_tab <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/annotation_20170710.filtered.summary_KOG_class.tab", sep="\t", header=TRUE)
tmp_names <- class_tab$KOG_class
t <- length(colnames(class_tab))-2
class_tab <- class_tab[,2:t]
class_tab <- mapply(class_tab, FUN=as.numeric)
class_tab <- class_tab[complete.cases(class_tab),]
rownames(class_tab) <- tmp_names

class_tab <- class_tab[2:length(rownames(class_tab)),]
class_tab <- decostand(class_tab,"total",2,na.rm=TRUE)

png(filename=paste0("NU_07_KOG_Class_Heatmap.png"), width = 1000, height = 1000)
plot.new()
test <- heatmap.2(class_tab)
heatmap.2(class_tab,margins=c(22,28),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ unlist(nu_order) ],cexRow = 1.1)
par(new=TRUE)
legend("bottomleft",c(unique(samp_assign)),border="white", bty = "n",pch=22,pt.bg=c(colorRamps::matlab.like2(length(unique(samp_assign)))))
dev.off()




#################################################################################################################################################################################################
# AMR ###########################################################################################################################################################################################
#################################################################################################################################################################################################

amr_tab <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/amr_count_report.tsv", sep="\t", header=TRUE)
tmp_names <- amr_tab$orf_id
t <- length(colnames(amr_tab))
amr_tab <- amr_tab[,2:t]
amr_tab <- mapply(amr_tab, FUN=as.numeric)
amr_tab <- amr_tab[complete.cases(amr_tab),]
rownames(amr_tab) <- tmp_names

amr_tab <- amr_tab[2:length(rownames(amr_tab)),]
amr_tab <- decostand(amr_tab,"total",2,na.rm=TRUE)
amr_tab <- -log(amr_tab+0.0000001)

png(filename=paste0("NU_08a_Only_AMR_Transcripts_Heatmap_SQUARE.png"), width = 1000, height = 1000)
plot.new()
test <- heatmap.2(amr_tab)
heatmap.2(amr_tab,margins=c(22,5),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ as.factor(samp_assign) ],cexRow = 0.1)
par(new=TRUE)
legend("bottomleft",c(unique(samp_assign)),border="white", bty = "n",pch=22,pt.bg=c(colorRamps::matlab.like2(length(unique(samp_assign)))))
dev.off()

png(filename=paste0("08b_Only_AMR_Transcripts_Heatmap_LONG.png"), width = 1000, height = 4200)
plot.new()
heatmap.2(amr_tab,margins=c(22,28),ColSideColors=c(colorRamps::matlab.like2(length(unique(samp_assign))))[ as.factor(samp_assign) ],cexRow = 0.6,main="AMR Transcripts Only")
par(new=TRUE)
legend("bottomleft",c(unique(samp_assign)),border="white", bty = "n",pch=22,pt.bg=c(colorRamps::matlab.like2(length(unique(samp_assign)))))
dev.off()

amr.pca <- prcomp(t(amr_tab))

covar.pvals <- data.frame(matrix(ncol=3))
for(i in rownames(amr_tab)){
  reg <- lm( as.numeric(amr_tab[i,]) ~ as.numeric(amr.pca$x[,1]) + as.numeric(amr.pca$x[,2]), na.action=na.omit )
  
  if (t(coef(summary(reg)))[2,1] > 0 && is.nan(t(coef(summary(reg)))[2,1]) == 0){
    cpus.lm2 <- stepAIC(reg, trace = FALSE)
    
    fstat <- summary(reg)$fstatistic 
    #covar.model_pvals[i] <- unname(p f(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
    covar.pvals[i,] <- unname(t(as.matrix(coef(summary(reg))[,4])))
    
    covar.deviance <- data.frame() #data.frame(matrix(ncol=dim(t(cpus.lm2$anova))[2],nrow=length(rownames(qc))))
    covar.resid_deviance <- data.frame() #data.frame(matrix(ncol=dim(t(cpus.lm2$anova))[2],nrow=length(rownames(qc))))
    
    L <- dim(t(cpus.lm2$anova))[2]
    if (L > 1){
      for (r in 2:L){
        rr <- t(cpus.lm2$anova)[,r][1]
        covar.deviance[i,rr] <- as.numeric(t(cpus.lm2$anova)[,r][3])
        covar.resid_deviance[i,rr] <- as.numeric(t(cpus.lm2$anova)[,r][5])
      }
    }else{
      covar.deviance[i,rr] <- 0
      covar.resid_deviance[i,rr] <- 0
    }
    
  }else{
    
    #covar.model_pvals[i] <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
    covar.pvals[i,] <- matrix(0,1,3)
    
  }
}
covar.pvals.sorted <- covar.pvals[with(covar.pvals,order(X1)),]
top25_signifIDs <- rownames(covar.pvals.sorted)[1:25]

par(mfrow = c(3,1))
boxplot(t(amr_tab[top25_signifIDs,prep_assign=="BloodAgar"]),col="red",main="Blood Agar (Top 25 Significant AMR Genes)")
boxplot(t(amr_tab[top25_signifIDs,prep_assign=="META"]),col="green",main="META")
boxplot(t(amr_tab[top25_signifIDs,prep_assign=="ORI"]),col="blue",main="ORI")

dfx <- data.frame(matrix(nrow=5,ncol=75))
amr_labels <- list()
for(p in 0:24){
  dfx[,(3*p)+1] <- c(t(amr_tab[top25_signifIDs,prep_assign=="BloodAgar"])[,(p+1)],median(t(amr_tab[top25_signifIDs,prep_assign=="BloodAgar"])[,(p+1)]))
  dfx[,(3*p)+2] <- c(t(amr_tab[top25_signifIDs,prep_assign=="META"])[,(p+1)],median(t(amr_tab[top25_signifIDs,prep_assign=="META"])[,(p+1)]))
  dfx[,(3*p)+3] <- t(amr_tab[top25_signifIDs,prep_assign=="ORI"])[,(p+1)]
  amr_labels[(3*p)+1] <- colnames(t(amr_tab[top25_signifIDs,prep_assign=="BloodAgar"]))[p+1]
  amr_labels[(3*p)+2] <- ""
  amr_labels[(3*p)+3] <- ""
}

par(mar = c(17, 4, 4, 2))
boxplot(dfx,col=c("red","green","lightblue"),names=c(unlist(amr_labels)),las=2,main="Top 25 Significant AMR Genes\n-log(RelAbs)")

par(mar = c(17, 4, 4, 2))
boxplot(dfx,col=c("red","green","lightblue"),names=c(unlist(amr_labels)),las=2,main="Top 25 Significant AMR Genes\n-log(RelAbs)", ylim=c(6.5,9.5))
legend("topleft",c("Blood Agar","Meta","Ori"),pch=22,pt.bg=c("red","green","lightblue"))




dfx.colSums <- colSums(dfx)
prep_dfx <- data.frame(matrix(nrow=25,ncol=3))
c <- 1
for(r in 0:24){
  rr <- r + 1
  prep_dfx[rr,1] <- dfx.colSums[c]
  c <- c + 1
  prep_dfx[rr,2] <- dfx.colSums[c]
  c <- c + 1
  prep_dfx[rr,3] <- dfx.colSums[c]
  c <- c + 1
}
prep_dfx <- prep_dfx[-35,]
colnames(prep_dfx) <- c("BloodAgar","Meta","Ori")
rownames(prep_dfx) <- colnames(t(amr_tab[top25_signifIDs,]))

png(filename=paste0("NU_F1_AMRbyPrep_neglogRelAbundance.png"), width = 1400, height = 900)
par(mfrow = c(1,1))
par(mar = c(15, 2, 4, 2))
barplot(as.matrix(t(-log(prep_dfx))),las=2,main="Cloned Sequence\n-Log(Relative Abundance)",col=c("red","green","lightblue"))
legend("topright",legend=c("Blood Agar","Meta","Ori"),pch=22,pt.bg=c("red","green","lightblue"))
dev.off()

par(mfrow = c(1,1))
barplot(as.matrix(-log(prep_dfx)),las=2,main="Cloned Sequence\n-Log(Relative Abundance)",col=c(rainbow(25)))
legend("topright",legend=colnames(t(amr_tab[top25_signifIDs,])),pch=22,pt.bg=c(rainbow(25)))


#library(ggplot2)    
#ggplot(students) + geom_boxplot(aes(x = success, y = WAM))
#t(amr_tab[top25_signifIDs,prep_assign=="BloodAgar"]
  
par(mfrow = c(1,1))
biplot(prcomp(amr_tab[rownames(covar.pvals.sorted)[1:25],]), cex=0.8, pch=20, expand=7, xlim=c(-2.6,0.5), ylim=c(-1,1))

#################################################################################################################################################################################################
# ALIGN FREQS ###################################################################################################################################################################################
#################################################################################################################################################################################################

align_tab <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/aligncounts.txt", sep="\t", header=TRUE)
rownames(align_tab) <- align_tab$X
align_tab <- align_tab[,-1]

sampleID <- c(29,29,29,33,33,43,43,59,68,68,68,68,73)
prepID <-   c("Blood Agar","LB_agar","Ori_biofilm","Blood Agar","LB_agar","Blood Agar","LB_agar","Ori_biofilm","Blood Agar","LB_agar","Ori_biofilm","Ori_biofilm","Ori_biofilm")

plot.new()
par(mfrow = c(3,1))
boxplot(-log(align_tab[,prepID=="Blood Agar"]),main="Alignment to Reference Genomes (-log(relAbs)) \nBlood Agar",col="red",las=2)
boxplot(-log(align_tab[,prepID=="LB_agar"]),main="LB_agar",col="green",las=2)
boxplot(-log(align_tab[,prepID=="Ori_biofilm"]),main="Ori_biofilm",col="lightblue",las=2)

dfx <- data.frame(matrix(nrow=5,ncol=dim(align_tab_rel)[1]*3))
amr_labels <- list()
align_tab_rel <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/aligncounts_rel.txt", sep="\t", header=TRUE)
align_tab_rel <- align_tab_rel[-35,]
rownames(align_tab_rel) <- align_tab_rel$X
align_tab_rel <- align_tab_rel[,-1]
for(p in 0:(dim(align_tab_rel)[1]-1)){
  pp <- p + 1
  dfx[,(3*p)+1] <- unlist(c(align_tab_rel[,prep_assign=="BloodAgar"][pp,],median(as.numeric(unlist(align_tab_rel[,prep_assign=="BloodAgar"][pp,])))))
  dfx[,(3*p)+2] <- unlist(c(align_tab_rel[,prep_assign=="META"][pp,],median(as.numeric(unlist(align_tab_rel[,prep_assign=="META"][pp,])))))
  # c(t(amr_tab[top25_signifIDs,prep_assign=="META"])[,(p+1)],median(t(amr_tab[top25_signifIDs,prep_assign=="META"])[,(p+1)]))
  dfx[,(3*p)+3] <- t(align_tab_rel[,prep_assign=="ORI"][pp,])
  # t(amr_tab[top25_signifIDs,prep_assign=="ORI"])[,(p+1)]
  amr_labels[(3*p)+1] <- rownames(align_tab_rel[,prep_assign=="BloodAgar"])[pp]
  amr_labels[(3*p)+2] <- ""
  amr_labels[(3*p)+3] <- ""
}
png(filename=paste0("NU_F2B_CloneFreq_logAbsAbundance.png"), width = 1400, height = 900)
par(mfrow = c(1,1))
par(mar = c(15, 2, 4, 2))
boxplot(log(dfx),col=c("red","green","lightblue"),names=c(unlist(amr_labels)),las=2,main="Alignment Frequencies (log(AbsoluteAbundance))")
legend("topright",legend=c("Blood Agar","Meta","Ori"),pch=22,pt.bg=c("red","green","lightblue"))
dev.off()


###############
# sum by prep #
###############

dfx.colSums <- colSums(dfx)
prep_dfx <- data.frame(matrix(nrow=34,ncol=3))
c <- 1
for(r in 0:34){
  rr <- r + 1
  prep_dfx[rr,1] <- dfx.colSums[c]
  c <- c + 1
  prep_dfx[rr,2] <- dfx.colSums[c]
  c <- c + 1
  prep_dfx[rr,3] <- dfx.colSums[c]
  c <- c + 1
}
prep_dfx <- prep_dfx[-35,]
colnames(prep_dfx) <- c("BloodAgar","Meta","Ori")
rownames(prep_dfx) <- rownames(align_tab_rel)

png(filename=paste0("F2_CloneByPrep_logAbsAbundance.png"), width = 1400, height = 900)
par(mfrow = c(1,1))
par(mar = c(15, 2, 4, 2))
barplot(as.matrix(t(log(prep_dfx))),las=2,main="Cloned Sequence\nAbsolute Abundance",col=c("red","green","lightblue"))
legend("topright",legend=c("Blood Agar","Meta","Ori"),pch=22,pt.bg=c("red","green","lightblue"))
dev.off()

###############
###############

align_tab_rel <- align_tab_rel/colSums(align_tab_rel)
dfx <- data.frame(matrix(nrow=5,ncol=dim(align_tab_rel)[1]*3))
amr_labels <- list()
#align_tab_rel <- t(align_tab_rel)
for(p in 0:(dim(align_tab_rel)[1]-1)){
  pp <- p + 1
  dfx[,(3*p)+1] <- unlist(c(align_tab_rel[,prep_assign=="BloodAgar"][pp,],median(as.numeric(unlist(align_tab_rel[,prep_assign=="BloodAgar"][pp,])))))
  dfx[,(3*p)+2] <- unlist(c(align_tab_rel[,prep_assign=="META"][pp,],median(as.numeric(unlist(align_tab_rel[,prep_assign=="META"][pp,])))))
    # c(t(amr_tab[top25_signifIDs,prep_assign=="META"])[,(p+1)],median(t(amr_tab[top25_signifIDs,prep_assign=="META"])[,(p+1)]))
  dfx[,(3*p)+3] <- t(align_tab_rel[,prep_assign=="ORI"][pp,])
    # t(amr_tab[top25_signifIDs,prep_assign=="ORI"])[,(p+1)]
  amr_labels[(3*p)+1] <- rownames(align_tab_rel[,prep_assign=="BloodAgar"])[pp]
  amr_labels[(3*p)+2] <- ""
  amr_labels[(3*p)+3] <- ""
}
png(filename=paste0("F2A_CloneFreq_neglogRelAbundance.png"), width = 1400, height = 900)
par(mfrow = c(1,1))
par(mar = c(15, 2, 4, 2))
boxplot(-log(dfx),col=c("red","green","lightblue"),names=c(unlist(amr_labels)),las=2,main="Alignment Frequencies (-log(RelativeAbundance))")
legend("topright",legend=c("Blood Agar","Meta","Ori"),pch=22,pt.bg=c("red","green","lightblue"))
dev.off()

###############
# sum by prep #
###############

dfx.colSums <- colSums(dfx)
prep_dfx <- data.frame(matrix(nrow=34,ncol=3))
c <- 1
for(r in 0:34){
  rr <- r + 1
  prep_dfx[rr,1] <- dfx.colSums[c]
  c <- c + 1
  prep_dfx[rr,2] <- dfx.colSums[c]
  c <- c + 1
  prep_dfx[rr,3] <- dfx.colSums[c]
  c <- c + 1
}
prep_dfx <- prep_dfx[-35,]
colnames(prep_dfx) <- c("BloodAgar","Meta","Ori")
rownames(prep_dfx) <- rownames(align_tab_rel)

png(filename=paste0("F2_CloneByPrep_neglogRelAbundance.png"), width = 1400, height = 900)
par(mfrow = c(1,1))
par(mar = c(15, 2, 4, 2))
barplot(as.matrix(t(-log(prep_dfx))),las=2,main="Cloned Sequence\n-Log(Relative Abundance)",col=c("red","green","lightblue"))
legend("topright",legend=c("Blood Agar","Meta","Ori"),pch=22,pt.bg=c("red","green","lightblue"))
dev.off()
