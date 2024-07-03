library(gplots)
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
library(broom)
library(MASS)

setwd ("/Users/apple/Desktop/DEV.gcid")

##############################################
# 0. curate output table ########################
##############################################

orfhits.t <- read.table("format2.assembly.orf.strict.txt", sep="\t", header=TRUE)

delivery_csv.t <- orfhits.t[,c("Best_Hit_evalue","Best_Hit_ARO","ARO","ARO_name","SNP","Best_Hit_ARO_category","ARO_category")]

rownames(delivery_csv.t) <-orfhits.t$ORF_ID

write.table(delivery_csv.t, file=paste0("ARO_summary.tsv"), quote=FALSE, sep='\t', col.names=TRUE)

########################################################################################################################################################################################################
# 1. Prepare input tables
########################################################################################################################################################################################################

best_hit_species <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/annotation_20170710.filtered.summary_best_hit_species.tab", sep="\t", header=TRUE)
rownames(best_hit_species) <- best_hit_species$best_hit_species

count_ABS_table <- best_hit_species[,4:16]
rownames(count_ABS_table) <- best_hit_species$best_hit_species

#count_REL_table <- best_hit_species[,4:16]/as.matrix(colSums(best_hit_species[,4:16]))
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

samp_assign <- c("x29UC",                    "x68UC",               "x33UC",                     "x43UC",                     "x68UC",                     "x59UC",             "x68UC",                               "x68UC",             "x73UC",             "x29UC",              "x29UC",                "x33UC",               "x43UC");
prep_assign <- c("BloodAgar",               "META",               "BloodAgar",                "BloodAgar",                "BloodAgar",                "ORI",              "ORI",                        "ORI",              "ORI",              "ORI",                "META",               "META",               "META") 

align_tab <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/aligncounts.txt", sep="\t", header=TRUE)
rownames(align_tab) <- align_tab$X
align_tab <- align_tab[,-1]


# ##############################################
# # CLC-ASSOCIATED AMR GENES #######################
# ##############################################
# 
# amr_tab <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/amr_count_report.tsv", sep="\t", header=TRUE)
# tmp_names <- amr_tab$orf_id
# t <- length(colnames(amr_tab))
# amr_tab <- amr_tab[,2:t]
# amr_tab <- mapply(amr_tab, FUN=as.numeric)
# amr_tab <- amr_tab[complete.cases(amr_tab),]
# rownames(amr_tab) <- tmp_names
# 
# library(colorRamps)
# 
# # devel viz
# amr_tab.pca <- decostand(amr_tab,"total",2,na.rm=TRUE)
# #plot(t(log(amr_tab.pca)), main="amr counts \nall ORFs : log(abs abundance)")
# 
# png(filename=paste0("CLCaligned_all.png"), width = 1840, height = 3680)
# heatmap.2(log(amr_tab+0.0000001), col=colorRamps::green2red, main="\namr counts \nORFs found with prior CLC alignments to AMR genes\nlog(abs abundance)", margins=c(12,16),cexCol = 0.6, trace='n')
# dev.off()
# 
# # save.image(file="LOST.delivery_csv.df.sort.ARO_name.Rdata")
# #load("LOST.delivery_csv.df.sort.ARO_name.Rdata")
# # 
# # # count
# # delivery.df <- as.data.frame(delivery_csv.t)
# # ARO_name <- delivery.df[order(delivery.df$ARO_name),]
# # 
# # sum(rownames(amr_tab) %in% rownames(delivery_csv.df.sort.ARO_name))
# # 
# # # subset to 12 defined amr genes
# # defined <- delivery_csv.df.sort.ARO_name[rownames(delivery_csv.df.sort.ARO_name) %in% rownames(amr_tab),]
# # 
# # png(filename=paste0("CLCaligned_and_RGIamr_defined.ARO_group_sorted.png"), width = 1840, height = 1840)
# # heatmap.2(log(amr_tab[rownames(defined),]+0.0000001), col=colorRamps::green2red, RowSideColors = colorRamps::matlab.like(length(unique(defined$ARO_name)))[factor(defined$ARO_name)], Rowv=colorRamps::matlab.like(length(unique(defined$ARO_name)))[factor(defined$ARO_name)], main="\n\nARO_name\namr counts \12 AMR-associated ORFs\namong all ORFs found with prior CLC alignments\nlog(abs abundance)", margins=c(12,16),cexCol = 0.6, trace='n')
# # legend("left",legend=unique(defined$ARO_name),pch=22,pt.bg=colorRamps::matlab.like(length(unique(defined$ARO_name))),border="white", bty = "n", fill=NULL, bg="white")
# # dev.off()
# # 
# # png(filename=paste0("CLCaligned_and_RGIamr_defined.png"), width = 1840, height = 1840)
# # heatmap.2(log(amr_tab[rownames(defined),]+0.0000001), col=colorRamps::green2red, RowSideColors = colorRamps::matlab.like(length(unique(defined$ARO_name)))[factor(defined$ARO_name)], main="\n\nARO_name\namr counts \12 AMR-associated ORFs\namong all ORFs found with prior CLC alignments\nlog(abs abundance)", margins=c(12,16),cexCol = 0.6, trace='n')
# # legend("left",legend=unique(defined$ARO_name),pch=22,pt.bg=colorRamps::matlab.like(length(unique(defined$ARO_name))),border="white", bty = "n", fill=NULL, bg="white")
# # dev.off()

##############################################
# OVER ALL GENES WITH RGI ASSOCIATIONS # a. ARO_NAME
##############################################

allcounts <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP/cath/annot3.tab", sep="\t", header=TRUE, fill=TRUE)
rownames(allcounts) <- allcounts$orf_id
allcounts <- allcounts[,colnames(allcounts)!="X67.68UC.ORI.NEBNEXT.PE.P51.1_S13_001"]

delivery.df <- as.data.frame(delivery_csv.t)
delivery_csv.df.sort.ARO_name <- delivery.df[order(delivery.df$ARO_name),]3

allcounts.rgi <- allcounts[rownames(delivery_csv.df.sort.ARO_name)[rownames(delivery_csv.df.sort.ARO_name) %in% rownames(allcounts)],34:46]
defined.rgi <- delivery_csv.df.sort.ARO_name[rownames(delivery_csv.df.sort.ARO_name) %in% rownames(allcounts.rgi),]

png(filename=paste0("RGIamr_defined.ARO_name.png"), width = 1840, height = 3680)
heatmap.2(log(as.matrix(allcounts.rgi)+0.0000001), RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi$ARO_name)))[factor(defined.rgi$ARO_name)], main="\n\nARO_name\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\nlog(abs abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
dev.off()

png(filename=paste0("RGIamr_defined.ARO_name.sorted.png"), width = 1840, height = 3680)
heatmap.2(log(as.matrix(allcounts.rgi)+0.0000001), Rowv = colorRamps::matlab.like(length(unique(defined.rgi$ARO_name)))[factor(defined.rgi$ARO_name)], RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi$ARO_name)))[factor(defined.rgi$ARO_name)], main="\n\nARO_name\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\nlog(abs abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
dev.off()

png(filename=paste0("RGIamr_defined.ARO_name.legend.png"), width = 1200, height = 2400)
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("left",legend=unique(defined.rgi$ARO_name),pch=22,pt.bg=colorRamps::matlab.like(length(unique(defined.rgi$ARO_name))),border="white", bty = "n", fill=NULL, bg="white",cex=0.9)
dev.off()

allcounts.rgi <- allcounts[,34:46]
allcounts.rgi.rel <- as.matrix(decostand(allcounts.rgi,"total",2,na.rm=TRUE))
#allcounts.rgi.rel <- allcounts.rgi.rel[) %in% rownames(allcounts)],]
allcounts.rgi.rel <- allcounts.rgi.rel[rownames(delivery_csv.df.sort.ARO_name)[rownames(delivery_csv.df.sort.ARO_name) %in% rownames(allcounts)],]
defined.rgi.rel <- delivery_csv.df.sort.ARO_name[rownames(delivery_csv.df.sort.ARO_name) %in% rownames(allcounts.rgi.rel),]

png(filename=paste0("RGIamr_defined.RelAbs.ARO_name.png"), width = 1840, height = 3680)
heatmap.2(-log(as.matrix(allcounts.rgi.rel)+0.0000001), RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi.rel$ARO_name)))[factor(defined.rgi.rel$ARO_name)], main="\n\nARO_name\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\n-log(rel abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
dev.off()

png(filename=paste0("RGIamr_defined.RelAbs.ARO_name.sorted.png"), width = 1840, height = 3680)
heatmap.2(-log(as.matrix(allcounts.rgi.rel)+0.0000001), Rowv = colorRamps::matlab.like(length(unique(defined.rgi.rel$ARO_name)))[factor(defined.rgi.rel$ARO_name)], RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi.rel$ARO_name)))[factor(defined.rgi.rel$ARO_name)], main="\n\nARO_name\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\n-log(rel abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
dev.off()

















# 
##############################################
# OVER ALL GENES WITH RGI ASSOCIATIONS # b. ARO_name
##############################################

#delivery_csv.df.sort.ARO_name

delivery.df <- as.data.frame(delivery_csv.t)
delivery_csv.df.sort.ARO_name <- delivery.df[order(delivery.df$ARO_name),]

#allcounts.rgi <- allcounts[rownames(delivery_csv.df.sort.ARO_name)[rownames(delivery_csv.df.sort.ARO_name) %in% rownames(allcounts)],31:42]
allcounts.rgi <- allcounts[rownames(delivery_csv.df.sort.ARO_name)[rownames(delivery_csv.df.sort.ARO_name) %in% rownames(allcounts)],34:45]
samp_assign <- c("x29UC",                    "x68UC",               "x33UC",                     "x43UC",                     "x68UC",                     "x59UC",                                    "x68UC",             "x73UC",             "x29UC",              "x29UC",                "x33UC",               "x43UC");
prep_assign <- c("BLOOD AGAR",                "LB AGAR",               "BLOOD AGAR",                "BLOOD AGAR",                "BLOOD AGAR",                     "ORI",                                      "ORI",              "ORI",              "ORI",                "LB AGAR",               "LB AGAR",               "LB AGAR") 


png(filename=paste0("legend.PrepType.png"), width = 1400, height = 1400)
plot.new()
#legend("center",pch=22,cex=4, legend=unique(prep_assign),pt.bg=c(primary.colors(length(unique(prep_assign)))),title="Prep Type")
legend("center",pch=22,cex=4, legend=c("BLOOD AGAR","LB AMP AGAR","ORI" ),pt.bg=c(primary.colors(length(unique(prep_assign)))),title="Prep Type")
dev.off()

png(filename=paste0("legend.SampleID.png"), width = 1400, height = 1400)
plot.new()
legend("center",pch=22,cex=4, legend=unique(samp_assign),pt.bg=c(matlab.like(length(unique(samp_assign))))[as.factor(unique(samp_assign))],title="Sample ID")
dev.off()


png(filename=paste0("legend.PrepType.2supp.png"), width = 1400, height = 1400)
plot.new()
#legend("center",pch=22,cex=4, legend=unique(prep_assign),pt.bg=c(primary.colors(length(unique(prep_assign)))),title="Prep Type")
legend("center",pch=22,cex=4, legend=c("BLOOD AGAR","LB AMP AGAR","ORI" ),pt.bg=c("red","green","lightblue"),title="Prep Type")
dev.off()

#remove NEBnext

defined.rgi <- delivery_csv.df.sort.ARO_name[rownames(delivery_csv.df.sort.ARO_name) %in% rownames(allcounts.rgi),]

png(filename=paste0("RGIamr_defined.ARO_name.png"), width = 1840, height = 3680)
heatmap.2(log(as.matrix(allcounts.rgi)+0.0000001), RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi$ARO_name)))[factor(defined.rgi$ARO_name)], main="\n\nARO_name\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\nlog(abs abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
dev.off()

png(filename=paste0("RGIamr_defined.ARO_name.sorted.png"), width = 1840, height = 3680)
heatmap.2(log(as.matrix(allcounts.rgi)+0.0000001), Rowv = colorRamps::matlab.like(length(unique(defined.rgi$ARO_name)))[factor(defined.rgi$ARO_name)], RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi$ARO_name)))[factor(defined.rgi$ARO_name)], main="\n\nARO_name\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\nlog(abs abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
dev.off()

png(filename=paste0("RGIamr_defined.ARO_name.legend.png"), width = 1200, height = 2400)
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("left",legend=unique(defined.rgi$ARO_name),pch=22,pt.bg=colorRamps::matlab.like(length(unique(defined.rgi$ARO_name))),border="white", bty = "n", fill=NULL, bg="white",cex=0.9)
dev.off()

write.table(allcounts.rgi, file=paste0("RGIamr_defined.RelAbs.ARO_name.tsv"), quote=FALSE, sep='\t', col.names=TRUE)


#allcounts.rgi <- allcounts[,34:46]
allcounts.rgi.rel <- as.matrix(decostand(allcounts.rgi,"total",2,na.rm=TRUE))
defined.rgi.rel <- delivery_csv.df.sort.ARO_name[rownames(delivery_csv.df.sort.ARO_name) %in% rownames(allcounts.rgi.rel),]

png(filename=paste0("carbapenem.LOOSE.negLogRelAbs.SampleID.ARO_name.png"), width = 1400, height = 1400)
heatmap.2(-log(allcounts.rgi.rel[rownames(defined.rgi.rel[grepl("sme|Nmc", defined.rgi.rel$ARO_name),]),]+0.0000000000000001),
          margins=c(15,10),trace='n',cexRow = 0.7, cexCol = 0.7,  
          col=colorRamps::green2red,
          labRow=paste0(rownames(defined.rgi.rel[grepl('sme', defined.rgi.rel$ARO_name),]),'\n',defined.rgi.rel[grepl("sme|Nmc", defined.rgi.rel$ARO_name),]$ARO_name),
          ColSideColors = c(matlab.like(length(unique(samp_assign))))[as.factor(samp_assign)]
)
dev.off()

png(filename=paste0("carbapenem.LOOSE.negLogRelAbs.PrepType.ARO_name.png"), width = 1400, height = 1400)
heatmap.2(-log(allcounts.rgi.rel[rownames(defined.rgi.rel[grepl('sme|Nmc', defined.rgi.rel$ARO_name),]),]+0.0000000000000001),
          margins=c(15,10),trace='n',cexRow = 0.7, cexCol = 0.7,  
          col=colorRamps::green2red,
          labRow=paste0(rownames(defined.rgi.rel[grepl('sme', defined.rgi.rel$ARO_name),]),'\n',defined.rgi.rel[grepl("sme|Nmc", defined.rgi.rel$ARO_name),]$ARO_name),
          ColSideColors = c(primary.colors(length(unique(prep_assign))))[as.factor(prep_assign)]
)
dev.off()


png(filename=paste0("legend.PrepType.png"), width = 1400, height = 1400)
plot.new()
legend("center",pch=22,cex=4, legend=unique(prep_assign),pt.bg=c(primary.colors(length(unique(prep_assign)))),title="Prep Type")
dev.off()

png(filename=paste0("legend.SampleID.png"), width = 1400, height = 1400)
plot.new()
legend("center",pch=22,cex=4, legend=unique(samp_assign),pt.bg=c(matlab.like(length(unique(samp_assign)))),title="Sample ID")
dev.off()

all_tmp <- list()
all_id <- list()
for(q in 1:6){
  all_id[q] <- rownames(allcounts.rgi.rel[rownames(defined.rgi.rel[grepl("sme|Nmc", as.character(plyr::ldply(strsplit(as.character(defined.rgi.rel$ARO_name), as.character(';')), rbind)[,1])),]),])[q]
  all_tmp[q] <- strsplit(as.character(defined.rgi.rel[rownames(allcounts.rgi.rel[rownames(defined.rgi.rel[grepl("sme|Nmc", as.character(plyr::ldply(strsplit(as.character(defined.rgi.rel$ARO_name), as.character(';')), rbind)[,1])),]),]),]$ARO_name),';')[[q]][1]
}

loose.df <- as.data.frame(matrix(ncol=length(colnames(allcounts.rgi.rel)),nrow=length(unique(all_tmp))))
rownames(loose.df) <- unique(all_tmp)
colnames(loose.df) <- colnames(allcounts.rgi.rel)

for(t in 1:6){
  if(is.null(loose.df[all_id[[t]][1],all_tmp[[t]][1]]) == TRUE){
    loose.df[ all_tmp[[t]][1], ] <- allcounts.rgi.rel[all_id[[t]][1],]
  }else{
    loose.df[ all_tmp[[t]][1], ] <- loose.df[all_tmp[[t]][1],]  + allcounts.rgi.rel[all_id[[t]][1],unlist(all_id)]
  }
}

png(filename=paste0("carbapenem.STRICT.negLogRelAbs.SampleID.ARO_name.png"), width = 1400, height = 1400)
heatmap.2(-log(as.matrix(loose.df)+0.0000000000000001),
          margins=c(25,15),trace='n',cexRow = 3.5, cexCol = 1.5,  
          col=colorRamps::green2red,
          #margins=c(10,10),
          #labRow=
            #paste0(rownames(defined.rgi.rel[grepl('sme', as.character(plyr::ldply(strsplit(as.character(defined.rgi.rel$ARO_name), as.character(';')), rbind)[,1])),]),'\n',defined.rgi.rel[grepl('sme', as.character(plyr::ldply(strsplit(as.character(defined.rgi.rel$ARO_name), as.character(';')), rbind)[,1])),]$ARO_name),
          ColSideColors = c(matlab.like(length(unique(samp_assign))))[as.factor(samp_assign)]
)
dev.off()
 
png(filename=paste0("carbapenem.STRICT.negLogRelAbs.PrepType.ARO_name.png"), width = 1400, height = 1400)
heatmap.2(-log(as.matrix(loose.df)+0.0000000000000001),
          margins=c(15,25),trace='n',cexRow = 3.5, cexCol = 1.5,  
          col=colorRamps::green2red,
          #labRow=
            #paste0(rownames(defined.rgi.rel[grepl('sme', as.character(plyr::ldply(strsplit(as.character(defined.rgi.rel$ARO_name), as.character(';')), rbind)[,1])),]),'\n',defined.rgi.rel[grepl('sme', as.character(plyr::ldply(strsplit(as.character(defined.rgi.rel$ARO_name), as.character(';')), rbind)[,1])),]$ARO_name),
            ColSideColors = c(primary.colors(length(unique(prep_assign))))[as.factor(prep_assign)]
)
dev.off()

tmpout <- allcounts.rgi.rel[rownames(defined.rgi.rel[grepl("sme|Nmc", ignore.case = TRUE, defined.rgi.rel$ARO_name),]),]
write.table(tmpout, file=paste0("carbapenem.RelAbs.tsv"), quote=FALSE, sep='\t', col.names=TRUE, 
            row.names=paste0(rownames(defined.rgi.rel[grepl('sme', defined.rgi.rel$ARO_name),]),' : ',defined.rgi.rel[grepl('sme', defined.rgi.rel$ARO_name),]$ARO_name))









png(filename=paste0("Blacatamase.LOOSE.negLogRelAbs.SampleID.ARO_name.png"), width = 1400, height = 1400)
heatmap.2(-log(allcounts.rgi.rel[rownames(defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]),]+0.0000000000000001), 
          margins=c(15,20),cexRow =3.5, cexCol = 1.5, col=colorRamps::green2red, trace='n',
          labRow=paste0(rownames(defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]),'\n',defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]$ARO_name),
          ColSideColors = c(matlab.like(length(unique(samp_assign))))[as.factor(samp_assign)]
)
dev.off()

png(filename=paste0("Blacatamase.LOOSE.negLogRelAbs.PrepType.ARO_name.png"), width = 1400, height = 1400)
heatmap.2(-log(allcounts.rgi.rel[rownames(defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]),]+0.0000000000000001), 
          margins=c(15,20),cexRow = 3.5, cexCol = 1.5, col=colorRamps::green2red, trace='n',
          labRow=paste0(rownames(defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]),'\n',defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]$ARO_name),
          ColSideColors = c(primary.colors(length(unique(prep_assign))))[as.factor(prep_assign)]
)
dev.off()

all_tmp <- list()
all_id <- list()
for(q in 1:4){
  all_id[q] <- rownames(allcounts.rgi.rel[rownames(defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", as.character(plyr::ldply(strsplit(as.character(defined.rgi.rel$ARO_name), as.character(';')), rbind)[,1])),]),])[q]
  all_tmp[q] <- strsplit(as.character(defined.rgi.rel[rownames(allcounts.rgi.rel[rownames(defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", as.character(plyr::ldply(strsplit(as.character(defined.rgi.rel$ARO_name), as.character(';')), rbind)[,1])),]),]),]$ARO_name),';')[[q]][1]
}

loose.df <- as.data.frame(matrix(ncol=length(colnames(allcounts.rgi.rel)),nrow=length(unique(all_tmp))))
rownames(loose.df) <- unique(all_tmp)
colnames(loose.df) <- colnames(allcounts.rgi.rel)

for(t in 1:4){
  if(is.null(loose.df[all_id[[t]][1],all_tmp[[t]][1]]) == TRUE){
    loose.df[ all_tmp[[t]][1], ] <- allcounts.rgi.rel[all_id[[t]][1],]
  }else{
    loose.df[ all_tmp[[t]][1], ] <- loose.df[all_tmp[[t]][1],]  + allcounts.rgi.rel[all_id[[t]][1],unlist(all_id)]
  }
}

png(filename=paste0("Blacatamase.STRICT.negLogRelAbs.SampleID.ARO_name.png"), width = 1400, height = 1400)
heatmap.2(-log(as.matrix(loose.df)+0.0000000000000001), 
          margins=c(15,40),cexRow =3.5, cexCol = 1.5, col=colorRamps::green2red, trace='n',
          #labRow=paste0(rownames(defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]),'\n',defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]$ARO_name),
          ColSideColors = c(matlab.like(length(unique(samp_assign))))[as.factor(samp_assign)]
)
dev.off()

png(filename=paste0("Blacatamase.STRICT.negLogRelAbs.PrepType.ARO_name.png"), width = 1400, height = 1400)
heatmap.2(-log(as.matrix(loose.df)+0.0000000000000001), 
          margins=c(15,40),cexRow = 3.5, cexCol = 1.5, col=colorRamps::green2red, trace='n',
          #labRow=paste0(rownames(defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]),'\n',defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]$ARO_name),
          ColSideColors = c(primary.colors(length(unique(prep_assign))))[as.factor(prep_assign)]
)
dev.off()





















tmpout <- allcounts.rgi.rel[rownames(defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]),]
write.table(tmpout, file=paste0("Blacatamase.LOOSE.RelAbs.tsv"), quote=FALSE, sep='\t', col.names=TRUE, 
            row.names=paste0(rownames(defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]),' : ',defined.rgi.rel[grepl("oxy-2-7|mecA|mecR1|mecI|cmy-14|sed1|pbp1a|pbp1x|oxa|shv-28|bla", ignore.case = TRUE, defined.rgi.rel$ARO_name),]$ARO_name))

png(filename=paste0("legend.PrepType.png"), width = 1400, height = 1400)
plot.new()
legend("center",pch=22,cex=4, legend=unique(prep_assign),pt.bg=c(primary.colors(length(unique(prep_assign)))),title="Prep Type")
dev.off()

png(filename=paste0("legend.SampleID.png"), width = 1400, height = 1400)
plot.new()
legend("center",pch=22,cex=4, legend=unique(samp_assign),pt.bg=c(matlab.like(length(unique(samp_assign)))),title="Sample ID")
dev.off()

# 
# png(filename=paste0("RGIamr_defined.RelAbs.Best_Hit_ARO.png"), width = 1840, height = 3680)
# heatmap.2(-log(allcounts.rgi.rel+0.0000001), RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi.rel$Best_Hit_ARO)))[factor(defined.rgi.rel$Best_Hit_ARO)], main="\n\nBest_Hit_ARO\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\n-log(rel abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
# dev.off()
# 
# png(filename=paste0("RGIamr_defined.RelAbs.Best_Hit_ARO.sorted.png"), width = 1840, height = 3680)
# heatmap.2(-log(allcounts.rgi.rel+0.0000001), Rowv = colorRamps::matlab.like(length(unique(defined.rgi.rel$Best_Hit_ARO)))[factor(defined.rgi.rel$Best_Hit_ARO)], RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi.rel$Best_Hit_ARO)))[factor(defined.rgi.rel$Best_Hit_ARO)], main="\n\nBest_Hit_ARO\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\n-log(rel abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
# dev.off()

# ##############################################
# # OVER ALL GENES WITH RGI ASSOCIATIONS # c. Best_Hit_ARO_category
# ##############################################
# 
# #delivery_csv.df.sort.ARO_name
# 
# delivery.df <- as.data.frame(delivery_csv.t)
# delivery_csv.df.sort.Best_Hit_ARO_category <- delivery.df[order(delivery.df$Best_Hit_ARO_category),]
# 
# allcounts.rgi <- allcounts[rownames(delivery_csv.df.sort.Best_Hit_ARO_category)[rownames(delivery_csv.df.sort.Best_Hit_ARO_category) %in% rownames(allcounts)],34:46]
# defined.rgi <- delivery_csv.df.sort.Best_Hit_ARO_category[rownames(delivery_csv.df.sort.Best_Hit_ARO_category) %in% rownames(allcounts.rgi),]
# 
# png(filename=paste0("RGIamr_defined.Best_Hit_ARO_category.png"), width = 1840, height = 3680)
# heatmap.2(log(as.matrix(allcounts.rgi)+0.0000001), RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi$Best_Hit_ARO_category)))[factor(defined.rgi$Best_Hit_ARO_category)], main="\n\nBest_Hit_ARO_category\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\nlog(abs abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
# dev.off()
# 
# png(filename=paste0("RGIamr_defined.Best_Hit_ARO_category.sorted.png"), width = 1840, height = 3680)
# heatmap.2(log(as.matrix(allcounts.rgi)+0.0000001), Rowv = colorRamps::matlab.like(length(unique(defined.rgi$Best_Hit_ARO_category)))[factor(defined.rgi$Best_Hit_ARO_category)], RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi$Best_Hit_ARO_category)))[factor(defined.rgi$Best_Hit_ARO_category)], main="\n\nBest_Hit_ARO_category\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\nlog(abs abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
# dev.off()
# 
# png(filename=paste0("RGIamr_defined.Best_Hit_ARO_category.legend.png"), width = 1200, height = 2400)
# plot(1, type="n", axes=FALSE, xlab="", ylab="")
# legend("left",legend=unique(defined.rgi$Best_Hit_ARO_category),pch=22,pt.bg=colorRamps::matlab.like(length(unique(defined.rgi$Best_Hit_ARO_category))),border="white", bty = "n", fill=NULL, bg="white",cex=0.9)
# dev.off()
# 
# allcounts.rgi <- allcounts[,34:46]
# allcounts.rgi.rel <- as.matrix(decostand(allcounts.rgi,"total",2,na.rm=TRUE))
# defined.rgi.rel <- delivery_csv.df.sort.Best_Hit_ARO_category[rownames(delivery_csv.df.sort.Best_Hit_ARO_category) %in% rownames(allcounts.rgi.rel),]
# 
# png(filename=paste0("RGIamr_defined.RelAbs.Best_Hit_ARO_category.png"), width = 1840, height = 3680)
# heatmap.2(-log(allcounts.rgi.rel+0.0000001), RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi.rel$Best_Hit_ARO_category)))[factor(defined.rgi.rel$Best_Hit_ARO_category)], main="\n\nBest_Hit_ARO_category\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\n-log(rel abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
# dev.off()
# 
# png(filename=paste0("RGIamr_defined.RelAbs.Best_Hit_ARO_category.sorted.png"), width = 1840, height = 3680)
# heatmap.2(-log(allcounts.rgi.rel+0.0000001), Rowv = colorRamps::matlab.like(length(unique(defined.rgi.rel$Best_Hit_ARO_category)))[factor(defined.rgi.rel$Best_Hit_ARO_category)], RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi.rel$Best_Hit_ARO_category)))[factor(defined.rgi.rel$Best_Hit_ARO_category)], main="\n\nBest_Hit_ARO_category\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\n-log(rel abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
# dev.off()

##############################################
# OVER ALL GENES WITH RGI ASSOCIATIONS # d. ARO_category
##############################################

#delivery_csv.df.sort.ARO_name

delivery.df <- as.data.frame(delivery_csv.t)
delivery_csv.df.sort.ARO_category <- delivery.df[order(delivery.df$ARO_category),]

allcounts.rgi <- allcounts[rownames(delivery_csv.df.sort.ARO_category)[rownames(delivery_csv.df.sort.ARO_category) %in% rownames(allcounts)],34:46]
defined.rgi <- delivery_csv.df.sort.ARO_category[rownames(delivery_csv.df.sort.ARO_category) %in% rownames(allcounts.rgi),]

png(filename=paste0("RGIamr_defined.ARO_category.png"), width = 1840, height = 3680)
heatmap.2(log(as.matrix(allcounts.rgi)+0.0000001), RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi$ARO_category)))[factor(defined.rgi$ARO_category)], main="\n\nARO_category\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\nlog(abs abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
dev.off()

png(filename=paste0("RGIamr_defined.ARO_category.sorted.png"), width = 1840, height = 3680)
heatmap.2(log(as.matrix(allcounts.rgi)+0.0000001), Rowv = colorRamps::matlab.like(length(unique(defined.rgi$ARO_category)))[factor(defined.rgi$ARO_category)], RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi$ARO_category)))[factor(defined.rgi$ARO_category)], main="\n\nARO_category\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\nlog(abs abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
dev.off()

png(filename=paste0("RGIamr_defined.ARO_category.legend.png"), width = 1200, height = 2400)
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("left",legend=unique(defined.rgi$ARO_category),pch=22,pt.bg=colorRamps::matlab.like(length(unique(defined.rgi$ARO_category))),border="white", bty = "n", fill=NULL, bg="white",cex=0.9)
dev.off()

allcounts.rgi <- allcounts[,34:46]
allcounts.rgi.rel <- as.matrix(decostand(allcounts.rgi,"total",2,na.rm=TRUE))[rownames(delivery_csv.df.sort.ARO_category)[rownames(delivery_csv.df.sort.ARO_category) %in% rownames(allcounts)],]
defined.rgi.rel <- delivery_csv.df.sort.ARO_category[rownames(delivery_csv.df.sort.ARO_category) %in% rownames(allcounts.rgi.rel),]

png(filename=paste0("RGIamr_defined.RelAbs.ARO_category.png"), width = 1840, height = 3680)
heatmap.2(-log(allcounts.rgi.rel+0.0000001), RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi.rel$ARO_category)))[factor(defined.rgi.rel$ARO_category)], main="\n\nARO_category\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\n-log(rel abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
dev.off()

png(filename=paste0("RGIamr_defined.RelAbs.ARO_category.sorted.png"), width = 1840, height = 3680)
heatmap.2(-log(allcounts.rgi.rel+0.0000001), Rowv = colorRamps::matlab.like(length(unique(defined.rgi.rel$ARO_category)))[factor(defined.rgi.rel$ARO_category)], RowSideColors = colorRamps::matlab.like(length(unique(defined.rgi.rel$ARO_category)))[factor(defined.rgi.rel$ARO_category)], main="\n\nARO_category\namr counts \364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\n-log(rel abundance)", col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
dev.off()

c <- 1
list_of_lists <- list()
for (category in unique(defined.rgi$ARO_category)){
  
  print(category)
  
  cat_ids <- rownames(delivery_csv.df.sort.ARO_name[delivery_csv.df.sort.ARO_name$ARO_category == category,])
  defined.rgi.rel <- delivery_csv.df.sort.ARO_name[rownames(delivery_csv.df.sort.ARO_name) %in% cat_ids,]
  
  #allcounts.rgi <- allcounts[,34:46]
  #allcounts.rgi.rel <- as.matrix(decostand(allcounts.rgi,"total",2,na.rm=TRUE))
  # rownames(allcounts) %in% cat_ids
  allcounts.rgi.rel.cat <- allcounts.rgi.rel[rownames(allcounts.rgi.rel) %in% cat_ids,]
  
  print(dim(defined.rgi.rel))
  print(dim(allcounts.rgi.rel.cat))
  
  #if(dim(defined.rgi.rel)[1] > 1){
  png(filename=paste0("RGIamr_defined.RelAbs.ARO_category.", c, ".", category,".png"), width = 1840, height = 1840)
  #if(length(cat_ids) > 1){
  if(sum(rownames(allcounts.rgi.rel) %in% cat_ids) > 1){
    heatmap.2(-log(as.matrix(allcounts.rgi.rel.cat)+0.0000001), main=paste0("\n\n\n\n", c," : ARO_category = \n", category,"\namr counts \n364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\n-log(rel abundance)"), col=colorRamps::green2red, margins=c(12,16),cexCol = 0.6, trace='n')
    summed.cat <- rowSums(allcounts.rgi.rel.cat, na.rm = TRUE)
  }else{
    par(mar = c(15, 2, 4, 2))
    barplot(-log(allcounts.rgi.rel.cat+0.0000001),las=2, cex.names = 0.8, main=paste0("\n\n\n\n", c," : ARO_category = \n", category,"\namr counts \n364 (with ID match) of 433 (found by RGI)\nAMR-associated ORFs\n-log(rel abundance)"))
    summed.cat <- as.list(allcounts.rgi.rel.cat)
  }
  dev.off()
  
  list_of_lists[c] <- list(summed.cat)
  
  c<-c+1
}

names(list_of_lists) <- unique(defined.rgi$ARO_category)


png(filename=paste0("RelAbs_of_each_ARO_category.png"), width = 900, height = 1680)
par(mar = c(55, 4, 2, 1))
boxplot(
  unlist(list_of_lists[1]),unlist(list_of_lists[2]),unlist(list_of_lists[3]),unlist(list_of_lists[4]),unlist(list_of_lists[5]),unlist(list_of_lists[6]),unlist(list_of_lists[7]),unlist(list_of_lists[8]),unlist(list_of_lists[9]),unlist(list_of_lists[10]),unlist(list_of_lists[11]),unlist(list_of_lists[12]),unlist(list_of_lists[13]),unlist(list_of_lists[14]),unlist(list_of_lists[15]),unlist(list_of_lists[16]),unlist(list_of_lists[17]),unlist(list_of_lists[18]),unlist(list_of_lists[19]),unlist(list_of_lists[20]),unlist(list_of_lists[21]),unlist(list_of_lists[22]),unlist(list_of_lists[23]),unlist(list_of_lists[24]),unlist(list_of_lists[25]),unlist(list_of_lists[26]),unlist(list_of_lists[27]),unlist(list_of_lists[28]),unlist(list_of_lists[29]),unlist(list_of_lists[30]),unlist(list_of_lists[31]),unlist(list_of_lists[32]),
  xlab="", ylab="Relative Abundance", margins = c(0, 0),las=3,
  #col=colorRamps::matlab.like(length(list_of_lists))[factor(unique(defined.rgi$ARO_category))],
  main="Abundance of all Card-Associated ORFs within 38 ARO Category Types", 
  #xlab="ARO_category Abundance\n", 
  ylab="Relative Abundance",
  las=3, 
  names=unique(defined.rgi$ARO_category), 
  cex=0.5,
  col=colorRamps::matlab.like(length(list_of_lists))[factor(unique(defined.rgi$ARO_category))],
  cex.axis=0.7
)
dev.off()

png(filename=paste0("negLogRelAbs_of_each_ARO_category.png"), width = 900, height = 1680)
par(mar = c(55, 4, 2, 1))
boxplot(
  -log(unlist(list_of_lists[1])+0.0000001),
  -log(unlist(list_of_lists[2])+0.0000001),
  -log(unlist(list_of_lists[3])+0.0000001),
  -log(unlist(list_of_lists[4])+0.0000001),
  -log(unlist(list_of_lists[5])+0.0000001),
  -log(unlist(list_of_lists[6])+0.0000001),
  -log(unlist(list_of_lists[7])+0.0000001),
  -log(unlist(list_of_lists[8])+0.0000001),
  -log(unlist(list_of_lists[9])+0.0000001),
  -log(unlist(list_of_lists[10])+0.0000001),
  -log(unlist(list_of_lists[11])+0.0000001),
  -log(unlist(list_of_lists[12])+0.0000001),
  -log(unlist(list_of_lists[13]) +0.0000001),
  -log(unlist(list_of_lists[14])+0.0000001),
  -log(unlist(list_of_lists[15])+0.0000001),
  -log(unlist(list_of_lists[16])+0.0000001),
  -log(unlist(list_of_lists[17])+0.0000001),
  -log(unlist(list_of_lists[18])+0.0000001),
  -log(unlist(list_of_lists[19])+0.0000001),
  -log(unlist(list_of_lists[20])+0.0000001),
  -log(unlist(list_of_lists[21])+0.0000001),
  -log(unlist(list_of_lists[22])+0.0000001),
  -log(unlist(list_of_lists[23])+0.0000001),
  -log(unlist(list_of_lists[24])+0.0000001),
  -log(unlist(list_of_lists[25])+0.0000001),
  -log(unlist(list_of_lists[26])+0.0000001),
  -log(unlist(list_of_lists[27])+0.0000001),
  -log(unlist(list_of_lists[28])+0.0000001),
  -log(unlist(list_of_lists[29])+0.0000001),
  -log(unlist(list_of_lists[30])+0.0000001),
  -log(unlist(list_of_lists[31])+0.0000001),
  -log(unlist(list_of_lists[32])+0.0000001),
  main="-log(Relative Abundance) of all Card-Associated ORFs within 38 ARO Category Types", 
  #xlab="ARO_category Abundance\n", 
  ylab="-log(Relative Abundance)",
  las=3, 
  names=unique(defined.rgi$ARO_category), 
  cex=0.5,
  col=colorRamps::matlab.like(length(list_of_lists))[factor(unique(defined.rgi$ARO_category))],
  cex.axis=0.7
)
legend("topleft","less abundant",cex=0.8,bty='n')
legend("bottomleft","more abundant",cex=0.8,bty='n')
dev.off()



##############################################
# MAX ORF REP PER SAMPLE ID
##############################################

# PREP INPUT
rgi.rel <- count_REL_table
# allcounts.rgi.rel

# ITERATE THROUGH ALL SAMP IDS
for(samp in unique(samp_assign)){
  orf.rel.samp <- allcounts.rgi.rel[,samp_assign == samp]
  rgi.rel.samp <- rgi.rel[,samp_assign == samp]
  if(sum(samp_assign == samp) > 1){
    
    orf.sorted_names <- names(rev(sort(rowSums(orf.rel.samp))))
    rgi.sorted_names <- names(rev(sort(rowSums(rgi.rel.samp))))
    
    metadata.top25 <- delivery_csv.df.sort.ARO_name[names(rowSums(orf.rel.samp[orf.sorted_names,])[1:25]),]
    bxplt <- t(orf.rel.samp[orf.sorted_names,])[,1:25]
    
    png(filename=paste0(samp,".relAbs.png"), width = 900, height = 1680)
    par(mfrow = c(2,1))
    par(mar = c(16, 4, 2, 1))
    colnames(bxplt) <- metadata.top25[names(rowSums(orf.rel.samp[orf.sorted_names,])[1:25]),"Best_Hit_ARO"]
    boxplot(bxplt,las=3,cex.axis=0.8)
    par(mar = c(13, 4, 24, 1))
    colnames(bxplt) <- names(rowSums(orf.rel.samp[orf.sorted_names,])[1:25])
    barplot(bxplt,las=3,cex.names=0.8,main=paste0(samp,"\nRelative Abundance : Top 25 ORF IDs"))
    dev.off()
    
  }else{
    
    png(filename=paste0(samp,".relAbs.png"), width = 900, height = 1000)
    par(mar = c(25, 4, 2, 1))
    par(mfrow = c(1,1))
    print(paste0("1sample: ",samp))
    metadata <- delivery_csv.df.sort.ARO_name[names(rev(sort(orf.rel.samp)))[1:25],]
    barplot(rev(sort(orf.rel.samp))[1:25],names=delivery_csv.df.sort.ARO_name[names(rev(sort(orf.rel.samp)))[1:25],"Best_Hit_ARO"],las=3,cex.names=0.8,main=paste0(samp,"\nRelative Abundance : Top 25 ORF IDs"))
    dev.off()
    
  }
}

uniq_aro_sum.df <- as.data.frame(matrix(nrow=length(unique(delivery_csv.df.sort.ARO_name$ARO_name)),ncol=length(unique(samp_assign))))
rownames(uniq_aro_sum.df) <- unique(delivery_csv.df.sort.ARO_name$ARO_name)
colnames(uniq_aro_sum.df) <- unique(samp_assign)
for(samp in unique(samp_assign)){
  orf.rel.samp <- allcounts.rgi.rel[,samp_assign == samp]
  for(uniq_aro_name in unique(delivery_csv.df.sort.ARO_name$ARO_name)){
    print(paste0(samp,' : ',uniq_aro_name))
    print(paste0(samp,' : ',uniq_aro_name))
    uniq_aro_sum.df[uniq_aro_name,samp]<- 
      sum(allcounts.rgi.rel[rownames(allcounts.rgi.rel) %in% rownames(delivery_csv.df.sort.ARO_name[delivery_csv.df.sort.ARO_name$ARO_name == uniq_aro_name,])
                            ,samp_assign == samp])
  }
}
png(filename="ARO_name_sum.png", width = 1000, height = 1000)
#par(mar = c(25, 4, 2, 15))
heatmap.2(-log(as.matrix(uniq_aro_sum.df)+0.0000001),margins=c(8,50),cex.axis=0.6,trace="none",main="\n\nARO_name\nSum Rel\nAbundance\nby Sample")
dev.off()

2
#for(samp in unique(samp_assign)){
#  uniq_aro_sum.df <- uniq_aro_sum.df[,samp_assign == samp]
#
#}