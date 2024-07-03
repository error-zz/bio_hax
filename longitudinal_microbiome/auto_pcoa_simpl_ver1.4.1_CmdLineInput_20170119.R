#!/usr/bin/env Rscript
library(vegan) 
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
library(heatmap.plus)


setwd ("/Users/apple/Desktop/exe_q_1.4_infantGlobal")
species <- read.table( "pcoa_otu_tax_vs_frequency.reduced.txt", sep="\t", header=TRUE )
#groups <- read.table( "global_analysis_prep.infants_A.txt", sep="\t", header=TRUE )
groups <- read.table( "infant_manifest_formatted.txt", sep="\t", header=TRUE )
tax_summary <- read.table( "otus.qfilter.EE0.15.curated_SILVA_123_SSURef_Nr99_tax_silva.wang.tax.summary", sep="\t", header=TRUE )
prefix <- "JMdebug_20180119"
cutoff <- 10^-4

# setwd ("/Users/apple/Desktop/exe_q_1.4_RITM054")
# groups <- read.table( "RITM054cmdLineRapManifest.SEP.tsv", sep="\t", header=TRUE )

# setwd ("/Users/apple/Desktop/exe_q_1.3prod")
# species <- read.table( "4.0_PCOA_OUTPUT/pcoa_otu_tax_vs_frequency.reduced.txt", sep="\t", header=TRUE )
# groups <- read.table( "manifest.tsv", sep="\t", header=TRUE )
# tax_summary <- read.table( "1.0_UPARSE_OUTPUT/otus.qfilter.EE0.15.curated_SILVA_123_SSURef_Nr99_tax_silva.wang.tax.summary", sep="\t", header=TRUE )
# prefix <- "JMDEBUG"
# cutoff <- 10^-4
 
#pcoa_manifest.tsv $sample_tsv PCOA

#input
#options(echo=TRUE) # if you want see commands in output file
#args <- commandArgs(trailingOnly = TRUE)
#print(args)
#species <- read.table( args[1], sep="\t", header=TRUE )
#groups <- read.table( args[2], sep="\t", header=TRUE )
#tax_summary <- read.table( args[3], sep="\t", header=TRUE )
#prefix <- args[4]
#print(species)

# parse out rownames and mice_meta
species.orig <- species
rownames(species) <- species$rankID
species <- species[,7:dim(species)[2]-1]
species <- t(species)
species <- as.data.frame(species)
species <- data.matrix(species)
species.pop <- species[!apply(is.na(species) | species == "", 1, all),]
species_relab<-decostand(species.pop, method="total")*100 # transfrom the species data.frame to OTU/species relative abundances

all_col <- c("yellowgreen","red","turquoise4","goldenrod1","slateblue4","lawngreen","royalblue1","orangered4","purple1","navyblue","mediumblue","lemonchiffon3","indianred1","lightblue","sienna1","steelblue4","lightcyan",'lightgreen','pink','yellow',"khaki","violetred2","honeydew","orangered2","navy","snow3")

#######################
# INCLUDING UNDEFINED #
#######################

dist.bray<-vegdist(species_relab, method="bray") # creates a Bray-Curtis dissimilarity matrix from the species data.frame. adding “pa” , would make the matrix unweighted (presence absence). Other
pcoa.dist<-pcoa(dist.bray) 

# PCOA viz
for(col_name in names(groups)[4:dim(groups)[2]]){
if(!is.null(col_name)){ 
if(sum(is.na(unlist(groups[col_name]))) != length(unlist(groups[col_name]))){
  
  lll <- length(unique(unlist(groups.single))[!is.na(unique(unlist(groups.single)))])
  cexsize <- 2
  if(lll > 20){
    cexsize <- 1.4
  }
  if(lll > 50){
    cexsize <- 1
  }
  if(lll > 100){
    cexsize <- 0.6
  }
  if(lll > 150){
    cexsize <- 0.4
  }
  if(lll > 150){
    cexsize <- 0.3
  }
  
  print(col_name)
  if( sum( grepl("\\|", unlist(groups[col_name])) ) == 0){
    # NO SPLIT ENTRIES
    groups.single            <- unlist(groups[col_name])
    groups.multi.super_group <- "NOPE"
    groups.multi.sub_group   <- "NOPE"
    pch_type <- 22
    if ( sum(grepl('[A-Za-z]', as.character(unlist(groups.single)) )) > 0 ){
      print("SINGLE F")
      col_type <- all_col[factor(as.character(unlist(groups.single)))]  
      png(filename=paste0(col_name,".pcoa.col_legend.png"), width = 800, height = 800 )
      plot.new()
      legend("center",
             legend=unique(unlist(unname(as.list(groups.single)))),
             pch=22,
             cex=cexsize,
             border="white", bty = "n", fill=NULL, bg="white",
             pt.bg=all_col[unique(unlist(groups.single))])
      dev.off()
    }else{
      print("SINGLE N")
      groups.single[as.character(unlist(groups[col_name])) == ""] <- 0
      pallette <- matlab.like2(max(as.numeric(groups.single[groups.single != "_"]), na.rm = TRUE))
      #pallette[-1] <- "#FFFFFF"
      col_type <- pallette[as.numeric(as.matrix(as.data.frame(groups.single)))]
      col_type[as.character(unlist(groups[col_name])) == ""] <- "#FFFFFF"
      png(filename=paste0(col_name,".pcoa.col_legend.png"), width = 800, height = 800 )
      plot.new()
      legend("center",
             legend=sort(unique(unlist(unname(as.list(groups.single))))),
             pch=22,
             cex=cexsize,
             border="white", bty = "n", fill=NULL, bg="white",
             pt.bg=pallette[sort(unique(unlist(unname(as.list(groups.single)))))])
      dev.off()
    }
    
  }else{
    # SPLIT ENTRIES INCLUDED
    groups.single            <- "NOPE"
    assignment_vector <- as.character(unlist(groups[col_name]))
    assignment_vector[assignment_vector == ""] <- "_|_"
    groups.multi.super_group <- as.data.frame(strsplit(assignment_vector,'\\|'))[1,]
    pch_type <- c(21:25)[factor(as.character(unlist(groups.multi.super_group)))]
    png(filename=paste0(col_name,".pcoa.pch_legend.png"), width = 800, height = 800 )
    plot.new()
    legend("center",
           legend=sort(unique(unlist(unname(as.list(groups.multi.super_group))))),
           pch=unique(c(21:25)[factor(as.character(unlist(groups.multi.super_group)))]),
           #pch=22,
           cex=cexsize,
           border="white", bty = "n", fill=NULL, bg="white",
           pt.bg="black")
    dev.off()
    
    groups.multi.sub_group   <- as.data.frame(strsplit(assignment_vector,'\\|'))[2,]
    png(filename=paste0(col_name,"pcoa.col_legend.png"), width = 800, height = 800 )
    if ( sum(grepl('[A-Za-z]', as.character(unlist(groups.multi.sub_group)) )) > 0 ){
      # as factor
      print("SUB F")
      col_type <- all_col[factor(as.character(unlist(groups.multi.sub_group)))]
      col_type[groups.multi.sub_group == "_"] <- "#FFFFFF"
      png(filename=paste0(col_name,".pcoa_LEGENDsub.all_samples.png"), width = 400, height = 700 )
      plot.new()
      legend("center",
             legend=levels(factor(as.character(unlist(groups.multi.sub_group)))),
             pch=22, 
             cex=cexsize,
             border="white", bty = "n", fill=NULL, bg="white",
             pt.bg=all_col [ unique(factor(as.character(unlist(groups.multi.sub_group)))) ] 
             )
      dev.off()
    }else{
      #as numeric
      print("SUB N")
      groups.multi.sub_group[as.character(unlist(groups[col_name])) == ""] <- 1
      pallette <- matlab.like2(max(as.numeric(groups.multi.sub_group[groups.multi.sub_group != "_"])))
      #pallette[-1] <- "#FFFFFF"
      col_type <- pallette[as.numeric(as.matrix(as.data.frame(groups.multi.sub_group)))]
      col_type[as.character(unlist(groups[col_name])) == ""] <- "#FFFFFF"
      groups.multi.sub_group[as.character(unlist(groups[col_name])) == ""] <- ""
      png(filename=paste0(col_name,".pcoa_LEGENDsub.all_samples.png"), width = 400, height = 700 )
      plot.new()
      legend("center",
             legend=sort(unique(as.numeric(as.matrix(as.data.frame(groups.multi.sub_group)))))[sort(unique(as.numeric(as.matrix(as.data.frame(groups.multi.sub_group))))) > 0],
             pch=22, cex=cexsize,
             border="white", bty = "n", fill=NULL, bg="white",
             pt.bg=pallette[ sort(unique(as.numeric(as.matrix(as.data.frame(groups.multi.sub_group)))))[sort(unique(as.numeric(as.matrix(as.data.frame(groups.multi.sub_group))))) > 0] ])
      dev.off()
    }
    
  }
   
  png(filename=paste0(col_name,".pcoa.all_samples.png"), width = 800, height = 800 )
  plot(pcoa.dist$vectors[,1:2],
       #bg=col_type,
       bg=alpha(col_type, 0.6),
       col="white",
       pch=pch_type,
       cex=2,
       alpha = 0.6,
       xlab=paste("PCoA.1(",round(100*pcoa.dist$values[1,3], digits=2),"%)"), 
       ylab=paste("PCoA.2(",round(100*pcoa.dist$values[2,3], digits=2),"%)"),
       main=paste0("All Samplse, Including Undefined, Labelled by Metadata Column\n",col_name)
  )
  
  if(dim(pcoa.dist$vectors)[1] <= 100){
    text( pcoa.dist$vectors, rownames( species ), pos=3, col="darkgray", cex=0.6 )
    text( pcoa.dist$vectors, as.character(unlist(groups[col_name])), pos=1, col="darkgray", cex=0.7 )
  }
  dev.off()
  
}
}
}
  
#######################
# EXCLUDING UNDEFINED #
#######################

for(col_name in names(groups)[6:dim(groups)[2]]){
if(sum(is.na(unlist(groups[col_name]))) != length(unlist(groups[col_name]))){
#for(col_name in names(groups)[4:dim(groups)[2]]){
  #species_relab.col <- species_relab[groups[groups[col_name]!="",]$SampleID,]
  species_relab<-decostand(species.pop, method="total")*100 # transfrom the species data.frame to OTU/species relative abundances
  #groups <- read.table( "infant_manifest_formatted.txt", sep="\t", header=TRUE )
  groups <- read.table( "RITM054cmdLineRapManifest.SEP.tsv", sep="\t", header=TRUE )
  
  species_relab.col <- species_relab[rownames(species_relab) %in% groups[groups[col_name]!="",]$SampleID,]
  print(paste(col_name,": calculating bray-curtis distance matrix"))
  dist.bray.col<-vegdist(species_relab.col, method="bray") # creates a Bray-Curtis dissimilarity matrix from the species data.frame. adding “pa” , would make the matrix unweighted (presence absence). Other
  print(paste(col_name,": PCoA analysis"))
  pcoa.dist.col<-pcoa(dist.bray.col) 
  groups.col <- groups[groups[col_name]!="",]
  
  if(!is.null(col_name)){ 
    
    groups <- groups.col
    species_relab <- species_relab.col
    dist.bray <- dist.bray.col
    pcoa.dist <- pcoa.dist.col
    groups.single <- groups.col[[col_name]]
    
    lll <- length(unique(unlist(groups.single))[!is.na(unique(unlist(groups.single)))])
    cexsize <- 2
    if(lll > 20){
      cexsize <- 1.4
    }
    if(lll > 50){
      cexsize <- 1
    }
    if(lll > 100){
      cexsize <- 0.6
    }
    if(lll > 150){
      cexsize <- 0.4
    }
    if(lll > 150){
      cexsize <- 0.3
    }
    
    print(col_name)
    if( sum( grepl("\\|", unlist(groups[col_name])) ) == 0){
      # NO SPLIT ENTRIES
      groups.single            <- unlist(groups[col_name])
      groups.multi.super_group <- "NOPE"
      groups.multi.sub_group   <- "NOPE"
      pch_type <- 22
      if ( sum(grepl('[A-Za-z]', as.character(unlist(groups.single)) )) > 0 ){
        print("SINGLE F")
        col_type <- all_col[factor(as.character(unlist(groups.single)))]  
        png(filename=paste0(col_name,".pcoa.col_legend.DEFINED_samples.png"), width = 800, height = 800 )
        plot.new()
        legend("center",
               legend=unique(unlist(unname(as.list(groups.single)))),
               pch=22,
               cex=cexsize,
               border="white", bty = "n", fill=NULL, bg="white",
               #pt.bg=all_col[unique(unlist(unname(as.list(groups.single))))])
               pt.bg=unique(all_col[factor(as.character(unlist(groups.single)))]  )
        )
        dev.off()
      }else{
        print("SINGLE N")
        groups.single[as.character(unlist(groups[col_name])) == ""] <- 0
        pallette <- matlab.like2(max(as.numeric(groups.single[groups.single != "_"]), na.rm = TRUE))
        #pallette[-1] <- "#FFFFFF"
        col_type <- pallette[as.numeric(as.matrix(as.data.frame(groups.single)))]
        col_type[as.character(unlist(groups[col_name])) == ""] <- "#FFFFFF"
        png(filename=paste0(col_name,".pcoa.col_legend.DEFINED_samples.png"), width = 800, height = 800 )
        plot.new()
        legend("center",
               legend=sort(unique(unlist(unname(as.list(groups.single))))),
               pch=22,
               cex=cexsize,
               border="white", bty = "n", fill=NULL, bg="white",
               pt.bg=pallette[sort(unique(unlist(unname(as.list(groups.single)))))])
        dev.off()
      }
      
    }else{
      # SPLIT ENTRIES INCLUDED
      groups.single            <- "NOPE"
      assignment_vector <- as.character(unlist(groups[col_name]))
      assignment_vector[assignment_vector == ""] <- "_|_"
      groups.multi.super_group <- as.data.frame(strsplit(assignment_vector,'\\|'))[1,]
      pch_type <- c(21:25)[factor(as.character(unlist(groups.multi.super_group)))]
      png(filename=paste0(col_name,".pcoa.pch_legend.DEFINED_samples.png"), width = 800, height = 800 )
      plot.new()
      legend("center",
             legend=sort(unique(unlist(unname(as.list(groups.multi.super_group))))),
             pch=unique(c(21:25)[factor(as.character(unlist(groups.multi.super_group)))]),
             #pch=22,
             cex=cexsize,
             border="white", bty = "n", fill=NULL, bg="white",
             pt.bg="black")
      dev.off()
      
      groups.multi.sub_group   <- as.data.frame(strsplit(assignment_vector,'\\|'))[2,]
      png(filename=paste0(col_name,"pcoa.col_legend.png"), width = 800, height = 800 )
      if ( sum(grepl('[A-Za-z]', as.character(unlist(groups.multi.sub_group)) )) > 0 ){
        # as factor
        print("SUB F")
        col_type <- all_col[factor(as.character(unlist(groups.multi.sub_group)))]
        col_type[groups.multi.sub_group == "_"] <- "#FFFFFF"
        png(filename=paste0(col_name,".pcoa_LEGENDsub.DEFINED_samples.png"), width = 400, height = 700 )
        plot.new()
        legend("center",
               legend=levels(factor(as.character(unlist(groups.multi.sub_group)))),
               pch=22, 
               cex=cexsize,
               border="white", bty = "n", fill=NULL, bg="white",
               pt.bg=all_col [ unique(factor(as.character(unlist(groups.multi.sub_group)))) ] 
        )
        dev.off()
      }else{
        #as numeric
        print("SUB N")
        groups.multi.sub_group[as.character(unlist(groups[col_name])) == ""] <- 1
        pallette <- matlab.like2(max(as.numeric(groups.multi.sub_group[groups.multi.sub_group != "_"])))
        #pallette[-1] <- "#FFFFFF"
        col_type <- pallette[as.numeric(as.matrix(as.data.frame(groups.multi.sub_group)))]
        col_type[as.character(unlist(groups[col_name])) == ""] <- "#FFFFFF"
        groups.multi.sub_group[as.character(unlist(groups[col_name])) == ""] <- ""
        png(filename=paste0(col_name,".pcoa_LEGENDsub.DEFINED_samples.png"), width = 400, height = 700 )
        plot.new()
        legend("center",
               legend=sort(unique(as.numeric(as.matrix(as.data.frame(groups.multi.sub_group)))))[sort(unique(as.numeric(as.matrix(as.data.frame(groups.multi.sub_group))))) > 0],
               pch=22, cex=cexsize,
               border="white", bty = "n", fill=NULL, bg="white",
               pt.bg=pallette[ sort(unique(as.numeric(as.matrix(as.data.frame(groups.multi.sub_group)))))[sort(unique(as.numeric(as.matrix(as.data.frame(groups.multi.sub_group))))) > 0] ])
        dev.off()
      }
      
    }
    
    png(filename=paste0(col_name,".pcoa.DEFINED_samples.png"), width = 800, height = 800 )
    plot(pcoa.dist$vectors[,1:2],
         #bg=col_type,
         bg=alpha(col_type, 0.6),
         col="white",
         pch=pch_type,
         cex=2,
         alpha = 0.6,
         xlab=paste("PCoA.1(",round(100*pcoa.dist$values[1,3], digits=2),"%)"), 
         ylab=paste("PCoA.2(",round(100*pcoa.dist$values[2,3], digits=2),"%)"),
         #main="query_tab.infant.genus_20overZero\nPCoA, Unclustered\nNP, Infant-only, No Unclassified Bacteria"
         main=paste0("Defined Samples within Metadata Column\n",col_name)
    )
    if(dim(pcoa.dist$vectors)[1] <= 100){
      text( pcoa.dist$vectors, rownames( species ), pos=3, col="darkgray", cex=0.6 )
      text( pcoa.dist$vectors, as.character(unlist(groups[col_name])), pos=1, col="darkgray", cex=0.7 )
    }
    dev.off()
    
  }
  
  
  # Cluster and PCOA viz
  
  max_clust <- 12
  if( max_clust > (dim(pcoa.dist.col$vectors)[1]-1) ){
    max_clust <- dim(pcoa.dist.col$vectors)[1]-1
  }
  asw <- numeric(max_clust)
  c_asw <- numeric(max_clust)
  png(filename=paste0(col_name,".clust_count.png"), width = 800, height = 2400)
  par(mfrow = c(max_clust-1,3))
  for (num_clust in 2:max_clust){
    print(num_clust)
    pk <- pam(as.data.frame(pcoa.dist.col$vectors[,1:2]),num_clust)
    asw[num_clust] <- pk $ silinfo $ avg.width
    c_asw[num_clust] <- c(pk $ silinfo $clus.avg.widths)
    clusplot(pk,main=paste("\nClusters (PAM,PC1+2), N=",num_clust),col.p="grey",color=T,metric="euclidean",shade=T,labels=0,cex=0.05) #`pam` does you partitioning
    barplot(pk$clusinfo[,1],col="darkgrey",names.arg=c(1:dim(pk$clusinfo)[1]),main="Cluster Size")
    barplot(pk$silinfo$clus.avg.widths,col="darkgrey",names.arg=c(1:dim(pk$clusinfo)[1]),main="Cluster Avg Width")
  }
  dev.off()
  
  png(filename=paste0(col_name,".clust_countBest.png"), width = 1200, height = 800)
  par(mfrow = c(1,1))
  k.best <- which.max(asw)
  barplot(asw, type= "h", main = paste("Avg.Sil.Width\noptimal number of clusters:", k.best), xlab= "k  (# clusters)", ylab = "average silhouette width")
  axis(1, k.best, paste("best",k.best,sep="\n"), col = "red", col.axis = "red")
  dev.off()
  
  png(filename=paste0(col_name,".clustered_pcoa.png"), width = 800, height = 800 )
  par(mfrow = c(1,1))
  pk <- pam(as.data.frame(-pcoa.dist.col$vectors[,1:2]),k.best)
  asw[k.best] <- pk $ silinfo $ avg.width
  c_asw[k.best] <- c(pk $ silinfo $clus.avg.widths)
  clusplot(pk,
           color=T,
           col.p = col_type,
           #col.p=c(matlab.like(30)[simple.groups[rownames(query_tab.infant.genus_20overZero),]$Infant_Age_in_Months_Pass_NP],alpha=0.8),
           metric="euclidean",
           shade=F,labels=5,
           cex=2,
           pch_type=pch_type,
           xlab=paste("PCoA.1(",round(100*pcoa.dist.col$values[1,3], digits=2),"%)"), 
           ylab=paste("PCoA.2(",round(100*pcoa.dist.col$values[2,3], digits=2),"%)"),
           #main="query_tab.infant.genus_20overZero\nPCoA, 3 Clusters\nNP, Infant-only, No Unclassified Bacteria",
           main=paste0("Defined Samples within Metadata Column\n",col_name),
           lines=0) #`pam` does you partitioning
  legend("bottomright",legend = c(1:k.best),
         #col=c("blue","red","hotpink3"), 
         border="white", 
         #bty = "n", 
         pch = c(1:k.best),
         fill="white", bg="white", cex=2)
  dev.off()
  
  library(gplots)
  #if( groups.col.single == "NOPE" ){
  if( sum( grepl("\\|", unlist(groups[col_name])) ) != 0){
    #par(mfrow = c(k.best,2))
    counter_sub   <- na.replace( as.data.frame(matrix(ncol=length(unique(unname(unlist(groups.col.multi.sub_group))))  , nrow=k.best)) , 0 )
    counter_super <- na.replace( as.data.frame(matrix(ncol=length(unique(unname(unlist(groups.col.multi.super_group)))), nrow=k.best)) , 0 )
    rownames(counter_super) <- rownames(counter_sub) <- c(1:k.best)
    if( sum(grepl('[A-Za-z]', as.character(unlist(groups.col.multi.sub_group)) )) == 0 ){
      colnames(counter_sub) <- unique(sort(as.numeric(unlist(as.matrix(groups.col.multi.sub_group)))))
    }else{
      colnames(counter_sub) <- unique(unname(unlist(groups.col.multi.sub_group)))
    }
    colnames(counter_super) <- unique(unname(unlist(groups.col.multi.super_group)))
    for(r in 1:k.best){
      for(x in names(table(as.character(unlist(groups.col.multi.sub_group[ groups.col$SampleID %in% names(pk$clustering)[pk$clustering==r] ]))))){
        counter_sub[r,x] <- table(as.character(unlist(groups.col.multi.sub_group[ groups.col$SampleID %in% names(pk$clustering)[pk$clustering==r] ])))[x]
      }
      for(x in names(table(as.character(unlist(groups.col.multi.super_group[ groups.col$SampleID %in% names(pk$clustering)[pk$clustering==r] ]))))){
        counter_super[r,x] <- table(as.character(unlist(groups.col.multi.super_group[ groups.col$SampleID %in% names(pk$clustering)[pk$clustering==r] ])))[x]
      }
    }
    
    # absolute
    png(filename=paste0(col_name,".ClusteredFactors_SUB_ABS.CLUSTcount.",k.best,".png"), width = 800, height = 800)
    heatmap.2(t(as.matrix(counter_sub))  ,col=redgreen(max(counter_sub)+1),   Colv = NA, Rowv = NA, xlab="Sub Groups",    ylab="Cluster ID", trace="none", 
              main=paste(col_name,"\nAbsolute Abundance)"), cexRow = cexsize )
    dev.off()
    png(filename=paste0(col_name,".ClusteredFactors_SUPER_ABS.CLUSTcount.",k.best,".png"), width = 800, height = 800)
    heatmap.2(t(as.matrix(counter_super)),col=redgreen(max(counter_super)+1), Colv = NA, Rowv = NA, xlab="Super Groups",  ylab="Cluster ID", trace="none", 
              main=paste(col_name,"\nAbsolute Abundance)"), cexRow = cexsize )
    dev.off()
    
    # relative
    png(filename=paste0(col_name,".ClusteredFactors_SUB_REL.CLUSTcount.",k.best,".png"), width = 800, height = 800)
    heatmap.2(t(decostand(as.matrix(counter_sub), method="total"))  ,col=redgreen(100),   Colv = NA, Rowv = NA, xlab="Sub Groups",  ylab="Cluster ID", trace="none", 
              main=paste(col_name,"\nRelative Abundance)"), cexRow = cexsize )
    dev.off()
    png(filename=paste0(col_name,".ClusteredFactors_SUPER_REL.CLUSTcount.",k.best,".png"), width = 800, height = 800)
    heatmap.2(t(decostand(as.matrix(counter_super), method="total")),col=redgreen(100), Colv = NA, Rowv = NA, xlab="Super Groups",  ylab="Cluster ID", trace="none", 
              main=paste(col_name,"\nRelative Abundance)"), cexRow = cexsize )
    dev.off()
    
  }else{
    counter_single <- na.replace( as.data.frame(matrix(ncol=length(unique(unname(unlist(unlist(groups[col_name])))))  , nrow=k.best)) , 0 )
    if( length(unique(unname(unlist(unlist(groups[col_name]))))) > 1 ){
      counter_single   <- na.replace( as.data.frame(matrix(ncol=length(unique(unname(unlist(unlist(groups[col_name])))))  , nrow=k.best)) , 0 )
      if( sum(grepl('[A-Za-z]', as.character(unlist(unlist(groups[col_name]))) )) == 0 ){
        colnames(counter_single) <- unique(sort(as.numeric(unlist(as.matrix(unlist(groups[col_name]))))))
      }else{
        colnames(counter_single) <- unique(unname(unlist(unlist(groups[col_name]))))
      }
      for(r in 1:k.best){
        for(x in names(table(as.character(unlist(unlist(groups[col_name])[ groups.col$SampleID %in% names(pk$clustering)[pk$clustering==r] ]))))){
          counter_single[r,x] <- table(as.character(unlist(unlist(groups[col_name])[ groups.col$SampleID %in% names(pk$clustering)[pk$clustering==r] ])))[x]
        }
      }
      png(filename=paste0(col_name,".ClusteredFactors_SINGLE_ABS.CLUSTcount.",k.best,".png"), width = 800, height = 800)
      heatmap.2(t(as.matrix(counter_single)),col=redgreen(max(counter_single)+1), 
                Colv = NA, Rowv = NA, ylab="Groups",  xlab="Cluster ID", trace="none", main=paste(col_name,"\nAbsolute Abundance"),
                margins=c(5,30), keysize=0.8, cexRow = cexsize )
      dev.off()
      png(filename=paste0(col_name,".ClusteredFactors_SINGLE_REL.CLUSTcount.",k.best,".png"), width = 800, height = 800)
      heatmap.2(t(decostand(as.matrix(counter_single), method="total")),col=redgreen(max(counter_single)+1), 
                Colv = NA, Rowv = NA, ylab="Groups",  xlab="Cluster ID", trace="none", main=paste(col_name,"\nRelative Abundance"),
                margins=c(5,30), keysize=0.8, cexRow = cexsize )
      dev.off()
    }else{
      # 1 entry - no need for heatmap
    }
  }
   
  # VERIFY CLUSTER COUNTS WITH PLS-DA
  
  library(mixOmics)
  plsda_overall_err <- list()
  pi <- 1
  if( dim(pcoa.dist.col$vectors)[1] > 5 ){
    #DEBUG: REMOVED DUE TO SLOW PROCESSING
    #for (num_clust in 2:max_clust){
      num_clust <- k.best
      
      print(num_clust)
      if(sum(table(pam(as.data.frame(pcoa.dist.col$vectors[,1:2]),num_clust)$clustering) == 1) == 0){
        pk <- pam(as.data.frame(pcoa.dist.col$vectors[,1:2]),num_clust)
        X <- pcoa.dist.col$vectors
        Y <- as.factor(pk$clustering)
        if(grepl("Error", try( plsda(X, Y, ncomp = 5) , TRUE )) == FALSE){
          plsda.res <- plsda(X, Y, ncomp = 5)
          perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = TRUE, auc = TRUE, nrepeat = 10) 
          png(filename=paste0("PLSDA.",col_name,".CLUSTcount.",num_clust,".png"), width = 800, height = 800)
          plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal", main=paste0("x=DISTANCE MATRIX\ny=CLUSTER ASSIGNMENTS AT N=",num_clust))
          write.table(perf.plsda$error.rate, file=paste0(col_name,".PCoA_PLSDAerrorRate_CLUSTcount.",num_clust,".tsv"))
          dev.off()
          plsda_overall_err[pi] <- min(perf.plsda$error.rate$overall[,1])
          pi <- pi + 1
        }
      }else{
        plsda_overall_err[pi] <- 0
        pi <- pi + 1
        print("skip")
      }
    #DEBUG: REMOVED DUE TO SLOW PROCESSING
    #}
    #names(plsda_overall_err) <- c(2:max_clust)
    #png(filename=paste0(col_name,"plsda_clustered_error.png"), width = 800, height = 800 )
    #barplot(unlist(plsda_overall_err),main=paste0(col_name,"\nPLSA-DA Re-sampling Error Rate At N=2 to ",max_clust))# , names.arg = c(2:max_clust))
    #dev.off()
    #plsda.k.best <- as.numeric(names(sort(unlist(plsda_overall_err))[ sort(unlist(plsda_overall_err)) != 0][1]))
    # 
    # if(k.best != plsda.k.best){
    #   
    #   png(filename=paste0(col_name,".sclustered_pcoa.AtPLSDAClusterCountEq",plsa.k.best,".png"), width = 800, height = 800 )
    #   par(mfrow = c(1,1))
    #   pk <- pam(as.data.frame(-pcoa.dist.col$vectors[,1:2]),plsda.k.best)
    #   asw[k.best] <- pk $ silinfo $ avg.width
    #   c_asw[k.best] <- c(pk $ silinfo $clus.avg.widths)
    #   clusplot(pk,
    #            color=T,
    #            col.p = col_type,
    #            #col.p=c(matlab.like(30)[simple.groups[rownames(query_tab.infant.genus_20overZero),]$Infant_Age_in_Months_Pass_NP],alpha=0.8),
    #            metric="euclidean",
    #            shade=F,labels=5,
    #            cex=2,
    #            pch_type=pch_type,
    #            xlab=paste("PCoA.1(",round(100*pcoa.dist.col$values[1,3], digits=2),"%)"), 
    #            ylab=paste("PCoA.2(",round(100*pcoa.dist.col$values[2,3], digits=2),"%)"),
    #            main=paste0("Samples Defined by Metadata Column\n",col_name) #`pam` does you partitioning
    #            )
    #   legend("bottomright",legend = c(1:plsda.k.best),
    #          #col=c("blue","red","hotpink3"), 
    #          border="white", 
    #          #bty = "n", 
    #          pch = c(1:plsda.k.best),
    #          fill="white", bg="white", 
    #          cex=2
    #          )
    #   dev.off()
    #   
    #   
    # }else{
    #   print("they match")
    # }
    
  }else{
    print("too few samples")
  }
  
} 
}

###########################################
# SPECIES ABUNDANCE AND DISTANCE HEATMAPS #
###########################################

labels <- list()

li <- 1
for(n in colnames(species_relab)){
  n.orig <- n
  #labels[li] <- tax_summary
  while(dim( tax_summary[ grepl( n, tax_summary$rankID ), ] )[1] == 0 || tax_summary[ grepl( n, tax_summary$rankID ), ]$taxon[1] == "unclassified"){
    n <- substr(n, 1, nchar(n)-1)
  }
  labels[li] <- paste(as.character(tax_summary[ grepl( n, tax_summary$rankID ), ][1,]$taxon),"(",n," <-",n.orig,")")
  li <- li + 1
}

cc <- 1.5
if(dim(species_relab)[2] < 150){
  cc = 1.5
  wc = 3600
}else if(dim(species_relab)[2] < 300){
  cc = 1
  wc = 4200
}else if(dim(species_relab)[2] < 600){
  cc = 0.7
  wc = 5600
}else{
  cc = 0.5
  wc = 5600
}
if(dim(species_relab)[1] < 20){
  hc = 3600
  cr = 4
}else if(dim(species_relab)[1] < 40){
  hc = 3600
  cr = 3
}else if(dim(species_relab)[1] < 100){
  hc = 4000
  cr = 2
}else if(dim(species_relab)[1] < 250){
  hc = 4600
  cr = 1
}else{
  hc = 4600
  cr = 0.5
}
png(filename="heat_relAbsGenus.png", width = wc, height = hc)
pkam <- pam(-log(as.matrix(species_relab)+0.00001),k.best)
h.tmp <- heatmap.2(-log(as.matrix(species_relab)+0.00001))
par(cex.main=4)
heatmap.2(-log(as.matrix(species_relab)+0.00001),
          col=c(rev(gray.colors(150)),matlab.like2(150)),
          trace="none",
          margins = c(60,30),
          keysize=0.7, key.par = list(cex=0.5),
          cexCol = cc,
          cexRow = cr,
          #cex.main=4,
          labCol = labels,
          main="\n\nEuclidean Abundance Heatmap\n-log( Genus Level Abundance )\n",
          #RowSideColors = all_col[pkam$clustering]
          colRow = c("darkred","darkblue","purple","darkgreen","gray")[cutree( as.hclust(h.tmp$colDendrogram),2)],
          RowSideColors = all_col[pkam$clustering],
          ColSideColors = primary.colors(10)[cutree( as.hclust(h.tmp$colDendrogram),10)] 
#legend("bottomright",legend=c("Samples:","left = Euclidean Hierarchical Clustering of Relative Abundances","right = PAM-derived PCoA Clusters of Bray-Curtis Distances","Genus IDs:","top = Euclidean Hierarchical Clustering of Sample Representation"))
)
dev.off()
   
  


#install.packages("heatmap.plus",repos='http://cran.us.r-project.org')


countem <- 0
coln <- list()
cn <- 1
for(c in colnames(groups)[4:dim(groups)[2]]){
  if( sum( grepl("\\|", unlist(groups[,c]) ) ) == 0){
    coln[cn] <- c
    cn <- cn + 1
  }else{
    coln[cn] <- paste0(c,".super")
    cn <- cn + 1
    coln[cn] <- paste0(c,".sub")
    cn <- cn + 1
    countem <- countem + 1
  }
}
color_matrix <- matrix(nrow=dim(groups)[1],ncol=(dim(groups)[2]-3+countem))
rownames(color_matrix) <- groups$SampleID
colnames(color_matrix) <- coln
  #colnames(groups)[4:dim(groups)[2]]
 
for(c in colnames(groups)[4:dim(groups)[2]]){ #colnames(color_matrix)){
  if( sum( grepl("\\|", unlist(groups[,c]) )) == 0){
    if ( sum(grepl('[A-Za-z]',  unlist(groups[,c]) )) > 0 ){
      color_matrix[,c] <- all_col[groups[,c]]
      color_matrix[is.na(color_matrix)] <- "#FFFFFF"
    }else{
      tmpgroup <- all_col[groups[,c]]
      tmpgroup[tmpgroup == "_"] <- ""
      col.ini <- matlab.like(max(na.replace(as.numeric(as.matrix(tmpgroup[tmpgroup!="_"])),0))+1)[as.numeric(as.matrix(tmpgroup))+1]
      col.ini[ is.na(col.ini) ] <- "#FFFFFF"
      color_matrix[,c] <- col.ini
    }
  }else{
    #groups[as.character(unlist(groups[,c])) == "",c] <- "_|_"
    tx <- as.character(unlist(groups[,c]))
    tx[tx == ""] <-  "_|_"
    tx[is.na(tx)] <-  "_|_"
    tmpsuper <- unlist(as.data.frame(strsplit(as.character(tx),'\\|'))[1,])
    if ( sum(grepl('[A-Za-z]', tmpsuper)) > 0 ){
      color_matrix[,paste0(c,".super")] <- all_col[tmpsuper]
      color_matrix[is.na(color_matrix)] <- "#FFFFFF"
    }else{
      col.ini <- matlab.like(max(na.replace(as.numeric(as.matrix(tmpsub[tmpsub!="_"])),0)))[as.numeric(as.matrix(tmpsuper))]
      col.ini[ is.na(col.ini) ] <- "#FFFFFF"
      color_matrix[,paste0(c,".super")] <- col.ini
    }
    tmpsub <- unlist(as.data.frame(strsplit(as.character(tx),'\\|'))[2,])
    if ( sum(grepl('[A-Za-z]', tmpsub)) > 0 ){
      color_matrix[,paste0(c,".sub")]   <- all_col[tmpsub]
      color_matrix[is.na(color_matrix)] <- "#FFFFFF"
    }else{
      tmpsub[tmpsub == "_"] <- ""
      col.ini <- matlab.like(max(na.replace(as.numeric(as.matrix(tmpsub[tmpsub!="_"])),0))+1)[as.numeric(as.matrix(tmpsub))+1]
      col.ini[ is.na(col.ini) ] <- "#FFFFFF"
      color_matrix[,paste0(c,".sub")] <- col.ini
    }
  }
} 
color_matrix <- color_matrix[rownames(species_relab),]
png(filename="heatPlus_relAbsGenus.png", width = (wc *4), height = (hc*4))
heatmap.plus(-log(as.matrix(species_relab)+0.00001),
          #dist="bray",
          col=c(rev(matlab.like2(150)),gray.colors(150)),
          margins = c(60,30),
          #keysize=0.7, key.par = list(cex=0.5),
          cexCol = cc,
          cexRow = cr,
          labCol = labels,
          main="\n\nEuclidean Abundance Heatmap\n-log( Genus Level Abundance )\n",
          RowSideColors = color_matrix#)#,
          #ColSideColors = as.matrix(primary.colors(10)[cutree( as.hclust(h.tmp$colDendrogram),10)] )
)
dev.off()
  


png(filename="heatPlus_relAbsGenus.small.png", width = 2400, height = 800)
heatmap.plus(-log(as.matrix(species_relab)+0.00001),
             #dist="bray",
             col=c(rev(matlab.like2(150)),gray.colors(150)),
             margins = c(15,9),
             keysize=0.7, key.par = list(cex=0.5),
             cexCol = cc,
             cexRow = 2,
             labCol = "", #labels,
             main="\n\n\n\nEuclidean Abundance Heatmap\n-log( Genus Level Abundance )\n",
             RowSideColors = color_matrix#)#,
             #ColSideColors = as.matrix(primary.colors(10)[cutree( as.hclust(h.tmp$colDendrogram),10)] )
)
dev.off()


png(filename="heat_BayCurtisDistGenus.png", width = 2400, height = 2400)
par(cex.main=4)
heatmap.2(as.matrix(dist.bray),
          col=c(rev(gray.colors(150)),matlab.like2(150)),
          trace="none",
          margins = c(20,20),
          keysize=0.7, key.par = list(cex=0.5),
          cexCol = 4,
          cexRow = 4,
          main="\n\nBray-Curtis Distance Heatmap\nGenus Level Abundance"
)
dev.off()

png(filename="heatPlus_BrayDist.small.png", width = 1600, height = 800)
heatmap.plus(as.matrix(dist.bray),
             #dist="bray",
             col=c(rev(gray.colors(150)),matlab.like2(150)),
             margins = c(15,9),
             keysize=0.7, key.par = list(cex=0.5),
             cexCol = 2,
             cexRow = 2,
             #labCol = "", #labels,
             main="\n\n\n\nBray Curtis Distance Heatmap",
             RowSideColors = color_matrix#)#,
             #ColSideColors = as.matrix(primary.colors(10)[cutree( as.hclust(h.tmp$colDendrogram),10)] )
)
dev.off()

png(filename="heatPlus_BayCurtisDistGenus.png", width = 2460, height = 2400)
par(cex.main=4)
heatmap.plus(as.matrix(dist.bray),
          col=c(rev(gray.colors(150)),matlab.like2(150)),
          trace="none",
          margins = c(20,20),
          keysize=0.7, key.par = list(cex=0.5),
          cexCol = 4,
          cexRow = 4,
          main="\n\nBray-Curtis Distance Heatmap\nGenus Level Abundance",
          RowSideColors = color_matrix#)#,
)
dev.off()
 

#######################################################
# SIGNIFICANT AND MAXIMUM SPECIES ABUNDANCE BARGRAPHS #
#######################################################

pc.pval <- list()
pc.pval.fdr <- list()
pi <- 1
for(i in colnames(species_relab)){
  reg <- lm( as.matrix(species_relab[,i]) ~ pcoa.dist$vectors[,1:2] + 1)
    if (t(coef(summary(reg)))[2,1] > 0 & is.nan(t(coef(summary(reg)))[2,1]) == 0){
      cpus.lm2 <- stepAIC(reg, trace = FALSE)
      fstat <- summary(reg)$fstatistic
      pc.pval[pi] <- unname(t(as.matrix(coef(summary(reg))[,4])))[1]
      pc.pval.fdr[pi] <- p.adjust( pc.pval[pi] , method="fdr" )
    }else{
      pc.pval[pi] <- "" #matrix(0,1,length(levels(factor(unlist(factor_sample_super_cluster_type)))))
      pc.pval.fdr[pi] <- "" #p.adjust( pc.pval[i] , method="fdr" )
    }
  pi <- pi + 1
}
names(pc.pval) <- colnames(species_relab)
names(pc.pval.fdr) <- colnames(species_relab)
  
nml <- as.numeric(as.list(as.matrix(pc.pval)))
names(nml) <- rownames(as.list(as.matrix(pc.pval)))
labels <- list()
li <- 1
for(n in names(nml)){
  n.orig <- n
  #labels[li] <- tax_summary
  while(dim( tax_summary[ grepl( n, tax_summary$rankID ), ] )[1] == 0 || tax_summary[ grepl( n, tax_summary$rankID ), ]$taxon[1] == "unclassified"){
    n <- substr(n, 1, nchar(n)-1)
  }
  labels[li] <- paste(as.character(tax_summary[ grepl( n, tax_summary$rankID ), ][1,]$taxon),"(",n.orig,")\n(",n,")")
  li <- li + 1
}
labels <- unlist(labels)

png(filename="pval_PC_lm_all.png", width = 800, height = 1600)
plot(-log(sort(nml)), xlim = c(-(length(nml)*0.2),length(nml)+(length(nml)*0.2)), main="\nALL\n-log(Pvalue from lm)\nSpecies Most Significantly Predicting All Samples'\nTaxa Abundance-derived Orientation",
     xlab=paste( sum(-log(sort(nml)) > -log(cutoff)), "significant genus level taxa"))
text(-log(sort(nml)), labels = c(labels), col=c("darkblue","darkred"),cex=0.8)
par(new=TRUE)
abline( h=-log(cutoff) , col="darkgreen", lwd=3 )
dev.off()

png(filename="pval_PC_lm_signif.png", width = 800, height = 1600)
plot(-log(sort(nml[names(nml[nml < cutoff])])), xlim = c(-(length(nml[nml < cutoff])*0.2),length(nml[nml < cutoff])+(length(nml[nml < cutoff])*0.2)), main="\nSIGNFICANT\n-log(Pvalue from lm)\nSpecies Most Significantly Predicting All Samples'\nTaxa Abundance-derived Orientation")
text(-log(sort(nml[names(nml[nml < cutoff])])), labels = c(labels), col=c("darkblue","darkred"),cex=0.8)
#par(new=TRUE)
#abline( h=-log(cutoff) , col="darkgreen", lwd=3 )
dev.off()

png(filename="pval_PC_lm_top25.png", width = 800, height = 1600)
plot(-log(sort(nml)[1:25]+0.000001), xlim = c(-10,35), main=paste("\nTOP 25(min p =",min(sort(nml)[1:25]),")\n-log(Pvalue from lm)\nSpecies Most Significantly Predicting All Samples'\nTaxa Abundance-derived Orientation"))
text(-log(sort(nml)[1:25]+0.000001), labels = c(labels), col=c("darkblue","darkred"),cex=0.8)
par(new=TRUE)
abline( h=-log(cutoff) , col="darkgreen", lwd=3 )
dev.off()
 
png(filename="pval_PC_lm_25maxAbundance.png", width = 800, height = 1600)
plot(-log(sort(nml[names(rev(sort(colSums(species_relab))[1:25]))])), xlim = c(-10,35), main=paste("\nTOP 25(min p =",min(sort(nml)[1:25]),")\n-log(Pvalue from lm)\nSpecies Most Significantly Predicting All Samples'\nTaxa Abundance-derived Orientation"))
text(-log(sort(nml[names(rev(sort(colSums(species_relab))[1:25]))])), labels = c(labels), col=c("darkblue","darkred"),cex=0.8)
par(new=TRUE)
abline( h=-log(cutoff) , col="darkgreen", lwd=3 )
dev.off()
 
L <- length(colSums(species_relab))
if(L > 25){ L <- 25 }
max_all <- decostand( t(species_relab) ,method="total")[names(sort(nml)),]
max_abs  <- decostand( t(species_relab) , method = "total" )[names(rev(sort(colSums(species_relab))[1:L])),]
max_signifP <- decostand( t(species_relab) , method = "total" )[ names(sort( nml[names(nml[nml < cutoff])])) , ] 
max_pval <- decostand( t(species_relab) , method = "total" )[ names(sort(nml)[1:L]) , ] 
 
QL <- length(names(rev(sort(t(decostand(t(max_abs), method="total"))[1,]))))
wd <- 1200
if(QL > 50){ wd <- 1800}
if(QL > 100){ wd <- 2000}
if(QL > 200){ wd <- 2200}
if(QL > 350){ wd <- 2600}
if(QL > 700){ wd <- 3200}
if(QL > 1400){ wd <- 4200}
png(filename="species_relAbs_25maxAbundance.png", width = wd, height = 2400)
par(mfrow = c(2,1))
par(mar = c(0, 8, 8, 30))
barplot(max_abs[,names(rev(sort(t(decostand(t(max_abs), method="total"))[1,])))],col=all_col,main="Maximum Species Abundance\n(Absolute)",las=3,cex.names = 0.01) 
par(new=TRUE)
par(mar = c(0, 0, 0, 0))
nml <- rownames(max_abs)
labels <- list()
li <- 1
for(n in nml){
  n.orig <- n
  #labels[li] <- tax_summary
  while(dim( tax_summary[ grepl( n, tax_summary$rankID ), ] )[1] == 0 || tax_summary[ grepl( n, tax_summary$rankID ), ]$taxon[1] == "unclassified"){
    n <- substr(n, 1, nchar(n)-1)
  }
  labels[li] <- paste(as.character(tax_summary[ grepl( n, tax_summary$rankID ), ][1,]$taxon),"(",n.orig,")\n(",n,")")
  li <- li + 1
}
labels <- unlist(labels)
if(wd == 1200){
  legend("right", #
         inset=c(-0.51,-0.05),
         legend=labels,pch=22,
         pt.bg=replicate(all_col,n = 10)[1:length(labels)],
         pt.cex=3,
         col=replicate(all_col,n = 10)[1:length(labels)],
         cex=0.8,border="white", bty = "n", fill=NULL, bg="white",
        y.intersp = 3   ) 
}else{
  legend("right", #
         legend=labels,pch=22,
         pt.bg=replicate(all_col,n = 10)[1:length(labels)],
         pt.cex=3,
         col=replicate(all_col,n = 10)[1:length(labels)],
         cex=0.8,border="white", bty = "n", fill=NULL, bg="white",
         y.intersp = 3   ) 
}
par(mar = c(11, 8, 1, 30))
#barplot(t(decostand(t(max_abs), method="total"))[,names(rev(sort(t(decostand(t(max_abs), method="total"))[1,])))],col=all_col,main="(Relative)",las=3) 
plotme <- t(decostand(t(max_abs), method="total")[ rowSums(decostand(t(max_abs), method="total")) > 0, ])
barplot(plotme,col=all_col,main="(Relative)",las=3)
dev.off()
    
png(filename="species_relAbs_25maxSignif.png", width = wd*5, height = 2400)
par(mfrow = c(2,1))
par(mar = c(0, 8, 8, 30))
barplot(max_pval[,names(rev(sort(t(decostand(t(max_pval), method="total"))[1,])))], border = FALSE, col=all_col,main="Maximum Species Abundance\n(Absolute)",las=3,cex.names = 0.01) 
par(new=TRUE)
par(mar = c(0, 0, 0, 0))
nml <- rownames(max_pval)
labels <- list()
li <- 1
for(n in nml){
  n.orig <- n
  #labels[li] <- tax_summary
  while(dim( tax_summary[ grepl( n, tax_summary$rankID ), ] )[1] == 0 || tax_summary[ grepl( n, tax_summary$rankID ), ]$taxon[1] == "unclassified"){
    n <- substr(n, 1, nchar(n)-1)
  }
  labels[li] <- paste(as.character(tax_summary[ grepl( n, tax_summary$rankID ), ][1,]$taxon),"(",n.orig,")\n(",n,")")
  li <- li + 1
}
labels <- unlist(labels)
legend("right", #
       inset=c(-0.51,-0.05),
       legend=labels,pch=22,
       pt.bg=replicate(all_col,n = 10)[1:length(labels)],
       pt.cex=3,
       col=replicate(all_col,n = 10)[1:length(labels)],
       cex=0.8,border="white", bty = "n", fill=NULL, bg="white",
       y.intersp = 3   ) 

par(mar = c(11, 8, 1, 30))
barplot(t(decostand(t(max_pval), method="total"))[,names(rev(sort(t(decostand(t(max_pval), method="total"))[1,])))], border = FALSE,col=all_col,main="(Relative)",las=3) 
dev.off()


  