setwd("C:\\Users\\jamis\\Desktop\\")
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
 
load("oops.Rdata") # <- figures_20200314.R

topDiffGenes <- function(allScore) {  return(allScore > 1)   }
#bm <- useMart("ensembl")
#bm <- useDataset("hsapiens_gene_ensembl", mart=bm)
bm <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")

###################################################################################################################################
####################################### TABLE 1 : NUMBER OF SIGNIF AT EACH CUTOFF (+SHARED) #######################################
###################################################################################################################################
  
  # save.image("20200320.Rdata")
  #load("20200320.Rdata")
  
############ 1.1 : Identify Significant OGs ############ 
  
v9_v10_OGs_map     <- read.table("v9_v10_OGs_map.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
odb10v1_OG_xrefs     <- read.table("odb10v1_OG_xrefs.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
odb10v1_OG2genes     <- read.table("odb10v1_OG2genes.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
  
# pull colnames
# pp <-  c( colnames(model.p.dos)[ grepl('p_val$', colnames(model.p.dos)) ],  colnames(model.p.dos)[ grepl('pv$', colnames(model.p.dos)) ] )

significant_og_counter_matrix.f_01 <- matrix(nrow=24,ncol=12)
significant_og_counter_matrix.f_02 <- matrix(nrow=24,ncol=12)
significant_og_counter_matrix.f_04 <- matrix(nrow=24,ncol=12)
significant_og_counter_matrix.f_05 <- matrix(nrow=24,ncol=12)
significant_og_counter_matrix.shared <- matrix(nrow=424,ncol=12)

#cutoff <- 0.05
for(cutoff in c(0.05,0.005)){
  # LM
  for(mod in 1:12){
        mtc <- 1
        for(met in c("f_01","f_02","f_04","f_05")){   # THEN RE_RUN TO SUM OVER 2 and 5!
          
            c.align.pv.lm <-  paste0( mod,".",met,".LM.p_val.FDR" )
            c.align.pv.caper <-  paste0( mod,".",met,".CAPER.p_val.FDR" )
            c.align.pv.mdmr_mash <-  paste0( mod,".",met,".MDMR_mash.","pv.FDR" )
            c.align.pv.mdmr_phylo <-  paste0( mod,".",met,".MDMR_phylo.","pv.FDR" )
            c.align.slp.lm <- paste0( mod,".",met,".LM.slope") 
            c.align.pos.lm <- paste0( mod,".",met,".LM.is_pos") 
            c.align.slp.caper <- paste0( mod,".",met,".CAPER.slope") 
            c.align.pos.caper <- paste0( mod,".",met,".CAPER.is_pos") 
            
            c.denov.pv.lm <-  paste0( mod+12,".",met,".LM.p_val.FDR" )
            c.denov.pv.caper <-  paste0( mod+12,".",met,".CAPER.p_val.FDR" )
            c.denov.pv.mdmr_mash <-  paste0( mod+12,".",met,".MDMR_mash.","pv.FDR" )
            c.denov.pv.mdmr_phylo <-  paste0( mod+12,".",met,".MDMR_phylo.","pv.FDR" )
            c.denov.slp.lm <- paste0( mod+12,".",met,".LM.slope") 
            c.denov.pos.lm <- paste0( mod+12,".",met,".LM.is_pos") 
            c.denov.slp.caper <- paste0( mod+12,".",met,".CAPER.slope") 
            c.denov.pos.caper <- paste0( mod+12,".",met,".CAPER.is_pos") 
            
            cmpr <- model.p.dos[ , c.align.pv.lm ] < cutoff
            cmpr[is.na(cmpr)] <- FALSE
            ltCutoff <- rownames(model.p.dos)[  cmpr ]
            ltCutoff.neg.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.lm] < 0]
            ltCutoff.pos.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.lm] > 0]
            ltCutoff.neg.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.lm] < 0]
            ltCutoff.pos.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.lm] > 0]
            if(met == "f_01"){
              significant_og_counter_matrix.f_01 [  1 ,  mod ] <- length(ltCutoff.pos.denov)
              significant_og_counter_matrix.f_01[  2 ,  mod ] <- length(ltCutoff.neg.denov)
              significant_og_counter_matrix.f_01[  3 ,  mod ] <- length(ltCutoff.pos.align)
              significant_og_counter_matrix.f_01[  4 ,  mod ] <- length(ltCutoff.neg.align)
              significant_og_counter_matrix.f_01[  5 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
              significant_og_counter_matrix.f_01[  6 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
              ltCutoff.neg.align.f_01 <- ltCutoff.neg.align
              ltCutoff.pos.align.f_01 <- ltCutoff.pos.align
              ltCutoff.neg.denov.f_01 <- ltCutoff.neg.denov
              ltCutoff.pos.denov.f_01 <- ltCutoff.pos.denov
              ltCutoff.neg.shared.f_01 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
              ltCutoff.pos.shared.f_01 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
            }
            if(met == "f_02"){
              significant_og_counter_matrix.f_02 [  1 ,  mod ] <- length(ltCutoff.pos.denov)
              significant_og_counter_matrix.f_02 [  2 ,  mod ] <- length(ltCutoff.neg.denov)
              significant_og_counter_matrix.f_02 [  3 ,  mod ] <- length(ltCutoff.pos.align)
              significant_og_counter_matrix.f_02 [  4 ,  mod ] <- length(ltCutoff.neg.align)
              significant_og_counter_matrix.f_02 [  5 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
              significant_og_counter_matrix.f_02 [  6 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
              ltCutoff.neg.align.f_02 <- ltCutoff.neg.align
              ltCutoff.pos.align.f_02 <- ltCutoff.pos.align
              ltCutoff.neg.denov.f_02 <- ltCutoff.neg.denov
              ltCutoff.pos.denov.f_02 <- ltCutoff.pos.denov
              ltCutoff.neg.shared.f_02 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
              ltCutoff.pos.shared.f_02 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
            }
            if(met == "f_04"){
              significant_og_counter_matrix.f_04 [  1 ,  mod ] <- length(ltCutoff.pos.denov)
              significant_og_counter_matrix.f_04 [  2 ,  mod ] <- length(ltCutoff.neg.denov)
              significant_og_counter_matrix.f_04 [  3 ,  mod ] <- length(ltCutoff.pos.align)
              significant_og_counter_matrix.f_04 [  4 ,  mod ] <- length(ltCutoff.neg.align)
              significant_og_counter_matrix.f_04 [  5 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
              significant_og_counter_matrix.f_04 [  6 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
              ltCutoff.neg.align.f_04 <- ltCutoff.neg.align
              ltCutoff.pos.align.f_04 <- ltCutoff.pos.align
              ltCutoff.neg.denov.f_04 <- ltCutoff.neg.denov
              ltCutoff.pos.denov.f_04 <- ltCutoff.pos.denov
              ltCutoff.neg.shared.f_04 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
              ltCutoff.pos.shared.f_04 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
            }
            if(met == "f_05"){
              significant_og_counter_matrix.f_05 [  1 ,  mod ] <- length(ltCutoff.pos.denov)
              significant_og_counter_matrix.f_05 [  2 ,  mod ] <- length(ltCutoff.neg.denov)
              significant_og_counter_matrix.f_05 [  3 ,  mod ] <- length(ltCutoff.pos.align)
              significant_og_counter_matrix.f_05 [  4 ,  mod ] <- length(ltCutoff.neg.align)
              significant_og_counter_matrix.f_05 [  5 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
              significant_og_counter_matrix.f_05 [  6 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
              ltCutoff.neg.align.f_05 <- ltCutoff.neg.align
              ltCutoff.pos.align.f_05 <- ltCutoff.pos.align
              ltCutoff.neg.denov.f_05 <- ltCutoff.neg.denov
              ltCutoff.pos.denov.f_05 <- ltCutoff.pos.denov
              ltCutoff.neg.shared.f_05 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
              ltCutoff.pos.shared.f_05 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
            }
            mtc <- mtc + 1
            
        }
        
        a <- ltCutoff.neg.denov.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
        a <- a[a %in% ltCutoff.neg.align.f_02]
        ltCutoff.neg.denov.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
        significant_og_counter_matrix.shared[ 1 , mod ] <- length( ltCutoff.neg.denov.f_shared )
          
        a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
        a <- a[a %in% ltCutoff.pos.align.f_02]
        ltCutoff.pos.denov.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
        significant_og_counter_matrix.shared[ 2 , mod ] <- length( ltCutoff.pos.denov.f_shared )
        
        a <- ltCutoff.neg.align.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
        a <- a[a %in% ltCutoff.neg.align.f_02]
        ltCutoff.neg.align.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
        significant_og_counter_matrix.shared[ 3 , mod ] <- length( ltCutoff.neg.align.f_shared )
        
        a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
        a <- a[a %in% ltCutoff.pos.align.f_02]
        ltCutoff.pos.align.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
        significant_og_counter_matrix.shared[ 4 , mod ] <- length( ltCutoff.pos.align.f_shared )
        
        a <- ltCutoff.neg.shared.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
        a <- a[a %in% ltCutoff.neg.align.f_02]
        ltCutoff.neg.shared.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
        significant_og_counter_matrix.shared[ 5 , mod ] <- length( ltCutoff.neg.shared.f_shared )
        
        a <- ltCutoff.pos.shared.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
        a <- a[a %in% ltCutoff.pos.align.f_02]
        ltCutoff.pos.shared.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
        significant_og_counter_matrix.shared[ 6 , mod ] <- length( ltCutoff.pos.shared.f_shared )
        
        #write.table(significant_og_counter_matrix.shared, paste0("significant_og_counter_matrix.",mod,".f_05.tsv"), sep='\t')
        
        ltCutoff.neg.align.f_shared.lm <- ltCutoff.neg.align
        ltCutoff.pos.align.f_shared.lm <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_shared.lm <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_shared.lm <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_shared.lm <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_shared.lm <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
        
        #for(met in c("f_01","f_02","f_04","f_05")){ 
        
        write.table(ltCutoff.neg.align.f_shared.lm,   paste0("outConvert/ltCutoff.neg.align.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.pos.align.f_shared.lm,   paste0("outConvert/ltCutoff.pos.align.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.neg.denov.f_shared.lm,   paste0("outConvert/ltCutoff.neg.denov.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.pos.denov.f_shared.lm,   paste0("outConvert/ltCutoff.pos.denov.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.neg.shared.f_shared.lm,   paste0("outConvert/ltCutoff.neg.shared.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.pos.shared.f_shared.lm,   paste0("outConvert/ltCutoff.pos.shared.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        met <- "f_01"
        write.table( ltCutoff.neg.align.f_01,   paste0("outConvert/ltCutoff.neg.align.f_01.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.pos.align.f_01,   paste0("outConvert/ltCutoff.pos.align.f_01.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.neg.denov.f_01,   paste0("outConvert/ltCutoff.neg.denov.f_01.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.pos.denov.f_01,   paste0("outConvert/ltCutoff.pos.denov.f_01.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.neg.shared.f_01,   paste0("outConvert/ltCutoff.neg.shared.f_01.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.pos.shared.f_01,   paste0("outConvert/ltCutoff.pos.shared.f_01.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        met <- "f_02"
        write.table( ltCutoff.neg.align.f_02,   paste0("outConvert/ltCutoff.neg.align.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.pos.align.f_02,   paste0("outConvert/ltCutoff.pos.align.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.neg.denov.f_02,   paste0("outConvert/ltCutoff.neg.denov.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.pos.denov.f_02,   paste0("outConvert/ltCutoff.pos.denov.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.neg.shared.f_02,   paste0("outConvert/ltCutoff.neg.shared.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.pos.shared.f_02,   paste0("outConvert/ltCutoff.pos.shared.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        met <- "f_04"
        write.table( ltCutoff.neg.align.f_04,   paste0("outConvert/ltCutoff.neg.align.f_04.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.pos.align.f_04,   paste0("outConvert/ltCutoff.pos.align.f_04.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.neg.denov.f_04,   paste0("outConvert/ltCutoff.neg.denov.f_04.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.pos.denov.f_04,   paste0("outConvert/ltCutoff.pos.denov.f_04.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.neg.shared.f_04,   paste0("outConvert/ltCutoff.neg.shared.f_04.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.pos.shared.f_04,   paste0("outConvert/ltCutoff.pos.shared.f_04.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        met <- "f_05"
        write.table( ltCutoff.neg.align.f_05,   paste0("outConvert/ltCutoff.neg.align.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.pos.align.f_05,   paste0("outConvert/ltCutoff.pos.align.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.neg.denov.f_05,   paste0("outConvert/ltCutoff.neg.denov.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table( ltCutoff.pos.denov.f_05,   paste0("outConvert/ltCutoff.pos.denov.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.neg.shared.f_05,   paste0("outConvert/ltCutoff.neg.shared.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(ltCutoff.pos.shared.f_05,   paste0("outConvert/ltCutoff.pos.shared.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        
  }
   
  # CAPER
  for(mod in 1:12){
    mtc <- 1
    for(met in c("f_01","f_02","f_04","f_05")){   # THEN RE_RUN TO SUM OVER 2 and 5!
      
      c.align.pv.lm <-  paste0( mod,".",met,".LM.p_val.FDR" )
      c.align.pv.caper <-  paste0( mod,".",met,".CAPER.p_val.FDR" )
      c.align.pv.mdmr_mash <-  paste0( mod,".",met,".MDMR_mash.","pv.FDR" )
      c.align.pv.mdmr_phylo <-  paste0( mod,".",met,".MDMR_phylo.","pv.FDR" )
      c.align.slp.lm <- paste0( mod,".",met,".LM.slope") 
      c.align.pos.lm <- paste0( mod,".",met,".LM.is_pos") 
      c.align.slp.caper <- paste0( mod,".",met,".CAPER.slope") 
      c.align.pos.caper <- paste0( mod,".",met,".CAPER.is_pos") 
      
      c.denov.pv.lm <-  paste0( mod+12,".",met,".LM.p_val.FDR" )
      c.denov.pv.caper <-  paste0( mod+12,".",met,".CAPER.p_val.FDR" )
      c.denov.pv.mdmr_mash <-  paste0( mod+12,".",met,".MDMR_mash.","pv.FDR" )
      c.denov.pv.mdmr_phylo <-  paste0( mod+12,".",met,".MDMR_phylo.","pv.FDR" )
      c.denov.slp.lm <- paste0( mod+12,".",met,".LM.slope") 
      c.denov.pos.lm <- paste0( mod+12,".",met,".LM.is_pos") 
      c.denov.slp.caper <- paste0( mod+12,".",met,".CAPER.slope") 
      c.denov.pos.caper <- paste0( mod+12,".",met,".CAPER.is_pos") 
      
      cmpr <- model.p.dos[ , c.align.pv.caper ] < cutoff
      cmpr[is.na(cmpr)] <- FALSE
      ltCutoff <- rownames(model.p.dos)[  cmpr ]
      ltCutoff.neg.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] < 0]
      ltCutoff.pos.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] > 0]
      ltCutoff.neg.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] < 0]
      ltCutoff.pos.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] > 0]
      if(met == "f_01"){
        significant_og_counter_matrix.f_01[  7 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_01[  8 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_01[  9 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  10 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_01[  11 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  12 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_01 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_01 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_01 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_01 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_01 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_01 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_02"){
        significant_og_counter_matrix.f_02 [  7 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_02 [  8 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_02 [  9 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_02 [  10 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_02 [  11 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_02 [  12 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_02 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_02 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_02 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_02 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_02 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_02 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_04"){
        significant_og_counter_matrix.f_04 [  7 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_04 [  8 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_04 [  9 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_04 [  10 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_04 [  11 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_04 [  12 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_04 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_04 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_04 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_04 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_04 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_04 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_05"){
        significant_og_counter_matrix.f_05 [  7 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_05 [  8 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_05 [  9 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_05 [  10 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_05 [  11 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_05 [  12 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_05 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_05 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_05 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_05 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_05 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_05 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      mtc <- mtc + 1
      
    }
    
    a <- ltCutoff.neg.denov.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    a <- a[a %in% ltCutoff.neg.align.f_02]
    ltCutoff.neg.denov.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    significant_og_counter_matrix.shared[ 7 , mod ] <- length( ltCutoff.neg.denov.f_shared )
    
    a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    a <- a[a %in% ltCutoff.pos.align.f_02]
    ltCutoff.pos.denov.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    significant_og_counter_matrix.shared[ 8 , mod ] <- length( ltCutoff.pos.denov.f_shared )
    
    a <- ltCutoff.neg.align.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    a <- a[a %in% ltCutoff.neg.align.f_02]
    ltCutoff.neg.align.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    significant_og_counter_matrix.shared[ 9 , mod ] <- length( ltCutoff.neg.align.f_shared )
    
    a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    a <- a[a %in% ltCutoff.pos.align.f_02]
    ltCutoff.pos.align.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    significant_og_counter_matrix.shared[ 10 , mod ] <- length( ltCutoff.pos.align.f_shared )
    
    a <- ltCutoff.neg.shared.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    a <- a[a %in% ltCutoff.neg.align.f_02]
    ltCutoff.neg.shared.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    significant_og_counter_matrix.shared[ 11 , mod ] <- length( ltCutoff.neg.shared.f_shared )
    
    a <- ltCutoff.pos.shared.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    a <- a[a %in% ltCutoff.pos.align.f_02]
    ltCutoff.pos.shared.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    significant_og_counter_matrix.shared[ 12 , mod ] <- length( ltCutoff.pos.shared.f_shared )
    
    #write.table(significant_og_counter_matrix.shared, paste0("significant_og_counter_matrix.",mod,".f_05.tsv"), sep='\t')
    
    ltCutoff.neg.align.f_shared.caper  <- ltCutoff.neg.align
    ltCutoff.pos.align.f_shared.caper  <- ltCutoff.pos.align
    ltCutoff.neg.denov.f_shared.caper  <- ltCutoff.neg.denov
    ltCutoff.pos.denov.f_shared.caper  <- ltCutoff.pos.denov
    ltCutoff.neg.shared.f_shared.caper <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
    ltCutoff.pos.shared.f_shared.caper <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
    
    write.table(ltCutoff.neg.align.f_shared.caper,   paste0("outConvert/ltCutoff.neg.align.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.align.f_shared.caper,   paste0("outConvert/ltCutoff.pos.align.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.denov.f_shared.caper,   paste0("outConvert/ltCutoff.neg.denov.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.denov.f_shared.caper,   paste0("outConvert/ltCutoff.pos.denov.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_shared.caper,   paste0("outConvert/ltCutoff.neg.shared.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_shared.caper,   paste0("outConvert/ltCutoff.pos.shared.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_01"
    write.table( ltCutoff.neg.align.f_01,   paste0("outConvert/ltCutoff.neg.align.f_01.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_01,   paste0("outConvert/ltCutoff.pos.align.f_01.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_01,   paste0("outConvert/ltCutoff.neg.denov.f_01.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_01,   paste0("outConvert/ltCutoff.pos.denov.f_01.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_01,   paste0("outConvert/ltCutoff.neg.shared.f_01.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_01,   paste0("outConvert/ltCutoff.pos.shared.f_01.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_02"
    write.table( ltCutoff.neg.align.f_02,   paste0("outConvert/ltCutoff.neg.align.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_02,   paste0("outConvert/ltCutoff.pos.align.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_02,   paste0("outConvert/ltCutoff.neg.denov.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_02,   paste0("outConvert/ltCutoff.pos.denov.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_02,   paste0("outConvert/ltCutoff.neg.shared.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_02,   paste0("outConvert/ltCutoff.pos.shared.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_04"
    write.table( ltCutoff.neg.align.f_04,   paste0("outConvert/ltCutoff.neg.align.f_04.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_04,   paste0("outConvert/ltCutoff.pos.align.f_04.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_04,   paste0("outConvert/ltCutoff.neg.denov.f_04.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_04,   paste0("outConvert/ltCutoff.pos.denov.f_04.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_04,   paste0("outConvert/ltCutoff.neg.shared.f_04.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_04,   paste0("outConvert/ltCutoff.pos.shared.f_04.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_05"
    write.table( ltCutoff.neg.align.f_05,   paste0("outConvert/ltCutoff.neg.align.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_05,   paste0("outConvert/ltCutoff.pos.align.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_05,   paste0("outConvert/ltCutoff.neg.denov.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_05,   paste0("outConvert/ltCutoff.pos.denov.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_05,   paste0("outConvert/ltCutoff.neg.shared.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_05,   paste0("outConvert/ltCutoff.pos.shared.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  }
  
  # MDMR_mash
  mod_list.mdmr_mash <- list()
  for(mod in 1:12){
    mtc <- 1
    for(met in c("f_01","f_02","f_04","f_05")){   # THEN RE_RUN TO SUM OVER 2 and 5!
      
      c.align.pv.lm <-  paste0( mod,".",met,".LM.p_val.FDR" )
      c.align.pv.caper <-  paste0( mod,".",met,".CAPER.p_val.FDR" )
      c.align.pv.mdmr_mash <-  paste0( mod,".",met,".MDMR_mash.","pv.FDR" )
      c.align.pv.mdmr_phylo <-  paste0( mod,".",met,".MDMR_phylo.","pv.FDR" )
      c.align.slp.lm <- paste0( mod,".",met,".LM.slope") 
      c.align.pos.lm <- paste0( mod,".",met,".LM.is_pos") 
      c.align.slp.caper <- paste0( mod,".",met,".CAPER.slope") 
      c.align.pos.caper <- paste0( mod,".",met,".CAPER.is_pos") 
      
      c.denov.pv.lm <-  paste0( mod+12,".",met,".LM.p_val.FDR" )
      c.denov.pv.caper <-  paste0( mod+12,".",met,".CAPER.p_val.FDR" )
      c.denov.pv.mdmr_mash <-  paste0( mod+12,".",met,".MDMR_mash.","pv.FDR" )
      c.denov.pv.mdmr_phylo <-  paste0( mod+12,".",met,".MDMR_phylo.","pv.FDR" )
      c.denov.slp.lm <- paste0( mod+12,".",met,".LM.slope") 
      c.denov.pos.lm <- paste0( mod+12,".",met,".LM.is_pos") 
      c.denov.slp.caper <- paste0( mod+12,".",met,".CAPER.slope") 
      c.denov.pos.caper <- paste0( mod+12,".",met,".CAPER.is_pos") 
      
      cmpr <- model.p.dos[ , c.align.pv.mdmr_mash ] < cutoff
      cmpr[is.na(cmpr)] <- FALSE
      ltCutoff <- rownames(model.p.dos)[  cmpr ]
      ltCutoff.neg.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] < 0]
      ltCutoff.pos.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] > 0]
      ltCutoff.neg.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] < 0]
      ltCutoff.pos.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] > 0]
      if(met == "f_01"){
        significant_og_counter_matrix.f_01[  13 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_01[  14 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_01[  15 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  16 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_01[  17 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  18 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_01 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_01 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_01 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_01 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_01 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_01 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_02"){
        significant_og_counter_matrix.f_02[  13 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_02[  14 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_02[  15 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_02[  16 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_02[  17 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_02[  18 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_02 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_02 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_02 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_02 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_02 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_02 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_04"){
        significant_og_counter_matrix.f_04[  13 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_04[  14 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_04[  15 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_04[  16 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_04[  17 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_04[  18 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_04 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_04 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_04 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_04 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_04 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_04 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_05"){
        significant_og_counter_matrix.f_05[  13 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_05[  14 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_05[  15 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_05[  16 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_05[  17 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_05[  18 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_05 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_05 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_05 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_05 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_05 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_05 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      mtc <- mtc + 1
      
    }
    
    a <- ltCutoff.neg.denov.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    a <- a[a %in% ltCutoff.neg.align.f_02]
    ltCutoff.neg.denov.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    significant_og_counter_matrix.shared[ 13 , mod ] <- length( ltCutoff.neg.denov.f_shared )
    
    a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    a <- a[a %in% ltCutoff.pos.align.f_02]
    ltCutoff.pos.denov.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    significant_og_counter_matrix.shared[ 14 , mod ] <- length( ltCutoff.pos.denov.f_shared )
    
    a <- ltCutoff.neg.align.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    a <- a[a %in% ltCutoff.neg.align.f_02]
    ltCutoff.neg.align.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    significant_og_counter_matrix.shared[ 15 , mod ] <- length( ltCutoff.neg.align.f_shared )
    
    a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    a <- a[a %in% ltCutoff.pos.align.f_02]
    ltCutoff.pos.align.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    significant_og_counter_matrix.shared[ 16 , mod ] <- length( ltCutoff.pos.align.f_shared )
    
    a <- ltCutoff.neg.shared.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    a <- a[a %in% ltCutoff.neg.align.f_02]
    ltCutoff.neg.shared.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    significant_og_counter_matrix.shared[ 17 , mod ] <- length( ltCutoff.neg.shared.f_shared )
    
    a <- ltCutoff.pos.shared.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    a <- a[a %in% ltCutoff.pos.align.f_02]
    ltCutoff.pos.shared.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    significant_og_counter_matrix.shared[ 18 , mod ] <- length( ltCutoff.pos.shared.f_shared )
    
    #write.table(significant_og_counter_matrix.shared, paste0("significant_og_counter_matrix.",mod,".f_05.tsv"), sep='\t')
    
    ltCutoff.neg.align.f_shared.mdmr_mash  <- ltCutoff.neg.align
    ltCutoff.pos.align.f_shared.mdmr_mash  <- ltCutoff.pos.align
    ltCutoff.neg.denov.f_shared.mdmr_mash  <- ltCutoff.neg.denov
    ltCutoff.pos.denov.f_shared.mdmr_mash  <- ltCutoff.pos.denov
    ltCutoff.neg.shared.f_shared.mdmr_mash <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
    ltCutoff.pos.shared.f_shared.mdmr_mash <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
    
    write.table(ltCutoff.neg.align.f_shared.mdmr_mash,   paste0("outConvert/ltCutoff.neg.align.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.align.f_shared.mdmr_mash,   paste0("outConvert/ltCutoff.pos.align.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.denov.f_shared.mdmr_mash,   paste0("outConvert/ltCutoff.neg.denov.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.denov.f_shared.mdmr_mash,   paste0("outConvert/ltCutoff.pos.denov.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_shared.mdmr_mash,   paste0("outConvert/ltCutoff.neg.shared.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_shared.mdmr_mash,   paste0("outConvert/ltCutoff.pos.shared.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_01"
    write.table( ltCutoff.neg.align.f_01,   paste0("outConvert/ltCutoff.neg.align.f_01.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_01,   paste0("outConvert/ltCutoff.pos.align.f_01.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_01,   paste0("outConvert/ltCutoff.neg.denov.f_01.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_01,   paste0("outConvert/ltCutoff.pos.denov.f_01.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_01,   paste0("outConvert/ltCutoff.neg.shared.f_01.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_01,   paste0("outConvert/ltCutoff.pos.shared.f_01.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_02"
    write.table( ltCutoff.neg.align.f_02,   paste0("outConvert/ltCutoff.neg.align.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_02,   paste0("outConvert/ltCutoff.pos.align.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_02,   paste0("outConvert/ltCutoff.neg.denov.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_02,   paste0("outConvert/ltCutoff.pos.denov.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_02,   paste0("outConvert/ltCutoff.neg.shared.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_02,   paste0("outConvert/ltCutoff.pos.shared.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_04"
    write.table( ltCutoff.neg.align.f_04,   paste0("outConvert/ltCutoff.neg.align.f_04.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_04,   paste0("outConvert/ltCutoff.pos.align.f_04.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_04,   paste0("outConvert/ltCutoff.neg.denov.f_04.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_04,   paste0("outConvert/ltCutoff.pos.denov.f_04.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_04,   paste0("outConvert/ltCutoff.neg.shared.f_04.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_04,   paste0("outConvert/ltCutoff.pos.shared.f_04.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_05"
    write.table( ltCutoff.neg.align.f_05,   paste0("outConvert/ltCutoff.neg.align.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_05,   paste0("outConvert/ltCutoff.pos.align.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_05,   paste0("outConvert/ltCutoff.neg.denov.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_05,   paste0("outConvert/ltCutoff.pos.denov.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_05,   paste0("outConvert/ltCutoff.neg.shared.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_05,   paste0("outConvert/ltCutoff.pos.shared.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    
  }
  
  # MDMR_phylo
  for(mod in 1:12){
    mtc <- 1
    for(met in c("f_01","f_02","f_04","f_05")){   # THEN RE_RUN TO SUM OVER 2 and 5!
      
      c.align.pv.lm <-  paste0( mod,".",met,".LM.p_val.FDR" )
      c.align.pv.caper <-  paste0( mod,".",met,".CAPER.p_val.FDR" )
      c.align.pv.mdmr_mash <-  paste0( mod,".",met,".MDMR_mash.","pv.FDR" )
      c.align.pv.mdmr_phylo <-  paste0( mod,".",met,".MDMR_phylo.","pv.FDR" )
      c.align.slp.lm <- paste0( mod,".",met,".LM.slope") 
      c.align.pos.lm <- paste0( mod,".",met,".LM.is_pos") 
      c.align.slp.caper <- paste0( mod,".",met,".CAPER.slope") 
      c.align.pos.caper <- paste0( mod,".",met,".CAPER.is_pos") 
      
      c.denov.pv.lm <-  paste0( mod+12,".",met,".LM.p_val.FDR" )
      c.denov.pv.caper <-  paste0( mod+12,".",met,".CAPER.p_val.FDR" )
      c.denov.pv.mdmr_mash <-  paste0( mod+12,".",met,".MDMR_mash.","pv.FDR" )
      c.denov.pv.mdmr_phylo <-  paste0( mod+12,".",met,".MDMR_phylo.","pv.FDR" )
      c.denov.slp.lm <- paste0( mod+12,".",met,".LM.slope") 
      c.denov.pos.lm <- paste0( mod+12,".",met,".LM.is_pos") 
      c.denov.slp.caper <- paste0( mod+12,".",met,".CAPER.slope") 
      c.denov.pos.caper <- paste0( mod+12,".",met,".CAPER.is_pos") 
      
      cmpr <- model.p.dos[ , c.align.pv.mdmr_phylo ] < cutoff
      cmpr[is.na(cmpr)] <- FALSE
      ltCutoff <- rownames(model.p.dos)[  cmpr ]
      ltCutoff.neg.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] < 0]
      ltCutoff.pos.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] > 0]
      ltCutoff.neg.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] < 0]
      ltCutoff.pos.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] > 0]
      if(met == "f_01"){
        significant_og_counter_matrix.f_01[  19 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_01[  20 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_01[  21 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  22 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_01[  23 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  24 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_01 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_01 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_01 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_01 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_01 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_01 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_02"){
        significant_og_counter_matrix.f_02[  19 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_02[  20 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_02[  21 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_02[  22 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_02[  23 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_02[  24 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_02 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_02 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_02 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_02 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_02 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_02 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_04"){
        significant_og_counter_matrix.f_04[  19 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_04[  20 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_04[  21 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_04[  22 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_04[  23 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_04[  24 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_04 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_04 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_04 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_04 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_04 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_04 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_05"){
        significant_og_counter_matrix.f_05[  19 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_05[  20 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_05[  21 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_05[  22 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_05[  23 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_05[  24 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_05 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_05 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_05 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_05 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_05 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_05 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      mtc <- mtc + 1
      
    }
    
    a <- ltCutoff.neg.denov.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    a <- a[a %in% ltCutoff.neg.align.f_02]
    ltCutoff.neg.denov.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    significant_og_counter_matrix.shared[ 19 , mod ] <- length( ltCutoff.neg.denov.f_shared )
    
    a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    a <- a[a %in% ltCutoff.pos.align.f_02]
    ltCutoff.pos.denov.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    significant_og_counter_matrix.shared[ 20 , mod ] <- length( ltCutoff.pos.denov.f_shared )
    
    a <- ltCutoff.neg.align.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    a <- a[a %in% ltCutoff.neg.align.f_02]
    ltCutoff.neg.align.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    significant_og_counter_matrix.shared[ 21 , mod ] <- length( ltCutoff.neg.align.f_shared )
    
    a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    a <- a[a %in% ltCutoff.pos.align.f_02]
    ltCutoff.pos.align.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    significant_og_counter_matrix.shared[ 22 , mod ] <- length( ltCutoff.pos.align.f_shared )
    
    a <- ltCutoff.neg.shared.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    a <- a[a %in% ltCutoff.neg.align.f_02]
    ltCutoff.neg.shared.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    significant_og_counter_matrix.shared[ 23 , mod ] <- length( ltCutoff.neg.shared.f_shared )
    
    a <- ltCutoff.pos.shared.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    a <- a[a %in% ltCutoff.pos.align.f_02]
    ltCutoff.pos.shared.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    significant_og_counter_matrix.shared[ 24 , mod ] <- length( ltCutoff.pos.shared.f_shared )
    
    #write.table(significant_og_counter_matrix.shared, paste0("significant_og_counter_matrix.",mod,".f_05.tsv"), sep='\t')
    
    ltCutoff.neg.align.f_shared.mdmr_phylo  <- ltCutoff.neg.align
    ltCutoff.pos.align.f_shared.mdmr_phylo  <- ltCutoff.pos.align
    ltCutoff.neg.denov.f_shared.mdmr_phylo  <- ltCutoff.neg.denov
    ltCutoff.pos.denov.f_shared.mdmr_phylo  <- ltCutoff.pos.denov
    ltCutoff.neg.shared.f_shared.mdmr_phylo <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
    ltCutoff.pos.shared.f_shared.mdmr_phylo <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
    
    write.table(ltCutoff.neg.align.f_shared.mdmr_phylo,   paste0("outConvert/ltCutoff.neg.align.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.align.f_shared.mdmr_phylo,   paste0("outConvert/ltCutoff.pos.align.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.denov.f_shared.mdmr_phylo,   paste0("outConvert/ltCutoff.neg.denov.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.denov.f_shared.mdmr_phylo,   paste0("outConvert/ltCutoff.pos.denov.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_shared.mdmr_phylo,   paste0("outConvert/ltCutoff.neg.shared.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_shared.mdmr_phylo,   paste0("outConvert/ltCutoff.pos.shared.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_01"
    write.table( ltCutoff.neg.align.f_01,   paste0("outConvert/ltCutoff.neg.align.f_01.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_01,   paste0("outConvert/ltCutoff.pos.align.f_01.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_01,   paste0("outConvert/ltCutoff.neg.denov.f_01.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_01,   paste0("outConvert/ltCutoff.pos.denov.f_01.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_01,   paste0("outConvert/ltCutoff.neg.shared.f_01.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_01,   paste0("outConvert/ltCutoff.pos.shared.f_01.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_02"
    write.table( ltCutoff.neg.align.f_02,   paste0("outConvert/ltCutoff.neg.align.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_02,   paste0("outConvert/ltCutoff.pos.align.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_02,   paste0("outConvert/ltCutoff.neg.denov.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_02,   paste0("outConvert/ltCutoff.pos.denov.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_02,   paste0("outConvert/ltCutoff.neg.shared.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_02,   paste0("outConvert/ltCutoff.pos.shared.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_04"
    write.table( ltCutoff.neg.align.f_04,   paste0("outConvert/ltCutoff.neg.align.f_04.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_04,   paste0("outConvert/ltCutoff.pos.align.f_04.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_04,   paste0("outConvert/ltCutoff.neg.denov.f_04.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_04,   paste0("outConvert/ltCutoff.pos.denov.f_04.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_04,   paste0("outConvert/ltCutoff.neg.shared.f_04.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_04,   paste0("outConvert/ltCutoff.pos.shared.f_04.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_05"
    write.table( ltCutoff.neg.align.f_05,   paste0("outConvert/ltCutoff.neg.align.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_05,   paste0("outConvert/ltCutoff.pos.align.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_05,   paste0("outConvert/ltCutoff.neg.denov.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_05,   paste0("outConvert/ltCutoff.pos.denov.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_05,   paste0("outConvert/ltCutoff.neg.shared.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_05,   paste0("outConvert/ltCutoff.pos.shared.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    
    
  }
  
  write.table(significant_og_counter_matrix.f_01,   paste0("significant_og_counter_matrix.f_01.",cutoff,".tsv"), sep='\t')
  write.table(significant_og_counter_matrix.f_02,   paste0("significant_og_counter_matrix.f_02.",cutoff,".tsv"), sep='\t')
  write.table(significant_og_counter_matrix.f_04,   paste0("significant_og_counter_matrix.f_04.",cutoff,".tsv"), sep='\t')
  write.table(significant_og_counter_matrix.f_05,   paste0("significant_og_counter_matrix.f_05.",cutoff,".tsv"), sep='\t')
  write.table(significant_og_counter_matrix.shared, paste0("significant_og_counter_matrix.shared.",cutoff,".tsv"), sep='\t')
}









######################################################################################################################################################################################################
# ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ##################
######################################################################################################################################################################################################
# ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ##################
######################################################################################################################################################################################################
# ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ##################
######################################################################################################################################################################################################
# ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ##################
######################################################################################################################################################################################################
# ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ################### ONLY 1 AND 2 ##################
######################################################################################################################################################################################################









############ 1.1 : Identify Significant OGs ############ 

v9_v10_OGs_map     <- read.table("v9_v10_OGs_map.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
odb10v1_OG_xrefs     <- read.table("odb10v1_OG_xrefs.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
odb10v1_OG2genes     <- read.table("odb10v1_OG2genes.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)

# pull colnames
# pp <-  c( colnames(model.p.dos)[ grepl('p_val$', colnames(model.p.dos)) ],  colnames(model.p.dos)[ grepl('pv$', colnames(model.p.dos)) ] )

significant_og_counter_matrix.f_01 <- matrix(nrow=24,ncol=12)
significant_og_counter_matrix.f_02 <- matrix(nrow=24,ncol=12)
significant_og_counter_matrix.f_04 <- matrix(nrow=24,ncol=12)
significant_og_counter_matrix.f_05 <- matrix(nrow=24,ncol=12)
significant_og_counter_matrix.shared <- matrix(nrow=424,ncol=12)

#cutoff <- 0.05
for(cutoff in c(0.05,0.005)){
  # LM
  for(mod in 1:12){
    mtc <- 1
    #for(met in c("f_01","f_02","f_04","f_05")){   # THEN RE_RUN TO SUM OVER 2 and 5!
    for(met in c("f_02","f_05")){   # THEN RE_RUN TO SUM OVER 2 and 5!
      
      c.align.pv.lm <-  paste0( mod,".",met,".LM.p_val.FDR" )
      c.align.pv.caper <-  paste0( mod,".",met,".CAPER.p_val.FDR" )
      c.align.pv.mdmr_mash <-  paste0( mod,".",met,".MDMR_mash.","pv.FDR" )
      c.align.pv.mdmr_phylo <-  paste0( mod,".",met,".MDMR_phylo.","pv.FDR" )
      c.align.slp.lm <- paste0( mod,".",met,".LM.slope") 
      c.align.pos.lm <- paste0( mod,".",met,".LM.is_pos") 
      c.align.slp.caper <- paste0( mod,".",met,".CAPER.slope") 
      c.align.pos.caper <- paste0( mod,".",met,".CAPER.is_pos") 
      
      c.denov.pv.lm <-  paste0( mod+12,".",met,".LM.p_val.FDR" )
      c.denov.pv.caper <-  paste0( mod+12,".",met,".CAPER.p_val.FDR" )
      c.denov.pv.mdmr_mash <-  paste0( mod+12,".",met,".MDMR_mash.","pv.FDR" )
      c.denov.pv.mdmr_phylo <-  paste0( mod+12,".",met,".MDMR_phylo.","pv.FDR" )
      c.denov.slp.lm <- paste0( mod+12,".",met,".LM.slope") 
      c.denov.pos.lm <- paste0( mod+12,".",met,".LM.is_pos") 
      c.denov.slp.caper <- paste0( mod+12,".",met,".CAPER.slope") 
      c.denov.pos.caper <- paste0( mod+12,".",met,".CAPER.is_pos") 
      
      cmpr <- model.p.dos[ , c.align.pv.lm ] < cutoff
      cmpr[is.na(cmpr)] <- FALSE
      ltCutoff <- rownames(model.p.dos)[  cmpr ]
      ltCutoff.neg.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.lm] < 0]
      ltCutoff.pos.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.lm] > 0]
      ltCutoff.neg.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.lm] < 0]
      ltCutoff.pos.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.lm] > 0]

      if(met == "f_02"){
        significant_og_counter_matrix.f_02 [  1 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_02 [  2 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_02 [  3 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_02 [  4 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_02 [  5 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_02 [  6 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_02 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_02 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_02 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_02 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_02 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_02 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_05"){
        significant_og_counter_matrix.f_05 [  1 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_05 [  2 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_05 [  3 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_05 [  4 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_05 [  5 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_05 [  6 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_05 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_05 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_05 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_05 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_05 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_05 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      mtc <- mtc + 1
      
    }
    
    significant_og_counter_matrix.shared[ 1 , mod ] <- sum( ltCutoff.neg.denov.f_05 %in%  ltCutoff.neg.denov.f_02 )
    significant_og_counter_matrix.shared[ 2 , mod ] <- sum( ltCutoff.pos.denov.f_05 %in%  ltCutoff.pos.denov.f_02 )
    significant_og_counter_matrix.shared[ 3 , mod ] <- sum( ltCutoff.neg.align.f_05 %in%  ltCutoff.neg.align.f_02 )
    significant_og_counter_matrix.shared[ 4 , mod ] <- sum( ltCutoff.pos.denov.f_05 %in%  ltCutoff.pos.denov.f_02 )
    significant_og_counter_matrix.shared[ 5 , mod ] <- sum( ltCutoff.neg.shared.f_05 %in% ltCutoff.neg.shared.f_02 )
    significant_og_counter_matrix.shared[ 6 , mod ] <- sum( ltCutoff.pos.shared.f_05 %in% ltCutoff.pos.shared.f_02 )
    
    ltCutoff.neg.align.f_shared.lm <- ltCutoff.neg.align
    ltCutoff.pos.align.f_shared.lm <- ltCutoff.pos.align
    ltCutoff.neg.denov.f_shared.lm <- ltCutoff.neg.denov
    ltCutoff.pos.denov.f_shared.lm <- ltCutoff.pos.denov
    ltCutoff.neg.shared.f_shared.lm <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
    ltCutoff.pos.shared.f_shared.lm <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
    
    write.table(ltCutoff.neg.align.f_shared.lm,   paste0("outConvert2/ltCutoff.neg.align.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.align.f_shared.lm,   paste0("outConvert2/ltCutoff.pos.align.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.denov.f_shared.lm,   paste0("outConvert2/ltCutoff.neg.denov.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.denov.f_shared.lm,   paste0("outConvert2/ltCutoff.pos.denov.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_shared.lm,   paste0("outConvert2/ltCutoff.neg.shared.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_shared.lm,   paste0("outConvert2/ltCutoff.pos.shared.f_shared.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_02"
    write.table( ltCutoff.neg.align.f_02,   paste0("outConvert2/ltCutoff.neg.align.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_02,   paste0("outConvert2/ltCutoff.pos.align.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_02,   paste0("outConvert2/ltCutoff.neg.denov.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_02,   paste0("outConvert2/ltCutoff.pos.denov.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_02,   paste0("outConvert2/ltCutoff.neg.shared.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_02,   paste0("outConvert2/ltCutoff.pos.shared.f_02.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_05"
    write.table( ltCutoff.neg.align.f_05,   paste0("outConvert2/ltCutoff.neg.align.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_05,   paste0("outConvert2/ltCutoff.pos.align.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_05,   paste0("outConvert2/ltCutoff.neg.denov.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_05,   paste0("outConvert2/ltCutoff.pos.denov.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_05,   paste0("outConvert2/ltCutoff.neg.shared.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_05,   paste0("outConvert2/ltCutoff.pos.shared.f_05.lm.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  }
  
  # CAPER
  for(mod in 1:12){
    mtc <- 1
    for(met in c("f_01","f_02","f_04","f_05")){   # THEN RE_RUN TO SUM OVER 2 and 5!
      
      c.align.pv.lm <-  paste0( mod,".",met,".LM.p_val.FDR" )
      c.align.pv.caper <-  paste0( mod,".",met,".CAPER.p_val.FDR" )
      c.align.pv.mdmr_mash <-  paste0( mod,".",met,".MDMR_mash.","pv.FDR" )
      c.align.pv.mdmr_phylo <-  paste0( mod,".",met,".MDMR_phylo.","pv.FDR" )
      c.align.slp.lm <- paste0( mod,".",met,".LM.slope") 
      c.align.pos.lm <- paste0( mod,".",met,".LM.is_pos") 
      c.align.slp.caper <- paste0( mod,".",met,".CAPER.slope") 
      c.align.pos.caper <- paste0( mod,".",met,".CAPER.is_pos") 
      
      c.denov.pv.lm <-  paste0( mod+12,".",met,".LM.p_val.FDR" )
      c.denov.pv.caper <-  paste0( mod+12,".",met,".CAPER.p_val.FDR" )
      c.denov.pv.mdmr_mash <-  paste0( mod+12,".",met,".MDMR_mash.","pv.FDR" )
      c.denov.pv.mdmr_phylo <-  paste0( mod+12,".",met,".MDMR_phylo.","pv.FDR" )
      c.denov.slp.lm <- paste0( mod+12,".",met,".LM.slope") 
      c.denov.pos.lm <- paste0( mod+12,".",met,".LM.is_pos") 
      c.denov.slp.caper <- paste0( mod+12,".",met,".CAPER.slope") 
      c.denov.pos.caper <- paste0( mod+12,".",met,".CAPER.is_pos") 
      
      cmpr <- model.p.dos[ , c.align.pv.caper ] < cutoff
      cmpr[is.na(cmpr)] <- FALSE
      ltCutoff <- rownames(model.p.dos)[  cmpr ]
      ltCutoff.neg.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] < 0]
      ltCutoff.pos.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] > 0]
      ltCutoff.neg.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] < 0]
      ltCutoff.pos.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] > 0]
      if(met == "f_01"){
        significant_og_counter_matrix.f_01[  7 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_01[  8 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_01[  9 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  10 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_01[  11 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  12 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_01 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_01 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_01 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_01 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_01 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_01 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_02"){
        significant_og_counter_matrix.f_02 [  7 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_02 [  8 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_02 [  9 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_02 [  10 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_02 [  11 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_02 [  12 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_02 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_02 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_02 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_02 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_02 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_02 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_04"){
        significant_og_counter_matrix.f_04 [  7 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_04 [  8 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_04 [  9 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_04 [  10 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_04 [  11 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_04 [  12 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_04 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_04 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_04 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_04 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_04 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_04 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_05"){
        significant_og_counter_matrix.f_05 [  7 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_05 [  8 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_05 [  9 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_05 [  10 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_05 [  11 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_05 [  12 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_05 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_05 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_05 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_05 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_05 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_05 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      mtc <- mtc + 1
      
    }
    
    # a <- ltCutoff.neg.denov.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    # a <- a[a %in% ltCutoff.neg.align.f_02]
    # ltCutoff.neg.denov.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    # significant_og_counter_matrix.shared[ 1 , mod ] <- length( ltCutoff.neg.denov.f_shared )
    # 
    # a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    # a <- a[a %in% ltCutoff.pos.align.f_02]
    # ltCutoff.pos.denov.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    # significant_og_counter_matrix.shared[ 2 , mod ] <- length( ltCutoff.pos.denov.f_shared )
    # 
    # a <- ltCutoff.neg.align.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    # a <- a[a %in% ltCutoff.neg.align.f_02]
    # ltCutoff.neg.align.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    # significant_og_counter_matrix.shared[ 3 , mod ] <- length( ltCutoff.neg.align.f_shared )
    # 
    # a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    # a <- a[a %in% ltCutoff.pos.align.f_02]
    # ltCutoff.pos.align.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    # significant_og_counter_matrix.shared[ 4 , mod ] <- length( ltCutoff.pos.align.f_shared )
    # 
    # a <- ltCutoff.neg.shared.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    # a <- a[a %in% ltCutoff.neg.align.f_02]
    # ltCutoff.neg.shared.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    # significant_og_counter_matrix.shared[ 5 , mod ] <- length( ltCutoff.neg.shared.f_shared )
    # 
    # a <- ltCutoff.pos.shared.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    # a <- a[a %in% ltCutoff.pos.align.f_02]
    # ltCutoff.pos.shared.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    # significant_og_counter_matrix.shared[ 6 , mod ] <- length( ltCutoff.pos.shared.f_shared )
    
    
    significant_og_counter_matrix.shared[ 7 , mod ] <- sum( ltCutoff.neg.denov.f_05 %in%  ltCutoff.neg.denov.f_02 )
    significant_og_counter_matrix.shared[ 8 , mod ] <- sum( ltCutoff.pos.denov.f_05 %in%  ltCutoff.pos.denov.f_02 )
    significant_og_counter_matrix.shared[ 9 , mod ] <- sum( ltCutoff.neg.align.f_05 %in%  ltCutoff.neg.align.f_02 )
    significant_og_counter_matrix.shared[ 10 , mod ] <- sum( ltCutoff.pos.denov.f_05 %in%  ltCutoff.pos.denov.f_02 )
    significant_og_counter_matrix.shared[ 11 , mod ] <- sum( ltCutoff.neg.shared.f_05 %in% ltCutoff.neg.shared.f_02 )
    significant_og_counter_matrix.shared[ 12 , mod ] <- sum( ltCutoff.pos.shared.f_05 %in% ltCutoff.pos.shared.f_02 )
    
    #write.table(significant_og_counter_matrix.shared, paste0("significant_og_counter_matrix.",mod,".f_05.tsv"), sep='\t')
    
    ltCutoff.neg.align.f_shared.caper <- ltCutoff.neg.align
    ltCutoff.pos.align.f_shared.caper <- ltCutoff.pos.align
    ltCutoff.neg.denov.f_shared.caper <- ltCutoff.neg.denov
    ltCutoff.pos.denov.f_shared.caper <- ltCutoff.pos.denov
    ltCutoff.neg.shared.f_shared.caper <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
    ltCutoff.pos.shared.f_shared.caper <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
    
    write.table(ltCutoff.neg.align.f_shared.caper,   paste0("outConvert2/ltCutoff.neg.align.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.align.f_shared.caper,   paste0("outConvert2/ltCutoff.pos.align.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.denov.f_shared.caper,   paste0("outConvert2/ltCutoff.neg.denov.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.denov.f_shared.caper,   paste0("outConvert2/ltCutoff.pos.denov.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_shared.caper,   paste0("outConvert2/ltCutoff.neg.shared.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_shared.caper,   paste0("outConvert2/ltCutoff.pos.shared.f_shared.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_02"
    write.table( ltCutoff.neg.align.f_02,   paste0("outConvert2/ltCutoff.neg.align.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_02,   paste0("outConvert2/ltCutoff.pos.align.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_02,   paste0("outConvert2/ltCutoff.neg.denov.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_02,   paste0("outConvert2/ltCutoff.pos.denov.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_02,   paste0("outConvert2/ltCutoff.neg.shared.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_02,   paste0("outConvert2/ltCutoff.pos.shared.f_02.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_05"
    write.table( ltCutoff.neg.align.f_05,   paste0("outConvert2/ltCutoff.neg.align.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_05,   paste0("outConvert2/ltCutoff.pos.align.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_05,   paste0("outConvert2/ltCutoff.neg.denov.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_05,   paste0("outConvert2/ltCutoff.pos.denov.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_05,   paste0("outConvert2/ltCutoff.neg.shared.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_05,   paste0("outConvert2/ltCutoff.pos.shared.f_05.caper.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  }
  
  # MDMR_mash
  mod_list.mdmr_mash <- list()
  for(mod in 1:12){
    mtc <- 1
    for(met in c("f_01","f_02","f_04","f_05")){   # THEN RE_RUN TO SUM OVER 2 and 5!
      
      c.align.pv.lm <-  paste0( mod,".",met,".LM.p_val.FDR" )
      c.align.pv.caper <-  paste0( mod,".",met,".CAPER.p_val.FDR" )
      c.align.pv.mdmr_mash <-  paste0( mod,".",met,".MDMR_mash.","pv.FDR" )
      c.align.pv.mdmr_phylo <-  paste0( mod,".",met,".MDMR_phylo.","pv.FDR" )
      c.align.slp.lm <- paste0( mod,".",met,".LM.slope") 
      c.align.pos.lm <- paste0( mod,".",met,".LM.is_pos") 
      c.align.slp.caper <- paste0( mod,".",met,".CAPER.slope") 
      c.align.pos.caper <- paste0( mod,".",met,".CAPER.is_pos") 
      
      c.denov.pv.lm <-  paste0( mod+12,".",met,".LM.p_val.FDR" )
      c.denov.pv.caper <-  paste0( mod+12,".",met,".CAPER.p_val.FDR" )
      c.denov.pv.mdmr_mash <-  paste0( mod+12,".",met,".MDMR_mash.","pv.FDR" )
      c.denov.pv.mdmr_phylo <-  paste0( mod+12,".",met,".MDMR_phylo.","pv.FDR" )
      c.denov.slp.lm <- paste0( mod+12,".",met,".LM.slope") 
      c.denov.pos.lm <- paste0( mod+12,".",met,".LM.is_pos") 
      c.denov.slp.caper <- paste0( mod+12,".",met,".CAPER.slope") 
      c.denov.pos.caper <- paste0( mod+12,".",met,".CAPER.is_pos") 
      
      cmpr <- model.p.dos[ , c.align.pv.mdmr_mash ] < cutoff
      cmpr[is.na(cmpr)] <- FALSE
      ltCutoff <- rownames(model.p.dos)[  cmpr ]
      ltCutoff.neg.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] < 0]
      ltCutoff.pos.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] > 0]
      ltCutoff.neg.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] < 0]
      ltCutoff.pos.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] > 0]
      if(met == "f_01"){
        significant_og_counter_matrix.f_01[  13 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_01[  14 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_01[  15 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  16 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_01[  17 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  18 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_01 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_01 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_01 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_01 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_01 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_01 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_02"){
        significant_og_counter_matrix.f_02[  13 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_02[  14 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_02[  15 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_02[  16 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_02[  17 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_02[  18 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_02 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_02 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_02 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_02 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_02 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_02 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_04"){
        significant_og_counter_matrix.f_04[  13 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_04[  14 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_04[  15 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_04[  16 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_04[  17 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_04[  18 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_04 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_04 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_04 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_04 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_04 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_04 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_05"){
        significant_og_counter_matrix.f_05[  13 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_05[  14 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_05[  15 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_05[  16 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_05[  17 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_05[  18 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_05 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_05 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_05 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_05 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_05 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_05 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      mtc <- mtc + 1
      
    }
    
    # a <- ltCutoff.neg.denov.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    # a <- a[a %in% ltCutoff.neg.align.f_02]
    # ltCutoff.neg.denov.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    # significant_og_counter_matrix.shared[ 1 , mod ] <- length( ltCutoff.neg.denov.f_shared )
    # 
    # a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    # a <- a[a %in% ltCutoff.pos.align.f_02]
    # ltCutoff.pos.denov.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    # significant_og_counter_matrix.shared[ 2 , mod ] <- length( ltCutoff.pos.denov.f_shared )
    # 
    # a <- ltCutoff.neg.align.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    # a <- a[a %in% ltCutoff.neg.align.f_02]
    # ltCutoff.neg.align.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    # significant_og_counter_matrix.shared[ 3 , mod ] <- length( ltCutoff.neg.align.f_shared )
    # 
    # a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    # a <- a[a %in% ltCutoff.pos.align.f_02]
    # ltCutoff.pos.align.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    # significant_og_counter_matrix.shared[ 4 , mod ] <- length( ltCutoff.pos.align.f_shared )
    # 
    # a <- ltCutoff.neg.shared.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    # a <- a[a %in% ltCutoff.neg.align.f_02]
    # ltCutoff.neg.shared.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    # significant_og_counter_matrix.shared[ 5 , mod ] <- length( ltCutoff.neg.shared.f_shared )
    # 
    # a <- ltCutoff.pos.shared.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    # a <- a[a %in% ltCutoff.pos.align.f_02]
    # ltCutoff.pos.shared.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    # significant_og_counter_matrix.shared[ 6 , mod ] <- length( ltCutoff.pos.shared.f_shared )
    
    
    significant_og_counter_matrix.shared[ 13 , mod ] <- sum( ltCutoff.neg.denov.f_05 %in%  ltCutoff.neg.align.f_02 )
    significant_og_counter_matrix.shared[ 14 , mod ] <- sum( ltCutoff.pos.denov.f_05 %in%  ltCutoff.pos.denov.f_02 )
    significant_og_counter_matrix.shared[ 15 , mod ] <- sum( ltCutoff.neg.align.f_05 %in%  ltCutoff.neg.align.f_02 )
    significant_og_counter_matrix.shared[ 16 , mod ] <- sum( ltCutoff.pos.denov.f_05 %in%  ltCutoff.pos.denov.f_02 )
    significant_og_counter_matrix.shared[ 17 , mod ] <- sum( ltCutoff.neg.shared.f_05 %in% ltCutoff.neg.shared.f_02 )
    significant_og_counter_matrix.shared[ 18 , mod ] <- sum( ltCutoff.pos.shared.f_05 %in% ltCutoff.pos.shared.f_02 )
    
    #write.table(significant_og_counter_matrix.shared, paste0("significant_og_counter_matrix.",mod,".f_05.tsv"), sep='\t')
    
    ltCutoff.neg.align.f_shared.mdmr_mash <- ltCutoff.neg.align
    ltCutoff.pos.align.f_shared.mdmr_mash <- ltCutoff.pos.align
    ltCutoff.neg.denov.f_shared.mdmr_mash <- ltCutoff.neg.denov
    ltCutoff.pos.denov.f_shared.mdmr_mash <- ltCutoff.pos.denov
    ltCutoff.neg.shared.f_shared.mdmr_mash <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
    ltCutoff.pos.shared.f_shared.mdmr_mash <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
    
    write.table(ltCutoff.neg.align.f_shared.mdmr_mash,   paste0("outConvert2/ltCutoff.neg.align.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.align.f_shared.mdmr_mash,   paste0("outConvert2/ltCutoff.pos.align.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.denov.f_shared.mdmr_mash,   paste0("outConvert2/ltCutoff.neg.denov.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.denov.f_shared.mdmr_mash,   paste0("outConvert2/ltCutoff.pos.denov.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_shared.mdmr_mash,   paste0("outConvert2/ltCutoff.neg.shared.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_shared.mdmr_mash,   paste0("outConvert2/ltCutoff.pos.shared.f_shared.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_02"
    write.table( ltCutoff.neg.align.f_02,   paste0("outConvert2/ltCutoff.neg.align.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_02,   paste0("outConvert2/ltCutoff.pos.align.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_02,   paste0("outConvert2/ltCutoff.neg.denov.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_02,   paste0("outConvert2/ltCutoff.pos.denov.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_02,   paste0("outConvert2/ltCutoff.neg.shared.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_02,   paste0("outConvert2/ltCutoff.pos.shared.f_02.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_05"
    write.table( ltCutoff.neg.align.f_05,   paste0("outConvert2/ltCutoff.neg.align.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_05,   paste0("outConvert2/ltCutoff.pos.align.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_05,   paste0("outConvert2/ltCutoff.neg.denov.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_05,   paste0("outConvert2/ltCutoff.pos.denov.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_05,   paste0("outConvert2/ltCutoff.neg.shared.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_05,   paste0("outConvert2/ltCutoff.pos.shared.f_05.mdmr_mash.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    
  }
  
  # MDMR_phylo
  for(mod in 1:12){
    mtc <- 1
    for(met in c("f_01","f_02","f_04","f_05")){   # THEN RE_RUN TO SUM OVER 2 and 5!
      
      c.align.pv.lm <-  paste0( mod,".",met,".LM.p_val.FDR" )
      c.align.pv.caper <-  paste0( mod,".",met,".CAPER.p_val.FDR" )
      c.align.pv.mdmr_mash <-  paste0( mod,".",met,".MDMR_mash.","pv.FDR" )
      c.align.pv.mdmr_phylo <-  paste0( mod,".",met,".MDMR_phylo.","pv.FDR" )
      c.align.slp.lm <- paste0( mod,".",met,".LM.slope") 
      c.align.pos.lm <- paste0( mod,".",met,".LM.is_pos") 
      c.align.slp.caper <- paste0( mod,".",met,".CAPER.slope") 
      c.align.pos.caper <- paste0( mod,".",met,".CAPER.is_pos") 
      
      c.denov.pv.lm <-  paste0( mod+12,".",met,".LM.p_val.FDR" )
      c.denov.pv.caper <-  paste0( mod+12,".",met,".CAPER.p_val.FDR" )
      c.denov.pv.mdmr_mash <-  paste0( mod+12,".",met,".MDMR_mash.","pv.FDR" )
      c.denov.pv.mdmr_phylo <-  paste0( mod+12,".",met,".MDMR_phylo.","pv.FDR" )
      c.denov.slp.lm <- paste0( mod+12,".",met,".LM.slope") 
      c.denov.pos.lm <- paste0( mod+12,".",met,".LM.is_pos") 
      c.denov.slp.caper <- paste0( mod+12,".",met,".CAPER.slope") 
      c.denov.pos.caper <- paste0( mod+12,".",met,".CAPER.is_pos") 
      
      cmpr <- model.p.dos[ , c.align.pv.mdmr_phylo ] < cutoff
      cmpr[is.na(cmpr)] <- FALSE
      ltCutoff <- rownames(model.p.dos)[  cmpr ]
      ltCutoff.neg.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] < 0]
      ltCutoff.pos.align <- ltCutoff[model.p.dos[ltCutoff,c.align.slp.caper] > 0]
      ltCutoff.neg.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] < 0]
      ltCutoff.pos.denov <- ltCutoff[model.p.dos[ltCutoff,c.denov.slp.caper] > 0]
      if(met == "f_01"){
        significant_og_counter_matrix.f_01[  19 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_01[  20 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_01[  21 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  22 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_01[  23 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_01[  24 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_01 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_01 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_01 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_01 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_01 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_01 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_02"){
        significant_og_counter_matrix.f_02[  19 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_02[  20 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_02[  21 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_02[  22 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_02[  23 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_02[  24 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_02 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_02 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_02 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_02 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_02 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_02 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_04"){
        significant_og_counter_matrix.f_04[  19 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_04[  20 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_04[  21 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_04[  22 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_04[  23 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_04[  24 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_04 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_04 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_04 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_04 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_04 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_04 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      if(met == "f_05"){
        significant_og_counter_matrix.f_05[  19 ,  mod ] <- length(ltCutoff.pos.denov)
        significant_og_counter_matrix.f_05[  20 ,  mod ] <- length(ltCutoff.neg.denov)
        significant_og_counter_matrix.f_05[  21 ,  mod ] <- length(ltCutoff.pos.align)
        significant_og_counter_matrix.f_05[  22 ,  mod ] <- length(ltCutoff.neg.align)
        significant_og_counter_matrix.f_05[  23 ,  mod ] <- sum(ltCutoff.pos.denov %in% ltCutoff.pos.align)
        significant_og_counter_matrix.f_05[  24 ,  mod ] <- sum(ltCutoff.neg.denov %in% ltCutoff.neg.align)
        ltCutoff.neg.align.f_05 <- ltCutoff.neg.align
        ltCutoff.pos.align.f_05 <- ltCutoff.pos.align
        ltCutoff.neg.denov.f_05 <- ltCutoff.neg.denov
        ltCutoff.pos.denov.f_05 <- ltCutoff.pos.denov
        ltCutoff.neg.shared.f_05 <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
        ltCutoff.pos.shared.f_05 <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
      }
      mtc <- mtc + 1
      
    }
    
    # a <- ltCutoff.neg.denov.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    # a <- a[a %in% ltCutoff.neg.align.f_02]
    # ltCutoff.neg.denov.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    # significant_og_counter_matrix.shared[ 1 , mod ] <- length( ltCutoff.neg.denov.f_shared )
    # 
    # a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    # a <- a[a %in% ltCutoff.pos.align.f_02]
    # ltCutoff.pos.denov.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    # significant_og_counter_matrix.shared[ 2 , mod ] <- length( ltCutoff.pos.denov.f_shared )
    # 
    # a <- ltCutoff.neg.align.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    # a <- a[a %in% ltCutoff.neg.align.f_02]
    # ltCutoff.neg.align.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    # significant_og_counter_matrix.shared[ 3 , mod ] <- length( ltCutoff.neg.align.f_shared )
    # 
    # a <- ltCutoff.pos.denov.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    # a <- a[a %in% ltCutoff.pos.align.f_02]
    # ltCutoff.pos.align.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    # significant_og_counter_matrix.shared[ 4 , mod ] <- length( ltCutoff.pos.align.f_shared )
    # 
    # a <- ltCutoff.neg.shared.f_05[ltCutoff.neg.align.f_05 %in% ltCutoff.neg.align.f_04]
    # a <- a[a %in% ltCutoff.neg.align.f_02]
    # ltCutoff.neg.shared.f_shared <- a[a %in% ltCutoff.neg.align.f_01]
    # significant_og_counter_matrix.shared[ 5 , mod ] <- length( ltCutoff.neg.shared.f_shared )
    # 
    # a <- ltCutoff.pos.shared.f_05[ltCutoff.pos.align.f_05 %in% ltCutoff.pos.align.f_04]
    # a <- a[a %in% ltCutoff.pos.align.f_02]
    # ltCutoff.pos.shared.f_shared <- a[a %in% ltCutoff.pos.align.f_01]
    # significant_og_counter_matrix.shared[ 6 , mod ] <- length( ltCutoff.pos.shared.f_shared )
    
    
    significant_og_counter_matrix.shared[ 19 , mod ] <- sum( ltCutoff.neg.denov.f_05 %in%  ltCutoff.neg.denov.f_02 )
    significant_og_counter_matrix.shared[ 20 , mod ] <- sum( ltCutoff.pos.denov.f_05 %in%  ltCutoff.pos.denov.f_02 )
    significant_og_counter_matrix.shared[ 21 , mod ] <- sum( ltCutoff.neg.align.f_05 %in%  ltCutoff.neg.align.f_02 )
    significant_og_counter_matrix.shared[ 22 , mod ] <- sum( ltCutoff.pos.denov.f_05 %in%  ltCutoff.pos.denov.f_02 )
    significant_og_counter_matrix.shared[ 23 , mod ] <- sum( ltCutoff.neg.shared.f_05 %in% ltCutoff.neg.shared.f_02 )
    significant_og_counter_matrix.shared[ 24 , mod ] <- sum( ltCutoff.pos.shared.f_05 %in% ltCutoff.pos.shared.f_02 )
    
    #write.table(significant_og_counter_matrix.shared, paste0("significant_og_counter_matrix.",mod,".f_05.tsv"), sep='\t')
    
    ltCutoff.neg.align.f_shared.mdmr_phylo <- ltCutoff.neg.align
    ltCutoff.pos.align.f_shared.mdmr_phylo <- ltCutoff.pos.align
    ltCutoff.neg.denov.f_shared.mdmr_phylo <- ltCutoff.neg.denov
    ltCutoff.pos.denov.f_shared.mdmr_phylo <- ltCutoff.pos.denov
    ltCutoff.neg.shared.f_shared.mdmr_phylo <- ltCutoff.neg.denov[ltCutoff.neg.denov %in% ltCutoff.neg.align]
    ltCutoff.pos.shared.f_shared.mdmr_phylo <- ltCutoff.pos.denov[ltCutoff.pos.denov %in% ltCutoff.pos.align]
    
    write.table(ltCutoff.neg.align.f_shared.mdmr_phylo,   paste0("outConvert2/ltCutoff.neg.align.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.align.f_shared.mdmr_phylo,   paste0("outConvert2/ltCutoff.pos.align.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.denov.f_shared.mdmr_phylo,   paste0("outConvert2/ltCutoff.neg.denov.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.denov.f_shared.mdmr_phylo,   paste0("outConvert2/ltCutoff.pos.denov.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_shared.mdmr_phylo,   paste0("outConvert2/ltCutoff.neg.shared.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_shared.mdmr_phylo,   paste0("outConvert2/ltCutoff.pos.shared.f_shared.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_02"
    write.table( ltCutoff.neg.align.f_02,   paste0("outConvert2/ltCutoff.neg.align.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_02,   paste0("outConvert2/ltCutoff.pos.align.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_02,   paste0("outConvert2/ltCutoff.neg.denov.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_02,   paste0("outConvert2/ltCutoff.pos.denov.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_02,   paste0("outConvert2/ltCutoff.neg.shared.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_02,   paste0("outConvert2/ltCutoff.pos.shared.f_02.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    met <- "f_05"
    write.table( ltCutoff.neg.align.f_05,   paste0("outConvert2/ltCutoff.neg.align.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.align.f_05,   paste0("outConvert2/ltCutoff.pos.align.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.neg.denov.f_05,   paste0("outConvert2/ltCutoff.neg.denov.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table( ltCutoff.pos.denov.f_05,   paste0("outConvert2/ltCutoff.pos.denov.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.neg.shared.f_05,   paste0("outConvert2/ltCutoff.neg.shared.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ltCutoff.pos.shared.f_05,   paste0("outConvert2/ltCutoff.pos.shared.f_05.mdmr_phylo.",cutoff,".mod_",mod,".tsv"), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    
    
  }
  
  write.table(significant_og_counter_matrix.f_01,   paste0("outConvert2/significant_og_counter_matrix.f_01.",cutoff,".tsv"), sep='\t')
  write.table(significant_og_counter_matrix.f_02,   paste0("outConvert2/significant_og_counter_matrix.f_02.",cutoff,".tsv"), sep='\t')
  write.table(significant_og_counter_matrix.f_04,   paste0("outConvert2/significant_og_counter_matrix.f_04.",cutoff,".tsv"), sep='\t')
  write.table(significant_og_counter_matrix.f_05,   paste0("outConvert2/significant_og_counter_matrix.f_05.",cutoff,".tsv"), sep='\t')
  write.table(significant_og_counter_matrix.shared, paste0("outConvert2/significant_og_counter_matrix.shared.",cutoff,".tsv"), sep='\t')
}



#######################################

cutoff <- 0.005

outtt02 <- matrix(ncol=9,nrow=24)
outtt02[is.na(outtt02)] <- 0


outtt05 <- matrix(ncol=9,nrow=24)
outtt05[is.na(outtt05)] <- 0

outtt <- matrix(ncol=9,nrow=24)
outtt[is.na(outtt)] <- 0
mm <- 1
for(mod in c(1,2,3,5,6,7,9,10,11)){
  for(roww in c(1:24)){
    if(roww %in% c(1,2,7,8,13,14,19,20)){ abund <- "denov" }
    else if(roww %in% c(3,4,9,10,15,16,21,22)){ abund <- "align" }
    else{ abund <- "shared" }
    if(roww %% 2 == 0){ dir <- "pos" }
    else{ dir <- "neg"}
    if(roww %in% c(1:6)){ typee <- "lm"  }
    if(roww %in% c(7:12)){ typee <- "caper"  }
    if(roww %in% c(13:18)){ typee <- "mdmr_mash"  }
    if(roww %in% c(19:24)){ typee <- "mdmr_phylo"  }
    outqqq <- NA
    
    bmg_a <- NA
    #bmg_a <- read.table(paste0("outConvert2/ltCutoff.",dir,".",abund,".f_02.",typee,".",cutoff,".mod_",mod,".tsv"), sep='\t')[,1]
    tryCatch( bmg_a <- read.table(paste0("outConvert2/ltCutoff.",dir,".",abund,".f_02.",typee,".",cutoff,".mod_",mod,".tsv"), sep='\t')[,1] ,
                        error = function(e) message("No bueno2."))
    aw_a <- NA
    #aw_a <- read.table(paste0("outConvert2/ltCutoff.",dir,".",abund,".f_05.",typee,".",cutoff,".mod_",mod,".tsv"), sep='\t')[,1]
    tryCatch( aw_a <- read.table(paste0("outConvert2/ltCutoff.",dir,".",abund,".f_05.",typee,".",cutoff,".mod_",mod,".tsv"), sep='\t')[,1] , 
              error = function(e) message("No bueno5."))
    
    if(!is.na(bmg_a)){
      outtt02[roww,mm] <- length( bmg_a )
    }  
    if(!is.na(aw_a)){
      outtt05[roww,mm] <- length( aw_a  )
    }
    if(!is.na(bmg_a) && !is.na(aw_a)){
      outtt[roww,mm] <- length( bmg_a[ bmg_a %in% aw_a ] )
    }
    
    #tryCatch( outqqq   <- sum( read.table(paste0("outConvert2/ltCutoff.",dir,".",abund,".f_02.",typee,".",cutoff,".mod_",mod,".tsv"), sep='\t')[,1] %in%  ) ,
    #          error = function(e) message("No bueno."))
    #if(!is.na(outqqq)){
    #  outtt[roww,mm] <- outqqq
    #}else{
    #  outtt[roww,mm] <- 0
    #}
  }
  mm <- mm + 1
}
 





############ 1.2 : Convert the identifiers to position information that can be compared to the massive gene reference database ############ 

# pull in match of indices to defined IDs
odb10v1_genes.2col    <- read.table("AvianLongevity_ManuscriptData/odb10v1_genes.2col.tsv", sep=" ", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
# 
odb10v1_genes.2col <- odb10v1_genes.2col[ grepl(":",odb10v1_genes.2col[,1]) , ]
odb10v1_genes  <- odb10v1_genes.2col[,2]
names(odb10v1_genes) <- odb10v1_genes.2col[,1]

#rm(odb10v1_genes.2col)

#cutoff <- 0.05 #for(cutoff in c(0.05,0.005)){
#fun <- "f_02" # f_05
#mod <- 11 # 11 #for(mod in c(3,7,11)){
#abund <- "shared" # "align" # 
#mname <- "caper" #for(mname in c("lm","caper","mdmr_mash","mdmr_phylo")){


cutoff <- 0.05
for(mname in c("caper","mdmr_mash","mdmr_phlyo")){
  #for(mod in c(11,7,3)){
  for(mod in c(6,2,9,5,1)){ #c(10,6,2,9,5,1)){
    for(abund in c("shared","align","denov")){
      #for(fun in c("f_02","f_05")){
      fun <- "f_02"
      if(! ( (abund == "shared") && (mod == 11) ) ){  
        print("!!!!!!")
        print("!!!!!!")
        print("!!!!!!")
        print("!!!!!!")
        print("!!!!!!")
        print("! NEW RUN !")
        print(paste(cutoff,mname,mod,abund,fun))
        print("!!!!!!")
        print("!!!!!!")
        print("!!!!!!")
        print("!!!!!!")
        print("!!!!!!")
        
        
        cur.pos <- list()
        cur.neg <- list()
        # ltCutoff.neg.align.f_01.caper.0.05.mod_7.tsv
        tryCatch( cur.pos   <- read.table( paste0("outConvert/ltCutoff.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") ,sep='\t',stringsAsFactors = FALSE)[,1] ,  error = function(e) message(paste("Empty input", paste0("ltCutoff.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv"))))
        tryCatch( cur.neg   <- read.table( paste0("outConvert/ltCutoff.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") ,sep='\t',stringsAsFactors = FALSE)[,1] ,  error = function(e) message(paste("Empty input", paste0("ltCutoff.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv"))))
        
        #tryCatch( cur.pos   <- read.table( paste0("outConvert/ltCutoff.pos.",abund,fun,".",mname,".",cutoff,".mod_",mod,".tsv") ,sep='\t',stringsAsFactors = FALSE)[,1] ,  error = function(e) message(paste("Empty input", paste0("ltCutoff.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv"))))
        #tryCatch( cur.neg   <- read.table( paste0("outConvert/ltCutoff.neg.",abund,fun,".",mname,".",cutoff,".mod_",mod,".tsv") ,sep='\t',stringsAsFactors = FALSE)[,1] ,  error = function(e) message(paste("Empty input", paste0("ltCutoff.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv"))))
        cur.all   <- c(cur.neg,cur.pos)
        
        pos.cur_og_id      <- v9_v10_OGs_map[ v9_v10_OGs_map[,2] %in% cur.pos , c(2,3,5) ]
        neg.cur_og_id      <- v9_v10_OGs_map[ v9_v10_OGs_map[,2] %in% cur.neg , c(2,3,5) ]
        all.cur_og_id      <- v9_v10_OGs_map[ v9_v10_OGs_map[,2] %in% cur.all , c(2,3,5) ]
        
        pos.cur_og_col      <- odb10v1_OG2genes[ odb10v1_OG2genes[,1] %in% pos.cur_og_id[,2] ,  ]
        neg.cur_og_col      <- odb10v1_OG2genes[ odb10v1_OG2genes[,1] %in% neg.cur_og_id[,2] ,  ]
        all.cur_og_col      <- odb10v1_OG2genes[ odb10v1_OG2genes[,1] %in% all.cur_og_id[,2] ,  ]
        
        # append negative rows
        
        neg.cur_og_col[,3] <- neg.cur_og_col[,2] %in% odb10v1_genes.2col[,1]
        
        neg.cur_og_col[,4] <- rep(NA,dim(neg.cur_og_col)[1])
        nn <- 1
        for(n in neg.cur_og_col[,2]){
          
          if(nn%%100==0){print((nn/dim(neg.cur_og_col)[1])*100)}
          if(neg.cur_og_col[nn,3] == TRUE){
            neg.cur_og_col[nn,4] <- odb10v1_genes.2col[ odb10v1_genes.2col[,1]==n , 2]
          }
          nn <- nn + 1
        }
        
        neg.cur_og_col[,5] <- rep(NA,dim(neg.cur_og_col)[1])
        nn <- 1
        for(d in neg.cur_og_col[,1]){
          if(nn%%100==0){print((nn/dim(neg.cur_og_col)[1])*100)}
          if( d %in% neg.cur_og_id[,2] ){
            neg.cur_og_col[nn,5] <-  neg.cur_og_id[ neg.cur_og_id[,2]==d , 1 ]
          }
          nn <- nn + 1
        }
        
        neg.cur_og_col[,6] <- rep(NA,dim(neg.cur_og_col)[1])
        nn<-1
        for(d in neg.cur_og_col[,5]){
          neg.cur_og_col[nn,6] <- model.p.dos[ d , "11.f_05.CAPER.p_val.FDR" ]
          nn <- nn + 1
        }
        
        neg.cur_og_col[,7] <- rep(NA,dim(neg.cur_og_col)[1])
        nn<-1
        for(d in neg.cur_og_col[,5]){
          neg.cur_og_col[nn,7] <- model.p.dos[ d , "11.f_05.CAPER.slope" ]
          nn <- nn + 1
        }
        
        neg.cur_og_col[,c(8:13)] <- rep(NA,dim(neg.cur_og_col)[1])
        nn<-1
        prev <- "x"
        for(d in neg.cur_og_col[,1]){
          if(d == prev){
            neg.cur_og_col[nn,c(8:13)]  <- neg.cur_og_col[nn-1,c(8:13)] 
          }else{
            print((nn/dim(neg.cur_og_col)[1])*100)
            if( "biological_process" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )       { neg.cur_og_col[nn,8]  <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
            if( "cellular_component" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )       { neg.cur_og_col[nn,9]  <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
            if( "molecular_function" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )       { neg.cur_og_col[nn,10] <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
            if( "unclassified GO terms" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )    { neg.cur_og_col[nn,11] <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
            if( "functional_category" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )      { neg.cur_og_col[nn,12] <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
            if( "interpro_domains" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )         { neg.cur_og_col[nn,13] <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
          }
          prev <- d
          nn <- nn + 1
        }
        
        colnames(neg.cur_og_col) <- c("Ver8ID","POS_ID","IsSymDef","Sym","OG_ID","P_FDR","SLOPE","BP","CC","MF","UNCLASSIFIED","FUNC_CAT","INTERPRO_DOMAINS")
        
        write.table(neg.cur_og_col, paste0("og_conversion_table.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") )
        
        write.table(unique(neg.cur_og_col$Sym[!is.na(neg.cur_og_col$Sym)]), 
                    file = paste0("og_conversion_UNIQsymbolsWITHloc.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv"),
                    row.names = FALSE, col.names = FALSE
        )
        
        write.table(unique(neg.cur_og_col$Sym[!is.na(neg.cur_og_col$Sym)])[ !grepl("LOC",unique(neg.cur_og_col$Sym[!is.na(neg.cur_og_col$Sym)])) ], 
                    file = paste0("og_conversion_UNIQsymbols.neg.a",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv"),
                    row.names = FALSE, col.names = FALSE
        )
        
        geneListX <- c(unique(neg.cur_og_col$Sym[!is.na(neg.cur_og_col$Sym)])[ !grepl("LOC",unique(neg.cur_og_col$Sym[!is.na(neg.cur_og_col$Sym)])) ],"FAKE")
        #c("AIG1","DYNLT1","EZR","GNS","IKBKE","LEMD3","MSRB3","NMBR","PEX3","RASSF3","RPL18A","RXYLT1","TBK1","WIF1","XPOT","FAKE")
        geneList <- factor(c(as.integer(rep(1,length(geneListX)-1)),0))
        names(geneList) <- geneListX
        
        #################### STOP ###############################################################################
        #################### STOP ###############################################################################
        #################### STOP ###############################################################################
        #################### ALWAYS STOP AND RUN THIS BY HAND UNTIL IT WORKS! WTF BIOMART! ######################
        results <- NA
        results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'go_id'), 
                         filter = "hgnc_symbol",
                         values = geneListX[1:length(geneListX)-1],
                         #c("AIG1","DYNLT1","EZR","GNS","IKBKE","LEMD3","MSRB3","NMBR","PEX3","RASSF3","RPL18A","RXYLT1","TBK1","WIF1","XPOT"), 
                         #geneList,
                         mart = bm)
        
        geneID2GO <- by(results$go_id, results$hgnc_symbol, function(x) as.character(x))
        
        queryGOdata.bp <- new( "topGOdata", 
                               description = "Simple session",
                               ontology = "BP",
                               allGenes = geneList,
                               geneSel = topDiffGenes,
                               nodeSize = 2, # <- limit to annotations held over > 10 genes
                               annot=annFUN.gene2GO,
                               gene2GO=geneID2GO
        )
        queryGOdata.mf <- new( "topGOdata", 
                               description = "Simple session",
                               ontology = "MF",
                               allGenes = geneList,
                               geneSel = topDiffGenes,
                               nodeSize = 2, # <- limit to annotations held over > 10 genes
                               annot=annFUN.gene2GO,
                               gene2GO=geneID2GO
        )
        queryGOdata.cc <- new( "topGOdata", 
                               description = "Simple session",
                               ontology = "CC",
                               allGenes = geneList,
                               geneSel = topDiffGenes,
                               nodeSize = 2, # <- limit to annotations held over > 10 genes
                               annot=annFUN.gene2GO,
                               gene2GO=geneID2GO
        )
        
        queryFisher.classic.fisher.bp <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "fisher")
        queryFisher.classic.ks.bp     <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "ks"    )
        queryFisher.elim.ks.bp     <- runTest(queryGOdata.bp, algorithm = "elim",    statistic = "ks"    )
        
        queryFisher.classic.fisher.mf <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "fisher")
        queryFisher.classic.ks.mf     <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "ks"    )
        queryFisher.elim.ks.mf     <- runTest(queryGOdata.mf, algorithm = "elim",    statistic = "ks"    )
        
        queryFisher.classic.fisher.cc <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "fisher")
        queryFisher.classic.ks.cc     <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "ks"    )
        queryFisher.elim.ks.cc     <- runTest(queryGOdata.cc, algorithm = "elim",    statistic = "ks"    )
        
        allRes.bp <- GenTable( queryGOdata.bp, classicFisher = queryFisher.classic.fisher.bp, classicKS = queryFisher.classic.ks.bp, elimKS = queryFisher.elim.ks.bp, 
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.bp)) )          
        
        allRes.mf <- GenTable( queryGOdata.mf, classicFisher = queryFisher.classic.fisher.mf, classicKS = queryFisher.classic.ks.mf, elimKS = queryFisher.elim.ks.mf, 
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.mf)) )          
        
        allRes.cc <- GenTable( queryGOdata.cc, classicFisher = queryFisher.classic.fisher.cc, classicKS = queryFisher.classic.ks.cc, elimKS = queryFisher.elim.ks.cc, 
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.cc)) )          
        
        
        #write.table(allRes.bp, file = "GO_table.BP.tsv", sep = '\t')
        #write.table(allRes.mf, file = "GO_table.MF.tsv", sep = '\t')
        #write.table(allRes.cc, file = "GO_table.CC.tsv", sep = '\t')
        write.table( allRes.bp , file = paste0("TopGO_BP.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") , sep = '\t', row.names = FALSE )
        write.table( allRes.mf , file = paste0("TopGO_MF.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") , sep = '\t', row.names = FALSE )
        write.table( allRes.cc , file = paste0("TopGO_CC.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") , sep = '\t', row.names = FALSE )
        
        pValue.classic.bp <- score(queryFisher.classic.ks.bp)
        pValue.elim.bp <- score(queryFisher.elim.ks.bp)[names(score(queryFisher.classic.ks.bp))]
        gstat.bp <- termStat(queryGOdata.bp, names(pValue.classic.bp))
        gSize.bp <- gstat.bp$Annotated / max(gstat.bp$Annotated) * 4
        
        png( paste0("TopGO_BP.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #"pValue.classic.BP.png" , height=1000,width=1000 )
        plot(pValue.classic.bp, pValue.elim.bp, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.bp, 
             col = redgreen(max(gstat.bp$Annotated))[ gstat.bp$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_BP_log.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "LOGpValue.classic.BP.png" , height=1000,width=1000 )
        plot(log(pValue.classic.bp), log(pValue.elim.bp), xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.bp, 
             col = redgreen(max(gstat.bp$Annotated))[ gstat.bp$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_BP_Tree_0.01.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.BP.0.01.png" , height=1000,width=1000 )
        #if(sum(allRes.bp$classicKS < 0.01) > 5){
        showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks.bp), firstSigNodes = sum(allRes.bp$classicKS < 0.01), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.01"))
        dev.off()
        #}else if(sum(allRes.bp$classicKS < 0.05) > 5 & sum(allRes.bp$classicKS < 0.05) < 35){
        
        png( paste0("TopGO_BP_Tree_0.05.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.BP.0.05.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks.bp), firstSigNodes = sum(allRes.bp$classicKS < 0.05), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.05"))
        dev.off()
        #}else{
        
        png( paste0("TopGO_BP_Tree_x10.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.BP.x10.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks.bp), firstSigNodes = 10, useInfo = 'all')
        title(paste0("Biological Processes\nCutoff count = 10"))
        dev.off() 
        #}
        
        
        
        pValue.classic.mf <- score(queryFisher.classic.ks.mf)
        pValue.elim.mf <- score(queryFisher.elim.ks.mf)[names(score(queryFisher.classic.ks.mf))]
        gstat.mf <- termStat(queryGOdata.mf, names(pValue.classic.mf))
        gSize.mf <- gstat.mf$Annotated / max(gstat.mf$Annotated) * 4
        
        png( paste0("TopGO_MF.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "pValue.classic.MF.png" , height=1000,width=1000 )
        plot(pValue.classic.mf, pValue.elim.mf, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.mf, 
             col = redgreen(max(gstat.mf$Annotated))[ gstat.mf$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_MF_log.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "LOGpValue.classic.MF.png" , height=1000,width=1000 )
        plot(log(pValue.classic.mf), log(pValue.elim.mf), xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.mf, 
             col = redgreen(max(gstat.mf$Annotated))[ gstat.mf$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_MF_Tree_0.01.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.MF.0.01.png" , height=1000,width=1000 )
        #if(sum(allRes.mf$classicKS < 0.01) > 5){
        showSigOfNodes(queryGOdata.mf, score(queryFisher.elim.ks.mf), firstSigNodes = sum(allRes.mf$classicKS < 0.01), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.01"))
        dev.off()
        #}else if(sum(allRes.mf$classicKS < 0.05) > 5 & sum(allRes.mf$classicKS < 0.05) < 35){
        
        png( paste0("TopGO_MF_Tree_0.05.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.MF.0.05.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.mf, score(queryFisher.elim.ks.mf), firstSigNodes = sum(allRes.mf$classicKS < 0.05), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.05"))
        dev.off()
        #}else{
        
        png( paste0("TopGO_MF_Tree_x10.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.MF.x10.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.mf, score(queryFisher.elim.ks.mf), firstSigNodes = 10, useInfo = 'all')
        title(paste0("Biological Processes\nCutoff count = 10"))
        #}
        dev.off() 
        
        
        
        pValue.classic.cc <- score(queryFisher.classic.ks.cc)
        pValue.elim.cc <- score(queryFisher.elim.ks.cc)[names(score(queryFisher.classic.ks.cc))]
        gstat.cc <- termStat(queryGOdata.cc, names(pValue.classic.cc))
        gSize.cc <- gstat.cc$Annotated / max(gstat.cc$Annotated) * 4
        
        png( paste0("TopGO_CC.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "pValue.classic.CC.png" , height=1000,width=1000 )
        plot(pValue.classic.cc, pValue.elim.cc, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.cc, 
             col = redgreen(max(gstat.cc$Annotated))[ gstat.cc$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_CC_log.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "LOGpValue.classic.CC.png" , height=1000,width=1000 )
        plot(log(pValue.classic.cc), log(pValue.elim.cc), xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.cc, 
             col = redgreen(max(gstat.cc$Annotated))[ gstat.cc$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_CC_log.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.CC.0.01.png" , height=1000,width=1000 )
        #if(sum(allRes.cc$classicKS < 0.01) > 5){
        showSigOfNodes(queryGOdata.cc, score(queryFisher.elim.ks.cc), firstSigNodes = sum(allRes.cc$classicKS < 0.01), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.01"))
        dev.off()
        #}else if(sum(allRes.cc$classicKS < 0.05) > 5 & sum(allRes.cc$classicKS < 0.05) < 35){
        
        png( paste0("TopGO_CC_Tree_0.05.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.CC.0.05.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.cc, score(queryFisher.elim.ks.cc), firstSigNodes = sum(allRes.cc$classicKS < 0.05), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.05"))
        dev.off()
        #}else{
        
        
        png( paste0("TopGO_CC_Tree_x10.neg.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.CC.x10.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.cc, score(queryFisher.elim.ks.cc), firstSigNodes = 10, useInfo = 'all')
        title(paste0("Biological Processes\nCutoff count = 10 (no GO terms with p < 0.05)"))
        #}
        dev.off() 
        
        #colnames(neg.cur_og_col) <- c("Ver8ID","POS_ID","IsSymDef","Sym","OG_ID","P_FDR","SLOPE")
        
        ####################################################################################################################################################################################################
        # append positive rows #############################################################################################################################################################################
        ####################################################################################################################################################################################################
        
        nn <- 1
        pos.cur_og_col[,3] <- pos.cur_og_col[,2] %in% odb10v1_genes.2col[,1]
        
        pos.cur_og_col[,4] <- rep(NA,dim(pos.cur_og_col)[1])
        nn <- yep <- nope <- 1
        for(n in pos.cur_og_col[,2]){
          
          if(nn%%100==0){print((nn/dim(pos.cur_og_col)[1])*100)}
          if(pos.cur_og_col[nn,3] == TRUE){
            pos.cur_og_col[nn,4] <- odb10v1_genes.2col[ odb10v1_genes.2col[,1]==n , 2]
          }
          nn <- nn + 1
        }
        
        pos.cur_og_col[,5] <- rep(NA,dim(pos.cur_og_col)[1])
        nn <- 1
        for(d in pos.cur_og_col[,1]){
          if(nn%%100==0){print((nn/dim(pos.cur_og_col)[1])*100)}
          if( d %in% pos.cur_og_id[,2] ){
            pos.cur_og_col[nn,5] <-  pos.cur_og_id[ pos.cur_og_id[,2]==d , 1 ]
          }
          nn <- nn + 1
        }
        
        
        #sort(table(pos.cur_og_col[,4]))
        
        # model.p.dos[ rownames( model.p.dos[ model.p.dos$`11.f_05.CAPER.p_val.FDR` < cutoff , ] ) , "11.f_05.CAPER.p_val.FDR" ]
        
        pos.cur_og_col[,6] <- rep(NA,dim(pos.cur_og_col)[1])
        nn<-1
        for(d in pos.cur_og_col[,5]){
          pos.cur_og_col[nn,6] <- model.p.dos[ d , "11.f_05.CAPER.p_val.FDR" ]
          nn <- nn + 1
        }
        
        pos.cur_og_col[,7] <- rep(NA,dim(pos.cur_og_col)[1])
        nn<-1
        for(d in pos.cur_og_col[,5]){
          pos.cur_og_col[nn,7] <- model.p.dos[ d , "11.f_05.CAPER.slope" ]
          nn <- nn + 1
        }
        
        pos.cur_og_col[,c(8:13)] <- rep(NA,dim(pos.cur_og_col)[1])
        nn<-1
        prev <- "x"
        for(d in pos.cur_og_col[,1]){
          if(nn%%100==0){}
          if(d == prev){
            if( "biological_process" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )       { pos.cur_og_col[nn,8]  <- pos.cur_og_col[nn-1,8] }
            if( "cellular_component" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )       { pos.cur_og_col[nn,9]  <- pos.cur_og_col[nn-1,9] }
            if( "molecular_function" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )       { pos.cur_og_col[nn,10] <- pos.cur_og_col[nn-1,10] }
            if( "unclassified GO terms" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )    { pos.cur_og_col[nn,11] <- pos.cur_og_col[nn-1,11] }
            if( "functional_category" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )      { pos.cur_og_col[nn,12] <- pos.cur_og_col[nn-1,12] }
            if( "interpro_domains" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )         { pos.cur_og_col[nn,13] <- pos.cur_og_col[nn-1,13] }
          }else{
            print((nn/dim(pos.cur_og_col)[1])*100)
            if( "biological_process" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )       { pos.cur_og_col[nn,8]  <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
            if( "cellular_component" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )       { pos.cur_og_col[nn,9]  <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
            if( "molecular_function" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )       { pos.cur_og_col[nn,10] <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
            if( "unclassified GO terms" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )    { pos.cur_og_col[nn,11] <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
            if( "functional_category" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )      { pos.cur_og_col[nn,12] <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
            if( "interpro_domains" %in% odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 2 ] )         { pos.cur_og_col[nn,13] <- paste(odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] == d , 3 ],collapse=",") }
          }
          prev <- d
          nn <- nn + 1
        }
        
        colnames(pos.cur_og_col) <- c("Ver8ID","POS_ID","IsSymDef","Sym","OG_ID","P_FDR","SLOPE","BP","CC","MF","UNCLASSIFIED","FUNC_CAT","INTERPRO_DOMAINS")
        
        write.table(pos.cur_og_col, paste0("og_conversion_table.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") )
        
        write.table(unique(pos.cur_og_col$Sym[!is.na(pos.cur_og_col$Sym)]), 
                    file = paste0("og_conversion_UNIQsymbolsWITHloc.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv"),
                    row.names = FALSE, col.names = FALSE
        )
        
        write.table(unique(pos.cur_og_col$Sym[!is.na(pos.cur_og_col$Sym)])[ !grepl("LOC",unique(pos.cur_og_col$Sym[!is.na(pos.cur_og_col$Sym)])) ], 
                    file = paste0("og_conversion_UNIQsymbols.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv"),
                    row.names = FALSE, col.names = FALSE
        )
        
        geneListX <- c(unique(pos.cur_og_col$Sym[!is.na(pos.cur_og_col$Sym)])[ !grepl("LOC",unique(pos.cur_og_col$Sym[!is.na(pos.cur_og_col$Sym)])) ],"FAKE")
        #c("AIG1","DYNLT1","EZR","GNS","IKBKE","LEMD3","MSRB3","NMBR","PEX3","RASSF3","RPL18A","RXYLT1","TBK1","WIF1","XPOT","FAKE")
        geneList <- factor(c(as.integer(rep(1,length(geneListX)-1)),0))
        names(geneList) <- geneListX
        
        #################### STOP ###############################################################################
        #################### STOP ###############################################################################
        #################### STOP ###############################################################################
        #################### ALWAYS STOP AND RUN THIS BY HAND UNTIL IT WORKS! WTF BIOMART! ######################
        results <- NA
        results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'go_id'), 
                         filter = "hgnc_symbol",
                         values = geneListX[1:length(geneListX)-1], #c("AIG1","DYNLT1","EZR","GNS","IKBKE","LEMD3","MSRB3","NMBR","PEX3","RASSF3","RPL18A","RXYLT1","TBK1","WIF1","XPOT"),
                         mart = bm)
        
        geneID2GO <- by(results$go_id, results$hgnc_symbol, function(x) as.character(x))
        
        queryGOdata.bp <- new( "topGOdata", 
                               description = "Simple session",
                               ontology = "BP",
                               allGenes = geneList,
                               geneSel = topDiffGenes,
                               nodeSize = 2, # <- limit to annotations held over > 10 genes
                               annot=annFUN.gene2GO,
                               gene2GO=geneID2GO
        )
        queryGOdata.mf <- new( "topGOdata", 
                               description = "Simple session",
                               ontology = "MF",
                               allGenes = geneList,
                               geneSel = topDiffGenes,
                               nodeSize = 2, # <- limit to annotations held over > 10 genes
                               annot=annFUN.gene2GO,
                               gene2GO=geneID2GO
        )
        queryGOdata.cc <- new( "topGOdata", 
                               description = "Simple session",
                               ontology = "CC",
                               allGenes = geneList,
                               geneSel = topDiffGenes,
                               nodeSize = 2, # <- limit to annotations held over > 10 genes
                               annot=annFUN.gene2GO,
                               gene2GO=geneID2GO
        )
        
        queryFisher.classic.fisher.bp <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "fisher")
        queryFisher.classic.ks.bp     <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "ks"    )
        queryFisher.elim.ks.bp     <- runTest(queryGOdata.bp, algorithm = "elim",    statistic = "ks"    )
        
        queryFisher.classic.fisher.mf <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "fisher")
        queryFisher.classic.ks.mf     <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "ks"    )
        queryFisher.elim.ks.mf     <- runTest(queryGOdata.mf, algorithm = "elim",    statistic = "ks"    )
        
        queryFisher.classic.fisher.cc <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "fisher")
        queryFisher.classic.ks.cc     <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "ks"    )
        queryFisher.elim.ks.cc     <- runTest(queryGOdata.cc, algorithm = "elim",    statistic = "ks"    )
        
        allRes.bp <- GenTable( queryGOdata.bp, classicFisher = queryFisher.classic.fisher.bp, classicKS = queryFisher.classic.ks.bp, elimKS = queryFisher.elim.ks.bp, 
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.bp)) )          
        
        allRes.mf <- GenTable( queryGOdata.mf, classicFisher = queryFisher.classic.fisher.mf, classicKS = queryFisher.classic.ks.mf, elimKS = queryFisher.elim.ks.mf, 
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.mf)) )          
        
        allRes.cc <- GenTable( queryGOdata.cc, classicFisher = queryFisher.classic.fisher.cc, classicKS = queryFisher.classic.ks.cc, elimKS = queryFisher.elim.ks.cc, 
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.cc)) )          
        
        
        #write.table(allRes.bp, file = "GO_table.BP.tsv", sep = '\t')
        #write.table(allRes.mf, file = "GO_table.MF.tsv", sep = '\t')
        #write.table(allRes.cc, file = "GO_table.CC.tsv", sep = '\t')
        write.table( allRes.bp , file = paste0("TopGO_BP.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") , sep = '\t', row.names = FALSE )
        write.table( allRes.mf , file = paste0("TopGO_MF.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") , sep = '\t', row.names = FALSE )
        write.table( allRes.cc , file = paste0("TopGO_CC.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") , sep = '\t', row.names = FALSE )
        
        pValue.classic.bp <- score(queryFisher.classic.ks.bp)
        pValue.elim.bp <- score(queryFisher.elim.ks.bp)[names(score(queryFisher.classic.ks.bp))]
        gstat.bp <- termStat(queryGOdata.bp, names(pValue.classic.bp))
        gSize.bp <- gstat.bp$Annotated / max(gstat.bp$Annotated) * 4
        
        png( paste0("TopGO_BP.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #"pValue.classic.BP.png" , height=1000,width=1000 )
        plot(pValue.classic.bp, pValue.elim.bp, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.bp, 
             col = redgreen(max(gstat.bp$Annotated))[ gstat.bp$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_BP_log.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "LOGpValue.classic.BP.png" , height=1000,width=1000 )
        plot(log(pValue.classic.bp), log(pValue.elim.bp), xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.bp, 
             col = redgreen(max(gstat.bp$Annotated))[ gstat.bp$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_BP_Tree_0.01.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.BP.0.01.png" , height=1000,width=1000 )
        #if(sum(allRes.bp$classicKS < 0.01) > 5){
        showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks.bp), firstSigNodes = sum(allRes.bp$classicKS < 0.01), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.01"))
        dev.off()
        #}else if(sum(allRes.bp$classicKS < 0.05) > 5 & sum(allRes.bp$classicKS < 0.05) < 35){
        
        png( paste0("TopGO_BP_Tree_0.05.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.BP.0.05.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks.bp), firstSigNodes = sum(allRes.bp$classicKS < 0.05), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.05"))
        dev.off()
        #}else{
        
        png( paste0("TopGO_BP_Tree_x10.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.BP.x10.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks.bp), firstSigNodes = 10, useInfo = 'all')
        title(paste0("Biological Processes\nCutoff count = 10"))
        dev.off() 
        #}
        
        
        
        pValue.classic.mf <- score(queryFisher.classic.ks.mf)
        pValue.elim.mf <- score(queryFisher.elim.ks.mf)[names(score(queryFisher.classic.ks.mf))]
        gstat.mf <- termStat(queryGOdata.mf, names(pValue.classic.mf))
        gSize.mf <- gstat.mf$Annotated / max(gstat.mf$Annotated) * 4
        
        png( paste0("TopGO_MF.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "pValue.classic.MF.png" , height=1000,width=1000 )
        plot(pValue.classic.mf, pValue.elim.mf, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.mf, 
             col = redgreen(max(gstat.mf$Annotated))[ gstat.mf$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_MF_log.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "LOGpValue.classic.MF.png" , height=1000,width=1000 )
        plot(log(pValue.classic.mf), log(pValue.elim.mf), xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.mf, 
             col = redgreen(max(gstat.mf$Annotated))[ gstat.mf$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_MF_Tree_0.01.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.MF.0.01.png" , height=1000,width=1000 )
        #if(sum(allRes.mf$classicKS < 0.01) > 5){
        showSigOfNodes(queryGOdata.mf, score(queryFisher.elim.ks.mf), firstSigNodes = sum(allRes.mf$classicKS < 0.01), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.01"))
        dev.off()
        #}else if(sum(allRes.mf$classicKS < 0.05) > 5 & sum(allRes.mf$classicKS < 0.05) < 35){
        
        png( paste0("TopGO_MF_Tree_0.05.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.MF.0.05.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.mf, score(queryFisher.elim.ks.mf), firstSigNodes = sum(allRes.mf$classicKS < 0.05), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.05"))
        dev.off()
        #}else{
        
        png( paste0("TopGO_MF_Tree_x10.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.MF.x10.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.mf, score(queryFisher.elim.ks.mf), firstSigNodes = 10, useInfo = 'all')
        title(paste0("Biological Processes\nCutoff count = 10"))
        #}
        dev.off() 
        
        
        
        pValue.classic.cc <- score(queryFisher.classic.ks.cc)
        pValue.elim.cc <- score(queryFisher.elim.ks.cc)[names(score(queryFisher.classic.ks.cc))]
        gstat.cc <- termStat(queryGOdata.cc, names(pValue.classic.cc))
        gSize.cc <- gstat.cc$Annotated / max(gstat.cc$Annotated) * 4
        
        png( paste0("TopGO_CC.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "pValue.classic.CC.png" , height=1000,width=1000 )
        plot(pValue.classic.cc, pValue.elim.cc, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.cc, 
             col = redgreen(max(gstat.cc$Annotated))[ gstat.cc$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_CC_log.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "LOGpValue.classic.CC.png" , height=1000,width=1000 )
        plot(log(pValue.classic.cc), log(pValue.elim.cc), xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.cc, 
             col = redgreen(max(gstat.cc$Annotated))[ gstat.cc$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_CC_log.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.CC.0.01.png" , height=1000,width=1000 )
        #if(sum(allRes.cc$classicKS < 0.01) > 5){
        showSigOfNodes(queryGOdata.cc, score(queryFisher.elim.ks.cc), firstSigNodes = sum(allRes.cc$classicKS < 0.01), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.01"))
        dev.off()
        #}else if(sum(allRes.cc$classicKS < 0.05) > 5 & sum(allRes.cc$classicKS < 0.05) < 35){
        
        png( paste0("TopGO_CC_Tree_0.05.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.CC.0.05.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.cc, score(queryFisher.elim.ks.cc), firstSigNodes = sum(allRes.cc$classicKS < 0.05), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.05"))
        dev.off()
        #}else{
        
        
        png( paste0("TopGO_CC_Tree_x10.pos.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.CC.x10.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.cc, score(queryFisher.elim.ks.cc), firstSigNodes = 10, useInfo = 'all')
        title(paste0("Biological Processes\nCutoff count = 10 (no GO terms with p < 0.05)"))
        #}
        dev.off() 
        
        ####################################################################################################################################################################################################
        # run merged queries ###############################################################################################################################################################################
        ####################################################################################################################################################################################################
        
        print("Populating")
        all.cur_og_col[,c(3:13)] <- rep(NA,dim(all.cur_og_col)[1])
        print("Positive")
        for(r in pos.cur_og_col[,2]){
          all.cur_og_col[all.cur_og_col[,2]==r,c(3:13)] <- pos.cur_og_col[pos.cur_og_col[,2]==r,c(3:13)] 
        } 
        print("Negative")
        for(r in neg.cur_og_col[,2]){
          all.cur_og_col[all.cur_og_col[,2]==r,c(3:13)] <- neg.cur_og_col[neg.cur_og_col[,2]==r,c(3:13)] 
        } 
        
        colnames(all.cur_og_col) <- c("Ver8ID","POS_ID","IsSymDef","Sym","OG_ID","P_FDR","SLOPE","BP","CC","MF","UNCLASSIFIED","FUNC_CAT","INTERPRO_DOMAINS")
        
        write.table(all.cur_og_col, paste0("og_conversion_table.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") )
        
        write.table(unique(all.cur_og_col$Sym[!is.na(all.cur_og_col$Sym)]), 
                    file = paste0("og_conversion_UNIQsymbolsWITHloc.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv"),
                    row.names = FALSE, col.names = FALSE
        )
        
        write.table(unique(all.cur_og_col$Sym[!is.na(all.cur_og_col$Sym)])[ !grepl("LOC",unique(all.cur_og_col$Sym[!is.na(all.cur_og_col$Sym)])) ], 
                    file = paste0("og_conversion_UNIQsymbols.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv"),
                    row.names = FALSE, col.names = FALSE
        )
        
        geneListX <- c(unique(all.cur_og_col$Sym[!is.na(all.cur_og_col$Sym)])[ !grepl("LOC",unique(all.cur_og_col$Sym[!is.na(all.cur_og_col$Sym)])) ],"FAKE")
        #c("AIG1","DYNLT1","EZR","GNS","IKBKE","LEMD3","MSRB3","NMBR","PEX3","RASSF3","RPL18A","RXYLT1","TBK1","WIF1","XPOT","FAKE")
        geneList <- factor(c(as.integer(rep(1,length(geneListX)-1)),0))
        names(geneList) <- geneListX
        
        
        #################### STOP ###############################################################################
        #################### STOP ###############################################################################
        #################### STOP ###############################################################################
        #################### ALWAYS STOP AND RUN THIS BY HAND UNTIL IT WORKS! WTF BIOMART! ######################
        results <- NA
        results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'go_id'), 
                         filter = "hgnc_symbol",
                         values = geneListX[1:length(geneListX)-1], #c("AIG1","DYNLT1","EZR","GNS","IKBKE","LEMD3","MSRB3","NMBR","PEX3","RASSF3","RPL18A","RXYLT1","TBK1","WIF1","XPOT"),
                         mart = bm)
        
        geneID2GO <- by(results$go_id, results$hgnc_symbol, function(x) as.character(x))
        
        queryGOdata.bp <- new( "topGOdata", 
                               description = "Simple session",
                               ontology = "BP",
                               allGenes = geneList,
                               geneSel = topDiffGenes,
                               nodeSize = 2, # <- limit to annotations held over > 10 genes
                               annot=annFUN.gene2GO,
                               gene2GO=geneID2GO
        )
        queryGOdata.mf <- new( "topGOdata", 
                               description = "Simple session",
                               ontology = "MF",
                               allGenes = geneList,
                               geneSel = topDiffGenes,
                               nodeSize = 2, # <- limit to annotations held over > 10 genes
                               annot=annFUN.gene2GO,
                               gene2GO=geneID2GO
        )
        queryGOdata.cc <- new( "topGOdata", 
                               description = "Simple session",
                               ontology = "CC",
                               allGenes = geneList,
                               geneSel = topDiffGenes,
                               nodeSize = 2, # <- limit to annotations held over > 10 genes
                               annot=annFUN.gene2GO,
                               gene2GO=geneID2GO
        )
        
        queryFisher.classic.fisher.bp <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "fisher")
        queryFisher.classic.ks.bp     <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "ks"    )
        queryFisher.elim.ks.bp     <- runTest(queryGOdata.bp, algorithm = "elim",    statistic = "ks"    )
        
        queryFisher.classic.fisher.mf <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "fisher")
        queryFisher.classic.ks.mf     <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "ks"    )
        queryFisher.elim.ks.mf     <- runTest(queryGOdata.mf, algorithm = "elim",    statistic = "ks"    )
        
        queryFisher.classic.fisher.cc <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "fisher")
        queryFisher.classic.ks.cc     <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "ks"    )
        queryFisher.elim.ks.cc     <- runTest(queryGOdata.cc, algorithm = "elim",    statistic = "ks"    )
        
        allRes.bp <- GenTable( queryGOdata.bp, classicFisher = queryFisher.classic.fisher.bp, classicKS = queryFisher.classic.ks.bp, elimKS = queryFisher.elim.ks.bp, 
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.bp)) )          
        
        allRes.mf <- GenTable( queryGOdata.mf, classicFisher = queryFisher.classic.fisher.mf, classicKS = queryFisher.classic.ks.mf, elimKS = queryFisher.elim.ks.mf, 
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.mf)) )          
        
        allRes.cc <- GenTable( queryGOdata.cc, classicFisher = queryFisher.classic.fisher.cc, classicKS = queryFisher.classic.ks.cc, elimKS = queryFisher.elim.ks.cc, 
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.cc)) )          
        
        
        #write.table(allRes.bp, file = "GO_table.BP.tsv", sep = '\t')
        #write.table(allRes.mf, file = "GO_table.MF.tsv", sep = '\t')
        #write.table(allRes.cc, file = "GO_table.CC.tsv", sep = '\t')
        write.table( allRes.bp , file = paste0("TopGO_BP.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") , sep = '\t', row.names = FALSE )
        write.table( allRes.mf , file = paste0("TopGO_MF.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") , sep = '\t', row.names = FALSE )
        write.table( allRes.cc , file = paste0("TopGO_CC.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".tsv") , sep = '\t', row.names = FALSE )
        
        pValue.classic.bp <- score(queryFisher.classic.ks.bp)
        pValue.elim.bp <- score(queryFisher.elim.ks.bp)[names(score(queryFisher.classic.ks.bp))]
        gstat.bp <- termStat(queryGOdata.bp, names(pValue.classic.bp))
        gSize.bp <- gstat.bp$Annotated / max(gstat.bp$Annotated) * 4
        
        png( paste0("TopGO_BP.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #"pValue.classic.BP.png" , height=1000,width=1000 )
        plot(pValue.classic.bp, pValue.elim.bp, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.bp, 
             col = redgreen(max(gstat.bp$Annotated))[ gstat.bp$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_BP_log.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "LOGpValue.classic.BP.png" , height=1000,width=1000 )
        plot(log(pValue.classic.bp), log(pValue.elim.bp), xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.bp, 
             col = redgreen(max(gstat.bp$Annotated))[ gstat.bp$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_BP_Tree_0.01.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.BP.0.01.png" , height=1000,width=1000 )
        #if(sum(allRes.bp$classicKS < 0.01) > 5){
        showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks.bp), firstSigNodes = sum(allRes.bp$classicKS < 0.01), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.01"))
        dev.off()
        #}else if(sum(allRes.bp$classicKS < 0.05) > 5 & sum(allRes.bp$classicKS < 0.05) < 35){
        
        png( paste0("TopGO_BP_Tree_0.05.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.BP.0.05.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks.bp), firstSigNodes = sum(allRes.bp$classicKS < 0.05), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.05"))
        dev.off()
        #}else{
        
        png( paste0("TopGO_BP_Tree_x10.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.BP.x10.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks.bp), firstSigNodes = 10, useInfo = 'all')
        title(paste0("Biological Processes\nCutoff count = 10"))
        dev.off() 
        #}
        
        
        
        pValue.classic.mf <- score(queryFisher.classic.ks.mf)
        pValue.elim.mf <- score(queryFisher.elim.ks.mf)[names(score(queryFisher.classic.ks.mf))]
        gstat.mf <- termStat(queryGOdata.mf, names(pValue.classic.mf))
        gSize.mf <- gstat.mf$Annotated / max(gstat.mf$Annotated) * 4
        
        png( paste0("TopGO_MF.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "pValue.classic.MF.png" , height=1000,width=1000 )
        plot(pValue.classic.mf, pValue.elim.mf, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.mf, 
             col = redgreen(max(gstat.mf$Annotated))[ gstat.mf$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_MF_log.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "LOGpValue.classic.MF.png" , height=1000,width=1000 )
        plot(log(pValue.classic.mf), log(pValue.elim.mf), xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.mf, 
             col = redgreen(max(gstat.mf$Annotated))[ gstat.mf$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_MF_Tree_0.01.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.MF.0.01.png" , height=1000,width=1000 )
        #if(sum(allRes.mf$classicKS < 0.01) > 5){
        showSigOfNodes(queryGOdata.mf, score(queryFisher.elim.ks.mf), firstSigNodes = sum(allRes.mf$classicKS < 0.01), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.01"))
        dev.off()
        #}else if(sum(allRes.mf$classicKS < 0.05) > 5 & sum(allRes.mf$classicKS < 0.05) < 35){
        
        png( paste0("TopGO_MF_Tree_0.05.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.MF.0.05.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.mf, score(queryFisher.elim.ks.mf), firstSigNodes = sum(allRes.mf$classicKS < 0.05), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.05"))
        dev.off()
        #}else{
        
        png( paste0("TopGO_MF_Tree_x10.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.MF.x10.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.mf, score(queryFisher.elim.ks.mf), firstSigNodes = 10, useInfo = 'all')
        title(paste0("Biological Processes\nCutoff count = 10"))
        #}
        dev.off() 
        
        
        
        pValue.classic.cc <- score(queryFisher.classic.ks.cc)
        pValue.elim.cc <- score(queryFisher.elim.ks.cc)[names(score(queryFisher.classic.ks.cc))]
        gstat.cc <- termStat(queryGOdata.cc, names(pValue.classic.cc))
        gSize.cc <- gstat.cc$Annotated / max(gstat.cc$Annotated) * 4
        
        png( paste0("TopGO_CC.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "pValue.classic.CC.png" , height=1000,width=1000 )
        plot(pValue.classic.cc, pValue.elim.cc, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.cc, 
             col = redgreen(max(gstat.cc$Annotated))[ gstat.cc$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_CC_log.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "LOGpValue.classic.CC.png" , height=1000,width=1000 )
        plot(log(pValue.classic.cc), log(pValue.elim.cc), xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.cc, 
             col = redgreen(max(gstat.cc$Annotated))[ gstat.cc$Annotated ] )#gCol)
        dev.off()
        
        png( paste0("TopGO_CC_log.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.CC.0.01.png" , height=1000,width=1000 )
        #if(sum(allRes.cc$classicKS < 0.01) > 5){
        showSigOfNodes(queryGOdata.cc, score(queryFisher.elim.ks.cc), firstSigNodes = sum(allRes.cc$classicKS < 0.01), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.01"))
        dev.off()
        #}else if(sum(allRes.cc$classicKS < 0.05) > 5 & sum(allRes.cc$classicKS < 0.05) < 35){
        
        png( paste0("TopGO_CC_Tree_0.05.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.CC.0.05.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.cc, score(queryFisher.elim.ks.cc), firstSigNodes = sum(allRes.cc$classicKS < 0.05), useInfo = 'all')
        title(paste0("Biological Processes\nCutoff p=0.05"))
        dev.off()
        #}else{
        
        
        png( paste0("TopGO_CC_Tree_x10.all.",abund,".",fun,".",mname,".",cutoff,".mod_",mod,".png"), height=1000,width=1000 )
        #png( "GO_tree.CC.x10.png" , height=1000,width=1000 )
        showSigOfNodes(queryGOdata.cc, score(queryFisher.elim.ks.cc), firstSigNodes = 10, useInfo = 'all')
        title(paste0("Biological Processes\nCutoff count = 10 (no GO terms with p < 0.05)"))
        #}
        dev.off() 

      }else{ print("Skipping complete.") }  
      #}
    }
  }
}
  

####################################################################






# Table 2:(by hand)

r <- rownames(model.p.dos)[ model.p.dos[ , "11.f_02.CAPER.p_val.FDR" ] < 0.05 ]
r <- r[!is.na(r)]
rtmp <- model.p.dos[r,c( "11.f_02.CAPER.p_val", "11.f_02.CAPER.p_val.FDR", "11.f_02.CAPER.slope", 
                 "11.f_02.LM.p_val.FDR",      "23.f_02.LM.p_val.FDR", 
                 "11.f_02.CAPER.p_val.FDR",   "23.f_02.CAPER.p_val.FDR", 
                 "11.f_02.MDMR_mash.pv.FDR",  "23.f_02.MDMR_mash.pv.FDR",
                 "11.f_02.MDMR_phylo.pv.FDR", "23.f_02.MDMR_phylo.pv.FDR", 
                 "11.f_05.LM.p_val.FDR",      "23.f_05.LM.p_val.FDR", 
                 "11.f_05.CAPER.p_val.FDR",   "23.f_05.CAPER.p_val.FDR", 
                 "11.f_05.MDMR_mash.pv.FDR",  "23.f_05.MDMR_mash.pv.FDR",
                 "11.f_05.MDMR_phylo.pv.FDR", "23.f_05.MDMR_phylo.pv.FDR"
               )
           ]
#View(rtmp)

r1 <- model.p.dos[ model.p.dos[ , "11.f_02.CAPER.slope" ] > 0 , ]
r2 <- model.p.dos[ model.p.dos[ , "11.f_02.CAPER.slope" ] < 0 , ]
r1 <- r1[ r1[ , "11.f_02.CAPER.p_val.FDR" ] < 0.05 , ]
r2 <- r2[ r2[ , "11.f_02.CAPER.p_val.FDR" ] < 0.05 , ]
r1 <- r1[order(r1[,"11.f_02.CAPER.p_val.FDR"]),]
r2 <- r2[order(r2[,"11.f_02.CAPER.p_val.FDR"]),]

rr <- c(rownames(r1),rownames(r2))
rr <- rr[ rr != "NA"]
rtmp <- model.p.dos[rr,c( "11.f_02.CAPER.slope",  "11.f_02.CAPER.p_val",       "11.f_02.CAPER.p_val.FDR", 
                          "11.f_02.LM.p_val.FDR",      "23.f_02.LM.p_val.FDR", 
                          "11.f_02.CAPER.p_val.FDR",   "23.f_02.CAPER.p_val.FDR", 
                          "11.f_02.MDMR_mash.pv.FDR",  "23.f_02.MDMR_mash.pv.FDR",
                          "11.f_02.MDMR_phylo.pv.FDR", "23.f_02.MDMR_phylo.pv.FDR", 
                          "11.f_05.LM.p_val.FDR",      "23.f_05.LM.p_val.FDR", 
                          "11.f_05.CAPER.p_val.FDR",   "23.f_05.CAPER.p_val.FDR", 
                          "11.f_05.MDMR_mash.pv.FDR",  "23.f_05.MDMR_mash.pv.FDR",
                          "11.f_05.MDMR_phylo.pv.FDR", "23.f_05.MDMR_phylo.pv.FDR"
)
]


correlation_report[  rr[1:25]  ,  c("ccc11.02.lm_est.rho","ccc11.02.caper_est.rho") ]


##############################################################################

#11b
r1 <- model.p.dos[ model.p.dos[ , "11.f_05.CAPER.slope" ] > 0 , ]
r2 <- model.p.dos[ model.p.dos[ , "11.f_05.CAPER.slope" ] < 0 , ]
r1 <- r1[ as.numeric(r1[ , "11.f_05.CAPER.p_val.FDR" ]) < 0.05 , ]
r2 <- r2[ as.numeric(r2[ , "11.f_05.CAPER.p_val.FDR" ]) < 0.025 , ]
r1 <- r1[order(r1[,"11.f_05.CAPER.p_val.FDR"]),]
r2 <- r2[order(r2[,"11.f_05.CAPER.p_val.FDR"]),]

rr <- c(rownames(r1),rownames(r2))
rr <- rr[ rr != "NA"]
rtmp <- model.p.dos[rr,c( "11.f_05.CAPER.slope",  "11.f_05.CAPER.p_val",       "11.f_05.CAPER.p_val.FDR", 
                        "11.f_05.LM.p_val.FDR",      "23.f_05.LM.p_val.FDR", 
                        "11.f_05.CAPER.p_val.FDR",   "23.f_05.CAPER.p_val.FDR", 
                        "11.f_05.MDMR_mash.pv.FDR",  "23.f_05.MDMR_mash.pv.FDR",
                        "11.f_05.MDMR_phylo.pv.FDR", "23.f_05.MDMR_phylo.pv.FDR", 
                        "11.f_02.LM.p_val.FDR",      "23.f_02.LM.p_val.FDR", 
                        "11.f_02.CAPER.p_val.FDR",   "23.f_02.CAPER.p_val.FDR", 
                        "11.f_02.MDMR_mash.pv.FDR",  "23.f_02.MDMR_mash.pv.FDR",
                        "11.f_02.MDMR_phylo.pv.FDR", "23.f_02.MDMR_phylo.pv.FDR"
)
]


correlation_report[  rr  ,  c("ccc11.02.lm_est.rho","ccc11.02.caper_est.rho") ]




qq <- c("EOG090F089Q","EOG090F0A3O","EOG090F0CAU","EOG090F06BO","EOG090F022V","EOG090F04AY","EOG090F00T9","EOG090F01TE","EOG090F07XQ","EOG090F040J","EOG090F007L","EOG090F02LH","EOG090F01W5","EOG090F039M","EOG090F052R")
correlation_report[  qq  ,  c("ccc11.02.lm_est.rho","ccc11.02.caper_est.rho") ]


##############################################################################



# grep 'EOG090F099I' out*/og_con* | awk '{print $5, $6}' | grep '"' | tr '\"' ' ' | awk ' {print $1, $2}' | sort | uniq | grep -v '^NA' | grep -v '^LOC'





