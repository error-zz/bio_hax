
library(ggplot2)
library(MASS)
library(vegan)
library(ggplot2)
library(Hmisc)
library(gplots)
library(ape)
library(caper)
setwd ("~/Desktop/exe.Avian")
load( "20190104_2pm.Rdata" )

lmlm  <- lm(log(mls$Body.Mass.Grams) ~ mls$MLS.Years)
MLSLW <- lmlm$residuals
names(MLSLW) <- compare_me$D 

lmlm  <- lm(mls$Body.Mass.Grams ~ mls$MLS.Years)
MLSW <- lmlm$residuals
names(MLSW) <- compare_me$D 

                 # dependent         # independent
f <- as.formula("MLS ~ OG.Abundance") 


for(modelid in c(#"CAPER.logwMLSresid_v_relAbs",
                 #"CAPER.wMLSresid_v_relAbs","CAPER.MLS_v_relAbs",
                 #"CAPER.logwMLSresid_v_absAbs","CAPER.wMLSresid_v_absAbs","CAPER.MLS_v_absAbs",
                 "CAPER.logwMLSresid_v_binaryAbs","CAPER.wMLSresid_v_binaryAbs","CAPER.MLS_v_binaryAbs",
                 "LM.logwMLSresid_v_relAbs",   "LM.wMLSresid_v_relAbs",   "LM.MLS_v_relAbs",
                 "LM.logwMLSresid_v_absAbs",   "LM.wMLSresid_v_absAbs",   "LM.MLS_v_absAbs",
                 "LM.logwMLSresid_v_binaryAbs",   "LM.wMLSresid_v_binaryAbs",   "LM.MLS_v_binaryAbs"
)){
  
  for(runtype in c("denovo_highConfidence","alignment")){    #,"denovo_lowConfidence"
  
    if(runtype=="alignment"){             
      shared_abundances    <- unlist(column_shared.hum)
      abundance_matrix.abs <- ag_count.hum#[ shared_abundances <= 46, ]
      abundance_matrix.rel <- ag_count.hum.rel#[ shared_abundances <= 46, ]
      shared_abundances    <- shared_abundances#[ shared_abundances <= 46 ]
      rownames( abundance_matrix.abs ) <- rownames( abundance_matrix.rel ) <- paste0( "orthoMCL_",rownames(abundance_matrix.abs) )
    } 
    if(runtype=="denovo_highConfidence"){ 
      abundance_matrix.abs <- og_counts.high_conf.abs 
      abundance_matrix.rel <- og_counts.high_conf.rel
      shared_abundances    <- shared_abundances.high_conf 
    }
    
    sum_mat <- matrix(ncol=7,nrow=length(rownames(abundance_matrix.abs)))
    colnames(sum_mat) <- c("cutoff","model","raw_pvalues","fdr_pvalues","posneg","slopes","intercepts") # ,c(phylotree$node.label))
    rownames(sum_mat) <- rownames(abundance_matrix.abs)
    
    start <- 0
    if( grepl( "binaryAbs",modelid)==TRUE ){
      start <- 4
    }
    for(cutoff in c(start:46)){
      runlist <- rownames(abundance_matrix.abs)[shared_abundances == cutoff]
      ogi <- 1
      for(og in runlist){

        
        if(ogi == 1 | ogi %% 20 == 0){  
          print(paste(runtype,modelid,": cutoff",cutoff, paste0( round((ogi/length( runlist ))*100, digits=2),"%")))
        }
        
        if(grepl("binaryAbs",modelid)==TRUE){
          compare_me     <- as.data.frame(matrix(ncol=8,nrow=50))
          if( grepl( ".MLS_",modelid)==TRUE  )      {  compare_me[,1] <- mls$MLS.Years }
          if( grepl( ".wMLSresid_",modelid)== TRUE ){  compare_me[,1] <- as.numeric(MLSW) }
          if( grepl( ".logwMLSresid_",modelid)==TRUE )      {  compare_me[,1] <- as.numeric(MLSLW) }
          #if( grepl( "relAbs",modelid)==TRUE  ){  compare_me[,2] <- as.numeric(abundance_matrix.rel[og,][abundance_matrix.rel[og,] > 0]) }
          #if( grepl( "absAbs",modelid)== TRUE ){  compare_me[,2] <- as.numeric(abundance_matrix.abs[og,][abundance_matrix.abs[og,] > 0]) }
          if( grepl( "binaryAbs",modelid)==TRUE )      {  compare_me[,2] <- as.factor(abundance_matrix.abs[og,] > 0) }
          compare_me[,3] <- as.character(mls$TreeFormat)
          compare_me[,4] <- as.character(mls$AbundanceOrder)
          compare_me[,5] <- as.character(mls$Best.Reference.Common.Name)
          compare_me[,6] <- as.character(mls$Scientific.Name)
          compare_me[,7] <- as.character(mls$Common.Name)
          compare_me[,8] <- as.character(mls$AlignLabel)
        }else{
          compare_me     <- as.data.frame(matrix(ncol=8,nrow=sum(abundance_matrix.rel[og,] > 0)))
          if( grepl( ".MLS_",modelid)==TRUE  )      {  compare_me[,1] <- mls$MLS.Years[ abundance_matrix.rel[og,] > 0 ] }
          if( grepl( ".wMLSresid_",modelid)== TRUE ){  compare_me[,1] <- as.numeric(MLSW[ abundance_matrix.rel[og,] > 0 ]) }
          if( grepl( ".logwMLSresid_",modelid)==TRUE )      {  compare_me[,1] <- as.numeric(MLSLW[ abundance_matrix.rel[og,] > 0 ]) }
          if( grepl( "relAbs",modelid)==TRUE  ){  compare_me[,2] <- as.numeric(abundance_matrix.rel[og,][abundance_matrix.rel[og,] > 0]) }
          if( grepl( "absAbs",modelid)== TRUE ){  compare_me[,2] <- as.numeric(abundance_matrix.abs[og,][abundance_matrix.abs[og,] > 0]) }
          #if( grepl( "binaryAbs",modelid)==TRUE )      {  compare_me[,2] <- as.factor(abundance_matrix.abs[og,] > 0) }
          compare_me[,3] <- as.character(mls$TreeFormat[abundance_matrix.rel[og,] > 0])
          compare_me[,4] <- as.character(mls$AbundanceOrder[abundance_matrix.rel[og,] > 0])
          compare_me[,5] <- as.character(mls$Best.Reference.Common.Name[abundance_matrix.rel[og,] > 0])
          compare_me[,6] <- as.character(mls$Scientific.Name[abundance_matrix.rel[og,] > 0])
          compare_me[,7] <- as.character(mls$Common.Name[abundance_matrix.rel[og,] > 0])
          compare_me[,8] <- as.character(mls$AlignLabel[abundance_matrix.rel[og,] > 0])
        }
        colnames(compare_me) <- c("MLS",
                                  "OG.Abundance",
                                  "TreeFormat",
                                  "AbundanceOrder",
                                  "Best.Reference.Common.Name",
                                  "Scientific.Name",
                                  "Common.Name",
                                  "AlignLabel")
        
        compdata <- NA
        tryCatch( compdata <- comparative.data(phylotree, compare_me, TreeFormat) ,  error = function(e) message(paste("ROBUST FAIL",og)))
        if(!is.na(compdata)){
          
          if(grepl( "LM.",modelid)==TRUE){
            reg   <- lm(f, data=compare_me)
          }else{
            if(grepl("binaryAbs",modelid)==TRUE){
              reg <- NA
              tryCatch( reg <- brunch(f, data=compdata, robust = 4) ,  error = function(e) message(paste("FAIL",og)))
            }else{
              reg   <- crunch(f, data=compdata)
            }
          }
          
          if(!is.na(reg)){
            if(!is.na(coef(reg))){
              if (t(coef(summary(reg)))[2,1] > 0 & is.nan(t(coef(summary(reg)))[2,1]) == 0){
                cpus.lm2 <- as.matrix( anova(reg) )
                fstat    <- summary(reg)$fstatistic
                posneg   <- as.numeric(slope > 0)
                raw_pval <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
                if(grepl("binaryAbs",modelid)==TRUE){
                  # binary run (brunch) does not return slope + intercept
                  slope    <- NA
                  int      <- NA
                } else{
                  slope    <- coef(pgls(f,data=compdata))[2]
                  int      <- coef(pgls(f,data=compdata))[1]
                }
              }else{
                raw_pval <- NA
                slope    <- NA
                posneg   <- NA
                int      <- NA
              }
            }
          }else{
            raw_pval <- NA
            slope    <- NA
            posneg   <- NA
            int      <- NA
          }
          
          sum_mat[ og , 1:7] <- c(cutoff, modelid, raw_pval, NA, posneg, slope, int)
          
          #print(paste(cutoff, modelid, raw_pval, posneg, slope, int))
          
        } 
        
        ogi <- ogi + 1
        
      }#og
    }#cutoff
    
    
    sum_mat.pop <- sum_mat[ !is.na( sum_mat[,"raw_pvalues"] ) , ]
    sum_mat.pop[,"fdr_pvalues"] <- p.adjust( as.numeric(sum_mat.pop[,3] ) , method="fdr" )
    sum_mat.pop <- sum_mat.pop[ order(sum_mat.pop[,"raw_pvalues"]) , ]
    write.table(sum_mat.pop, file=paste0("RERUN.",runtype,".",modelid,".report.tsv"), sep = "\t")
    
  }#runtype
}#modelid
 