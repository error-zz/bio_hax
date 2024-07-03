# SOFTviz.pcoa.v2.R
# updated to support:
#  prefix|suffix input
#  factor vs. numeric coloring gradients
#  
# support: jmccorri@jcvi.org

# install.packages('vegan', repos='http://cran.us.r-project.org')
# install.packages('ape', repos='http://cran.us.r-project.org')
# install.packages('gtools', repos='http://cran.us.r-project.org')
# install.packages('colorRamps', repos='http://cran.us.r-project.org')
# install.packages('fpc', repos='http://cran.us.r-project.org')
# install.packages('dplyr', repos='http://cran.us.r-project.org')
# install.packages('mixOmics', repos='http://cran.us.r-project.org')
# install.packages('png', repos='http://cran.us.r-project.org')
# install.packages('gridExtra', repos='http://cran.us.r-project.org')
library(vegan) 
library(ape)
library(gtools)
library(colorRamps)
library(fpc)
library(cluster)
library(dplyr)
# library(mixOmics)
# library(png)
library(ade4)
#library(gridExtra)
library(tidyr)
#library(dplyr)
library(ggplot2)
#library(readr)
library(broom)
library(MASS)



# library(data.table)



infant_lm_wrapper <- function(qc_data, factors, numbers, grouped_by = NULL, unlisted = TRUE){
  # Factors is a list of qc metrics that need to be treated as factors
  # numbers is a list of qc metrics that need to be treated as numberics
  
  new_qc <- data.frame(row.names = row.names(qc_data))
  name_list <- c()
  
  for(col_name in names(qc_data)){
    if(!is.null(grouped_by) & grouped_by != col_name){
      if(col_name %in% factors){
        new_qc <- cbind(new_qc, as.factor(qc_data[,col_name]))
        name_list <- c(name_list, col_name)
      }
      else if( col_name %in% numbers){
        new_qc <- cbind(new_qc, as.numeric(qc_data[,col_name]))
        name_list <- c(name_list, col_name)
      }
      else{
        if(unlisted){
          new_qc <- cbind(new_qc, qc_data[,col_name])
          name_list <- c(name_list, col_name)
        }
      }
    }
  }
  
  names(new_qc) <- name_list
  return(new_qc)
}




### inputs required:
## species

# MOTHUR EXAMPLE / LOCAL #
#setwd ("/Users/apple/Desktop/yap_pcoa_test_mothur")
#species <- read.table("pcoa_otu_tax_vs_frequency.reduced.txt", sep="\t", header=TRUE)
#prefix <- "Mothur_Example.Species.Reduced"
#### groups <- read.table("MOTHUR_EXAMPLE.manifest.tsv", sep="\t", header=TRUE)
#groups <- read.table("20161120_multi_pch_devel_test.txt", sep="\t", header=TRUE)

# BETTER MOTHUR EXAMPLE (v7)
#setwd ("/Users/apple/Desktop/yap_pcoa_test_infant/01_24_v7/mothur/")
#species <- read.table("pcoa_manifest.tsv", sep="\t", header=TRUE)
#groups <- read.table("multi_pch_test.tsv", sep="\t", header=TRUE)
#tax_summary <- read.table("otus.qfilter.EE0.15.curated_SILVA_123_SSURef_Nr99_tax_silva.wang.tax.summary", sep="\t", header=TRUE)
#prefix <- "v7.01_24"

# INFANT EXAMPLE / LOCAL #
setwd ("/Users/apple/Desktop/yap_pcoa_test_infant")
species <- read.table("pcoa_otu_tax_vs_frequency.reduced.txt", sep="\t", header=TRUE)
#groups <- read.table("infant_manifest_formatted.txt", sep="\t", header=TRUE)
#groups <- read.table("20161120_multi_pch_tests.A.txt", sep="\t", header=TRUE)
# groups <- read.table("20161122_multi_pch_tests.BtoD.txt", sep="\t", header=TRUE)
# groups <- read.table("20170126_multi_pch_tests.C3plus.txt", sep="\t", header=TRUE)

groups <- read.table("20170427.NP_clinical_manifest_A.txt", sep="\t", header=TRUE)

tax_summary <- read.table("otus.qfilter.EE0.15.curated_SILVA_123_SSURef_Nr99_tax_silva.wang.tax.summary", sep="\t", header=TRUE)
prefix <- "NP_clinical_manifest_A_testX"

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
rownames(species) <- species$rankID
species <- species[,-1]
species <- species[,-1]
species <- species[,-1]
species <- species[,-1]
species <- species[,-1]
species <- t(species)

species <- species[,-1]

species <- as.data.frame(species)
species <- data.matrix(species)
species.pop <- species[!apply(is.na(species) | species == "", 1, all),]
species_relab<-decostand(species.pop, method="total")*100 # transfrom the species data.frame to OTU/species relative abundances
species_relab.orig<-decostand(species.pop, method="total")*100 # dupe as log

# remove columns without 5 observations theres an occurrence
# species.clean <- data.frame()
ts <- 1
sc <- 1
kc <- 1
r <- colnames(species_relab)
remove_list <- list()
keepme <- list()
for(spec_col_name in r){
  tot <- length(species_relab.orig[,spec_col_name])
  ct <- sum(species_relab.orig[,spec_col_name] > 0)
  pct <- (ct/tot)*100
  print(paste(spec_col_name, "COUNT: ",ct," (",pct,"%) of TOT:",tot))
  if(pct < 0.2){
    # remove
    remove_list[[sc]] <- spec_col_name
    sc <- sc+1
  }else{
    # keep
    keepme[[kc]] <- spec_col_name
    kc <- kc+1
  }
  ts <- ts+1
}

print("REMOVING...")
for (r in remove_list){
  print(r)
}
print("END OF REMOVE LIST.")

print(paste("removing",sc,"entries."))
#species_relab <- subset(species_relab.orig, ,colnames(species_relab.orig) == keepme)
species_relab.df <-as.data.frame(species_relab.orig)
species_relab <- species_relab.df[ , names(species_relab.df) %in% keepme]
#species_relab <- species_relab.df[, colnames(data) %in% keepme]
dim(species_relab.orig)
dim(species_relab)

mice_bray<-vegdist(species_relab, method="bray") # creates a Bray-Curtis dissimilarity matrix from the species data.frame. adding “pa” , would make the matrix unweighted (presence absence). Other

# calculates a principal coordinates analyses from the mice_bray Bray-Curtis dissimilarity matrix
pcoa_mice<-pcoa(mice_bray) 

# summarizes eigenvector and cummulative variance information of the principal coordinates analyses from the mice_bray Bray-Curtis dissimilarity matrix
pcoa_mice_values <- pcoa_mice$values 


################### save.image() ###################
tmpTITLE <- paste(prefix, ".Rdata",sep="") 
save.image(file=tmpTITLE)

# 
# # DEBUG LOAD # 
# library(vegan) 
# library(ape)
# library(gtools)
# library(colorRamps)
# library(fpc)
# library(cluster)
# library(dplyr)
# library(mixOmics)
# library(png)
# library(ade4)
### setwd ("/Users/apple/Desktop/yap_pcoa_test_infant")
### prefix <- "Infant.Species"
### tmpTITLE <- paste(prefix, ".Rdata",sep="")  
### load(file=tmpTITLE)
#  
# load(file="v7.bday_x2.Rdata")
# prefix <- "NP_clinical_manifest_A"



# # DEBUG COLUMN
# col_name <- "02.C1c.Infant_Caregiver_Simple_Pass_NP.RMPRU"
# col_count <- 14
#  

###   > col_count <- 16
###   > col_name <- "X05_NPSwabForm.AntibiotTakenANDname"

# PCOA viz
col_count <- 4
for(col_name in names(groups)){
  if(!is.null(col_name)){ 
    if((col_name != "SampleID")&&(col_name != "R1")&&(col_name!="R2")){
      
      
      ######START DEBUG######
      
      ######################
      # PARSE OUT SUBGROUP #
      ######################
      
      print(paste("Calculating group [col_count : ",col_name,"], pcoa..."))
      
      num_subgroups <- length(t(unique(groups[col_count])))
      tmp <- groups[col_count][,1]
      
      print(paste(col_count,":",col_name))
      subgroup_list <- list()
      species_relab.subgroup <- data.frame()
      sgc <- 1
      nc <- 1
      rownames(groups) <- groups$SampleID
      for (n in rownames(groups)){
        print(paste("N:",n))
        #if (!(is.na(groups[nc,col_name])) && (groups[nc,col_name]!="") ){
        if (!(is.na(groups[nc,col_count])) && (groups[nc,col_count]!="") ){
          print(paste("X - ",groups[nc,col_count]," N:",n))
          subgroup_list[[sgc]] <- n
          sgc <- sgc + 1
        }else{
          print(groups[n,col_count])
        }
        nc <- nc + 1
      }
      species_relab.subgroup <- species_relab.df[rownames(species_relab.df) %in% subgroup_list , ]
      
      write.table(species_relab.subgroup, file=paste(prefix,".",col_name,".subgroup.tsv", sep = ""), quote=FALSE, sep='\t')
      
      # visualize across "all sample" pcoa
      
      print(col_name)
      
      # V7 UPDATE : FORCE TO 2-SIDED AND SKIP EMPTY
      tmp2 <- list()
      if( grepl("\\|", levels(factor(tmp))[2]) == FALSE ){ 
        print("Auto-populating NULL supergroup...")
        tq_q <- 1
        for(tq in 1:length(tmp)){
          if (is.na(tmp[[tq_q]])){
            # skip
            print(paste(tmp[[tq]],":skip"))
            tmp2[[tq_q]] <- paste(as.character(tmp[[tq_q]]),"|~",sep="")
          }else if (tmp[[tq_q]]==""){
            # skip
            print(paste(tmp[[tq]],":skip"))
            tmp2[[tq_q]] <- paste(as.character(tmp[[tq_q]]),"",sep="")
          }else{
            print(paste(tmp[[tq]],":keep"))
            tmp2[[tq_q]] <- paste(as.character(tmp[[tq_q]]),"|~",sep="")
          }
          tq_q <- tq_q+1
        }
        tmp <- factor(unlist(tmp2))
      }
      
      #if((col_name != "SampleID")&&(col_name != "R1")&&(col_name!="R2")){
      #if( grepl("\\|", levels(tmp)[2]) == FALSE ){ 
      if( grepl("\\|", levels(factor(tmp))[2]) == FALSE ){ 
        print("False")
        #####################################################################################
        ######################## ENTRIES WITHOUT "|" DELIMITER ##############################          
        #####################################################################################
        
        # all samples is.numeric check
        
        #print(paste("ROWS WITH NUMS: ",sum(  as.data.frame(gregexpr("[0-9]+", unlist(tmp)))[1,] < 0 )))
        #print(paste("ROWS WITH CHARS: ",sum(as.data.frame(gregexpr("[a-z]+", unlist(tmp)))[1,] < 0)))
        #if(sum(  as.data.frame(gregexpr("[0-9]+", unlist(tmp)))[1,] < 0  ) == length(tmp)){
        
        numeric_check_tmp <- list()
        ls <- 1
        for(l in tmp){
          if (!is.na(l)){
            numeric_check_tmp[[ls]] <- l
          }
          ls <- ls + 1
        }
        #if(sum(  as.data.frame(gregexpr("[0-9]+", unlist(tmp)))[1,] < 0  ) > 0){  
        if((sum(  as.data.frame(gregexpr("[0-9]+", unlist(numeric_check_tmp)))[1,] < 0  ) > 0)||(grep(unlist(numeric_check_tmp)[1],"~")==1)){
          # as factors
          
          # UPDATE 11/26
          
          png(filename=paste(prefix,".all_samples.",col_name,".PCoA.png", sep = ""), width = 1400, height = 1000)
          plot(pcoa_mice$vectors, bg=c(rainbow(num_subgroups)[tmp],alpha=0.8), pch=21, cex=1.3, data=row.names(pcoa_mice$vectors), main=paste("[All Samples, Factors, 1-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values[2,3], digits=2),"%)"))
          abline(v=0,lty=5) 
          abline(h=0,lty=5) 
          if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
            text( pcoa_mice$vectors, rownames( species ), pos= 3, col="black", cex=0.7 )
            text( pcoa_mice$vectors, as.character(tmp), pos=1, col="black", cex=0.7 )
          }
          legend( "topleft", border="white", bty = "n", fill=NULL, bg="white", legend=unique(tmp), pch=22, pt.bg=c(rainbow(length(unique(tmp)))[unique(tmp)]), cex=0.75 ) 
          dev.off()
          
        }else{
          # as numeric
          print("treating all-numeric factors as numeric")
          
          #all numeric, populate spectrum across full range of colors
          # find min and max
          print("Treating as Numeric")
          min <- 9999
          max <- 0
          for (w in unlist(strtoi(unlist(groups.multi.sub_group)))){
            if(!(is.na(w))){
              if (w != "_"){
                print(w)
                if(w < min){
                  min <- w
                }
                if(w > max){
                  max <- w
                }
              }
            }
          }
          
          png(filename=paste(prefix,".",col_name,".all_samples.numeric_spectrum_debug.png", sep = ""), width = 1400, height = 1000)
          # FULL # spectrum debug
          qq <- as.numeric(max)-as.numeric(min)
          col_type <- c(colorRamps::matlab.like(qq+2))
          plot(rep(1.2,qq+2),ylim=c(0.75,1.25),col=col_type,pch=19,cex=3,main=paste("[",col_name," debug] Numeric column scaled from ",qq+2," colors"))
          legend("top","full scale assigned to nodes",cex=0.7)
          col_type_full <- col_type
          
          # PLOT # spectrum debug
          wc <- 1
          col_type_subset <- list()
          for (w in unlist(strtoi(unlist(groups.multi.sub_group)))){
            if(!is.na(w)){
              print(w)
              col_type_subset[[wc]] <- col_type_full[[w+2]]
              wc <- wc + 1
            }
          }
          col_type_subset[[wc]] <- "#D3D3D3"
          col_type <- as.character(col_type_subset)
          # correct plot order
          par(new=TRUE)
          plot(rep(1,length(col_type)),ylim=c(0.75,1.25),col=col_type,pch=19,cex=3)
          legend("center",paste(length(col_type)," nodes assigned to ",qq+2, "nodes"),cex=0.7)
          #text(rep(1,length(col_type)),groups.multi.sub_group,pos=4, col="black", cex=0.7)
          
          # LEGEND # spectrum debug
          wc <- 1
          col_type_subset <- list()
          for (w in unique(sort(unlist(strtoi(unlist(groups.multi.sub_group)))))){
            #for (w in unique(mixedsort(groups.multi.sub_group))){
            print(w)
            #if (w > -1){
            #if(!(is.na(w))){
            ww <- w+2
            col_type_subset[[wc]] <- col_type_full[[ww]]
            #}else{
            
            #}
            wc <- wc + 1
          }
          col_type_subset[[wc]] <- "#D3D3D3"
          col_type_legend <- as.character(col_type_subset)
          # correct plot order
          par(new=TRUE)
          plot(rep(0.8,length(col_type_legend)),ylim=c(0.75,1.25),col=col_type_legend,pch=19,cex=3)
          legend("bottom","assigned to legend entries before NA removal",cex=0.7)
          dev.off()
          
          # row contains only single entries
          
          print(paste("Visualizing group [",col_name," (1-column, numeric)], across all samples..."))
          
          png(filename=paste(prefix,".all_samples.",col_name,".PCoA.png", sep = ""), width = 1400, height = 1000)
          plot(pcoa_mice$vectors, bg=col_type, pch=21, data=rownames(pcoa_mice$vectors), cex=1.3, main=paste("[All Samples, Numeric, 1-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values[2,3], digits=2),"%)"))
          abline(v=0,lty=5) 
          abline(h=0,lty=5) 
          if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
            text( pcoa_mice$vectors, rownames( species ), pos=3, col="black", cex=0.6 )
            text( pcoa_mice$vectors, as.character(tmp), pos=1, col="black", cex=0.7 )
          }
          legend( "topleft", border="white", bty = "n", fill=NULL, bg="white", legend=mixedsort(unique(tmp)), pch=22, pt.bg=col_type_legend, cex=0.75 ) 
          dev.off()
        }
        
        
        
      }else{
        print("True")
        
        #####################################################################################
        ################## ENTRIES WITH "|" DELIMITER (SUB-GROUPS) ##########################          
        #####################################################################################
        
        
        # row contains multiple subgroups  
        
        print(paste("Populating NAs in [",col_name," (2-column)]..."))
        
        nu_tmp <- list()
        sub_tmp <- list()
        nt <- 1
        st <- 1
        for(t in tmp){
          print(paste(t,":",length(t)))
          if(nchar(t) > 1){
            nu_tmp[[nt]] <- t
            sub_tmp[[nt]] <- t
            st <- st + 1
            print(paste("FOUND:",nu_tmp[[nt]]))
          }else{
            nu_tmp[[nt]] <- "_|-1"
            print(paste("B:",nu_tmp[[nt]]))
            print ("EMPTY")
          }
          nt <- nt + 1
        }
        tmp <- nu_tmp
        
        print(paste("Parsing subgroups for column [",col_name," (2-column)]..."))
        
        groups.multi.super_group <- as.data.frame(strsplit(as.character(tmp),'\\|'))[1,]
        groups.multi.sub_group <- as.data.frame(strsplit(as.character(tmp),'\\|'))[2,]
        pch_type <- c(1:14)[factor(unlist(groups.multi.super_group))]
        
        #if(is.numeric(strtoi(unlist(groups.multi.sub_group))) == FALSE){
        #if(sum(  as.data.frame(gregexpr("[0-9]+", unlist(groups.multi.sub_group)))[1,] < 0  ) == length(groups.multi.sub_group)){
        if(sum(grepl('[A-Za-z]', unlist(groups.multi.sub_group))) > 0){
          
          # not all numeric, treat as factors
          print("Treating as Factors")
          if (is.na(unique(mixedsort(strtoi(unlist(groups.multi.sub_group))))[1])==TRUE){
            # input as 1-column with forced NA column
            print("coercing 1-column entry colors")
            col_type <- c(colorRamps::matlab.like2(length(unique(unlist(groups.multi.super_group)))))[factor(unlist(groups.multi.super_group))]
          }else{
            # input as 2-column  
            col_type <- c(colorRamps::matlab.like2(length(unique(unlist(groups.multi.sub_group)))))[factor(unlist(groups.multi.sub_group))]
          }
          
        }else{
          
          #all numeric, populate spectrum across full range of colors
          # find min and max
          print("Treating as Numeric")
          min <- 9999
          max <- 0
          for (w in unlist(strtoi(unlist(groups.multi.sub_group)))){
            if(!(is.na(w))){
              if (w != "_"){
                print(w)
                if(w < min){
                  min <- w
                }
                if(w > max){
                  max <- w
                }
              }
            }
          }
          
          
          png(filename=paste(prefix,".",col_name,".all_samples.numeric_spectrum_debug.png", sep = ""), width = 1400, height = 1000)
          
          # FULL # spectrum debug
          qq <- as.numeric(max)-as.numeric(min)
          col_type <- c(colorRamps::matlab.like(qq+2))
          plot(rep(1.2,qq+2),ylim=c(0.75,1.25),col=col_type,pch=19,cex=3,main=paste("[col_name] Numeric column scaled from ",qq+2," colors"))
          legend("top","full scale assigned to nodes",cex=0.7)
          col_type_full <- col_type
          
          # PLOT # spectrum debug
          wc <- 1
          col_type_subset <- list()
          for (w in unlist(strtoi(unlist(groups.multi.sub_group)))){
            #for (w in groups.multi.sub_group){
            print(w)
            if (is.na(w)){
              col_type_subset[[wc]] <- "#D3D3D3"
            }else if (w > -1){
              col_type_subset[[wc]] <- col_type_full[[w+2]]
            }else{
              col_type_subset[[wc]] <- "#D3D3D3"
            }
            wc <- wc + 1
          }
          col_type <- as.character(col_type_subset)
          # correct plot order
          par(new=TRUE)
          plot(rep(1,length(col_type)),ylim=c(0.75,1.25),col=col_type,pch=19,cex=3)
          #legend("center",paste(length(col_type_legend)," nodes assigned to ",qq+2, "nodes"),cex=0.7)
          legend("center",paste(length(col_type)," nodes assigned to ",qq+2, "nodes"),cex=0.7)
          #text(rep(1,length(col_type)),groups.multi.sub_group,pos=4, col="black", cex=0.7)
          
          # LEGEND # spectrum debug
          wc <- 1
          col_type_subset <- list()
          for (w in sort(unlist(strtoi(unlist(groups.multi.sub_group))))){
            #for (w in unique(mixedsort(groups.multi.sub_group))){
            print(w)
            if (is.na(w)){
              col_type_subset[[wc]] <- "#D3D3D3"
            }else if (w > -1){
              ww <- w+2
              col_type_subset[[wc]] <- col_type_full[[ww]]
            }else{
              col_type_subset[[wc]] <- "#D3D3D3"
            }
            wc <- wc + 1
          }
          col_type_legend <- as.character(col_type_subset)
          # correct plot order
          par(new=TRUE)
          plot(rep(0.8,length(col_type_legend)),ylim=c(0.75,1.25),col=col_type_legend,pch=19,cex=3)
          legend("bottom","assigned to legend entries before NA removal",cex=0.7)
          dev.off()
        }
        
        print(paste("Visualizing group [",col_name," (2-column)], across all samples..."))
        png(filename=paste(prefix,".all_samples.",col_name,".PCoA.png", sep = ""), width = 1400, height = 1000)
        if (length(rownames(species_relab))<=100){
          plot(pcoa_mice$vectors, col=col_type, pch=pch_type, data=rownames(pcoa_mice$vectors), cex=2.0, main=paste("[All Samples, 2-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values[2,3], digits=2),"%)"))
        }else{
          plot(pcoa_mice$vectors, col=col_type, pch=pch_type, data=rownames(pcoa_mice$vectors), cex=1.3, main=paste("[All Samples, 2-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values[2,3], digits=2),"%)"))
        }
        abline(v=0,lty=5) 
        abline(h=0,lty=5) 
        if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
          text( pcoa_mice$vectors, rownames( species ), pos=3, col="black", cex=0.6 )
          text( pcoa_mice$vectors, as.character(tmp), pos=1, col="black", cex=0.7 )
        }
        if(length(unique(unlist(groups.multi.sub_group))) < 30){
          legend( "topleft", legend=mixedsort(unique(factor(unlist(groups.multi.super_group)))), pch=c(1:14)[mixedsort(unique(factor(unlist(groups.multi.super_group))))], col=1, cex=1.0 ) 
          #legend( "topright", legend=mixedsort(levels(factor(unlist(groups.multi.sub_group)))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(unlist(groups.multi.sub_group))))[factor(unique(mixedsort(unlist(tmp))))]), cex=1.0 ) 
          
          #if(sum(  as.data.frame(gregexpr("[0-9]+", unlist(groups.multi.sub_group)))[1,] < 0  ) == length(groups.multi.sub_group)){
          if(sum(grepl('[A-Za-z]', unlist(groups.multi.sub_group))) > 0){
            
              # not all numeric, treat as factors
            print("Legend as factors")
            #legend( "topright", legend=unique(mixedsort(unlist(groups.multi.sub_group))), pch=22, pt.bg=col_type[factor(unique(mixedsort(unlist(groups.multi.sub_group))))])
            legend( "topright", legend=unique(mixedsort(unlist(groups.multi.sub_group))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
          }else{
            print("Legend as numeric")
            groups.multi.sub_group.int <- unlist(strtoi(unlist(groups.multi.sub_group)))
            groups.multi.sub2_group <- mixedsort(groups.multi.sub_group.int)
            #legend( "topright", legend=mixedsort(strtoi(unlist(groups.multi.sub_group))), pch=22, pt.bg=col_type_legend)
            legend( "topright", legend=unique(mixedsort(strtoi(unlist(groups.multi.sub_group)))), pch=22, pt.bg=unique(col_type_legend))
          }
          
        }else{
          
          #if(sum(  as.data.frame(gregexpr("[0-9]+", unlist(groups.multi.sub_group)))[1,] < 0  ) == length(groups.multi.sub_group)){
          if(sum(grepl('[A-Za-z]', unlist(groups.multi.sub_group))) > 0){
            
            # not all numeric, treat as factors
            print("Legend as factors")
            legend( "topright", legend=unique(mixedsort(unlist(groups.multi.sub_group))), pch=22, cex=0.7, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
          }else{
            print("Legend as numeric")
            groups.multi.sub_group.int <- unlist(strtoi(unlist(groups.multi.sub_group)))
            groups.multi.sub2_group <- mixedsort(groups.multi.sub_group.int)
            legend( "topright", legend=unique(mixedsort(strtoi(unlist(groups.multi.sub_group)))), pch=22, cex=0.7, pt.bg=unique(col_type_legend))
          }            
          
        }
        dev.off()
        
        
        
        ########################################################## SWITCH TO DEFINED SAMP<LES ONLY
        
        
        print(paste("Visualizing group [",col_name," (2-column)], across ONLY THE DEFINED samples..."))
        
        # select reduced tmp
        tmp <- sub_tmp
        
        # select reduced species_relab.subgroup
        
        print(paste(col_count,":",col_name))
        subgroup_list <- list()
        species_relab.subgroup <- data.frame()
        sgc <- 1
        nc <- 1
        rownames(groups) <- groups$SampleID
        for (n in rownames(groups)){
          print(paste("N:",n))
          #if (!(is.na(groups[nc,col_name])) && (groups[nc,col_name]!="") ){
          if (!(is.na(groups[nc,col_count])) && (groups[nc,col_count]!="") ){
            print(paste("X - ",groups[nc,col_count]," N:",n))
            subgroup_list[[sgc]] <- n
            sgc <- sgc + 1
          }else{
            print(groups[n,col_count])
          }
          nc <- nc + 1
        }
        species_relab.subgroup <- species_relab.df[rownames(species_relab.df) %in% subgroup_list , ]
        
        # within subgroup colnames, etc.
        
        groups.multi.super_group <- as.data.frame(strsplit(as.character(tmp),'\\|'))[1,]
        # ??? # if (is.na(unique(mixedsort(strtoi(unlist(groups.multi.sub_group))))[1])==TRUE){
        # ??? #   groups.multi.sub_group <- c(rep("~",dim(groups.multi.super_group)[2]))
        # ??? #   groups.multi.sub_group <- as.list(groups.multi.sub_group)
        # ??? #   names(groups.multi.sub_group) <- names(as.list(groups.multi.super_group))
        # ??? # }else{
          groups.multi.sub_group <- as.data.frame(strsplit(as.character(tmp),'\\|'))[2,]
        # ??? # }
        
        p_t <- 1
        pch_list <- list()
        for(q in factor(unlist(groups.multi.super_group))){
          if (!(q == "NULL")){
            # nothing
            pch_list[[p_t]] <- q
            p_t <- p_t + 1
          }
        }
        #pch_type <- c(1:14)[factor(unlist(groups.multi.super_group))]
        pch_type <- c(1:14)[factor(unlist(pch_list))]
        
        #if(is.numeric(strtoi(unlist(groups.multi.sub_group))) == FALSE){
        #if(sum(  as.data.frame(gregexpr("[0-9]+", unlist(groups.multi.sub_group)))[1,] < 0  ) == length(groups.multi.sub_group)){
        if(sum(grepl('[A-Za-z]', unlist(groups.multi.sub_group))) > 0){
          
          # not all numeric, treat as factors
          print("Treating as Factors")
          if (is.na(unique(mixedsort(strtoi(unlist(groups.multi.sub_group))))[1])==TRUE){
            # input as 1-column with forced NA column
            print("coercing 1-column entry colors")
            col_type <- c(colorRamps::matlab.like2(length(unique(unlist(groups.multi.super_group)))))[factor(unlist(groups.multi.super_group))]
          }else{
            # input as 2-column  
            col_type <- c(colorRamps::matlab.like2(length(unique(unlist(groups.multi.sub_group)))))[factor(unlist(groups.multi.sub_group))]
          }
          
        }else{
          
          
          png(filename=paste(prefix,".",col_name,".within_subgroup.numeric_spectrum_debug.png", sep = ""), width = 1400, height = 1000)
          # PLOT # spectrum debug
          wc <- 1
          col_type_subset <- list()
          for (w in unlist(strtoi(unlist(groups.multi.sub_group)))){
            #for (w in groups.multi.sub_group){
            if(!(is.na(w))){
              print(w)
              col_type_subset[[wc]] <- col_type_full[[w+2]]
              wc <- wc + 1
            }
          }
          col_type <- as.character(col_type_subset)
          # correct plot order
          plot(rep(0.8,length(col_type)),ylim=c(0.75,1.25),col=col_type,pch=19,cex=3,main=paste(length(col_type)," colors"))
          #text(rep(1,length(col_type)),groups.multi.sub_group,pos=4, col="black", cex=0.7)
          
          # LEGEND # spectrum debug
          wc <- 1
          col_type_subset <- list()
          for (w in sort(unlist(strtoi(unlist(groups.multi.sub_group))))){
            #for (w in unique(mixedsort(groups.multi.sub_group))){
            print(w)
            ww <- w+2
            col_type_subset[[wc]] <- col_type_full[[ww]]
            wc <- wc + 1
          }
          col_type_legend <- as.character(col_type_subset)
          # correct plot order
          plot(rep(1.2,length(col_type_legend)),ylim=c(0.75,1.25),col=col_type_legend,pch=19,cex=3,main=paste(length(col_type_legend)," colors"))
          dev.off()
        }
        
        print(paste("Calculating group [",col_name,"], pcoa..."))
        
        mice_bray.subgroup<-vegdist(species_relab.subgroup, method="bray") # creates a Bray-Curtis dissimilarity matrix from the species data.frame. adding “pa” , would make the matrix unweighted (presence absence). Other
        pcoa_mice.subgroup<-pcoa(mice_bray.subgroup)
        pcoa_mice_values.subgroup <- pcoa_mice.subgroup$values
        
        # update text entries
        topleft_text <- list() 
        mt <- 1
        #for (m in mixedsort(unique(factor(unlist(groups.multi.super_group))))){
        for(m in unique(unlist(groups.multi.super_group))){
          if (!(m == "NULL")){
            #if(is.null(m) == FALSE){
            topleft_text[[mt]] <- m
            mt <- mt + 1
          }
        }
        topleft_text <- unique(topleft_text)
        # plot
        
        png(filename=paste(prefix,".within_subgroup.",col_name,".PCoA.png", sep = ""), width = 1400, height = 1000)
        if (length(rownames(species_relab))<=100){
          plot(pcoa_mice.subgroup$vectors, col=col_type, pch=pch_type, data=rownames(pcoa_mice.subgroup$vectors), cex=2.0, main=paste("[Within Subgroup, 2-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[2,3], digits=2),"%)"))
        }else{
          plot(pcoa_mice.subgroup$vectors, col=col_type, pch=pch_type, data=rownames(pcoa_mice.subgroup$vectors), cex=1.3, main=paste("[Within Subgroup, 2-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[2,3], digits=2),"%)"))
        }
        abline(v=0,lty=5) 
        abline(h=0,lty=5) 
        if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
          text( pcoa_mice.subgroup$vectors, row.names(species_relab.subgroup), pos=3, col="black", cex=0.6 )
          text( pcoa_mice.subgroup$vectors, as.character(unlist(tmp)), pos=1, col="black", cex=0.7 )
        }
        #if(length(unique(unlist(groups.multi.sub_group))) < 30){
        #legend( "topleft", legend=mixedsort(unique(factor(unlist(groups.multi.super_group)))), pch=c(1:14)[mixedsort(unique(factor(unlist(groups.multi.super_group))))], col=1, cex=1.0 ) 
        #legend( "topleft", legend=topleft_text, pch=c(1:14)[mixedsort(unique(factor(unlist(groups.multi.super_group))))], col=1, cex=1.0 ) 
        legend( "topleft", legend=topleft_text, pch=unique(pch_type), col=1, cex=1.0 ) 
        #if(sum(  as.data.frame(gregexpr("[0-9]+", unlist(groups.multi.sub_group)))[1,] < 0  ) == length(groups.multi.sub_group)){
        if(sum(grepl('[A-Za-z]', unlist(groups.multi.sub_group))) > 0){
          
          # not all numeric, treat as factors
          print("Legend as factors")
          #legend( "topright", legend=unique(mixedsort(unlist(groups.multi.sub_group))), pch=22, pt.bg=col_type[factor(unique(mixedsort(unlist(groups.multi.sub_group))))])
          legend( "topright", legend=unique(mixedsort(unlist(groups.multi.sub_group))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
        }else{
          print("Legend as numeric")
          
          # update legend text to remove empty
          topright_text <- list()
          ts <- 1
          for(a in unique(mixedsort(strtoi(unlist(groups.multi.sub_group))))){
            if (!(is.na(a))){
              topright_text[[ts]] <- a
              ts <- ts + 1
            }
          }
          
          groups.multi.sub_group.int <- unlist(strtoi(unlist(groups.multi.sub_group)))
          groups.multi.sub2_group <- mixedsort(groups.multi.sub_group.int)
          #legend( "topright", legend=mixedsort(strtoi(unlist(groups.multi.sub_group))), pch=22, pt.bg=col_type_legend)
          #legend( "topright", legend=unique(mixedsort(strtoi(unlist(groups.multi.sub_group)))), pch=22, pt.bg=unique(col_type_legend))
          legend( "topright", legend=topright_text, pch=22, pt.bg=unique(col_type_legend))
        }
        #}else{
        #  legend( "topleft", legend=mixedsort(unique(factor(unlist(groups.multi.super_group)))), pch=c(1:14)[mixedsort(unique(factor(unlist(groups.multi.super_group))))], col=1, cex=0.7 ) 
        ##  legend( "topright", legend= mixedsort(levels(factor(unlist(groups.multi.sub_group)))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(unlist(groups.multi.sub_group))))[factor(unique(mixedsort(unlist(tmp))))]), cex=0.7 ) 
        #}
        dev.off()
        
      }
      #}else{
      # skip columns that delineate uparse run methods
      #}
      
      
      print(paste("Visualizing group [",col_name," (2-column)], across all samples, WITHIN CLUSTERS..."))
      
      
      mykde2d = function (dfxy, xax = 1, yax = 2, pch = 20, cpoint = 1, neig = NULL, 
                          cneig = 2, xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, 
                          cgrid = 1, include.origin = TRUE, origin = c(0, 0), sub = "", 
                          csub = 1.25, possub = "bottomleft", pixmap = NULL, contour = NULL, 
                          area = NULL, add.plot = FALSE, ... ) 
      {
        if (!require(MASS)) 
          stop("library MASS required for kde2d")
        par(mar = c(0.1, 0.1, 0.1, 0.1))
        #s.label(dfxy, xax = xax, yax = yax, clab = 0, pch = pch_type, 
        #        cpoint = cpoint, neig = neig, cneig = cneig, xlim = xlim, 
        #        ylim = ylim, grid = grid, addaxes = addaxes, cgrid = cgrid, 
        #        include.origin = include.origin, origin = origin, sub = sub, 
        #        csub = csub, possub = possub, pixmap = pixmap, contour = contour, 
        #        area = area, add.plot = add.plot)
        if(length(rownames(pcoa_mice.subgroup$vectors))<=50){
          plot(pcoa_mice.subgroup$vectors, col=col_type, pch=pch_type, data=rownames(pcoa_mice.subgroup$vectors), cex=3.0, main=paste("[Within Subgroup, 2-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[2,3], digits=2),"%)"))
        }else{
          plot(pcoa_mice.subgroup$vectors, col=col_type, pch=pch_type, data=rownames(pcoa_mice.subgroup$vectors), cex=0.7, main=paste("[Within Subgroup, 2-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[2,3], digits=2),"%)"))
        }
        x <- as.numeric(pcoa_mice.subgroup$vectors[, 1])
        y <- as.numeric(pcoa_mice.subgroup$vectors[, 2])
        xykde = MASS::kde2d(x, y, lims = par("usr"))
        zlim = range(xykde$z, finite = TRUE)
        lev = seq(zlim[1], zlim[2], le = 8)
        lev = lev[2:7]
        contour(xykde, add = TRUE, levels = lev, drawlabels = FALSE, ...)
        invisible(match.call())
      }
      
      
      
      #colors = factor(unique(col_type_legend))
      colors = unique(col_type)
      ids = rownames(pcoa_mice.subgroup)
      
      png(filename=paste(prefix,".all_samples.",col_name,".ClusterDensity.png", sep = ""), width = 1800, height = 1800)
      par(mfrow = c(1,1))
      mykde2d(pcoa_mice.subgroup$vectors, col="red", lwd=1.4, lty=1,main=NULL)
      if((length(rownames(species_relab))>100)==FALSE){
        # s.class(pcoa_mice.subgroup$vectors, factor(groups[,col_count]), cpoint = 1, csub=1, col=colors, add.p=TRUE, label="")
      }
      #s.corcircle(as.numeric(pcoa_mice.subgroup$vectors[,1]), full = TRUE, box = TRUE, csub=1, possub="bottom")
      legend( "topleft", legend=topleft_text, pch=unique(pch_type), col=1, cex=1.0 ) 
      #legend( "topright", legend=mixedsort(as.numeric(unique(unlist(groups.multi.sub_group)))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group))))))) 
      legend( "topright", legend=mixedsort(unique(unlist(groups.multi.sub_group))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group))))))) 
      dev.off()
      
      ####################################### IDEAL SILHOUETTE SELECTION ROUTINE ########################################
      
      ## Use the silhouette widths for assessing the best number of clusters,
      ## following a one-dimensional example from Christian Hennig :
      ##
      #x <- c(rnorm(50), rnorm(50,mean=5), rnorm(30,mean=15))
      
      max_clust <- 10
      if(max_clust > (sum(tmp!="NULL")-2)){
        max_clust <- sum(tmp!="NULL")-2
      }
      asw <- numeric(max_clust)
      c_asw <- numeric(max_clust)
      ## Note that "k=1" won't work!pk
      
      png(filename=paste(prefix,".all_samples.",col_name,".cluster_count_PCoA.png", sep = ""), width = 800, height = 2400)
      par(mfrow = c(max_clust-1,3))
      for (num_clust in 2:max_clust){
        pk <- pam(as.data.frame(pcoa_mice.subgroup$vectors[,1:2]),num_clust)
        asw[num_clust] <- pk $ silinfo $ avg.width
        c_asw[num_clust] <- c(pk $ silinfo $clus.avg.widths)
        clusplot(pk,main=paste("\nClusters (PAM,PC1+2), N=",num_clust),col.p="grey",color=T,metric="euclidean",shade=T,labels=0,cex=0.05) #`pam` does you partitioning
        barplot(pk$clusinfo[,1],col="darkgrey",names.arg=c(1:dim(pk$clusinfo)[1]),main="Cluster Size")
        barplot(pk$silinfo$clus.avg.widths,col="darkgrey",names.arg=c(1:dim(pk$clusinfo)[1]),main="Cluster Avg Width")
      }
      dev.off()
      
      
      png(filename=paste(prefix,".all_samples.",col_name,".cluster_count_histo.png", sep = ""), width = 1200, height = 800)
      par(mfrow = c(1,1))
      k.best <- which.max(asw)
      barplot(asw, type= "h", main = paste("Avg.Sil.Width\noptimal number of clusters:", k.best), xlab= "k  (# clusters)", ylab = "average silhouette width")
      axis(1, k.best, paste("best",k.best,sep="\n"), col = "red", col.axis = "red")
      dev.off()
      
      # tmpTITLE <- paste(prefix, ".debug.Rdata",sep="") 
      # save.image(file=tmpTITLE)
      
      # initiate cluster type matrix
      group_classification <- groups[col_count][,1]
      rownames(group_classification) <- rownames(groups[col_count][,1])
      
      num_clust <- k.best
      super_cluster_type <- matrix(nrow=num_clust,ncol=length(topleft_text))
      colnames(super_cluster_type) <- topleft_text
      rownames(super_cluster_type) <- c(1:num_clust)
      
      
      
      #if (is.na(unique(mixedsort(strtoi(unlist(groups.multi.sub_group))))[1])!=TRUE){
      #if (is.na(sort(unique(mixedsort(strtoi(unlist(groups.multi.sub_group))))))==TRUE){
      
      topright_text <-  unique(mixedsort(unlist(groups.multi.sub_group)))
      sub_cluster_type <- matrix(nrow=num_clust,ncol=length(topright_text))
      print(paste("TR:",topright_text))
      colnames(sub_cluster_type) <- topright_text
      rownames(sub_cluster_type) <- c(1:num_clust)
      super_cluster_type <- as.data.frame(super_cluster_type)
      sub_cluster_type <- as.data.frame(sub_cluster_type)
      
      sample_super_cluster_type <- matrix(nrow=length(rownames(species_relab.subgroup)),ncol=length(topleft_text))
      colnames(sample_super_cluster_type) <- topleft_text
      rownames(sample_super_cluster_type) <- rownames(species_relab.subgroup)
      
      sample_sub_cluster_type <- matrix(nrow=length(rownames(species_relab.subgroup)),ncol=length(topright_text))
      colnames(sample_sub_cluster_type) <- topright_text
      rownames(sample_sub_cluster_type) <- rownames(species_relab.subgroup)
      
      factor_sample_super_cluster_type <- matrix(nrow=length(rownames(species_relab.subgroup)),ncol=1)
      factor_sample_super_cluster_type <- as.data.frame(factor_sample_sub_cluster_type)
      factor_sample_sub_cluster_type <- matrix(nrow=length(rownames(species_relab.subgroup)),ncol=1)
      factor_sample_sub_cluster_type <- as.data.frame(factor_sample_sub_cluster_type)
      
      colnames(factor_sample_super_cluster_type) <- col_name
      rownames(factor_sample_super_cluster_type) <- rownames(species_relab.subgroup)
      colnames(factor_sample_sub_cluster_type) <- col_name
      rownames(factor_sample_sub_cluster_type) <- rownames(species_relab.subgroup)
      
      #}else{
      #  sub_cluster_type <- super_cluster_type
      #}
      
      for (tt in rownames(super_cluster_type)){
        for (rr in colnames(super_cluster_type)){
          super_cluster_type[tt,rr] <- 0
        }
      }
      for (tt in rownames(sub_cluster_type)){
        for (rr in colnames(sub_cluster_type)){
          sub_cluster_type[tt,rr] <- 0
        }
      }
      
      for (tt in rownames(sample_super_cluster_type)){
        for (rr in colnames(sample_super_cluster_type)){
          sample_super_cluster_type[tt,rr] <- 0
        }
      }
      for (tt in rownames(sample_sub_cluster_type)){
        for (rr in colnames(sample_sub_cluster_type)){
          sample_sub_cluster_type[tt,rr] <- 0
        }
      }
      
      for (tt in rownames(factor_sample_super_cluster_type)){
          factor_sample_super_cluster_type[tt,col_name] <- 0
      }
      for (tt in rownames(factor_sample_sub_cluster_type)){
          factor_sample_sub_cluster_type[tt,col_name] <- 0
      }
      
      # make sure pk matches k.best
      pk <- pam(as.data.frame(pcoa_mice.subgroup$vectors[,1:2]),k.best)
      # assign cluster type and super/sub groups
      for (r in rownames(groups)){
        #print(r)
        if(r %in% names(pk$clustering)){
          #print(r)
          rc <- 1
          for (m in groups[col_count][,1]){
            #print(r)
            if(!is.na(m)){
              if ( (m =="")==FALSE ){
                if (rownames(groups)[rc] == r){
                  
                  super_type <- as.character(as.data.frame(strsplit(as.character(m),'\\|'))[1,])
                  sub_type <- as.character(as.data.frame(strsplit(as.character(m),'\\|'))[2,])
                  # ??? # if (is.na(unique(mixedsort(strtoi(unlist(groups.multi.sub_group))))[1])==TRUE){
                  # ??? #   sub_type <- super_type
                  # ??? #   sub_cluster_type <- super_cluster_type
                  # ??? # }
                  
                  cluster_id <- pk$clustering[[r]]
                  print(paste(r,super_type,cluster_id,super_cluster_type[cluster_id,super_type]," to ",super_cluster_type[cluster_id,super_type]))
                  print(paste(r,super_type,cluster_id,sub_cluster_type[cluster_id,sub_type]," to ",sub_cluster_type[cluster_id,sub_type]))
                  super_cluster_type[cluster_id,super_type] <- super_cluster_type[cluster_id,super_type] + 1
                  sub_cluster_type[cluster_id,sub_type] <- sub_cluster_type[cluster_id,sub_type] + 1
                  
                  sample_super_cluster_type[r,super_type] <- sample_super_cluster_type[r,super_type] + 1
                  sample_sub_cluster_type[r,sub_type] <- sample_sub_cluster_type[r,sub_type] + 1
                  
                  factor_sample_super_cluster_type[r,col_name] <- super_type
                  factor_sample_sub_cluster_type[r,col_name] <- sub_type
                  
                  
                }
              }
            }
            rc <- rc + 1
          }
        }
      }
      
      # super_cluster_type
      write.table(super_cluster_type, file=paste(prefix,".",col_name,".clusterassigned_supergroup.tsv", sep = ""), quote=FALSE, sep='\t')
      write.table(sample_super_cluster_type, file=paste(prefix,".",col_name,".clusterassigned_factoredsubset_supergroup.tsv", sep = ""), quote=FALSE, sep='\t')
      
      # sub_cluster_type
      write.table(sub_cluster_type, file=paste(prefix,".",col_name,".clusterassigned_subgroup.tsv", sep = ""), quote=FALSE, sep='\t')
      write.table(sample_sub_cluster_type, file=paste(prefix,".",col_name,".clusterassigned_factoredsubset_subgroup.tsv", sep = ""), quote=FALSE, sep='\t')
      
      # initialize matrix with species-level taxa
      col_x8 <- list()
      ty <- 1
      for (m in colnames(species)){
        col_x8[[ty]] <- paste(unlist(strsplit(m, '[.]'))[1:8],collapse=".")
        ty <- ty+1
      }
      species_type <- matrix(nrow=k.best,ncol=length(col_x8))
      colnames(species_type) <- col_x8
      rownames(species_type) <- c(1:k.best)
      for (tt in rownames(species_type)){
        for (rr in colnames(species_type)){
          species_type[tt,rr] <- 0
        }
      }
      
      #######################################
      # determine most common taxa by species
      species_sum <- list()
      sc <- 1
      for (r in rownames(species)){ # sample id w species association
        sum <- 0
        if(r %in% names(pk$clustering)){ # sample id in clustering
          # sample with species association has cluster
          rc <- 1
          cluster_id <- pk$clustering[[r]]
          for (m in colnames(species)){
            species_id <- paste(unlist(strsplit(colnames(species)[rc], '[.]'))[1:8],collapse=".")
            species_type[cluster_id,species_id] <- species_type[cluster_id,species_id] + species[r,m]
            sum <- sum + species[r,m]
            rc <- rc + 1
          }
        }
        species_sum[[sc]] <- sum
        sc <- sc + 1
      }
      names(species_sum) <- rownames(species)
      
      # SPECIES
      write.table(species_type, file=paste(prefix,".",col_name,".allSPECIESassigned_subgroup.tsv", sep = ""), quote=FALSE, sep='\t')
      
      # get the top 10 for each cluster
      top <- list()      
      for (k in 1:k.best){
        top <- append(top, names(rev(sort(species_type[k,]))[1:10]))
      }
      # and the unique set amongst all of those top 10s
      top <- unique(top)
      
      species_type_topSpecies <- matrix(nrow=k.best,ncol=length(top))
      colnames(species_type_topSpecies) <- top
      rownames(species_type_topSpecies) <- c(1:k.best)
      for (tt in rownames(species_type_topSpecies)){
        for (rr in colnames(species_type_topSpecies)){
          species_type_topSpecies[tt,rr] <- species_type[tt,rr]
        }
      }
      
      # SPECIES
      write.table(species_type_topSpecies, file=paste(prefix,".",col_name,".topSPECIESassigned_subgroup.tsv", sep = ""), quote=FALSE, sep='\t')
      
      
      ############################################################################################################
      ############################################################################################################
      ################# REBUILD species_relab.subgroup3 STATISTICAL ANALYSIS ################################################
      ############################################################################################################
      ############################################################################################################
      
      
      
      
      ############################################################################################################
      ############################################################################################################
      ################# SPECIES v SUPERGROUP STATISTICAL ANALYSIS ################################################
      ############################################################################################################
      ############################################################################################################
      
      list_of_lms <- list()
      reg <-data.frame()
      
      #lm.factor_sample_super_cluster_type <- infant_lm_wrapper(factor_sample_super_cluster_type, factors = colnames(factor_sample_super_cluster_type), numbers = c(), grouped_by = 'none')
      
      N <- ncol(species_relab.subgroup)
      report.table <- data.frame(matrix(ncol=8,nrow=length(rownames(factor_sample_super_cluster_type))))
      rownames(report.table) <- rownames(factor_sample_super_cluster_type)
      covar.pvals <- data.frame(matrix(ncol=length(levels(factor(factor_sample_super_cluster_type)))))
      covar.model_pvals <- data.frame(ncol=1)
      #covar.pvals.i <- data.frame()
      
      plot.new()
      #par(mfrow = c(4,4))
      for(i in 1:N){
        
        print(i)
        reg <- lm( as.numeric(species_relab.subgroup[,i]) ~ as.factor(factor_sample_super_cluster_type[,col_name]), na.action=na.omit )
        
        png(filename=paste(prefix,".subset.",col_name,".species_vs_supergroup_factor_residuals.png", sep = ""), width = 500, height = 500)
        plot.new()
        ggplot(data=reg,aes(x=as.factor(factor_sample_super_cluster_type[,col_name]),y=species_relab.subgroup[,i]))+geom_point()
        dev.off()
        
        if (t(coef(summary(reg)))[2,1] > 0 & is.nan(t(coef(summary(reg)))[2,1]) == 0){
          cpus.lm2 <- stepAIC(reg, trace = FALSE)
        
          fstat <- summary(reg)$fstatistic
          covar.model_pvals[i] <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
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
          
          covar.model_pvals[i] <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
          covar.pvals[i,] <- matrix(0,1,length(levels(factor(factor_sample_super_cluster_type))))
          
        }
 
      }
      
      ###### SIGNIFICANT TAXA, P-VALUE TABLE ########
      
      png(filename=paste(prefix,".subset.",col_name,".significant_species_by_factor.png", sep = ""), width = 500, height = 500)
      plot.new()
      color_q <- c(colorRamps::matlab.like2(length(colnames(covar.pvals))))
      legend_q <- c()
      plot(sort(as.numeric(covar.model_pvals)),col="black",pch=15,cex=1.0,ylim=c(0,1),ylab="P-value",xlab="Taxa Index (Independently Sorted by Column)")
      par(new=TRUE)
      for(q in 1:length(colnames(covar.pvals))){
        plot(sort(covar.pvals[,q]),col=color_q[q],pch=22,cex=0.5,ylim=c(0,1),main=paste(col_name,"\nFactor : Signif. Taxa Count"),ylab="P-value",xlab="Taxa Index (Independently Sorted by Column)")
        legend_q[q] <- paste(unlist(as.list(unique(factor_sample_super_cluster_type)))[q]," : ",sum(sort(covar.pvals[,q])<0.05)," (",round(sum(sort(covar.pvals[,q])<0.05)/length(factor_sample_super_cluster_type),digits=2)*100,"%)")
        par(new=TRUE)
      }
      color_q[q+1] <- "black"
      legend_q[q+1] <- paste("Model : ",sum(sort(as.numeric(covar.model_pvals))<0.05)," (",round(sum(sort(as.numeric(covar.model_pvals))<0.05)/length(covar.model_pvals),digits=2)*100,"%)")
      abline(h=0.05,lty=4, col="red")
      legend("topleft",legend=legend_q,pch=22,pt.bg=color_q,cex=0.7,bty='n')
      #colnames(covar.model_pvals) <- colnames(species_relab.subgroup)
      colnames(covar.model_pvals) <- colnames(sub_cluster_type)
      species_signif_pvals <- covar.model_pvals[,covar.model_pvals<0.05]
      colnames(species_signif_pvals) <- colnames(covar.model_pvals)[covar.model_pvals<0.05]
      species_signif_pvals <- sort(species_signif_pvals)
      if(length(species_signif_pvals)>=5){
        legend("top",legend=paste(colnames(species_signif_pvals)[1:5]," p=",sort(species_signif_pvals)[1:5]),bty='n',cex=0.7)
      }else{
        legend("top",legend=paste(colnames(species_signif_pvals)[1:length(species_signif_pvals)]," p=",sort(species_signif_pvals)[1:length(species_signif_pvals)]),bty='n',cex=0.7)
      }
      dev.off()
      
      write.table(t(sort(covar.model_pvals)), file=paste(prefix,".",col_name,".signifSPECIESassigned_subgroup.tsv", sep = ""), quote=FALSE, sep='\t')
      rownames(covar.pvals) <- colnames(covar.model_pvals)
      for(q in 1:length(colnames(covar.pvals))){
        tmp <- as.data.frame(covar.pvals[,q])
        #rownames(tmp) <- rownames(covar.pvals)
        rownames(tmp) <- colnames(species_relab.subgroup)
        tmp <- tmp[order(tmp),,drop=FALSE]
        write.table(tmp, file=paste(prefix,".",col_name,".FactorEquals_",unlist(as.list(unique(factor_sample_super_cluster_type)))[q],".signifSPECIESassigned_subgroup.tsv", sep = ""), quote=FALSE, sep='\t')
      }
      SPECIESvSUPERGROUP.covar.pvals <- covar.model_pvals
      #colnames(SPECIESvSUPERGROUP.covar.pvals) <- colnames(covar.model_pvals)
      #SPECIESvSUPERGROUP.species_type_topSpecies <- species_type[,as.factor(colnames(sort(as.numeric(SPECIESvSUPERGROUP.covar.pvals)))[1:10])]
      
      SPECIESvSUPERGROUP.species_type_topSpecies <- species_type[,as.factor(colnames(sort(SPECIESvSUPERGROUP.covar.pvals)[1:10]))]
      
      
      ##### barplot
      
      tax_sum.nu <- tax_summary[3]
      rownames(tax_sum.nu) <- as.factor(unlist(tax_summary[2]))
      
      nu_colnames <- c()
      c <- 1
      for (spec_id in as.factor(colnames(SPECIESvSUPERGROUP.species_type_topSpecies))){
        old_spec_id <- unlist(strsplit(spec_id, '[.]'))
        nu_spec_id <- list()
        nu_spec_id[[1]] <- 0 # set root
        zt <- 2
        for(zzz in 2:length(old_spec_id)){
          if(old_spec_id[[zzz]] != 0){
            nu_spec_id[[zt]] <- old_spec_id[[zzz]]
            zt <- zt + 1
          }
        }
        query_id <- paste(as.list(nu_spec_id),collapse=".")
        print(paste(spec_id," to ",query_id))
        nu_colnames[c] <- query_id
        c <- c + 1
      }
      
      png(filename=paste(prefix,".all_samples.",col_name,".SPECIESvSUPERGROUP.influential_species_frequencies.png", sep = ""), width = 1400, height = 1400)
      
      plot.new()
      locs <- list()
      par(mfrow = c(1,1))
      par(fig=c(0, 1, 0.2, 1), new = TRUE)
      barp <- barplot(SPECIESvSUPERGROUP.species_type_topSpecies,xaxt='n',main=paste("Species Frequencies\n(Explaining Most Variation Amongst Super-Group)\n",col_name),col=c(colorRamps::primary.colors(k.best)),log="y",beside=T)
      d <- 1
      for (r in 1:length(barp)){
        if (r %% 3 == 0){ 
          locs[d] <- r+d-2
          d <- d + 1
        }
      }
      axis(1,at=locs,labels=paste(tax_sum.nu[as.character(as.factor(nu_colnames)),],"\n",colnames(SPECIESvSUPERGROUP.species_type_topSpecies),"\n(",nu_colnames,")"), las=3,cex.axis=1.2)
      legend("topright",legend=c(1:k.best),pch=22,pt.bg=c(colorRamps::primary.colors(k.best)),y.intersp=2)
      
      dev.off()
      
      
      
      
      #SPECIESvSUPERGROUP.species_type_topSpecies <- matrix(nrow=k.best,ncol=10)
      #colnames(SPECIESvSUPERGROUP.species_type_topSpecies) <- colnames(sort(SPECIESvSUPERGROUP.covar.pvals)[1:10])
      #rownames(SPECIESvSUPERGROUP.species_type_topSpecies) <- c(1:k.best)
      #for (tt in rownames(SPECIESvSUPERGROUP.species_type_topSpecies)){
      #  for (rr in colnames(SPECIESvSUPERGROUP.species_type_topSpecies)){
      #    SPECIESvSUPERGROUP.species_type_topSpecies[tt,rr] <- species_type[tt,rr]
      #  }
      #}
      
      # SPECIES
      #write.table(SPECIESvSUPERGROUP.species_type_topSpecies, file=paste(prefix,".",col_name,"SPECIESvSUPERGROUP.species_type_topSpecies.tsv", sep = ""), quote=FALSE, sep='\t')
      
      
      
      ############################################################################################################
      ############################################################################################################
      ################# SPECIES v PC VECTOR ######################################################################
      ############################################################################################################
      ############################################################################################################
      
      covar.pvals <- data.frame(matrix(ncol=2))
      covar.model_pvals <- data.frame(ncol=1)
      
      for(i in 1:N){
        
        print(i)
        reg <- lm( as.numeric(species_relab.subgroup[,i]) ~ as.numeric(pcoa_mice.subgroup$vectors[,1]) + as.numeric(pcoa_mice.subgroup$vectors[,2]), na.action=na.omit )

        #png(filename=paste(prefix,".subset.",col_name,".",i,".species_vs_pcoa_PC_residuals.png", sep = ""), width = 500, height = 500)
        #plot.new()
        #ggplot(data=reg,aes(x=as.factor(factor_sample_super_cluster_type[,col_name]),y=species_relab.subgroup[,i]))+geom_point()
        #dev.off()
        
        if (t(coef(summary(reg)))[2,1] > 0 & is.nan(t(coef(summary(reg)))[2,1]) == 0){
          cpus.lm2 <- stepAIC(reg, trace = FALSE)
          
          fstat <- summary(reg)$fstatistic
          covar.model_pvals[i] <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
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
          
          covar.model_pvals[i] <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
          covar.pvals[i,] <- matrix(0,1,2)
          
        }
        
      }
      
      ###### SIGNIFICANT TAXA, P-VALUE TABLE ########
      
      png(filename=paste(prefix,".subset.",col_name,".significant_species_by_PC.png", sep = ""), width = 500, height = 500)
      plot.new()
      color_q <- c(colorRamps::matlab.like2(length(colnames(covar.pvals))))
      legend_q <- c()
      plot(sort(as.numeric(covar.model_pvals)),col="black",pch=15,cex=1.0,ylim=c(0,1),ylab="P-value",xlab="Taxa Index (Independently Sorted by Column)")
      par(new=TRUE)
      for(q in 1:length(colnames(covar.pvals))){
        plot(sort(covar.pvals[,q]),col=color_q[q],pch=22,cex=0.5,ylim=c(0,1),main=paste(col_name,"\nFactor : Signif. Taxa Count"),ylab="P-value",xlab="Taxa Index (Independently Sorted by Column)")
        #legend_q[q] <- paste(unlist(as.list(unique(factor_sample_super_cluster_type)))[q]," : ",sum(sort(covar.pvals[,q])<0.05)," (",round(sum(sort(covar.pvals[,q])<0.05)/length(factor_sample_super_cluster_type),digits=2),"%)")
        legend_q[q] <- paste("PC",q," : ",sum(sort(covar.pvals[,q])<0.05)," (",round(sum(sort(covar.pvals[,q])<0.05)/length(factor_sample_super_cluster_type),digits=2)*100,"%)")
        par(new=TRUE)
      }
      color_q[q+1] <- "black"
      legend_q[q+1] <- paste("Model : ",sum(sort(as.numeric(covar.model_pvals))<0.05)," (",round(sum(sort(as.numeric(covar.model_pvals))<0.05)/length(covar.model_pvals),digits=2)*100,"%)")
      abline(h=0.05,lty=4, col="red")
      legend("topleft",legend=legend_q,pch=22,pt.bg=color_q,cex=0.7,bty='n')
      colnames(covar.model_pvals) <- colnames(species_relab.subgroup)
      species_signif_pvals <- covar.model_pvals[,covar.model_pvals<0.05]
      colnames(species_signif_pvals) <- colnames(covar.model_pvals)[covar.model_pvals<0.05]
      species_signif_pvals <- sort(species_signif_pvals)
      legend("top",legend=paste(colnames(species_signif_pvals)[1:5]," p=",sort(species_signif_pvals)[1:5]),bty='n',cex=0.7)
      dev.off()
      
      write.table(t(sort(covar.model_pvals)), file=paste(prefix,".",col_name,".signifSPECIESassigned_PC.tsv", sep = ""), quote=FALSE, sep='\t')
      rownames(covar.pvals) <- colnames(covar.model_pvals)
      for(q in 1:length(colnames(covar.pvals))){
        tmp <- as.data.frame(covar.pvals[,q])
        rownames(tmp) <- rownames(covar.pvals)
        tmp <- tmp[order(tmp),,drop=FALSE]
        write.table(tmp, file=paste(prefix,".",col_name,".",unlist(as.list(unique(factor_sample_super_cluster_type)))[q],".signifSPECIESassigned_PC.tsv", sep = ""), quote=FALSE, sep='\t')
      }
        
      SPECIESvPC.covar.pvals <- covar.model_pvals
      
      SPECIESvPC.covar.pvals <- species_type[,as.factor(colnames(sort(SPECIESvPC.covar.pvals))[1:10])]
        
      
      ##### barplot
      
      tax_sum.nu <- tax_summary[3]
      rownames(tax_sum.nu) <- as.factor(unlist(tax_summary[2]))
      
      nu_colnames <- c()
      c <- 1
      for (spec_id in as.factor(colnames(SPECIESvSUPERGROUP.species_type_topSpecies))){
        old_spec_id <- unlist(strsplit(spec_id, '[.]'))
        nu_spec_id <- list()
        nu_spec_id[[1]] <- 0 # set root
        zt <- 2
        for(zzz in 2:length(old_spec_id)){
          if(old_spec_id[[zzz]] != 0){
            nu_spec_id[[zt]] <- old_spec_id[[zzz]]
            zt <- zt + 1
          }
        }
        query_id <- paste(as.list(nu_spec_id),collapse=".")
        print(paste(spec_id," to ",query_id))
        nu_colnames[c] <- query_id
        c <- c + 1
      }
      
      png(filename=paste(prefix,".all_samples.",col_name,".SPECIESvSUPERGROUP.influential_species_frequencies.png", sep = ""), width = 1400, height = 1400)
      
      plot.new()
      locs <- list()
      par(mfrow = c(1,1))
      par(fig=c(0, 1, 0.2, 1), new = TRUE)
      barp <- barplot(SPECIESvSUPERGROUP.species_type_topSpecies,xaxt='n',main=paste("Species Frequencies\n(Explaining Most Variation Amongst Super-Group)\n",col_name),col=c(colorRamps::primary.colors(k.best)),log="y",beside=T)
      d <- 1
      for (r in 1:length(barp)){
        if (r %% 3 == 0){ 
          locs[d] <- r+d-2
          d <- d + 1
        }
      }
      axis(1,at=locs,labels=paste(tax_sum.nu[as.character(as.factor(nu_colnames)),],"\n",colnames(SPECIESvSUPERGROUP.species_type_topSpecies),"\n(",nu_colnames,")"), las=3,cex.axis=1.2)
      legend("topright",legend=c(1:k.best),pch=22,pt.bg=c(colorRamps::primary.colors(k.best)),y.intersp=2)
      
      dev.off()
      
        
        
        
        
        
      
      genus <- species
      
      #########################################
      # initialize matrix with genus-level taxa
      col_x8 <- list()
      ty <- 1
      for (m in colnames(genus)){
        col_x8[[ty]] <- paste(unlist(strsplit(m, '[.]'))[1:7],collapse=".")
        ty <- ty+1
      }
      genus_type <- matrix(nrow=k.best,ncol=length(col_x8))
      colnames(genus_type) <- col_x8
      rownames(genus_type) <- c(1:k.best)
      for (tt in rownames(genus_type)){
        for (rr in colnames(genus_type)){
          genus_type[tt,rr] <- 0
        }
      }
      
      # determine most common taxa
      genus_sum <- list()
      sc <- 1
      for (r in rownames(genus)){ # sample id w genus association
        sum <- 0
        if(r %in% names(pk$clustering)){ # sample id in clustering
          # sample with genus association has cluster
          rc <- 1
          cluster_id <- pk$clustering[[r]]
          for (m in colnames(genus)){
            genus_id <- paste(unlist(strsplit(colnames(genus)[rc], '[.]'))[1:7],collapse=".")
            genus_type[cluster_id,genus_id] <- genus_type[cluster_id,genus_id] + genus[r,m]
            sum <- sum + genus[r,m]
            rc <- rc + 1
          }
        }
        genus_sum[[sc]] <- sum
        sc <- sc + 1
      }
      names(genus_sum) <- rownames(genus)
      
      # SPECIES
      write.table(genus_type, file=paste(prefix,".",col_name,".allGENUSassigned_subgroup.tsv", sep = ""), quote=FALSE, sep='\t')
      
      # get the top 10 for each cluster
      top <- list()      
      for (k in 1:k.best){
        top <- append(top, names(rev(sort(genus_type[k,]))[1:10]))
      }
      # and the unique set amongst all of those top 10s
      top <- unique(top)
      
      genus_type_topgenus <- matrix(nrow=k.best,ncol=length(top))
      colnames(genus_type_topgenus) <- top
      rownames(genus_type_topgenus) <- c(1:k.best)
      for (tt in rownames(genus_type_topgenus)){
        for (rr in colnames(genus_type_topgenus)){
          genus_type_topgenus[tt,rr] <- genus_type[tt,rr]
        }
      }
      
      # SPECIES
      write.table(genus_type_topgenus, file=paste(prefix,".",col_name,".topGENUSassigned_subgroup.tsv", sep = ""), quote=FALSE, sep='\t')
      
      
      ###############################################
      # parse tax_summary to get labels
      #tax_summary <- read.table("otus.qfilter.EE0.15.curated_SILVA_123_SSURef_Nr99_tax_silva.wang.tax.summary", sep="\t", header=TRUE)
      
      ##### PARSE TAXONOMY LABELS ######
      #tax_summary <- read.table("otus.qfilter.EE0.15.curated_SILVA_123_SSURef_Nr99_tax_silva.wang.tax.summary", sep="\t", header=TRUE)
      species.taxon_id <- list()
      species.short_rankid <- list()
      si <- 1
      for(spec_id in colnames(species_type_topSpecies)){
        
        old_spec_id <- unlist(strsplit(spec_id, '[.]'))
        nu_spec_id <- list()
        nu_spec_id[[1]] <- 0 # set root
        zt <- 2
        for(zzz in 2:length(old_spec_id)){
          if(old_spec_id[[zzz]] != 0){
            nu_spec_id[[zt]] <- old_spec_id[[zzz]]
            zt <- zt + 1
          }
        }
        query_id <- paste(as.list(nu_spec_id),collapse=".")
        tax_assign <- filter(tax_summary,rankID==query_id)$taxon
        print(paste("tax:",tax_assign))
        if(as.character(tax_assign) == "unclassified"){ # if unclassified, iterate down until defined)
          print(paste("!",tax_assign))
          while((as.character(tax_assign) == "unclassified")==TRUE){
            cat(paste("       ",spec_id,"->",query_id,"~>",tax_assign,"\n"))
            old_spec_id <- unlist(strsplit(query_id, '[.]'))
            nu_spec_id <- old_spec_id[1:(length(old_spec_id)-1)]
            query_id <- paste(as.list(nu_spec_id),collapse=".")
            tax_assign <- filter(tax_summary,rankID==query_id)$taxon
            print(paste("new taxassign : ",tax_assign))
          }
        }
        species.taxon_id[[si]] <- as.character(tax_assign)
        species.short_rankid[[si]] <- query_id
        si <- si+1
        cat(paste(spec_id,"->",query_id,"~>",tax_assign,"\n"))
      }
      names(species.taxon_id) <- spec_id
      
      genus.taxon_id <- list()
      genus.short_rankid <- list()
      si <- 1
      for(spec_id in colnames(genus_type_topgenus)){
        
        old_spec_id <- unlist(strsplit(spec_id, '[.]'))
        nu_spec_id <- list()
        nu_spec_id[[1]] <- 0 # set root
        zt <- 2
        for(zzz in 2:length(old_spec_id)){
          if(old_spec_id[[zzz]] != 0){
            nu_spec_id[[zt]] <- old_spec_id[[zzz]]
            zt <- zt + 1
          }
        }
        query_id <- paste(as.list(nu_spec_id),collapse=".")
        tax_assign <- filter(tax_summary,rankID==query_id)$taxon
        print(paste("tax:",tax_assign))
        if(as.character(tax_assign) == "unclassified"){ # if unclassified, iterate down until defined)
          print(paste("!",tax_assign))
          while((as.character(tax_assign) == "unclassified")==TRUE){
            cat(paste("       ",spec_id,"->",query_id,"~>",tax_assign,"\n"))
            old_spec_id <- unlist(strsplit(query_id, '[.]'))
            nu_spec_id <- old_spec_id[1:(length(old_spec_id)-1)]
            query_id <- paste(as.list(nu_spec_id),collapse=".")
            tax_assign <- filter(tax_summary,rankID==query_id)$taxon
            print(paste("new taxassign : ",tax_assign))
          }
        }
        genus.taxon_id[[si]] <- as.character(tax_assign)
        genus.short_rankid[[si]] <- query_id
        si <- si+1
        cat(paste(spec_id,"->",query_id,"~>",tax_assign,"\n"))
      }
      names(genus.taxon_id) <- spec_id
      
      
      
      ################### BUILD IDEAL SUMMARY PLOT ###################
      
      png(filename=paste(prefix,".all_samples.",col_name,".auto_clustered_PCoA.png", sep = ""), width = 2000, height = 2000)
      
      plot.new()
      par(mfrow = c(1,1))
      
      par(fig=c(0, 0.4, 0.04, 0.2), new = TRUE)
      barplot(species_type_topSpecies)
      barp <- barplot(species_type_topSpecies,xaxt='n',main="Species Frequencies (Max Count)",col=c(colorRamps::primary.colors(k.best)),new=TRUE)
      par(fig=c(0, 0.4, 0.04, 0.2), new = TRUE)
      axis(1,at=barp,labels=paste(species.taxon_id,"\n",colnames(species_type_topSpecies),"\n(",species.short_rankid,")"), las=3,cex.axis=1.2)
      legend("topright",legend=c(1:k.best),pch=22,pt.bg=c(colorRamps::primary.colors(k.best)),y.intersp=2)
      
      par(fig=c(0.4, 0.8, 0.04, 0.2), new = TRUE)
      barp <- barplot(genus_type_topgenus,xaxt='n',main="Genus Frequencies (Max Count)",col=c(colorRamps::primary.colors(k.best)),new=TRUE)
      par(fig=c(0.4, 0.8, 0.04, 0.2), new = TRUE)
      axis(1,at=barp,labels=paste(genus.taxon_id,"\n",colnames(genus_type_topgenus),"\n(",genus.short_rankid,")"), las=3,cex.axis=1.2)
      legend("topright",legend=c(1:k.best),pch=22,pt.bg=c(colorRamps::primary.colors(k.best)),y.intersp=2)
      
      pk <- pam(as.data.frame(pcoa_mice.subgroup$vectors[,1:2]),k.best)
      par(fig=c(0, 0.8, 0.2, 1), new = TRUE)
      if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
        clusplot(pk,main=paste("\n",prefix," : ", col_name,"\nClusters (PAM,PC1+2), Ideal Cluster Count =",k.best),col.p=col_type,color=T,metric="euclidean",shade=T,labels=5,cex=3.0)
      }else{
        clusplot(pk,main=paste("\n",prefix," : ", col_name,"\nClusters (PAM,PC1+2), Ideal Cluster Count =",k.best),col.p=col_type,color=T,metric="euclidean",shade=T,labels=5,cex=1.0)
      }
      legend( "topright", legend=mixedsort(unique(unlist(groups.multi.sub_group))), pch=22, bty='n', pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
      legend( "topleft", legend=c("Lowest Density","","","Highest Density"), pch=22, bty='n', pt.bg=c("blue","green","red","purple") )
      
      rowsize = 0.84/(k.best+2)
      row_loc <- 1-rowsize
      par(fig=c(.78, 1, row_loc, 1), new = TRUE)
      par(new=TRUE)
      barplot(pk$clusinfo[,1],col=c(colorRamps::primary.colors(k.best)),names.arg=c(1:dim(pk$clusinfo)[1]),main="Cluster Size")
      
      par(fig=c(.78, 1, row_loc-rowsize, row_loc), new = TRUE)
      par(new=TRUE)
      barplot(pk$silinfo$clus.avg.widths,col=c(colorRamps::primary.colors(k.best)),names.arg=c(1:dim(pk$clusinfo)[1]),main="Cluster Avg Width")
      row_loc <- row_loc-rowsize
      
      for(k in 1:k.best){
        
        if(k.best <= 4){
          par(fig=c(.75, .89, row_loc-rowsize, row_loc), new = TRUE)
        }else{
          par(fig=c(.75, .89, row_loc-(rowsize+0.03), (row_loc+0.03)), new = TRUE)
        }
        par(new=TRUE)
        if (k==1){
          #pie(t(as.data.frame(as.matrix(table(factor(unlist(super_cluster_type[k,])))))),col=c(colorRamps::matlab.like(dim(super_cluster_type)[2])),cex = 0.1,main="Super-Group")
          #par(mar=c(0,0,0,0))
          pie(unlist(super_cluster_type[k,]),col=c(colorRamps::matlab.like(dim(super_cluster_type)[2])),cex = 0.1,main=paste("\nSuper-Group\n(",k,")"))
        }else{
          pie(unlist(super_cluster_type[k,]),col=c(colorRamps::matlab.like(dim(super_cluster_type)[2])),cex = 0.1,main=paste("\n(",k,")"))
        }
        if(k.best <= 4){
          par(fig=c(.84, .98, row_loc-rowsize, row_loc), new = TRUE)
        }else{
          par(fig=c(.84, .98, row_loc-(rowsize+0.03), (row_loc+0.03)), new = TRUE)
        }
        par(new=TRUE)
        if (k==1){
          pie(unlist(sub_cluster_type[k,]),col=c(colorRamps::matlab.like(dim(sub_cluster_type)[2])),cex = 0.1,main=paste("\nSub-Group\n(",k,")"))
        }else{
          pie(unlist(sub_cluster_type[k,]),col=c(colorRamps::matlab.like(dim(sub_cluster_type)[2])),cex = 0.1,main=paste("\n(",k,")"))
        }
        #legend( "topleft", legend=levels(factor(unlist(q_supers))), pch=22, col="grey", pt.bg=c(colorRamps::matlab.like(length(unique(factor(unlist(q_supers)))))), cex=3 ) 
        #legend( "topright", legend=mixedsort(levels(factor(unlist(q_subs)))), pch=22, col="grey", pt.bg=c(colorRamps::matlab.like(length(unique(factor(unlist(q_subs)))))))
        row_loc <- row_loc-rowsize
        
      }
      
      par(fig=c(0.8, 1, 1-((k.best+2)*rowsize), 1-(2*rowsize)), new = TRUE)
      legend("right",
             legend = colnames(sub_cluster_type),
             pch=22,
             pt.bg=c(colorRamps::matlab.like(dim(sub_cluster_type)[2])),
             bty='y')
      
      par(fig=c(0.8, 1, 0, 1-1.5*rowsize), new = TRUE)
      legend("right",inset=0.75,
             pch=22,
             legend = colnames(super_cluster_type),
             pt.bg=c(colorRamps::matlab.like(dim(super_cluster_type)[2])),
             bty='y')
      
      par(fig=c(.8, 1, 0, 0.2), new = TRUE)
      mykde2d(pcoa_mice.subgroup$vectors, col="red", lwd=1.4, lty=1,main=NULL)
      if(length(rownames(species_relab))<100){
        #s.class(pcoa_mice.subgroup$vectors, factor(groups[,col_count]), cpoint = 1, sub = paste(main, "original", sep="\n"), csub=1, col=colors, add.p=TRUE, label="")
        # s.class(pcoa_mice.subgroup$vectors, factor(unlist(tmp)), cpoint = 1, sub = paste(main, "original", sep="\n"), csub=1, col=colors, add.p=TRUE, label="")
      }
      #s.corcircle(as.numeric(pcoa_mice.subgroup$vectors[, 1]), full = TRUE, box = TRUE, sub = paste(main, "Contribution of variables to axes.",sep="\n"), csub=1, possub="bottom")
      
      dev.off()
      
      
      if(k.best <= 4){
        png(filename=paste(prefix,".all_samples.",col_name,".auto_clustered_pies.png", sep = ""), width = 1400, height = 1500)
      }else{
        png(filename=paste(prefix,".all_samples.",col_name,".auto_clustered_pies.png", sep = ""), width = 1400, height = 3000)
      }
      plot.new()
      par(mfrow = c(1,1))
      
      row_loc <- 1
      rowsize <- 1/k.best
      for(k in 1:k.best){
        if((row_loc-rowsize) < 0){row_loc<-rowsize}
        par(fig=c(0, 0.5, row_loc-rowsize, row_loc), new = TRUE)
        par(new=TRUE)
        if (k==1){
          pie(unlist(super_cluster_type[k,]),col=c(colorRamps::matlab.like(dim(super_cluster_type)[2])),cex = 0.1,main=paste("\n\nSuper-Group\n(",k,")"))
        }else{
          pie(unlist(super_cluster_type[k,]),col=c(colorRamps::matlab.like(dim(super_cluster_type)[2])),cex = 0.1,main=paste("\n(",k,")"))
        }
        par(fig=c(.5, 1, row_loc-rowsize, row_loc), new = TRUE)
        par(new=TRUE)
        if (k==1){
          pie(unlist(sub_cluster_type[k,]),col=c(colorRamps::matlab.like(dim(sub_cluster_type)[2])),cex = 0.1,main=paste("\n\nSub-Group\n(",k,")"))
        }else{
          pie(unlist(sub_cluster_type[k,]),col=c(colorRamps::matlab.like(dim(sub_cluster_type)[2])),cex = 0.1,main=paste("\n(",k,")"))
        }
        row_loc <- row_loc-rowsize
        
      }
      par(fig=c(0, 1, 0, 1), new = TRUE)
      legend("topright",
             legend = colnames(sub_cluster_type),
             pch=22,
             pt.bg=c(colorRamps::matlab.like(dim(sub_cluster_type)[2])),
             bty='n')
      legend("topleft",
             legend = colnames(super_cluster_type),
             pch=22,
             pt.bg=c(colorRamps::matlab.like(dim(super_cluster_type)[2])),
             bty='n')
      
      dev.off()
      
      
      png(filename=paste(prefix,".all_samples.",col_name,".species_frequencies.png", sep = ""), width = 1400, height = 1400)
      
      plot.new()
      par(mfrow = c(1,1))
      par(fig=c(0, 1, 0.2, 1), new = TRUE)
      barp <- barplot(species_type_topSpecies,xaxt='n',main="Species Frequencies (Max Count)",col=c(colorRamps::primary.colors(k.best)))
      axis(1,at=barp,labels=paste(species.taxon_id,"\n",colnames(species_type_topSpecies),"\n(",species.short_rankid,")"), las=3,cex.axis=1.2)
      legend("topright",legend=c(1:k.best),pch=22,pt.bg=c(colorRamps::primary.colors(k.best)),y.intersp=2)
      
      dev.off()
      
      png(filename=paste(prefix,".all_samples.",col_name,".genus_frequencies.png", sep = ""), width = 1400, height = 1400)
      
      plot.new()
      par(mfrow = c(1,1))
      par(fig=c(0, 1, 0.2, 1), new = TRUE)
      barp <- barplot(genus_type_topgenus,xaxt='n',main="Genus Frequencies (Max Count)",col=c(colorRamps::primary.colors(k.best)))
      axis(1,at=barp,labels=paste(genus.taxon_id,"\n",colnames(genus_type_topgenus),"\n(",genus.short_rankid,")"), las=3,cex.axis=1.2)
      legend("topright",legend=c(1:k.best),pch=22,pt.bg=c(colorRamps::primary.colors(k.best)),y.intersp=2)
      
      dev.off()
      
      
      
      
      png(filename=paste(prefix,".all_samples.",col_name,".ClusterSummary.png", sep = ""), width = 1800, height = 1800)
      par(mfrow = c(2,2))
      if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
        plot(pcoa_mice.subgroup$vectors, col=col_type, pch=pch_type, data=rownames(pcoa_mice.subgroup$vectors), cex=3.0, main=paste("[Within Subgroup, 2-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[2,3], digits=2),"%)"))
      }else{
        plot(pcoa_mice.subgroup$vectors, col=col_type, pch=pch_type, data=rownames(pcoa_mice.subgroup$vectors), cex=0.7, main=paste("[Within Subgroup, 2-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[2,3], digits=2),"%)"))
      }
      legend( "topleft", legend=topleft_text, pch=unique(pch_type), col=1, cex=1.0 ) 
      legend( "topright", legend=mixedsort(unique(unlist(groups.multi.sub_group))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
      
      mykde2d(pcoa_mice.subgroup$vectors, col="red", lwd=1.4, lty=1)
      if(length(rownames(species_relab))<100){
        # s.class(pcoa_mice.subgroup$vectors, factor(groups[,col_count]), cpoint = 1, sub = paste(main, "original", sep="\n"), csub=1, col=colors, add.p=TRUE, label="")
      }
      #s.corcircle(as.numeric(pcoa_mice.subgroup$vectors[, 1]), full = TRUE, box = TRUE, sub = paste(main, "Contribution of variables to axes.",sep="\n"), csub=1, possub="bottom")
      
      num_clust <- k.best
      pk <- pam(as.data.frame(pcoa_mice.subgroup$vectors[,1:2]),num_clust)
      if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
        clusplot(pk,main="\nClusters (PAM) : First and Second Component",col.p=col_type,color=TRUE,labels=4,span=TRUE) #`pam` does you partitioning
      }else{
        clusplot(pk,main="\nClusters (PAM) : First and Second Component",col.p=col_type,color=TRUE,labels=4,span=TRUE) #`pam` does you partitioning
      }
      legend( "topleft", legend=c(1:num_clust), pch=c(1:num_clust), col=1, cex=1.0 ) 
      legend( "topright", legend=mixedsort(unique(unlist(groups.multi.sub_group))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
      
      #plot.new()
      
      # v1.3 update
      #img = readPNG(paste(prefix,".all_samples.",col_name,".auto_clustered_pies.png", sep = ""))
      #lim <- par()
      #rasterImage(img, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4],add=TRUE)
      
      dev.off()
      
      
      png(filename=paste(prefix,".all_samples.",col_name,".PAM_Cluster.png", sep = ""), width = 1400, height = 1400)
      
      if (is.na(unique(mixedsort(strtoi(unlist(groups.multi.sub_group))))[1])==TRUE){
        topright_text <- topleft_text  
      }
      
      par(mfrow = c(2,2))
      if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
        plot(pcoa_mice.subgroup$vectors, col=col_type, pch=pch_type, data=rownames(pcoa_mice.subgroup$vectors), cex=3.0, main=paste("[Within Subgroup, 2-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[2,3], digits=2),"%)"))
      }else{
        plot(pcoa_mice.subgroup$vectors, col=col_type, pch=pch_type, data=rownames(pcoa_mice.subgroup$vectors), cex=0.7, main=paste("[Within Subgroup, 2-column]\n",prefix, " : ", col_name), xlab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[1,3], digits=2),"%)"), ylab=paste("PCoA.1(",round(100*pcoa_mice_values.subgroup[2,3], digits=2),"%)"))
      }
      legend( "topleft", legend=topleft_text, pch=unique(pch_type), col=1, cex=1.0 ) 
      legend( "topright", legend=mixedsort(unique(unlist(groups.multi.sub_group))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
      
      pk <- pam(as.data.frame(pcoa_mice.subgroup$vectors[,1:2]),num_clust)
      if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
        clusplot(pk,main="\nClusters (PAM) : First and Second Component (Euclidean)",col.p=col_type,color=TRUE,labels=5,metric="euclidean",span=FALSE,cex=3.0) #`pam` does you partitioning
      }else{
        clusplot(pk,main="\nClusters (PAM) : First and Second Component (Euclidean)",col.p=col_type,color=TRUE,labels=5,metric="euclidean",span=FALSE) #`pam` does you partitioning
      }
      legend( "topleft", legend=c(1:num_clust), pch=c(1:num_clust), col=1, cex=1.0 ) 
      legend( "topright", legend=mixedsort(unique(unlist(groups.multi.sub_group))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
      
      pk <- pam(as.data.frame(pcoa_mice.subgroup$vectors[,1:2]),num_clust)
      if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
        clusplot(pk,main="\nClusters (PAM) : First and Second Component (Euclidean)\nZoomed on Tightest Cluster Bound",col.p=col_type,color=TRUE,labels=5,metric="euclidean",span=TRUE,shade=TRUE,cex=3.0) #`pam` does you partitioning
      }else{
        clusplot(pk,main="\nClusters (PAM) : First and Second Component (Euclidean)\nZoomed on Tightest Cluster Bound",col.p=col_type,color=TRUE,labels=5,metric="euclidean",span=TRUE,shade=TRUE) #`pam` does you partitioning
      }
      legend( "topleft", legend=c(1:num_clust), pch=c(1:num_clust), col=1, cex=1.0 ) 
      legend( "topright", legend=mixedsort(unique(unlist(groups.multi.sub_group))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
      
      #pk <- pam(as.data.frame(pcoa_mice.subgroup$vectors[,1:2]),num_clust)
      #clusplot(pk,main="\nClusters (PAM) : First and Second Component (Manhattan)",col.p=col_type,color=TRUE,metric="manhattan") #`pam` does you partitioning
      #legend( "topleft", legend=c(1:num_clust), pch=c(1:num_clust), col=1, cex=1.0 ) 
      #legend( "topright", legend=unique(mixedsort(unlist(groups.multi.sub_group))), pch=22, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
      
      pk <- pam(as.data.frame(pcoa_mice.subgroup$vectors[,1:3]),num_clust)
      if((num_subgroups<2)||(length(rownames(species_relab))<=50)){
        clusplot(pk,main="\nClusters (PAM) : First 3 Components (Euclidean)",col.p=col_type,labels=5,color=TRUE,cex=3.0) #`pam` does you partitioning
      }else{
        clusplot(pk,main="\nClusters (PAM) : First 3 Components (Euclidean)",col.p=col_type,labels=5,color=TRUE) #`pam` does you partitioning
      }
      legend( "topleft", legend=c(1:num_clust), pch=c(1:num_clust), col=1, cex=1.0 ) 
      #if (is.na(unique(mixedsort(strtoi(unlist(groups.multi.sub_group))))[1])==FALSE){
      #  legend( "topright", legend=mixedsort(as.numeric(unique(unlist(groups.multi.sub_group)))), pch=22, labels=4, pt.bg=c(colorRamps::matlab.like2(length(unique(mixedsort(unlist(groups.multi.sub_group)))))) )
      #}
      dev.off()
      
      
      ######END DEBUG######
      
      col_count <- col_count + 1    
    }
  }
}
