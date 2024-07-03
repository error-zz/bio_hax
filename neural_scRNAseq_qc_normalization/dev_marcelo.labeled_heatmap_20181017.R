#!/usr/bin/env Rscript
library(gplots)
library(colorRamps)
library(ape)
library(vegan)
library(caper)
library(plot3D)
library(rgl)
library(heatmap.plus)
#install.packages("beanplot",repos='http://cran.us.r-project.org')
library(beanplot)

setwd ("/Users/apple/Desktop/develop.Marcelo")
 
CC <-  read.csv("sumCC.tsv")
MP <-  read.csv("sumMP.tsv")
BP <-  read.csv("sumBP.tsv")
KD <-  read.csv("KEGG_Diseases_pval.csv")
KB <-  read.csv("KEGG_Biological_Pathways.csv")
GBP <- read.csv("GO_Biological_Processes.csv")
GMP <- read.csv("GO_Molecular_Functions.csv")
GCC <- read.csv("GO_Cellular_Components.csv")

HvsT2D.all.out.annotated <- read.table("HvsT2D.all.out.annotated.txt", sep="\t", header=TRUE, fill=TRUE)
HvsT2D.all.out.annotated[ , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ] <- t(decostand( t(HvsT2D.all.out.annotated[ , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ]) , method="total" ))
 
########################################################################################################################
########### (1) CC : STANDARD ##########################################################################################
######################################################################################################################## 
 
CC_matrix <- matrix(ncol=8,nrow=sum(as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1])))
rownames(CC_matrix) <- as.character( HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]) , ]$Gene_name )
colnames(CC_matrix) <- as.character(GCC[1:8,2])
CC_matrix[is.na(CC_matrix)] <- "#FFFFFF"
CC_color=primary.colors(8)
for( q in 1:8 ){
  filename <- paste0("CC",q,".csv")
  infile <-  read.csv(filename)
  CC_matrix[unique(as.character(infile[,1])),q] <- CC_color[q]
}
colnames(CC_matrix) <- c("1","2       ",
                         "3","4       ",
                         "5","6       ",
                         "7","8       ")

h<- heatmap.2(log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01))

h.orderedColumnIDs <- colnames(HvsT2D.all.out.annotated[,4:47])[h$colInd]
h.color_col <- c("blue","purple")[as.factor(grepl("Healthy", h.orderedColumnIDs))]
h.orderedRowIDs    <- rev(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]) , ]$Gene_name[h$rowInd])

htmp <- strsplit(h.orderedColumnIDs,"\\.")
tt <- list()
ti <- 1
for(t in htmp){
  ttt <- strsplit(t[2],"_")[[1]][1]
  tt[ti] <- ttt
  ti <- ti + 1
}
metadata          <- read.table("metadata.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$Sample
md <- metadata[ unlist(tt) , c(2,6,8,9,10) ]
md.col <- md

md.col[md[,5] == "healthy",5] <- c(green2red(5),"black")[1]
md.col[md[,5] == "Mild",5] <- c(green2red(5),"black")[2]
md.col[md[,5] == "Moderate",5] <- c(green2red(5),"black")[3]
md.col[md[,5] == "Severe",5] <- c(green2red(5),"black")[4]
md.col[md[,5] == "Sever",5] <- c(green2red(5),"black")[5]
md.col[md[,5] == "Gingivitis",5] <- c(green2red(5),"black")[6]
 
md.col[as.character(md[,4]) == "B",4] <- "turquoise3"
md.col[as.character(md[,4]) == "W",4] <- "blue4"
md.col[as.character(md[,4]) == "H",4] <- "slateblue"

md.col[as.character(md[,3]) == "S",3] <- "brown4"
md.col[as.character(md[,3]) == "FS",3] <- "sienna4"
md.col[as.character(md[,3]) == "N",3] <- "springgreen4"
md.col[as.character(md[,3]) == "NS",3] <- "greenyellow"

md.col[as.character(md[,2]) == "M",2] <- "navyblue"
md.col[as.character(md[,2]) == "F",2] <- "orchid2"

md.col[as.character(md[,1]) == "Healthy",1] <- "green"
md.col[as.character(md[,1]) == "T2D",1] <- "red"

colnames(md.col) <- c("a       ","b",
                      "c       ","c",
                      "e       ")

png(filename="CellularComponents.png", width = 1000, height = 1000)
heatmap.plus(  log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01) ,
               col=c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="Cellular Components",
               mar=c(20,5),
               labRow=as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]) , ]$Gene_name),
               #ColSideColors = t(as.matrix(c("blue","purple")[as.factor(grepl("Healthy", h.orderedColumnIDs))])),
               ColSideColors = as.matrix(md.col),
               RowSideColors = CC_matrix, #t(CC_matrix),
               cexCol=0.9
              )
legend("bottomleft",
        legend = paste(c(1:8),"=",as.character(GCC[,2])[1:8]),
        #inset=c(-0.05,-0.05),
        cex=0.6
        )
legend("topleft",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=0.6
)
dev.off()
 
########### (1) CC : SORTED ##########################################################################################

reOrd <- as.numeric(CC_matrix[,1] != "#FFFFFF")
names(reOrd) <- HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]) , ]$Gene_name

# 2
fill <- 2
for(z in c(2:8)){
  
  r <- reOrd[ names(CC_matrix[ CC_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF")[CC_matrix[ CC_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF"] ]
  reOrd[ names(CC_matrix[ CC_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF")[CC_matrix[ CC_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF"] ][ r >= (fill-1) ] <- fill
  tmp <- as.numeric(CC_matrix[,z][reOrd == 0] != "#FFFFFF" )
  fill <- fill + 1
  
  tmp[tmp==1] <- fill
  reOrd[ ! reOrd %in% c(1:as.numeric(z-1)) ][ reOrd[ ! reOrd %in% c(1:as.numeric(z-1)) ]==0 ] <- tmp
  fill <- fill + 1
  
}

# name it

subOrd <- list()
iter <- sum(as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]))
for(y in max(reOrd):1){
  if(sum(reOrd==y) == 1){
    subOrd[iter] <- names(reOrd)[ reOrd==y ]
    iter <- iter - 1
  }
  else if(sum(reOrd==y) > 0){
    tq <- log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
    pvals <- list()
    slopes <- list()
    p <- 1
    for( t in 1:dim( tq[reOrd==y,] )[1] ){
      lm_reg <- lm(   as.numeric(tq[reOrd==y,][t,]) ~ as.factor(grepl("Healthy", colnames(tq))))
      fstat  <- summary(lm_reg)$fstatistic
      pvals[p]   <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
      slopes[p]  <- as.numeric(as.matrix(coef(lm_reg))[2,1])
      p <- p + 1
    }
    slopes <- as.numeric(slopes)
    names(slopes) <-  HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]) , ]$Gene_name[reOrd==y]
    slopes.sort <- sort(slopes)
    for(s in names(slopes.sort)){
      subOrd[iter] <- s
      iter <- iter - 1
    }
  }
}
labrow <- as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]) , ]$Gene_name)
ordOrd <- list()
for(l in 1:length(labrow)){
  for(q in 1:length(unlist(subOrd))){
    if( labrow[l] == unlist(subOrd)[q] ){
      ordOrd[l] <- q
    }
  }
}
ordOrd <- unlist(ordOrd)

#reOrd <- reOrd[ ! reOrd %in% c(1) ][ as.numeric(CC_matrix[,2][ ! reOrd %in% c(1) ] != "#FFFFFF") ]
 
tq <- log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(CC[,1]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
rownames(tq) <- labrow

colOrder <- c("MF.17_T2D_M_70_N_W_norm","MF.53_T2D_M_70_NS_W_norm","MF.20_T2D_F_68_N_B_norm","MF.54_T2D_F_68_N_B_norm","MF.52_T2D_M_64_FS_W_norm",
              "MF.16_T2D_F_61_N_B_norm","MF.19_T2D_M_57_N_B_norm",
              "MF.15_T2D_M_52_NS_B_norm","MF.18_T2D_F_46_S_B_norm","MF.11_T2D_M_44_NS_B_norm","MF.13_T2D_F_37_FS_H_norm",               
              
              "MF.1_Healthy_M_60_NS_W_norm","MF.48_Healthy_M_59_NS_W_norm","MF.45_Healthy_M_45_NS_W_norm","MF.46_Healthy_M_40_NS_W_norm",   
              "MF.4_Healthy_F_37_N_B_norm", "MF.5_Healthy_F_35_S_W_norm", "MF.3_Healthy_F_35_S_W_norm",
              "MF.7_Healthy_F_35_S_W_norm","MF.47_Healthy_F_33_N_B_norm","MF.2_Healthy_F_29_NS_B_norm", "MF.9_Healthy_F_28_N_B_norm")

md.order <- c(17,53,20,54,52,16,19,
              15,18,11,13,
              1,48,45,46,  
              4,5,3,
              7,47,2,9)

CC_matrix.sort <- CC_matrix[  unlist(subOrd), ] 
tq.sort <- tq[ unlist(subOrd) , colOrder ] 

png(filename="CellularComponents_SORT.png", width = 800, height = 800)
heatmap.plus(  scale( tq.sort ) ,
               col=redblue(150),#c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="\nCellular Components",
               mar=c(20,5),
               ColSideColors = as.matrix(as.data.frame(md.col)[as.character(md.order),]), #as.matrix(md.col),
               RowSideColors = CC_matrix.sort, # CC_matrix, #t(CC_matrix),
               cexCol=0.9,
              Rowv = NA,
              Colv = NA
)

dev.off()

png(filename="CellularComponents_LEGEND.png", width = 500, height = 500)
plot.new()
legend("top",
       legend = paste(c(1:8),"=",as.character(GCC[,2])[1:8]),
       #inset=c(-0.05,-0.05),
       cex=1
)
legend("bottom",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=1
)

dev.off()



##################################
# BP #

BP_matrix <- matrix(ncol=16,nrow=sum(as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(BP[,1])))
rownames(BP_matrix) <- as.character( HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(BP[,1]) , ]$Gene_name )
colnames(BP_matrix) <- as.character(GBP[1:16,2])
BP_matrix[is.na(BP_matrix)] <- "#FFFFFF"
BP_color=c("red","blue","darkgreen","purple","red","blue","darkgreen","purple","red","blue","darkgreen","purple","red","blue","darkgreen","purple")
for( q in 1:16 ){
  if(q<10){
    filename <- paste0("BP0",q,".csv")
  }else{
    filename <- paste0("BP",q,".csv")
  }
  infile <-  read.csv(filename)
  BP_matrix[unique(as.character(infile[,1])),q] <- BP_color[q]
}
colnames(BP_matrix) <- c("1","2       ",
                         "3","4       ",
                         "5","6       ",
                         "7","8       ",
                         "9","10      ",
                         "11","12      ",
                         "13","14      ",
                         "15","16      ")

png(filename="BiologicalProcesses.png", width = 1600, height = 1000)
heatmap.plus(  log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(BP[,1]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01) ,
               col=c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="Biological Pathways",
               mar=c(20,5),
               labRow=as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(BP[,1]) , ]$Gene_name),
               #ColSideColors = t(as.matrix(c("blue","purple")[as.factor(grepl("Healthy", h.orderedColumnIDs))])),
               ColSideColors = as.matrix(md.col),
               RowSideColors = BP_matrix, #t(CC_matrix),
               cexCol=0.9
)
legend("bottomleft",
       legend = paste(c(1:16),"=",as.character(GBP[,2])[1:16]),
       #inset=c(-0.05,-0.05),
       cex=0.6
)
legend("topleft",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=0.6
)
dev.off()
 

# SORT IT (BP) #


reOrd <- as.numeric(BP_matrix[,1] != "#FFFFFF")
names(reOrd) <- HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(BP[,1]) , ]$Gene_name

# 2
fill <- 2
for(z in c(2:16)){
  
  r <- reOrd[ names(BP_matrix[ BP_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF")[BP_matrix[ BP_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF"] ]
  reOrd[ names(BP_matrix[ BP_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF")[BP_matrix[ BP_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF"] ][ r >= (fill-1) ] <- fill
  tmp <- as.numeric(BP_matrix[,z][reOrd == 0] != "#FFFFFF" )
  fill <- fill + 1
  
  tmp[tmp==1] <- fill
  reOrd[ ! reOrd %in% c(1:as.numeric(z-1)) ][ reOrd[ ! reOrd %in% c(1:as.numeric(z-1)) ]==0 ] <- tmp
  fill <- fill + 1
  
}

# name it

subOrd <- list()
iter <- sum(as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(BP[,1]))
for(y in max(reOrd):1){
  if(sum(reOrd==y) == 1){
    subOrd[iter] <- names(reOrd)[ reOrd==y ]
    iter <- iter - 1
  }
  else if(sum(reOrd==y) > 0){
    tq <- log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(BP[,1]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
    pvals <- list()
    slopes <- list()
    p <- 1
    for( t in 1:dim( tq[reOrd==y,] )[1] ){
      lm_reg <- lm(   as.numeric(tq[reOrd==y,][t,]) ~ as.factor(grepl("Healthy", colnames(tq))))
      fstat  <- summary(lm_reg)$fstatistic
      pvals[p]   <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
      slopes[p]  <- as.numeric(as.matrix(coef(lm_reg))[2,1])
      p <- p + 1
    }
    slopes <- as.numeric(slopes)
    names(slopes) <-  HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(BP[,1]) , ]$Gene_name[reOrd==y]
    slopes.sort <- sort(slopes)
    for(s in names(slopes.sort)){
      subOrd[iter] <- s
      iter <- iter - 1
    }
  }
}
labrow <- as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(BP[,1]) , ]$Gene_name)
ordOrd <- list()
for(l in 1:length(labrow)){
  for(q in 1:length(unlist(subOrd))){
    if( labrow[l] == unlist(subOrd)[q] ){
      ordOrd[l] <- q
    }
  }
}
ordOrd <- unlist(ordOrd)

tq <- log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(BP[,1]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
rownames(tq) <- labrow

BP_matrix.sort <- BP_matrix[  unlist(subOrd), ] 
tq.sort <- tq[ unlist(subOrd) , colOrder ] 

png(filename="BiologicalProcesses_SORT.png", width = 1400, height = 800)
heatmap.plus(  scale( tq.sort ) ,
               col=redblue(150),#c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="\nBiological Pathways",
               mar=c(20,5),
               ColSideColors = as.matrix(as.data.frame(md.col)[as.character(md.order),]), #as.matrix(md.col),
               RowSideColors = BP_matrix.sort, # CC_matrix, #t(CC_matrix),
               cexCol=0.9,
               Rowv = NA,
               Colv = NA
)

dev.off()

png(filename="BiologicalProcesses_LEGEND.png", width = 500, height = 500)
plot.new()
legend("top",
       legend = paste(c(1:16),"=",as.character(GBP[,2])[1:16]),
       #inset=c(-0.05,-0.05),
       cex=1
)
legend("bottom",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=1
)

dev.off()



#####################################################
# MP #

MP_matrix <- matrix(ncol=6,nrow=sum(as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1])))
rownames(MP_matrix) <- as.character( HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , ]$Gene_name )
colnames(MP_matrix) <- as.character(GMP[1:6,2])
MP_matrix[is.na(MP_matrix)] <- "#FFFFFF"
MP_color=primary.colors(6)
for( q in 1:6 ){
  if(q<4){
    filename <- paste0("MF",q,".csv")
  }else{
    filename <- paste0("MP",q,".csv")
  }
  infile <-  read.csv(filename)
  MP_matrix[unique(as.character(infile[,1])),q] <- MP_color[q]
}
colnames(MP_matrix) <- c("1","2       ",
                         "3","4       ",
                         "5","6       ")

png(filename="MolecularFunctions.png", width = 1000, height = 1000)
heatmap.plus(  log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01) ,
               col=c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="Molecular Functions",
               mar=c(20,5),
               labRow=as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , ]$Gene_name),
               #ColSideColors = t(as.matrix(c("blue","purple")[as.factor(grepl("Healthy", h.orderedColumnIDs))])),
               ColSideColors = as.matrix(md.col),
               RowSideColors = MP_matrix, #t(CC_matrix),
               cexCol=0.9
)
legend("bottomleft",
       legend = paste(c(1:16),"=",as.character(GMP[,2])[1:16]),
       #inset=c(-0.05,-0.05),
       cex=0.6
)
legend("topleft",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=0.6
)
dev.off()

# and sort (MP) # 


reOrd <- as.numeric(MP_matrix[,1] != "#FFFFFF")
names(reOrd) <- HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , ]$Gene_name

# 2
fill <- 2
for(z in c(2:6)){
  
  r <- reOrd[ names(MP_matrix[ MP_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF")[MP_matrix[ MP_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF"] ]
  reOrd[ names(MP_matrix[ MP_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF")[MP_matrix[ MP_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF"] ][ r >= (fill-1) ] <- fill
  tmp <- as.numeric(MP_matrix[,z][reOrd == 0] != "#FFFFFF" )
  fill <- fill + 1
  
  tmp[tmp==1] <- fill
  reOrd[ ! reOrd %in% c(1:as.numeric(z-1)) ][ reOrd[ ! reOrd %in% c(1:as.numeric(z-1)) ]==0 ] <- tmp
  fill <- fill + 1
  
}

# name it

subOrd <- list()
iter <- sum(as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]))
for(y in max(reOrd):1){
  if(sum(reOrd==y) == 1){
    subOrd[iter] <- names(reOrd)[ reOrd==y ]
    iter <- iter - 1
  }
  else if(sum(reOrd==y) > 0){
    tq <- log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
    pvals <- list()
    slopes <- list()
    p <- 1
    for( t in 1:dim( tq[reOrd==y,] )[1] ){
      lm_reg <- lm(   as.numeric(tq[reOrd==y,][t,]) ~ as.factor(grepl("Healthy", colnames(tq))))
      fstat  <- summary(lm_reg)$fstatistic
      pvals[p]   <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
      slopes[p]  <- as.numeric(as.matrix(coef(lm_reg))[2,1])
      p <- p + 1
    }
    slopes <- as.numeric(slopes)
    names(slopes) <-  HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , ]$Gene_name[reOrd==y]
    slopes.sort <- sort(slopes)
    for(s in names(slopes.sort)){
      subOrd[iter] <- s
      iter <- iter - 1
    }
  }
}
labrow <- as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , ]$Gene_name)
ordOrd <- list()
for(l in 1:length(labrow)){
  for(q in 1:length(unlist(subOrd))){
    if( labrow[l] == unlist(subOrd)[q] ){
      ordOrd[l] <- q
    }
  }
}
ordOrd <- unlist(ordOrd)

tq <- log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
rownames(tq) <- labrow

MP_matrix.sort <- MP_matrix[  unlist(subOrd), ] 
tq.sort <- tq[ unlist(subOrd) , colOrder ] 

png(filename="MolecularFunctions_SORT.png", width = 800, height = 800)
heatmap.plus(  scale(tq.sort) ,
               col=redblue(150),#c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="\nMolecular Functions",
               mar=c(20,5),
               ColSideColors = as.matrix(as.data.frame(md.col)[as.character(md.order),]), #as.matrix(md.col),
               RowSideColors = MP_matrix.sort, # CC_matrix, #t(CC_matrix),
               cexCol=0.9,
               Rowv = NA,
               Colv = NA
)

dev.off()

png(filename="MolecularFunctions_LEGEND.png", width = 500, height = 500)
plot.new()
legend("top",
       legend = paste(c(1:6),"=",as.character(GMP[,2])[1:6]),
       #inset=c(-0.05,-0.05),
       cex=1
)
legend("bottom",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=1
)

dev.off()



######################################################################
# KD # 

  uniqKD <- as.character(unique(KD[1:12,1]))
  KD_matrix <- matrix(ncol=10,nrow=11)
  shared <- unlist(as.character(HvsT2D.all.out.annotated$Gene_name [ HvsT2D.all.out.annotated$Gene_name %in% as.character(KD[,5][1:12]) ]))
  rownames(KD_matrix) <- shared
  colnames(KD_matrix) <- uniqKD
  KD_matrix[is.na(KD_matrix)] <- "#FFFFFF"
  KD_color=primary.colors(10)

for(k in 1:12){
  if(as.character(KD[k,5]) == "ABCCA1"){  KD[k,5] <- "ABCA1" }
  KD_matrix[ as.character(KD[k,5]) , uniqKD == KD[k,1] ] <- KD_color[uniqKD == KD[k,1]]
  #KD_matrix[ as.character(KD[k,5]) , k ] <- KD_color[k]
} 
  
colnames(KD_matrix) <- c("1","2       ",
                         "3","4       ",
                         "5","6       ",
                         "7","8       ",
                         "9","10      ")

#png(filename="MolecularFunctions.png", width = 1000, height = 1000)

tq <- log(as.matrix(HvsT2D.all.out.annotated[ HvsT2D.all.out.annotated$Gene_name %in% shared , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
rownames(tq) <- as.character(HvsT2D.all.out.annotated$Gene_name [ HvsT2D.all.out.annotated$Gene_name %in% shared ])
  
png(filename="Diseases.png", width = 1000, height = 1000)

heatmap.plus(  scale(tq),
               col=redblue(150),#c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="Diseases",
               mar=c(20,8),
               #labRow=as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , ]$Gene_name),
               #ColSideColors = t(as.matrix(c("blue","purple")[as.factor(grepl("Healthy", h.orderedColumnIDs))])),
               ColSideColors = as.matrix(md.col),
               RowSideColors = KD_matrix, #t(CC_matrix),
               #Rowv = NA,
               #Colv = NA,
               cexCol=0.9
)
 
legend("bottomleft",
       legend = paste(c(1:6),"=",as.character(GMP[,2])[1:6]),
       #inset=c(-0.05,-0.05),
       cex=0.6
)
legend("topleft",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=0.6
)
dev.off()

# and sort 

#MP_matrix.sort <- KD_matrix[  unlist(subOrd), ] 
tq.sort <- tq[ , colOrder ]  
 
png(filename="Diseases_SORT.png", width = 1000, height = 1000)

heatmap.plus(  scale(tq.sort),
               col=redblue(150),#c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="Diseases",
               mar=c(20,8),
               #labRow=as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , ]$Gene_name),
               #ColSideColors = t(as.matrix(c("blue","purple")[as.factor(grepl("Healthy", h.orderedColumnIDs))])),
               ColSideColors = as.matrix(as.data.frame(md.col)[as.character(md.order),]),
               RowSideColors = KD_matrix, #t(CC_matrix),
               Rowv = NA,
               Colv = NA,
               cexCol=0.9
)


dev.off()

png(filename="MolecularFunctions_LEGEND.png", width = 700, height = 500)
plot.new()
legend("top",
       legend = paste(c(1:10),"=",as.character(uniqKD)),
       #cex=0.7,
       #inset=c(-0.05,-0.05),
       cex=1
)
legend("bottom",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=1
)

dev.off()


png(filename="Diseases_LEGEND.png", width = 600, height = 500)
plot.new()
legend("top",
       legend = paste(c(1:10),"=",as.character(uniqKD)),
       #inset=c(-0.05,-0.05),
       cex=0.7
)
legend("bottom",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=1
)

dev.off()



######################################################################
# BIOLOGICAL PATHWAYS # 



uniqKB <- as.character(unique(KB[1:53,1]))
shared <- unlist(as.character(HvsT2D.all.out.annotated$Gene_name [ HvsT2D.all.out.annotated$Gene_name %in% as.character(KB[,3][1:53]) ]))

KB_matrix <- matrix(ncol=length(uniqKB),nrow=length(shared))
rownames(KB_matrix) <- shared
colnames(KB_matrix) <- uniqKB
KB_matrix[is.na(KB_matrix)] <- "#FFFFFF"
KB_color=primary.colors(8)

for(k in 1:53){
  #if(as.character(KB[k,5]) == "ABCCA1"){  KB[k,3] <- "ABCA1" }
  KB_matrix[ as.character(KB[k,3]) , uniqKB == KB[k,1] ] <- KB_color[uniqKB == KB[k,1]]
} 

colnames(KB_matrix) <- c("1","2       ",
                         "3","4       ",
                         "5","6       ",
                         "7","8       ")

tq <- log(as.matrix(HvsT2D.all.out.annotated[ HvsT2D.all.out.annotated$Gene_name %in% shared , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
rownames(tq) <- as.character(HvsT2D.all.out.annotated$Gene_name [ HvsT2D.all.out.annotated$Gene_name %in% shared ])

png(filename="BiologicalPathways.png", width = 1400, height = 1000)

heatmap.plus(  scale(tq),
               col=redblue(150),#c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="Biological Pathways",
               mar=c(20,8),
               #labRow=as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , ]$Gene_name),
               #ColSideColors = t(as.matrix(c("blue","purple")[as.factor(grepl("Healthy", h.orderedColumnIDs))])),
               ColSideColors = as.matrix(md.col),
               RowSideColors = KB_matrix, #t(CC_matrix),
               #Rowv = NA,
               #Colv = NA,
               cexCol=0.9
)

legend("bottomleft",
       legend = paste(c(1:8),"=",uniqKB),
       #inset=c(-0.05,-0.05),
       cex=0.6
)
legend("topleft",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=0.6
)
dev.off()


# and sort 


reOrd <- as.numeric(KB_matrix[,1] != "#FFFFFF")
names(reOrd) <- HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(KB[,3]) , ]$Gene_name

# 2
fill <- 2
for(z in c(2:6)){
  
  r <- reOrd[ names(KB_matrix[ KB_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF")[KB_matrix[ KB_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF"] ]
  reOrd[ names(KB_matrix[ KB_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF")[KB_matrix[ KB_matrix[,1] != "#FFFFFF" ,z] != "#FFFFFF"] ][ r >= (fill-1) ] <- fill
  tmp <- as.numeric(KB_matrix[,z][reOrd == 0] != "#FFFFFF" )
  fill <- fill + 1
  
  tmp[tmp==1] <- fill
  reOrd[ ! reOrd %in% c(1:as.numeric(z-1)) ][ reOrd[ ! reOrd %in% c(1:as.numeric(z-1)) ]==0 ] <- tmp
  fill <- fill + 1
  
}

# name it

subOrd <- list()
iter <- sum(as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(KB[,3]))
for(y in max(reOrd):1){
  if(sum(reOrd==y) == 1){
    subOrd[iter] <- names(reOrd)[ reOrd==y ]
    iter <- iter - 1
  }
  else if(sum(reOrd==y) > 0){
    tq <- log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(KB[,3]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
    pvals <- list()
    slopes <- list()
    p <- 1
    for( t in 1:dim( tq[reOrd==y,] )[1] ){
      lm_reg <- lm(   as.numeric(tq[reOrd==y,][t,]) ~ as.factor(grepl("Healthy", colnames(tq))))
      fstat  <- summary(lm_reg)$fstatistic
      pvals[p]   <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
      slopes[p]  <- as.numeric(as.matrix(coef(lm_reg))[2,1])
      p <- p + 1
    }
    slopes <- as.numeric(slopes)
    names(slopes) <-  HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(KB[,3]) , ]$Gene_name[reOrd==y]
    slopes.sort <- sort(slopes)
    for(s in names(slopes.sort)){
      subOrd[iter] <- s
      iter <- iter - 1
    }
  }
}
labrow <- as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(KB[,3]) , ]$Gene_name)
ordOrd <- list()
for(l in 1:length(labrow)){
  for(q in 1:length(unlist(subOrd))){
    if( labrow[l] == unlist(subOrd)[q] ){
      ordOrd[l] <- q
    }
  }
}
ordOrd <- unlist(ordOrd)

tq <- log(as.matrix(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(KB[,3]) , colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
rownames(tq) <- labrow

KB_matrix.sort <- KB_matrix[  unlist(subOrd), ] 
tq.sort <- tq[ unlist(subOrd) , colOrder ] 

png(filename="BiologicalPathways_SORT.png", width = 800, height = 800)
heatmap.plus(  scale(tq.sort) ,
               col=redblue(150),#c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="\nBiological Pathways",
               mar=c(20,5),
               ColSideColors = as.matrix(as.data.frame(md.col)[as.character(md.order),]), #as.matrix(md.col),
               RowSideColors = KB_matrix.sort, # CC_matrix, #t(CC_matrix),
               cexCol=0.9,
               Rowv = NA,
               Colv = NA
) 

dev.off()

png(filename="BiologicalPathways_LEGEND.png", width = 500, height = 500)
plot.new()
legend("top",
       legend = paste(c(1:8),"=",uniqKB),
       #inset=c(-0.05,-0.05),
       cex=1
)
legend("bottom",
       legend = paste(c("a","b","c","d","e"),"=",colnames(metadata[,c(2,6,8,9,10)])),
       #inset=c(-0.05,-0.05),
       cex=1
)

dev.off()























########################################################################################################################
########### (1) VOLCANO PLOTS ##########################################################################################
######################################################################################################################## 
library(vioplot)
 
#l <- log(as.matrix(HvsT2D.all.out.annotated[1:25,4:47])+0.0001)
#l <- as.matrix(HvsT2D.all.out.annotated[1:25,4:47])
l <- as.matrix(HvsT2D.all.out.annotated[,4:47])
#rownames(l) <- as.character(HvsT2D.all.out.annotated[1:25,]$Gene_name)
rownames(l) <- as.character(HvsT2D.all.out.annotated$Gene_name)

# png(filename="VioPlot_HealthyVsT2D_WIDE.png", width = 1200, height = 240)
# 
# plot.new()
# par(mfrow=c(1,25))
# max <- 0
# min <- 0
# for(r in 1:25){
#   par(mar=c(0,0,0,0))
#   vioplot( t(l)[grepl("Healthy",colnames(l)),r], t(l)[!grepl("Healthy",colnames(l)),r],
#           col = c("cyan3","darkred"),
#           names=c("Healthy","T2D"),
#           ylim=c(min(l),max(l)),
#           drawRect = TRUE,
#           border = FALSE
#           #axis=FALSE
#           )
#   title(main=paste0("\n",rownames(l)[r]))
# }
#    
# dev.off()



# for(u in unique(perturbation$key)){
# 
#   if( length(unique(perturbation[ perturbation$key==u , ]$FDR)) > 1 ){
#     
#     nam <- as.character(perturbation[ perturbation$key==u ,][1:25,]$Gene_name)
#     l <- log(as.matrix(HvsT2D.all.out.annotated[as.character(HvsT2D.all.out.annotated$Gene_name) %in% nam,4:47])+0.0001)
#     rownames(l) <- as.character(HvsT2D.all.out.annotated$Gene_name)[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% nam ]
#       
#     png(filename=paste0("VioPlot_",u,"_HealthyVsT2D.png"), width = 1200, height = 240)
#     
#     plot.new()
#     par(mfrow=c(1,25))
#     for(r in 1:25){
#       par(mar=c(0,0,0,0))
#       vioplot( t(l)[grepl("Healthy",colnames(l)),r], t(l)[!grepl("Healthy",colnames(l)),r],
#                col = "gray",
#                names=c("Healthy","T2D"),1
#                ylim=c(min(l),max(l)),
#                drawRect = TRUE,
#                border = FALSE
#       )
#       title(main=paste0("\n",rownames(l)[r]))
#     }
#     
#     dev.off()
#     
#   }
#   
# }

####
# l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , ]
####
# 
# 
# png(filename="VioPlot_HealthyVsT2D_WIDE.png", width = 1200, height = 240)
# 
# plot.new()
# par(mfrow=c(1,25))
# max <- 0
# min <- 0
# for(r in 1:25){
#   par(mar=c(0,0,0,0))
#   vioplot( t(l)[grepl("Healthy",colnames(l)),r], t(l)[!grepl("Healthy",colnames(l)),r],
#            col = c("cyan3","darkred"),
#            names=c("Healthy","T2D"),
#            ylim=c(min(l),max(l)),
#            drawRect = TRUE,
#            border = FALSE
#            #axis=FALSE
#   )
#   title(main=paste0("\n",rownames(l)[r]))
# }
# 
# dev.off()
# 
# png(filename="spreadsheetorder25.groupscaled.png", width = 320, height = 1205)
# plot.new()
# par( mfrow=c(25,1) )
# par( mar = c(0,0,0,0) )
# miny1 <- 9999
# miny2 <- -9999
# for(q in 1:25){
#   y1 <- min( t( l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , grepl("Healthy",colnames(l)) ] )[,q] ,t( l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , !grepl("Healthy",colnames(l)) ] )[,1] )
#   if(y1 < miny1){ miny1 <- y1 }
#   y2 <- max( t( l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , grepl("Healthy",colnames(l)) ] )[,q] ,t( l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , !grepl("Healthy",colnames(l)) ] )[,1] )
#   if(y2 > miny2){ miny2 <- y2 }
# }
# for(q in 1:25){
#   par(mar=c(0,14,0,0))
#   vioplot( t( l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , grepl("Healthy",colnames(l)) ] )[,q] ,
#            t( l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , !grepl("Healthy",colnames(l)) ] )[,q],
#            col = "darkslategray4", border = FALSE, names = NA, ylim = c(miny1,miny2)
#           )
#   legend("left",inset=c(-0.3,0),
#          legend = rownames(l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , grepl("Healthy",colnames(l)) ])[q],
#          border="white", bty = "n", fill=NULL, bg="white",
#          cex=1.2)
# }
# dev.off()
#  
# png(filename="spreadsheetorder25.png", width = 320, height = 1205)
# plot.new()
# par( mfrow=c(25,1) )
# par( mar = c(0,14,0,0) )
# par(xpd=TRUE) 
# for(q in 1:25){
#   set.seed(1) # just to get the same random numbers
#   vioplot( t( l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , grepl("Healthy",colnames(l)) ] )[,q] ,
#            t( l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , !grepl("Healthy",colnames(l)) ] )[,q],
#            col = "darkslategray4", border = FALSE, names = NA#, #ylim = c(miny1,miny2)
#            #, add=TRUE
#   )
#   legend("left",inset=c(-0.7,0),
#          legend = rownames(l[ as.character(HvsT2D.all.out.annotated$Gene_name[1:1000][as.character(HvsT2D.all.out.annotated$Gene_name[1:1000]) %in% rownames(l)]) , grepl("Healthy",colnames(l)) ])[q],
#          border="white", bty = "n", fill=NULL, bg="white",
#          cex=1.2)
# }
# dev.off()
#  

HvsT2D.all.out.annotated.rel <- HvsT2D.all.out.annotated
HvsT2D.all.out.annotated.rel[26:47] <- t(decostand(t(HvsT2D.all.out.annotated.rel[26:47]), method="total"))
l <- HvsT2D.all.out.annotated.rel[HvsT2D.all.out.annotated.rel$Gene_name %in% names(sort(pvals.fdr))[1:25],c(1:3,26:47)]
uniq.l_names <- as.character(unique(l$Gene_name))
l.uniq <- matrix(ncol=22,nrow=length(as.character(unique(l$Gene_name))))
colnames(l.uniq)
rownames(l.uniq) <- uniq.l_names
for(u in uniq.l_names){
  l.uniq[ u , ] <-  colSums(l[l$Gene_name == u,4:25])
}
colnames(l.uniq) <- colnames(HvsT2D.all.out.annotated.rel[26:47])

# foresight's ranks
# top 25 in T2DvsHealthy.M.W.corrected.selected

derek_top25 <- c("NECTIN2","LILRB4","DST","INPP4B","VWF","CRYBG3","EMP1","CAMK2D","SERPINB2","ZNF204P","KBTBD8","CD86","CTSL","GREM2","GPRIN3","IGHG2","EPB41L3","PID1","ANKRD36BP2","CD300E","VCAN","MALT1","FAM153B","CYP1B1","LRP1")
l <- HvsT2D.all.out.annotated.rel[HvsT2D.all.out.annotated.rel$Gene_name %in% derek_top25,c(26:47)]
rownames(l) <- as.character(HvsT2D.all.out.annotated.rel[HvsT2D.all.out.annotated.rel$Gene_name %in% derek_top25,"Gene_name"]) 

plot.new()
par( mfrow=c(25,1) )
par( mar = c(0,8,0,0) )
par(xpd=TRUE) 
for(q in 1:25){
  set.seed(1) # just to get the same random numbers
  #vioplot( t( l[ q , grepl("Healthy",colnames(l)) ] )  ,
  #         t( l[ q , !grepl("Healthy",colnames(l)) ] ),
  #         col = "darkslategray4", border = FALSE, names = NA#, #ylim = c(miny1,miny2)
  #         #, add=TRUE
  #)
  #legend("left",inset=c(-0.4,0),
  #       legend = rownames(l)[q],
  #      border="white", bty = "n", fill=NULL, bg="white",
  #      cex=1)
  png(filename = paste0("beanplot.foresite.",q,".",rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q],".png"))
  beanplot( as.numeric(t( l[ q , grepl("Healthy",colnames(l)) ] ))[ !is.na(as.numeric(t( l[ q , grepl("Healthy",colnames(l)) ] ))) ]  ,
            as.numeric(t( l[ q , !grepl("Healthy",colnames(l)) ]))[ !is.na(as.numeric(t( l[ q , !grepl("Healthy",colnames(l)) ]))) ],
            bw="nrd0" ,
            col = list(c("darkgreen","green","black","gray3"),c("red","black","black","gray3")), names = c("Healthy","T2D") #,
            #main=rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q]
  )
  dev.off()
}

miny1 <- 9999
miny2 <- -9999
for(q in 1:25){
  y1 <- min( l[ q               , grepl("Healthy",colnames(l.uniq)) ] )
  if(y1 < miny1){ miny1 <- y1 }
  y2 <- max( l[ q               , grepl("Healthy",colnames(l.uniq)) ] )
  if(y2 > miny2){ miny2 <- y2 }
}
plot.new()
par( mfrow=c(25,1) )
par( mar = c(0,8,0,0) )
par(xpd=TRUE) 
for(q in 1:25){
  # set.seed(1) # just to get the same random numbers
  # vioplot( t( l[ q , grepl("Healthy",colnames(l)) ] )  ,
  #          t( l[ q , !grepl("Healthy",colnames(l)) ] ),
  #          col = "darkslategray4", border = FALSE, names = NA, ylim = c(miny1,miny2)
  #          #, add=TRUE
  # )
  # legend("left",inset=c(-0.4,0),
  #        legend = rownames(l)[q],
  #        border="white", bty = "n", fill=NULL, bg="white",
  #        cex=1)
  
  png(filename = paste0("beanplot.foresite.ylim.",q,".",rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q],".png"))
  beanplot( as.numeric(t( l[ q , grepl("Healthy",colnames(l)) ] ))[ !is.na(as.numeric(t( l[ q , grepl("Healthy",colnames(l)) ] ))) ]  ,
            as.numeric(t( l[ q , !grepl("Healthy",colnames(l)) ]))[ !is.na(as.numeric(t( l[ q , !grepl("Healthy",colnames(l)) ]))) ],
            bw="nrd0" ,
            col = list(c("darkgreen","green","black","gray3"),c("red","black","black","gray3")), names = c("Healthy","T2D"), ylim = c(miny1,miny2)  #,
            #main=rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q]
  )
  dev.off()

}


# compare to jamisons

load("pvals.Rdata")

#top43 <- HvsT2D.all.out.annotated[ HvsT2D.all.out.annotated$Ensembl_GeneID %in% names(sort(pvals.fdr))[1:44], ]
# names(sort(pvals.fdr))[1:44] %in% as.character(HvsT2D.all.out.annotated$Gene_name)

HvsT2D.all.out.annotated.rel <- HvsT2D.all.out.annotated
HvsT2D.all.out.annotated.rel[26:47] <- t(decostand(t(HvsT2D.all.out.annotated.rel[26:47]), method="total"))
l <- HvsT2D.all.out.annotated.rel[HvsT2D.all.out.annotated.rel$Gene_name %in% names(sort(pvals.fdr))[1:25],c(1:3,26:47)]
uniq.l_names <- as.character(unique(l$Gene_name))
l.uniq <- matrix(ncol=22,nrow=length(as.character(unique(l$Gene_name))))
colnames(l.uniq)
rownames(l.uniq) <- uniq.l_names
for(u in uniq.l_names){
  l.uniq[ u , ] <-  colSums(l[l$Gene_name == u,4:25])
}
colnames(l.uniq) <- colnames(HvsT2D.all.out.annotated.rel[26:47])
plot.new()
#par( mfrow=c(5,2) )
par( mfrow=c(1,1) )
par( mar = c(2,8,0,0) )
par(xpd=TRUE) 
for(q in 1:25){
  #set.seed(1) # just to get the same random numbers
  #vioplot( l.uniq[ names(sort(pvals.fdr))[q]               , grepl("Healthy",colnames(l.uniq)) ], 
  #         l.uniq[ names(sort(pvals.fdr))[q]               , !grepl("Healthy",colnames(l.uniq)) ],
  #         col = "darkslategray4", border = FALSE, names = NA#, #ylim = c(miny1,miny2)
  #         #, add=TRUE
  #)
  #legend("left",inset=c(-0.4,0),
  #       legend = rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q],
  #       border="white", bty = "n", fill=NULL, bg="white",
  #       cex=1)
  
  #par( mfrow=c(3,1) )
  #beanplot(c(1:20))
  png(filename = paste0("beanplot.",q,".",rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q],".png"))
  beanplot( l.uniq[ names(sort(pvals.fdr))[q]               , grepl("Healthy",colnames(l.uniq)) ] , 
            l.uniq[ names(sort(pvals.fdr))[q]               , !grepl("Healthy",colnames(l.uniq))] ,
            bw="nrd0",
            col = list(c("darkgreen","green","black","gray3"),c("red","black","black","gray3")), names = c("Healthy","T2D"),
            main=rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q]
          )
  #legend("left",inset=c(-0.25,0),
  #       legend = rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q],
  #       border="white", bty = "n", fill=NULL, bg="white",
  #       cex=1)
  dev.off()
  
}


miny1 <- 9999
miny2 <- -9999
for(q in 1:25){
  y1 <- min( l.uniq[ names(sort(pvals.fdr))[q]               , grepl("Healthy",colnames(l.uniq)) ] )
  if(y1 < miny1){ miny1 <- y1 }
  y2 <- max( l.uniq[ names(sort(pvals.fdr))[q]               , grepl("Healthy",colnames(l.uniq)) ] )
  if(y2 > miny2){ miny2 <- y2 }
}
plot.new()
par( mfrow=c(25,1) )
par( mar = c(0,8,0,0) )
par(xpd=TRUE) 
for(q in 1:25){
  #set.seed(1) # just to get the same random numbers
  #vioplot( l.uniq[ names(sort(pvals.fdr))[q]               , grepl("Healthy",colnames(l.uniq)) ], 
  #         l.uniq[ names(sort(pvals.fdr))[q]               , !grepl("Healthy",colnames(l.uniq)) ],
  #         col = "darkslategray4", border = FALSE, names = NA, ylim = c(miny1,miny2)
  #         #, add=TRUE
  #)
  #legend("left",inset=c(-0.4,0),
  #       legend = rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q],
  #       border="white", bty = "n", fill=NULL, bg="white",
  #       cex=1)
  png(filename = paste0("beanplot.ylim.",q,".",rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q],".png"))
  beanplot( l.uniq[ names(sort(pvals.fdr))[q]               , grepl("Healthy",colnames(l.uniq)) ] , 
            l.uniq[ names(sort(pvals.fdr))[q]               , !grepl("Healthy",colnames(l.uniq))] ,
            bw="nrd0",
            col = list(c("darkgreen","green","black","gray3"),c("red","black","black","gray3")), names = c("Healthy","T2D"),
            main=rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q], ylim=c(miny1,miny2)
  )
  #legend("left",inset=c(-0.25,0),
  #       legend = rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q],
  #       border="white", bty = "n", fill=NULL, bg="white",
  #       cex=1)
  dev.off()
}




# ! # 
library(beanplot)
plot.new()
par( mfrow=c(25,1) )
par( mar = c(0,8,0,0) )
par(xpd=TRUE) 
for(q in 1:25){
  #set.seed(1) # just to get the same random numbers
  #beanplot( l.uniq[ names(sort(pvals.fdr))[q]               , grepl("Healthy",colnames(l.uniq)) ], 
  #         l.uniq[ names(sort(pvals.fdr))[q]               , !grepl("Healthy",colnames(l.uniq)) ],
  #         col = "darkslategray4", border = FALSE, names = NA#, #ylim = c(miny1,miny2)
  #         #, add=TRUE
  #)
  #legend("left",inset=c(-0.4,0),
  #       legend = rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q],
  #       border="white", bty = "n", fill=NULL, bg="white",
  #       cex=1)
  
  png(filename = paste0("beanplot.",rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q],".png"))
  beanplot( l.uniq[ names(sort(pvals.fdr))[q]               , grepl("Healthy",colnames(l.uniq)) ] , 
            l.uniq[ names(sort(pvals.fdr))[q]               , !grepl("Healthy",colnames(l.uniq))] ,
            bw="nrd0",
            col = list(c("darkgreen","green","black","gray3"),c("red","black","black","gray3")), names = c("Healthy","T2D"),
            main=rownames( l.uniq[ names(sort(pvals.fdr)[1:25]) %in% rownames(l.uniq) , grepl("Healthy",colnames(l.uniq)) ])[q]
  )
  dev.off()
  
}




# CellMobility
# 
CellMobility <- read.table("CellMobility.tsv", sep="\t", header=FALSE, fill=TRUE)
CellMobility[,2]
 
tq <- log(as.matrix(HvsT2D.all.out.annotated[ as.character( HvsT2D.all.out.annotated$Gene_name) %in% unique(toupper(unlist(as.character(CellMobility[,2])))) , colOrder ])+0.001)  #colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
rownames(tq) <- as.character(HvsT2D.all.out.annotated$Gene_name [ as.character( HvsT2D.all.out.annotated$Gene_name) %in% unique(toupper(unlist(as.character(CellMobility[,2])))) ])
tq <- tq[ rowSums(tq) > -130 , ] 

png(filename="CellMobilityGT130.png", width = 1000, height = 1400)

#heatmap.plus(  scale(tq),
heatmap.2( scale(tq),
               col=redblue(150),#c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="Cell Mobility",
               mar=c(20,8),
               #labRow=as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , ]$Gene_name),
               ColSideColors = t(as.matrix(c("red","green")[as.factor(grepl("Healthy", colOrder))])),
#               ColSideColors = as.matrix(md.col),
#               RowSideColors = KB_matrix, #t(CC_matrix),
               #Rowv = NA,
               Colv = FALSE, dendrogram='row', 
               cexCol=0.9
)
 
legend("topright",
       legend = c("T2D","Healthy"),
       inset=c(0.08,0.08),
       pch=24,
       col=c("red","green"),
       cex=0.6
)
dev.off()

tq <- log(as.matrix(HvsT2D.all.out.annotated[ as.character( HvsT2D.all.out.annotated$Gene_name) %in% unique(toupper(unlist(as.character(CellMobility[,2])))) , colOrder ])+0.001)  #colnames(HvsT2D.all.out.annotated)[4:47][23:44] ])+0.01)
rownames(tq) <- as.character(HvsT2D.all.out.annotated$Gene_name [ as.character( HvsT2D.all.out.annotated$Gene_name) %in% unique(toupper(unlist(as.character(CellMobility[,2])))) ])
tq <- tq[ rowSums(tq) > -140 , ] 

png(filename="CellMobilityGT140.png", width = 1000, height = 1400)

#heatmap.plus(  scale(tq),
heatmap.2( scale(tq),
           col=redblue(150),#c(rev(matlab.like2(150))),#,gray.colors(150)),
           main="Cell Mobility",
           mar=c(20,8),
           #labRow=as.character(HvsT2D.all.out.annotated[ as.character(HvsT2D.all.out.annotated$Gene_name) %in% as.character(MP[,1]) , ]$Gene_name),
           ColSideColors = t(as.matrix(c("red","green")[as.factor(grepl("Healthy", colOrder))])),
           #               ColSideColors = as.matrix(md.col),
           #               RowSideColors = KB_matrix, #t(CC_matrix),
           #Rowv = NA,
           Colv = FALSE, dendrogram='row', 
           cexCol=0.9
)

legend("topright",
       legend = c("T2D","Healthy"),
       inset=c(0.08,0.08),
       pch=24,
       col=c("red","green"),
       cex=0.6
)
dev.off()

# 













save.image(file="20181012star.Rdata")













derekPvals <- read.table("T2DvsHealthy.M.W.corrected.selected.txt", sep="\t", header=TRUE, fill=TRUE)
l <- HvsT2D.all.out.annotated.rel[HvsT2D.all.out.annotated.rel$Gene_name %in% derekPvals$Gene_name[1:50],c(1:3,26:47)]
uniq.l_names <- as.character(unique(l$Gene_name))
l.uniq <- matrix(ncol=22,nrow=length(as.character(unique(l$Gene_name))))
colnames(l.uniq) <- colnames(HvsT2D.all.out.annotated.rel)[26:47]
rownames(l.uniq) <- uniq.l_names
for(u in uniq.l_names){
  l.uniq[ u , ] <-  colSums(l[l$Gene_name == u,4:25])
}
colnames(l.uniq) <- colnames(HvsT2D.all.out.annotated.rel[26:47])

h <- heatmap.plus( scale(log(l.uniq+0.0001)))
heatmap.2(     scale(log(l.uniq+0.0001)) ,
               col=redblue(200),#c(rev(matlab.like2(150))),#,gray.colors(150)),
               main="\nDerek Significant",
               mar=c(20,5),
               ColSideColors = c("blue","purple")[as.factor(grepl("Healthy", h.orderedColumnIDs))],
               #ColSideColors = as.matrix(as.data.frame(md.col)[as.character(md.order),]), #as.matrix(md.col),
               #RowSideColors = KB_matrix.sort, # CC_matrix, #t(CC_matrix),
               cexCol=0.9,
               trace = "none"
               #Rowv = NA,
               #Colv = NA
) 

l <- HvsT2D.all.out.annotated.rel[HvsT2D.all.out.annotated.rel$Gene_name %in% derekPvals$Gene_name,c(1:3,26:47)]
uniq.l_names <- as.character(unique(l$Gene_name))
l.uniq <- matrix(ncol=22,nrow=length(as.character(unique(l$Gene_name))))
colnames(l.uniq) <- colnames(HvsT2D.all.out.annotated.rel)[26:47]
rownames(l.uniq) <- uniq.l_names
for(u in uniq.l_names){
  l.uniq[ u , ] <-  colSums(l[l$Gene_name == u,4:25])
}
plot.new()
par(mfrow=c(1,1))
pca_box_hash <- prcomp(t(l.uniq), center = TRUE) 
biplot(pca_box_hash, main="PCA + Eigenvectors (All Highlighted IP Gene IDs)", main.cex = 0.5, 
       xlim=c(-0.3,1.1),
       cex=0.5,
       col=c("gray","red"),
       
)


PCA <- prcomp(t(l.uniq), scale=T, center = T) 
write.table ( as.table( sort(PCA$center) ), file="derek.scaled_residual_lengths.txt" )
write.table ( as.table( sort(PCA$center) ), file="derek.unscaled_residual_lengths.txt" )
       








l <- HvsT2D.all.out.annotated.rel[ , c(1:3,26:47) ]
uniq.l_names <- as.character(unique(l$Gene_name))
l.uniq <- matrix(ncol=22,nrow=length(as.character(unique(l$Gene_name))))
colnames(l.uniq) <- colnames(HvsT2D.all.out.annotated.rel)[26:47]
rownames(l.uniq) <- uniq.l_names
for(u in uniq.l_names){
  l.uniq[ u , ] <-  colSums(l[l$Gene_name == u,4:25])
}
plot.new()
par(mfrow=c(1,1))
pca_box_hash <- prcomp(t(l.uniq), center = TRUE) 
biplot(pca_box_hash, main="PCA + Eigenvectors (All Highlighted IP Gene IDs)", main.cex = 0.5, 
       xlim=c(-1,0.6),
       cex=0.5,
       col=c("gray","red")
)

PCA <- prcomp(t(l.uniq), scale=T, center = T) 
write.table ( as.table( sort(pca_box_hash$center) ), file="all.scaled_residual_lengths.txt" )

PCA <- prcomp(t(l.uniq), center = T) 
write.table ( as.table( sort(PCA$center) ), file="all.unscaled_residual_lengths.txt" )










## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ##

## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ##

## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ##

## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ## !! ##







plot.new()
par(mfrow=c(6,3))
par(mar=c(2,2,2,2))
for(u in unique(perturbation$key)){
  plot( perturbation[ perturbation$key==u , ]$logFC , -log(perturbation[ perturbation$key==u , ]$PValue),
        main=u,
        xlim = c(-max(abs(perturbation[ perturbation$key==u , ]$logFC)),max(abs(perturbation[ perturbation$key==u , ]$logFC))),
        #ylim = c(0,20)
  )
  if(sum(-log(perturbation[ perturbation$key==u , ]$PValue) >= 6) > 0){
    text( perturbation[ -log(perturbation[ perturbation$key==u , ]$PValue) >= 6 , ][ perturbation[ -log(perturbation[ perturbation$key==u , ]$PValue) >= 6 , ]$key==u , ]$logFC , -log(perturbation[ -log(perturbation[ perturbation$key==u , ]$PValue) >= 6 , ][ perturbation[ -log(perturbation[ perturbation$key==u , ]$PValue) >= 6 , ]$key==u , ]$PValue),
          perturbation[ -log(perturbation[ perturbation$key==u , ]$PValue) >= 6 , ]$Gene_name,
          col="red", cex = 0.5,
          pos=2
    )
  }
  abline(h=6,col="red")
}





# FDR(PVALUE) -- some of these arnet logged!
plot.new()
par(mfrow=c(6,3))
par(mar=c(2,2,2,2))
for(u in unique(perturbation$key)){
  plot( perturbation[ perturbation$key==u , ]$logFC , -log(perturbation[ perturbation$key==u , ]$FDR),
        main=u,
        xlim = c(-max(abs(perturbation[ perturbation$key==u , ]$logFC)),max(abs(perturbation[ perturbation$key==u , ]$logFC))),
        #ylim = c(0,20)
  )
  if(sum(-log(perturbation[ perturbation$key==u , ]$FDR) >= 6) > 0){
    text( perturbation[ -log(perturbation[ perturbation$key==u , ]$FDR) >= 6 , ][ perturbation[ -log(perturbation[ perturbation$key==u , ]$FDR) >= 6 , ]$key==u , ]$logFC , -log(perturbation[ -log(perturbation[ perturbation$key==u , ]$FDR) >= 6 , ][ perturbation[ -log(perturbation[ perturbation$key==u , ]$FDR) >= 6 , ]$key==u , ]$FDR),
          perturbation[ -log(perturbation[ perturbation$key==u , ]$FDR) >= 6 , ]$Gene_name,
          col="red", cex = 0.5,
          pos=2
    )
  }
  abline(h=6,col="red")
  
  yo <- perturbation[ perturbation$key==u , ][ order( perturbation[ perturbation$key==u , ]$FDR ) , ]
  yo[1:25,]
   
}








