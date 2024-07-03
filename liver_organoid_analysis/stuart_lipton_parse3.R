setwd("C:\\Users\\jamis\\Desktop\\minibrains")

all <- read.table("ALLSAMPLES_expression.csv", sep = ",", stringsAsFactors = FALSE)
abund <- all[2:dim(all)[1],2:dim(all)[2]]
rownames(abund) <- all[2:dim(all)[1],1]
colnames(abund) <- as.character(unlist(all[1,2:dim(all)[2]])) 
abund <- sapply( abund , as.numeric )
rownames(abund) <- all[2:dim(all)[1],1]
#remove(all)

abund.rel <- abund
for(y in 1:dim(abund)[2]){
  denom <- sum(abund[,y])
  for(x in 1:dim(abund)[1]){
    abund.rel[x,y] <- abund[x,y]/denom
  }
}
 
labels <- read.table("outs/filtered_feature_bc_matrix/features.tsv", sep = "\t")

barcodes <- as.character(labels[labels[,3] != "Gene Expression",1])
 
dim(abund[barcodes,])

#install.packages("beanplot")
library(beanplot)

abund.rel.viz <- abund.rel
abund.rel.viz[abund.rel.viz == 0]
beanplot(log(abund.rel[barcodes[1],][abund.rel[barcodes[1],] != 0]), 
         log(abund.rel[barcodes[2],][abund.rel[barcodes[2],] != 0]), 
         log(abund.rel[barcodes[3],][abund.rel[barcodes[3],] != 0]), 
         names = barcodes)

sum(log(abund.rel[barcodes[1],][abund.rel[barcodes[1],] != 0]) > -4)
sum(log(abund.rel[barcodes[2],][abund.rel[barcodes[2],] != 0]) > -3.47)
sum(log(abund.rel[barcodes[3],][abund.rel[barcodes[3],] != 0]) > -4)

  
wtype <- names(log(abund.rel[barcodes[1],][abund.rel[barcodes[1],] != 0])[ log(abund.rel[barcodes[1],][abund.rel[barcodes[1],] != 0]) > -4 ])
ps1 <- names(log(abund.rel[barcodes[2],][abund.rel[barcodes[2],] != 0])[ log(abund.rel[barcodes[2],][abund.rel[barcodes[2],] != 0]) > -3.47 ])
app <- names(log(abund.rel[barcodes[3],][abund.rel[barcodes[3],] != 0])[ log(abund.rel[barcodes[3],][abund.rel[barcodes[3],] != 0]) > -4 ])
sum(wtype %in% ps1)
sum(wtype %in% app)
sum(ps1 %in% app)
length( wtype[!wtype %in% ps1][ ! wtype[!wtype %in% ps1] %in% app ] )
length( ps1[!ps1 %in% wtype][ ! ps1[!ps1 %in% wtype] %in% app ] )
length( app[!app %in% wtype][ ! app[!app %in% wtype] %in% ps1 ] )

length(wtype[wtype %in% ps1][ wtype[wtype %in% ps1] %in% app ])

sum_mat <- matrix(nrow=length(c(as.character(labels[labels[,3] == "Gene Expression",1]))),ncol=3)
sum_mat_multi_m <- matrix(nrow=length(c(as.character(labels[labels[,3] == "Gene Expression",1]))),ncol=3)
for(m in c(1:3)){
  qi <- 1
  for(q in c(as.character(labels[labels[,3] == "Gene Expression",1]))){
    reg_1    <- lm( abund.rel[q,] ~ abund.rel[barcodes[m],] + 1 )
    pvalue_1 <- as.matrix(anova(reg_1))[1,5]
    sum_mat[qi, m] <- pvalue_1
    if(m == 1){
      reg_2    <- lm( abund.rel[q,] ~ abund.rel[barcodes[1],] + abund.rel[barcodes[2],] + abund.rel[barcodes[3],] + 1 )
      pvalue_2 <- unlist(as.matrix(anova(reg_2))[1:3,5])
      sum_mat_multi_m[qi,] <- pvalue_2
    }
    qi <- qi + 1
  }
}

# save.image("20190911.am.Rdata")
log(abund.rel[barcodes[1],][abund.rel[barcodes[1],] != 0])
rownames(sum_mat)   <- rownames(sum_mat_multi_m) <- c(as.character(labels[labels[,3] == "Gene Expression",1]))
sum_mat.pop         <- sum_mat[rowSums(is.na(sum_mat)) < 3,]
sum_mat_multi_m.pop <- sum_mat_multi_m[rowSums(is.na(sum_mat_multi_m)) < 3,]

dim(sum_mat.pop) 
dim(sum_mat_multi_m.pop) 

sum_mat.viz <- sum_mat.pop
sum_mat.viz[sum_mat.pop == 0] <- min(sum_mat.pop[sum_mat.pop != 0])

png("pvalues_iterateGandM.png", width=900, height=900)
plot.new()
plot(-log(sort(sum_mat.viz[,1])),col="darkred",ylim=c(0,750), main="-log(pvalues) : relAbund(geneG) ~ relAbund(markerM)")
par(new=TRUE)
plot(-log(sort(sum_mat.viz[,2])),col="darkviolet",ylim=c(0,750))
par(new=TRUE)
plot(-log(sort(sum_mat.viz[,3])),col="blue",ylim=c(0,750))
legend("topright", barcodes, col=c("darkred","darkviolet","blue"), pch=15)
dev.off()

sum_mat.fdr <- sum_mat
sum_mat.fdr[,1] <- p.adjust( sum_mat[,1] , method="fdr" )
sum_mat.fdr[,2] <- p.adjust( sum_mat[,2] , method="fdr" )
sum_mat.fdr[,3] <- p.adjust( sum_mat[,3] , method="fdr" )
sum_mat.fdr.pop <- sum_mat[rowSums(is.na(sum_mat.fdr)) < 3,]
sum_mat.fdr.viz <- sum_mat.fdr.pop
sum_mat.fdr.viz[sum_mat.fdr.pop == 0] <- min(sum_mat.fdr.pop[sum_mat.fdr.pop != 0])

png("pvaluesFDR_iterateGandM.png", width=900, height=900)
plot.new()
plot(-log(sort(sum_mat.fdr.viz[,1])),col="darkred",ylim=c(0,750), main="-log(pvalues), FDR-adjusted : relAbund(geneG) ~ relAbund(markerM)")
par(new=TRUE)
plot(-log(sort(sum_mat.fdr.viz[,2])),col="darkviolet",ylim=c(0,750))
par(new=TRUE)
plot(-log(sort(sum_mat.fdr.viz[,3])),col="blue",ylim=c(0,750))
legend("topright", barcodes, col=c("darkred","darkviolet","blue"), pch=15)
dev.off()



sum_mat_multi_m.viz <- sum_mat_multi_m.pop
sum_mat_multi_m.viz[sum_mat_multi_m.pop == 0] <- min(sum_mat_multi_m.pop[sum_mat_multi_m.pop != 0])

png("pvalues_iterateGmultiM.png", width=900, height=900)
plot.new()
plot(-log(sort(sum_mat_multi_m.viz[,1])),col="darkred",ylim=c(0,750), main="-log(pvalues) : relAbund(geneG) ~ relAbund(marker1) + relAbund(marker2) + relAbund(marker3)")
par(new=TRUE)
plot(-log(sort(sum_mat_multi_m.viz[,2])),col="darkviolet",ylim=c(0,750))
par(new=TRUE)
plot(-log(sort(sum_mat_multi_m.viz[,3])),col="blue",ylim=c(0,750))
legend("topright", barcodes, col=c("darkred","darkviolet","blue"), pch=15)
dev.off()

sum_mat_multi_m.fdr <- sum_mat_multi_m
sum_mat_multi_m.fdr[,1] <- p.adjust( sum_mat_multi_m[,1] , method="fdr" )
sum_mat_multi_m.fdr[,2] <- p.adjust( sum_mat_multi_m[,2] , method="fdr" )
sum_mat_multi_m.fdr[,3] <- p.adjust( sum_mat_multi_m[,3] , method="fdr" )
sum_mat_multi_m.fdr.pop <- sum_mat_multi_m[rowSums(is.na(sum_mat_multi_m.fdr)) < 3,]
sum_mat_multi_m.fdr.viz <- sum_mat_multi_m.fdr.pop
sum_mat_multi_m.fdr.viz[sum_mat_multi_m.fdr.pop == 0] <- min(sum_mat_multi_m.fdr.pop[sum_mat_multi_m.fdr.pop != 0])

png("pvaluesFDR_iterateGmultiM.png", width=900, height=900)
plot.new()
plot(-log(sort(sum_mat_multi_m.fdr.viz[,1])),col="darkred",ylim=c(0,750), main="-log(pvalues), FDR-adjusted : relAbund(geneG) ~ relAbund(marker1) + relAbund(marker2) + relAbund(marker3)")
par(new=TRUE)
plot(-log(sort(sum_mat_multi_m.fdr.viz[,2])),col="darkviolet",ylim=c(0,750))
par(new=TRUE)
plot(-log(sort(sum_mat_multi_m.fdr.viz[,3])),col="blue",ylim=c(0,750))
legend("topright", barcodes, col=c("darkred","darkviolet","blue"), pch=15)
dev.off()


sum( sum_mat_multi_m[,1][ !is.na(sum_mat_multi_m[,1]) ] < 10e-50 )
sum( sum_mat_multi_m[,2][ !is.na(sum_mat_multi_m[,2]) ] < 10e-50 )
sum( sum_mat_multi_m[,3][ !is.na(sum_mat_multi_m[,3]) ] < 10e-50 )
sum( sum_mat_multi_m.fdr[,1][ !is.na(sum_mat_multi_m.fdr[,1]) ] < 10e-50 )
sum( sum_mat_multi_m.fdr[,2][ !is.na(sum_mat_multi_m.fdr[,2]) ] < 10e-50 )
sum( sum_mat_multi_m.fdr[,3][ !is.na(sum_mat_multi_m.fdr[,3]) ] < 10e-50 )


sum( sum_mat[,1][ !is.na(sum_mat[,1]) ] < 10e-50 )
sum( sum_mat[,2][ !is.na(sum_mat[,2]) ] < 10e-50 )
sum( sum_mat[,3][ !is.na(sum_mat[,3]) ] < 10e-50 )
sum( sum_mat.fdr[,1][ !is.na(sum_mat.fdr[,1]) ] < 10e-50 )
sum( sum_mat.fdr[,2][ !is.na(sum_mat.fdr[,2]) ] < 10e-50 )
sum( sum_mat.fdr[,3][ !is.na(sum_mat.fdr[,3]) ] < 10e-50 )

 


sum( sum_mat_multi_m[,1][ !is.na(sum_mat_multi_m[,1]) ] < 10e-100 )
sum( sum_mat_multi_m[,2][ !is.na(sum_mat_multi_m[,2]) ] < 10e-100 )
sum( sum_mat_multi_m[,3][ !is.na(sum_mat_multi_m[,3]) ] < 10e-100 )
sum( sum_mat_multi_m.fdr[,1][ !is.na(sum_mat_multi_m.fdr[,1]) ] < 10e-100 )
sum( sum_mat_multi_m.fdr[,2][ !is.na(sum_mat_multi_m.fdr[,2]) ] < 10e-100 )
sum( sum_mat_multi_m.fdr[,3][ !is.na(sum_mat_multi_m.fdr[,3]) ] < 10e-100 )


sum( sum_mat[,1][ !is.na(sum_mat[,1]) ] < 10e-100 )
sum( sum_mat[,2][ !is.na(sum_mat[,2]) ] < 10e-100 )
sum( sum_mat[,3][ !is.na(sum_mat[,3]) ] < 10e-100 )
sum( sum_mat.fdr[,1][ !is.na(sum_mat.fdr[,1]) ] < 10e-100 )
sum( sum_mat.fdr[,2][ !is.na(sum_mat.fdr[,2]) ] < 10e-100 )
sum( sum_mat.fdr[,3][ !is.na(sum_mat.fdr[,3]) ] < 10e-100 )
 

most_signif_50 <- c( names(sort(sum_mat_multi_m.fdr[,1])[1:25]), 
                     names(sort(sum_mat_multi_m.fdr[,2])[1:25]), 
                     names(sort(sum_mat_multi_m.fdr[,2])[1:25]) )
most_signif_50 <- unique(most_signif_50)

library(gplots)
library(beanplot)
abund.rel.viz[abund.rel.viz == 0] <- min(abund.rel.viz[abund.rel.viz != 0])
beanplot( log( abund.rel.viz[ most_signif_50[1] ,  ] ),
          log( abund.rel.viz[ most_signif_50[2] ,  ] ),
          #log( abund.rel.viz[ most_signif_50[3] ,  ] ),
          log( abund.rel.viz[ most_signif_50[4] ,  ] ),
          log( abund.rel.viz[ most_signif_50[5] ,  ] ),
          log( abund.rel.viz[ most_signif_50[6] ,  ] ),
          log( abund.rel.viz[ most_signif_50[7] ,  ] ),
          log( abund.rel.viz[ most_signif_50[8] ,  ] ),
          log( abund.rel.viz[ most_signif_50[9] ,  ] ),
          log( abund.rel.viz[ most_signif_50[10] ,  ] ),
          what=c(TRUE,TRUE,TRUE,FALSE)
          
          )
 
beanplot( log( abund.rel.viz[ most_signif_50[41] ,  ] ),
          log( abund.rel.viz[ most_signif_50[42] ,  ] ),
          log( abund.rel.viz[ most_signif_50[43] ,  ] ),
          log( abund.rel.viz[ most_signif_50[44] ,  ] ),
          log( abund.rel.viz[ most_signif_50[45] ,  ] ),
          log( abund.rel.viz[ most_signif_50[46] ,  ] ),
          log( abund.rel.viz[ most_signif_50[47] ,  ] ),
          log( abund.rel.viz[ most_signif_50[48] ,  ] ),
          log( abund.rel.viz[ most_signif_50[49] ,  ] ),
          log( abund.rel.viz[ most_signif_50[50] ,  ] ),
          what=c(TRUE,TRUE,TRUE,FALSE)
          
)

diffexp_orig <- read.table("ALLSAMPLES_analysis/diffexp/graphclust/differential_expression.csv", sep = ",", stringsAsFactors = FALSE)
diffexp <- diffexp_orig[2:dim(diffexp_orig)[1],2:dim(diffexp_orig)[2]]
rownames(diffexp) <- diffexp_orig[2:dim(diffexp_orig)[1],1]
colnames(diffexp) <- as.character(unlist(diffexp_orig[1,2:dim(diffexp_orig)[2]])) 

most_signif <- c( names(sum_mat_multi_m.fdr[,1][ sum_mat_multi_m.fdr[,1] < 10e-100 ])#, 
                  #names(sum_mat_multi_m.fdr[,2][ sum_mat_multi_m.fdr[,2] < 10e-100 ]), 
                  #names(sum_mat_multi_m.fdr[,3][ sum_mat_multi_m.fdr[,3] < 10e-100 ]) 
                 )
most_signif <- unique(most_signif)


signifnames <- names(sum_mat_multi_m[,1][ !is.na(sum_mat_multi_m[,1]) ][ sum_mat_multi_m[,1][ !is.na(sum_mat_multi_m[,1]) ] < 10e-100 ])
counts <- sapply(diffexp[ rownames(diffexp) %in% most_signif , grepl("p value", colnames(diffexp)) ], as.numeric )
rownames(counts) <- rownames(diffexp[ rownames(diffexp) %in% most_signif , ])

png("wildType_e100.png",height=500,width=800)
heatmap.2( -log(counts),
           margins = c(15,15), trace='none', #col=c(rep("blue",50),rev(redblue(150)),rep("red",50)) )#, scale='column')
           col=c(rev(redblue(50)),rep("red",100)))
dev.off()

most_signif <- c( names(sum_mat_multi_m.fdr[,2][ sum_mat_multi_m.fdr[,2] < 10e-100 ]))
most_signif <- unique(most_signif)

signifnames <- names(sum_mat_multi_m[,1][ !is.na(sum_mat_multi_m[,1]) ][ sum_mat_multi_m[,1][ !is.na(sum_mat_multi_m[,1]) ] < 10e-100 ])
counts <- sapply(diffexp[ rownames(diffexp) %in% most_signif , grepl("p value", colnames(diffexp)) ], as.numeric )
rownames(counts) <- rownames(diffexp[ rownames(diffexp) %in% most_signif , ])
png("PS1_e100.png",height=500,width=800)
heatmap.2( -log(counts),
           margins = c(15,10), trace='none', col=c(rev(redblue(50)),rep("red",550))) #, scale='column')
dev.off()

most_signif <- c( names(sum_mat_multi_m.fdr[,3][ sum_mat_multi_m.fdr[,3] < 10e-100 ]))
most_signif <- unique(most_signif)

signifnames <- names(sum_mat_multi_m[,1][ !is.na(sum_mat_multi_m[,1]) ][ sum_mat_multi_m[,1][ !is.na(sum_mat_multi_m[,1]) ] < 10e-100 ])
counts <- sapply(diffexp[ rownames(diffexp) %in% most_signif , grepl("p value", colnames(diffexp)) ], as.numeric )
rownames(counts) <- rownames(diffexp[ rownames(diffexp) %in% most_signif , ])
png("APP_e100.png",height=500,width=800)
heatmap.2( -log(counts),
           margins = c(15,10), trace='none', col=c(rev(redblue(50)),rep("red",550))) #, scale='column')
dev.off()
 




most_signif <- c( names(sum_mat_multi_m.fdr[,1][ sum_mat_multi_m.fdr[,1] < 10e-100 ]), 
                  names(sum_mat_multi_m.fdr[,2][ sum_mat_multi_m.fdr[,2] < 10e-100 ]), 
                  names(sum_mat_multi_m.fdr[,3][ sum_mat_multi_m.fdr[,3] < 10e-100 ]) 
)
most_signif <- unique(most_signif)

dd <- diffexp[ most_signif , grepl("p value", colnames(diffexp))  ]
dd <- dd[ -1 , ] 
ddd <- sapply(dd, as.numeric)
rownames(ddd) <- most_signif[!is.na(most_signif)]

colSums(ddd < 0.05)

# save.image("20190912.am.Rdata")





############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################


wt_enriched  <- names(log(abund.rel[barcodes[1],][abund.rel[barcodes[1],] != 0]))[ log(abund.rel[barcodes[1],][abund.rel[barcodes[1],] != 0]) > -4 ]
ps1_enriched <- names(log(abund.rel[barcodes[2],][abund.rel[barcodes[2],] != 0]))[ log(abund.rel[barcodes[2],][abund.rel[barcodes[2],] != 0]) > -3.47 ]
app_enriched <- names(log(abund.rel[barcodes[3],][abund.rel[barcodes[3],] != 0]))[ log(abund.rel[barcodes[3],][abund.rel[barcodes[3],] != 0]) > -4 ]

enriched <- matrix(nrow=3,ncol=length(colnames(abund.rel)))
colnames(enriched) <- colnames(abund.rel)
rownames(enriched) <- c("Control","PS1","APP")
enriched[is.na(enriched)] <- 0
enriched[1,wt_enriched]   <- 1
enriched[2,ps1_enriched]  <- 1
enriched[3,app_enriched]  <- 1
 
 
sum_mat_1 <- sum_mat_2 <- sum_mat_3 <- sum_mat_4 <- fold_changes <- matrix(nrow=length(c(as.character(labels[labels[,3] == "Gene Expression",1]))),ncol=3)
QQQ <- length( c(as.character(labels[labels[,3] == "Gene Expression",1])) )

for(m in c(1:3)){
  qi <- 1
  for(q in c(as.character(labels[labels[,3] == "Gene Expression",1]))){
    
    
    if( qi %% 50 == 0){ print(paste(m,(qi/QQQ)*100,"%")) }
    
    reg_1    <- lm( abund.rel[q,] ~ abund.rel[barcodes[m],] + 1 )
    sum_mat_1[qi, m] <- as.matrix(anova(reg_1))[1,5]
    reg_4    <- lm( abund.rel[q,] ~ as.factor(enriched[m,]) + 1 )
    sum_mat_4[qi, m] <- as.matrix(anova(reg_4))[1,5]
    
    
    if(m == 1){
      
      reg_2    <- lm( abund.rel[q,] ~ abund.rel[barcodes[1],] + abund.rel[barcodes[2],] + abund.rel[barcodes[3],] + 1 )
      sum_mat_2[qi, ] <- unlist(as.matrix(anova(reg_2))[1:3,5])
      reg_3    <- lm( abund.rel[q,] ~ as.factor(enriched[1,]) + as.factor(enriched[2,]) + as.factor(enriched[3,]) + 1 )
      sum_mat_3[qi, ] <- unlist(as.matrix(anova(reg_3))[1:3,5])
      
      if(mean(abund.rel[q,!colnames(abund.rel) %in% wt_enriched] ) > 0){
        FC_wt  <- mean(abund.rel[q,wt_enriched])  / mean(abund.rel[q,!colnames(abund.rel) %in% wt_enriched] )
      }else{ FC_wt <- 0 }
      
      if(mean(abund.rel[q,!colnames(abund.rel) %in% ps1_enriched] ) > 0){
       FC_ps1 <- mean(abund.rel[q,ps1_enriched]) / mean(abund.rel[q,!colnames(abund.rel) %in% ps1_enriched])
      }else{ FC_ps1 <- 0 }
      
      if(mean(abund.rel[q,!colnames(abund.rel) %in% app_enriched] ) > 0){
        FC_app <- mean(abund.rel[q,app_enriched]) / mean(abund.rel[q,!colnames(abund.rel) %in% app_enriched])
      }else{ FC_app <- 0 }
      
      if(FC_wt != 0){  FC_wt  <- log(FC_wt) }
      if(FC_ps1 != 0){ FC_ps1 <- log(FC_ps1) }
      if(FC_app != 0){ FC_app <- log(FC_app) }
      
      fold_changes[qi,] <- c(FC_wt, FC_ps1, FC_app)
      
    }
    qi <- qi + 1
    
  }
}

# save.image("q.Rdata")

rownames(sum_mat_1) <- rownames(sum_mat_2) <- rownames(sum_mat_3) <- rownames(sum_mat_4) <- rownames(fold_changes) <- c(as.character(labels[labels[,3] == "Gene Expression",1]))
colnames(sum_mat_1) <- colnames(sum_mat_2) <- colnames(sum_mat_3) <- colnames(sum_mat_4) <- colnames(fold_changes) <- c("Control", "PS1+", "APP+")
write.table(sum_mat_1, file = "pvalues_univariate_barcodes.tsv", sep="\t")
write.table(sum_mat_4, file = "pvalues_univariate_barcodesEnnriched.tsv", sep="\t")
write.table(sum_mat_2, file = "pvalues_multivariate_barcodes.tsv", sep="\t")
write.table(sum_mat_3, file = "pvalues_multivariate_barcodesEnnriched.tsv", sep="\t")
write.table(fold_changes, file = "fold_changes.tsv", sep="\t")
 
sum_mat_1.fdr <- sum_mat_1
sum_mat_1.fdr[,1] <- p.adjust( sum_mat_1[,1] , method="fdr" )
sum_mat_1.fdr[,2] <- p.adjust( sum_mat_1[,2] , method="fdr" )
sum_mat_1.fdr[,3] <- p.adjust( sum_mat_1[,3] , method="fdr" )

png(filename="pvalues_univariate_barcodes.control.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )[ !is.infinite(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )[ !is.infinite(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )) ],
      pch=15,col="darkblue",cex=0.7
    )  
dev.off()
png(filename="pvalues_univariate_barcodes.ps1.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )[ !is.infinite(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )[ !is.infinite(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()
png(filename="pvalues_univariate_barcodes.app.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )[ !is.infinite(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )[ !is.infinite(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()


sum_mat_1.fdr <- sum_mat_4
sum_mat_1.fdr[,1] <- p.adjust( sum_mat_4[,1] , method="fdr" )
sum_mat_1.fdr[,2] <- p.adjust( sum_mat_4[,2] , method="fdr" )
sum_mat_1.fdr[,3] <- p.adjust( sum_mat_4[,3] , method="fdr" )

png(filename="pvalues_univariate_enriched.control.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )[ !is.infinite(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )[ !is.infinite(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()
png(filename="pvalues_univariate_enriched.ps1.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )[ !is.infinite(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )[ !is.infinite(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()
png(filename="pvalues_univariate_enriched.app.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )[ !is.infinite(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )[ !is.infinite(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()



sum_mat_1.fdr <- sum_mat_2
sum_mat_1.fdr[,1] <- p.adjust( sum_mat_2[,1] , method="fdr" )
sum_mat_1.fdr[,2] <- p.adjust( sum_mat_2[,2] , method="fdr" )
sum_mat_1.fdr[,3] <- p.adjust( sum_mat_2[,3] , method="fdr" )

png(filename="pvalues_multivariate_barcodes.control.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )[ !is.infinite(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )[ !is.infinite(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()
png(filename="pvalues_multivariate_barcodes.ps1.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )[ !is.infinite(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )[ !is.infinite(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()
png(filename="pvalues_multivariate_barcodes.app.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )[ !is.infinite(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )[ !is.infinite(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()


sum_mat_1.fdr <- sum_mat_3
sum_mat_1.fdr[,1] <- p.adjust( sum_mat_3[,1] , method="fdr" )
sum_mat_1.fdr[,2] <- p.adjust( sum_mat_3[,2] , method="fdr" )
sum_mat_1.fdr[,3] <- p.adjust( sum_mat_3[,3] , method="fdr" )

png(filename="pvalues_multivariate_enriched.control.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )[ !is.infinite(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )[ !is.infinite(-log( sum_mat_1.fdr[,1][!is.na(sum_mat_1.fdr[,1])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()
png(filename="pvalues_multivariate_enriched.ps1.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )[ !is.infinite(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )[ !is.infinite(-log( sum_mat_1.fdr[,2][!is.na(sum_mat_1.fdr[,2])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()
png(filename="pvalues_multivariate_enriched.app.png",height=800,width=800)
plot( fold_changes[names(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )[ !is.infinite(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )) ]), 1]  , 
      -log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )[ !is.infinite(-log( sum_mat_1.fdr[,3][!is.na(sum_mat_1.fdr[,3])] )) ],
      pch=15,col="darkblue",cex=0.7
)  
dev.off()




# ONCE PER MODEL
# set cmpr to the summary matrix, fdr-adjust it, and print the number over the cutoff

cmpr  <- sum_mat_3

cmpr.fdr <- cmpr
cmpr.fdr[,1] <- p.adjust( cmpr[,1] , method="fdr" )
cmpr.fdr[,2] <- p.adjust( cmpr[,2] , method="fdr" )
cmpr.fdr[,3] <- p.adjust( cmpr[,3] , method="fdr" )

q_wt_pos  <- rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-50 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-50 , ]) ,1] > 0  ]
q_wt_neg  <- rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-50 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-50 , ]) ,1] < 0  ]

q_ps1_pos  <- rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-50 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-50 , ]) ,1] > 0  ]
q_ps1_neg  <- rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-50 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-50 , ]) ,1] < 0  ]

q_app_pos  <- rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-50 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-50 , ]) ,1] > 0  ]
q_app_neg  <- rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-50 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-50 , ]) ,1] < 0  ]

paste( length(q_wt_pos) , length(q_wt_neg) , length(q_ps1_pos) , length(q_ps1_neg) , length(q_app_pos) , length(q_app_neg) )


q_wt_pos  <- rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-100 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-100 , ]) ,1] > 0  ]
q_wt_neg  <- rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-100 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-100 , ]) ,1] < 0  ]

q_ps1_pos  <- rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-100 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-100 , ]) ,1] > 0  ]
q_ps1_neg  <- rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-100 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-100 , ]) ,1] < 0  ]

q_app_pos  <- rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-100 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-100 , ]) ,1] > 0  ]
q_app_neg  <- rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-100 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-100 , ]) ,1] < 0  ]

paste( length(q_wt_pos) , length(q_wt_neg) , length(q_ps1_pos) , length(q_ps1_neg) , length(q_app_pos) , length(q_app_neg) )





dd <- diffexp[ most_signif , grepl("p value", colnames(diffexp))  ]
dd <- dd[ -1 , ] 
ddd <- sapply(dd, as.numeric)
rownames(ddd) <- most_signif[!is.na(most_signif)]

colSums(ddd < 0.05)

# save.image("20190912.am.Rdata")








nureport <- read.table("all_pvalues.20190913.txt", sep = "\t", stringsAsFactors = FALSE)
nureportq <- nureport[2:dim(nureport)[1],2:dim(nureport)[2]]
rownames(nureportq) <- nureport[2:dim(nureport)[1],1]
colnames(nureportq) <- as.character(unlist(nureport[1,2:dim(nureport)[2]])) 

nureportq$labels <- rep(NA,dim(nureportq)[1])
for(r in rownames(nureportq)){
  nureportq[r,"labels"] <- as.character(labels[labels[,1]==r,2])
}

nureportq$Cluster1_logFC <- nureportq$Cluster1_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster2_logFC <- nureportq$Cluster2_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster3_logFC <- nureportq$Cluster3_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster4_logFC <- nureportq$Cluster4_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster5_logFC <- nureportq$Cluster5_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster6_logFC <- nureportq$Cluster6_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster7_logFC <- nureportq$Cluster7_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster8_logFC <- nureportq$Cluster8_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster9_logFC <- nureportq$Cluster9_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster10_logFC <- nureportq$Cluster10_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster11_logFC <- nureportq$Cluster11_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster12_logFC <- nureportq$Cluster12_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster13_logFC <- nureportq$Cluster13_P <- rep(NA,dim(nureportq)[1])
nureportq$Cluster14_logFC <- nureportq$Cluster14_P <- rep(NA,dim(nureportq)[1])

#c( colnames( diffexp[  , grepl("p value", colnames(diffexp)) ,   ] ) , colnames( diffexp[  , grepl("fold change", colnames(diffexp)) ,   ] ) )

for(r in rownames(nureportq)){
  if(r %in% rownames(diffexp)){
    nureportq[r,colnames(diffexp[ r, c(3,4,6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28,30,31,33,34,36,37,39,40,42,43) ])] <- diffexp[ r, c(3,4,6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28,30,31,33,34,36,37,39,40,42,43) ]
  }
}

write.table(nureportq, "comparable_pvalue_matrix.20190913.tsv", sep = "\t")








# ONCE PER MODEL
# set cmpr to the summary matrix, fdr-adjust it, et the names over a looser cutoff

cmpr  <- sum_mat_4

cmpr.fdr <- cmpr
cmpr.fdr[,1] <- p.adjust( cmpr[,1] , method="fdr" )
cmpr.fdr[,2] <- p.adjust( cmpr[,2] , method="fdr" )
cmpr.fdr[,3] <- p.adjust( cmpr[,3] , method="fdr" )

q_wt_pos  <- rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-5 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-5 , ]) ,1] > 0  ]
q_wt_neg  <- rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-5 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,1][  !is.na(cmpr.fdr[,1])  ] < 10e-5 , ]) ,1] < 0  ]

q_ps1_pos  <- rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-5 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-5 , ]) ,1] > 0  ]
q_ps1_neg  <- rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-5 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,2][  !is.na(cmpr.fdr[,2])  ] < 10e-5 , ]) ,1] < 0  ]

q_app_pos  <- rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-5 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-5 , ]) ,1] > 0  ]
q_app_neg  <- rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-5 , ])[  fold_changes[ rownames(cmpr.fdr[ cmpr.fdr[,3][  !is.na(cmpr.fdr[,3])  ] < 10e-5 , ]) ,1] < 0  ]

paste( length(q_wt_pos) , length(q_wt_neg) , length(q_ps1_pos) , length(q_ps1_neg) , length(q_app_pos) , length(q_app_neg) )

library(gplots)

# ctrl # 

most_signif <- q_wt_pos
counts <- sapply(diffexp[ rownames(diffexp) %in% most_signif , grepl("p value", colnames(diffexp)) ], as.numeric )
rownames(counts) <- rownames(diffexp[ rownames(diffexp) %in% most_signif , ])
nams <- unique(rownames(counts)[ rowSums(counts < 1e-1) ])
nams_orig <- nams
ri <- 1
for(r in nams){
  if(r %in% as.character(labels[,1])){    nams[ri] <- paste(r, labels[labels[,1] == r,2])   }
  ri <- ri + 1
}
result <- counts[ rownames(counts) %in% nams_orig , ]
colSums(  result < 0.01  )
png("q_wt_pos.cluster_specific_pvalues.png",height=800,width=1300)
heatmap.2( -log(counts[ nams_orig , ]),
           margins = c(15,20), trace='none',
           col=c(rev(redblue(50)),rep("red",100)), 
           main = paste0("Control : Multivariate Continuous Model (Positive)\n", length(nams), " genes"),
           #cexRow = 0.7, 
           keysize = 0.7,
           labRow = nams
           )
dev.off()

most_signif <- q_wt_neg
counts <- sapply(diffexp[ rownames(diffexp) %in% most_signif , grepl("p value", colnames(diffexp)) ], as.numeric )
rownames(counts) <- rownames(diffexp[ rownames(diffexp) %in% most_signif , ])
nams <- unique(rownames(counts)[ rowSums(counts < 1e-1) ])
nams_orig <- nams
ri <- 1
for(r in nams){
  if(r %in% as.character(labels[,1])){    nams[ri] <- paste(r, labels[labels[,1] == r,2])   }
  ri <- ri + 1
}
result <- counts[ rownames(counts) %in% nams_orig , ]
colSums(  result < 0.01  )
png("q_wt_neg.cluster_specific_pvalues.png",height=800,width=1300)
heatmap.2( -log(counts[ nams_orig , ]),
           margins = c(15,20), trace='none',
           col=c(rev(redblue(50)),rep("red",100)), 
           main = paste0("Control : Multivariate Continuous Model (Negative)\n", length(nams), " genes"),
           #cexRow = 0.7, 
           keysize = 0.7,
           labRow = nams
)
dev.off()

# ps1 #

most_signif <- q_ps1_pos
counts <- sapply(diffexp[ rownames(diffexp) %in% most_signif , grepl("p value", colnames(diffexp)) ], as.numeric )
rownames(counts) <- rownames(diffexp[ rownames(diffexp) %in% most_signif , ])
nams <- unique(rownames(counts)[ rowSums(counts < 1e-1) ])
nams_orig <- nams
ri <- 1
for(r in nams){
  if(r %in% as.character(labels[,1])){    nams[ri] <- paste(r, labels[labels[,1] == r,2])   }
  ri <- ri + 1
}
result <- counts[ rownames(counts) %in% nams_orig , ]
colSums(  result < 0.01  )
png("q_ps1_pos.cluster_specific_pvalues.png",height=800,width=1300)
heatmap.2( -log(counts[ nams_orig , ]),
           margins = c(15,20), trace='none',
           col=c(rev(redblue(50)),rep("red",100)), 
           main = paste0("PS1+ : Multivariate Continuous Model (Positive)\n", length(nams), " genes"),
           #cexRow = 0.7, 
           keysize = 0.7,
           labRow = nams
)
dev.off()

most_signif <- q_ps1_neg
counts <- sapply(diffexp[ rownames(diffexp) %in% most_signif , grepl("p value", colnames(diffexp)) ], as.numeric )
rownames(counts) <- rownames(diffexp[ rownames(diffexp) %in% most_signif , ])
nams <- unique(rownames(counts)[ rowSums(counts < 1e-1) ])
nams_orig <- nams
ri <- 1
for(r in nams){
  if(r %in% as.character(labels[,1])){    nams[ri] <- paste(r, labels[labels[,1] == r,2])   }
  ri <- ri + 1
}
result <- counts[ rownames(counts) %in% nams_orig , ]
colSums(  result < 0.01  )
png("q_ps1_neg.cluster_specific_pvalues.png",height=800,width=1300)
heatmap.2( -log(counts[ nams_orig , ]),
           margins = c(15,20), trace='none',
           col=c(rev(redblue(50)),rep("red",100)), 
           main = paste0("PS1+ : Multivariate Continuous Model (Negative)\n", length(nams), " genes"),
           #cexRow = 0.7, 
           keysize = 0.7,
           labRow = nams
)
dev.off()

# app #

most_signif <- q_app_pos
counts <- sapply(diffexp[ rownames(diffexp) %in% most_signif , grepl("p value", colnames(diffexp)) ], as.numeric )
rownames(counts) <- rownames(diffexp[ rownames(diffexp) %in% most_signif , ])
nams <- unique(rownames(counts)[ rowSums(counts < 1e-1) ])
nams_orig <- nams
ri <- 1
for(r in nams){
  if(r %in% as.character(labels[,1])){    nams[ri] <- paste(r, labels[labels[,1] == r,2])   }
  ri <- ri + 1
}
result <- counts[ rownames(counts) %in% nams_orig , ]
colSums(  result < 0.01  )
png("q_app_pos.cluster_specific_pvalues.png",height=800,width=1300)
heatmap.2( -log(counts[ nams_orig , ]),
           margins = c(15,20), trace='none',
           col=c(rev(redblue(50)),rep("red",100)), 
           main = paste0("APP+ : Multivariate Continuous Model (Positive)\n", length(nams), " genes"),
           #cexRow = 0.7, 
           keysize = 0.7,
           labRow = nams
)
dev.off()

most_signif <- q_app_neg
counts <- sapply(diffexp[ rownames(diffexp) %in% most_signif , grepl("p value", colnames(diffexp)) ], as.numeric )
rownames(counts) <- rownames(diffexp[ rownames(diffexp) %in% most_signif , ])
nams <- unique(rownames(counts)[ rowSums(counts < 1e-1) ])
nams_orig <- nams
ri <- 1
for(r in nams){
  if(r %in% as.character(labels[,1])){    nams[ri] <- paste(r, labels[labels[,1] == r,2])   }
  ri <- ri + 1
}
result <- counts[ rownames(counts) %in% nams_orig , ]
colSums(  result < 0.01  )
png("q_app_neg.cluster_specific_pvalues.png",height=800,width=1300)
heatmap.2( -log(counts[ nams_orig , ]),
           margins = c(15,20), trace='none',
           col=c(rev(redblue(50)),rep("red",100)), 
           main = paste0("APP+ : Multivariate Continuous Model (Negative)\n", length(nams), " genes"),
           #cexRow = 0.7, 
           keysize = 0.7,
           labRow = nams
)
dev.off()
 









# 

install.packages("ccmap")



















# compare "case" and "control" # 

# save.image("www.Rdata")

over_undefined   <- colnames(abund.rel)[ !colnames(abund.rel) %in% ps1_enriched ] 
over_undefined   <- over_undefined[ !over_undefined %in% app_enriched ]

fold_changes_only <- matrix(nrow=length(c(as.character(labels[labels[,3] == "Gene Expression",1]))),ncol=15)
qi <- 1
for(q in c(as.character(labels[labels[,3] == "Gene Expression",1]))){

      if(mean(abund.rel[q,colnames(abund.rel) %in% wt_enriched] ) > 0){
        FC_ps1_only          <- mean(abund.rel[q,ps1_enriched])  / mean(abund.rel[q,colnames(abund.rel) %in% wt_enriched] )
      }else{ 
        if(mean(abund.rel[q,ps1_enriched]) > 0){
          FC_ps1_only <- Inf
        }else{
          FC_ps1_only <- 0
        }
      }
  
      if(mean(abund.rel[q,colnames(abund.rel) %in% wt_enriched] ) > 0){
        FC_app_only  <- mean(abund.rel[q,app_enriched])  / mean(abund.rel[q,colnames(abund.rel) %in% wt_enriched] )
      }else{ 
        if(mean(abund.rel[q,app_enriched]) > 0){
          FC_app_only <- Inf
        }else{
          FC_app_only <- 0
        }
      }
      
      if(mean(abund.rel[q,over_undefined]) > 0){
        FC_ps1_and_undefined <- mean(abund.rel[q,ps1_enriched]) / mean(abund.rel[q,over_undefined])
        FC_app_and_undefined <- mean(abund.rel[q,app_enriched]) / mean(abund.rel[q,over_undefined])
      }else{
        if(mean(abund.rel[q,ps1_enriched]) > 0){
          FC_ps1_and_undefined <- Inf
        }else{
          FC_ps1_and_undefined <- 0
        }
        if(mean(abund.rel[q,app_enriched]) > 0){
          FC_app_and_undefined <- Inf
        }else{
          FC_app_and_undefined <- 0
        }
      }
      
      # replace previous #
      # # # # # # # # # # 
      # # # # # # # # # # 
        
        if(mean(abund.rel[q,!colnames(abund.rel) %in% wt_enriched] ) > 0){
          FC_wt  <- mean(abund.rel[q,wt_enriched])  / mean(abund.rel[q,!colnames(abund.rel) %in% wt_enriched] )
        }else{ 
          if( mean(abund.rel[q,wt_enriched]) > 0 ){
            FC_wt <- 0 
          }else{
            FC_wt <- Inf
          }
        }
        
        if(mean(abund.rel[q,!colnames(abund.rel) %in% ps1_enriched] ) > 0){
          FC_ps1 <- mean(abund.rel[q,ps1_enriched]) / mean(abund.rel[q,!colnames(abund.rel) %in% ps1_enriched])
        }else{ 
          if( mean(abund.rel[q,ps1_enriched]) > 0 ){
            FC_ps1 <- 0 
          }else{
            FC_ps1 <- Inf
          }
        }
        
        if(mean(abund.rel[q,!colnames(abund.rel) %in% app_enriched] ) > 0){
          FC_app <- mean(abund.rel[q,app_enriched]) / mean(abund.rel[q,!colnames(abund.rel) %in% app_enriched])
        }else{ 
          if( mean(abund.rel[q,app_enriched]) > 0 ){
            FC_app <- 0
          }else{
            FC_app <- Inf
          }
        }
  
      # # # # # # # # # # 
      # # # # # # # # # # 
      # # # # # # # # # # 
  
      if(FC_ps1_only != 0 && !is.infinite(FC_ps1_only)){ FC_ps1_only <- log(FC_ps1_only) }
      if(FC_app_only != 0 && !is.infinite(FC_app_only)){ FC_app_only <- log(FC_app_only) }
      if(FC_ps1_and_undefined != 0 && !is.infinite(FC_ps1_and_undefined)){ FC_ps1_and_undefined <- log(FC_ps1_and_undefined) }
      if(FC_app_and_undefined != 0 && !is.infinite(FC_app_and_undefined)){ FC_app_and_undefined <- log(FC_app_and_undefined) }
      
      #fold_changes_only[qi,] <- c(FC_ps1_only, FC_ps1_and_undefined, FC_app_only, FC_app_and_undefined, FC_wt, FC_ps1, FC_app, sum( abund.rel[q,] > 0 ))
    
      aa <- colnames(abund.rel) %in% ps1_enriched
      bb <- colnames(abund.rel) %in% app_enriched
      cc <- colnames(abund.rel) %in% over_undefined
      reg_q    <- lm( abund.rel[q,] ~ aa + bb + cc + 1 )
      reg_r    <- lm( abund.rel[q,] ~ aa + cc + 1 )
      reg_s    <- lm( abund.rel[q,] ~ bb + cc + 1 )
      p_q <- unlist(as.matrix(anova(reg_q))[1:3,5])
      p_r <- unlist(as.matrix(anova(reg_r))[1:2,5])
      p_s <- unlist(as.matrix(anova(reg_s))[1:2,5])
      
      fold_changes_only[qi,] <- unlist(c(FC_ps1_only, FC_ps1_and_undefined, FC_app_only, FC_app_and_undefined, FC_wt, FC_ps1, FC_app, 
                                  sum( abund.rel[q,] > 0 ),
                                  unname(p_q),
                                  unname(p_r),
                                  unname(p_s)))
      
      
      qi <- qi + 1
    
}
 
# save.image("zzz.Rdata")

load("zzz.Rdata")

#sum_mat.bonf <- sum_mat
#sum_mat.bonf[,1] <- p.adjust( sum_mat[,1] , method="bonferroni" )
#sum_mat.bonf[,2] <- p.adjust( sum_mat[,2] , method="bonferroni" )
#sum_mat.bonf[,3] <- p.adjust( sum_mat[,3] , method="bonferroni" )
#sum_mat.bonf.pop <- sum_mat[rowSums(is.na(sum_mat.bonf)) < 3,]
#sum_mat.bonf.viz <- sum_mat.bonf.pop
#sum_mat.bonf.viz[sum_mat.bonf.pop == 0] <- min(sum_mat.bonf.pop[sum_mat.bonf.pop != 0])

# Yes, you can do the overall model, but the important follow-up comparisons should be CONTROL vs. APP and CONTROL vs. PS1.

# You could do a separate model for each, or you could do an ANOVA with post-hoc comparisons (using, e.g., Scheffe, Tukey, or Bonferroni corrections):
  
fold_changes_only.bonf <- fold_changes_only
fold_changes_only.bonf[,9] <- p.adjust( fold_changes_only[,9] , method="bonferroni" )
fold_changes_only.bonf[,10] <- p.adjust( fold_changes_only[,9] , method="bonferroni" )
fold_changes_only.bonf[,11] <- p.adjust( fold_changes_only[,9] , method="bonferroni" )
fold_changes_only.bonf[,12] <- p.adjust( fold_changes_only[,9] , method="bonferroni" )
fold_changes_only.bonf[,13] <- p.adjust( fold_changes_only[,9] , method="bonferroni" )
fold_changes_only.bonf[,14] <- p.adjust( fold_changes_only[,9] , method="bonferroni" )
fold_changes_only.bonf[,15] <- p.adjust( fold_changes_only[,9] , method="bonferroni" )

rownames(fold_changes_only.bonf) <- as.character(labels[labels[,3] == "Gene Expression",1])

par(mfrow=c(2,1))
#png(filename="20190915.FC_ps1_only.png",height=800,width=800)
plot( fold_changes_only.bonf[names(-log( fold_changes_only.bonf[,9][!is.na(fold_changes_only.bonf[,9])] )[ !is.infinite(-log( fold_changes_only.bonf[,9][!is.na(fold_changes_only.bonf[,9])] )) ]), 1]  , 
      -log( fold_changes_only.bonf[,9][!is.na(fold_changes_only.bonf[,9])] )[ !is.infinite(-log( fold_changes_only.bonf[,9][!is.na(fold_changes_only.bonf[,9])] )) ],
      pch=15,col="darkblue",cex=0.7
)  

#dev.off()



# ALL OUT # 

remove(nu_mat)

nu_mat <- matrix(nrow = 33538, ncol = 53)

#nu_mat[,1:6]   <- lapply( nureportq[,1:6] ) # fold changes : test 1 (in/out)
nu_mat[,1]     <- nureportq[,1]
nu_mat[,2]     <- nureportq[,2]
nu_mat[,3]     <- nureportq[,3]
nu_mat[,4]     <- nureportq[,4]
nu_mat[,5]     <- nureportq[,5]
nu_mat[,6]     <- nureportq[,6]
#nu_mat[,7:9]   <- sum_mat_1[,1:3] # test 1 : univariate   : numeric
nu_mat[,7]   <- sum_mat_1[,1] # test 1 : univariate   : numeric
nu_mat[,8]   <- sum_mat_1[,2] # test 1 : univariate   : numeric
nu_mat[,9]   <- sum_mat_1[,3] # test 1 : univariate   : numeric
#nu_mat[,10:12] <- sum_mat_4 # test 1 : univariate   : enriched
nu_mat[,10] <- sum_mat_4[,1] # test 1 : univariate   : enriched
nu_mat[,11] <- sum_mat_4[,2] # test 1 : univariate   : enriched
nu_mat[,12] <- sum_mat_4[,3] # test 1 : univariate   : enriched
#nu_mat[,13:15] <- sum_mat_2 # test 1 : multivariate : numeric
nu_mat[,13] <- sum_mat_2[,1] # test 1 : multivariate : numeric
nu_mat[,14] <- sum_mat_2[,2] # test 1 : multivariate : numeric
nu_mat[,15] <- sum_mat_2[,3] # test 1 : multivariate : numeric
#nu_mat[,16:18] <- sum_mat_3 # test 1 : multivariate : enriched
nu_mat[,16] <- sum_mat_3[,1] # test 1 : multivariate : enriched
nu_mat[,17] <- sum_mat_3[,2] # test 1 : multivariate : enriched
nu_mat[,18] <- sum_mat_3[,3] # test 1 : multivariate : enriched

nu_mat[,19]   <- p.adjust( sum_mat_1[,1] , method="bonferroni" ) # test 1 : univariate   : numeric  : bonferroni
nu_mat[,20]   <- p.adjust( sum_mat_1[,2] , method="bonferroni" ) # test 1 : univariate   : numeric  : bonferroni
nu_mat[,21]   <- p.adjust( sum_mat_1[,3] , method="bonferroni" ) # test 1 : univariate   : numeric  : bonferroni

nu_mat[,22]   <- p.adjust( sum_mat_4[,1] , method="bonferroni" ) # test 1 : univariate   : enriched  : bonferroni
nu_mat[,23]   <- p.adjust( sum_mat_4[,2] , method="bonferroni" ) # test 1 : univariate   : enriched  : bonferroni
nu_mat[,24]   <- p.adjust( sum_mat_4[,3] , method="bonferroni" ) # test 1 : univariate   : enriched  : bonferroni

nu_mat[,25]   <- p.adjust( sum_mat_2[,1] , method="bonferroni" ) # test 1 : multivariate   : numeric  : bonferroni
nu_mat[,26]   <- p.adjust( sum_mat_2[,2] , method="bonferroni" ) # test 1 : multivariate   : numeric  : bonferroni
nu_mat[,27]   <- p.adjust( sum_mat_2[,3] , method="bonferroni" ) # test 1 : multivariate   : numeric  : bonferroni

nu_mat[,28]   <- p.adjust( sum_mat_3[,1] , method="bonferroni" ) # test 1 : multivariate   : enriched  : bonferroni
nu_mat[,29]   <- p.adjust( sum_mat_3[,2] , method="bonferroni" ) # test 1 : multivariate   : enriched  : bonferroni
nu_mat[,30]   <- p.adjust( sum_mat_3[,3] , method="bonferroni" ) # test 1 : multivariate   : enriched  : bonferroni

#nu_mat[,31:45] <- fold_changes_only
w <- 1
for(z in c(31:45)){
  if(z %in% c(35:37)){
    tmp <- log(fold_changes_only[,w])
    tmp[ tmp == -Inf ] <- 0
    nu_mat[,z] <- tmp
  }else{
    nu_mat[,z] <- fold_changes_only[,w]
  }
  w <- w + 1
}
#nu_mat[,46:52] <- fold_changes_only.bonf[,9:15]
w <- 9
for(z in c(46:52)){
  nu_mat[,z] <- fold_changes_only.bonf[,w]
  w <- w + 1
}

nu_mat[,53] <- nureportq[,"labels"]

rownames(nu_mat) <- c(as.character(labels[labels[,3] == "Gene Expression",1]))
write.table(nu_mat, file="summary.20190916.tsv", sep="\t")
 









zreport <- read.table("all_pvalues.20190916_reformatted_pm.txt", sep = "\t", stringsAsFactors = FALSE)
zreportz <- zreport[2:dim(zreport)[1],2:dim(zreport)[2]]
rownames(zreportz) <- zreport[2:dim(zreport)[1],1]
colnames(zreportz) <- as.character(unlist(zreport[1,2:dim(zreport)[2]])) 


plot.new()
par(mfrow=c(3,2))

fcs <- as.numeric(zreportz$T3.FC_ps1_only)
nodepch <- rep(15, length(fcs))
nodepch[is.infinite(fcs)] <- 3
fcs[is.infinite(fcs)] <- max(fcs[!is.infinite(fcs)])+1
#logfcs <- log(fcs)
plot(fcs  ,
     -log(as.numeric(zreportz$T3.multi.Enriched.PS1.Bonferroni)) ,
     pch=nodepch,col="darkblue",cex=0.7,
     main="T3.multi.Enriched.PS1.Bonferroni",
     xlab="log2(FoldChange)", ylab="-log(Bonferroni-adjusted P-value)"
)

fcs <- as.numeric(zreportz$T3.FC_app_only)
nodepch <- rep(15, length(fcs))
nodepch[is.infinite(fcs)] <- 3
#fcs[is.infinite(fcs)] <- max(fcs[!is.infinite(fcs)])+1
#logfcs <- log(fcs)
plot(fcs  ,
     -log(as.numeric(zreportz$T3.multi.Enriched.APP.Bonferroni)) ,
     pch=nodepch,col="darkblue",cex=0.7,
     main="T3.multi.Enriched.APP.Bonferroni",
     xlab="log2(FoldChange)", ylab="-log(Bonferroni-adjusted P-value)"
     )


fcs <- as.numeric(zreportz$T3.FC_ps1_only)
nodepch <- rep(15, length(fcs))
nodepch[is.infinite(fcs)] <- 3
fcs[is.infinite(fcs)] <- max(fcs[!is.infinite(fcs)])+1
#logfcs <- log(fcs)
plot(fcs  ,
     -log(as.numeric(zreportz$T3.multi.Enriched.PS1andUNDEF.Bonferroni)) ,
     pch=nodepch,col="darkblue",cex=0.7,
     main="T3.multi.Enriched.PS1andUNDEF.Bonferroni",
     xlab="log2(FoldChange)", ylab="-log(Bonferroni-adjusted P-value)"
)

fcs <- as.numeric(zreportz$T3.FC_app_only)
nodepch <- rep(15, length(fcs))
nodepch[is.infinite(fcs)] <- 3
#fcs[is.infinite(fcs)] <- max(fcs[!is.infinite(fcs)])+1
#logfcs <- log(fcs)
plot(fcs  ,
     -log(as.numeric(zreportz$T3.multi.Enriched.APPandUNDEF.Bonferroni)) ,
     pch=nodepch,col="darkblue",cex=0.7,
     main="T3.multi.Enriched.APPandUNDEF.Bonferroni",
     xlab="log2(FoldChange)", ylab="-log(Bonferroni-adjusted P-value)"
)

fcs <- as.numeric(zreportz$T3.FC_ps1_andUndefined)
nodepch <- rep(15, length(fcs))
nodepch[is.infinite(fcs)] <- 3
fcs[is.infinite(fcs)] <- max(fcs[!is.infinite(fcs)])+1
#logfcs <- log(fcs)
plot(fcs,
     -log(as.numeric(zreportz$T4.multi.Enriched.PS1.Bonferroni)) ,
     pch=nodepch,col="darkblue",cex=0.7,
     main="T4.multi.Enriched.PS1.Bonferroni",
     xlab="log2(FoldChange)", ylab="-log(Bonferroni-adjusted P-value)"
)

fcs <- as.numeric(zreportz$T3.FC_app_andUndefined)
nodepch <- rep(15, length(fcs))
nodepch[is.infinite(fcs)] <- 3
fcs[is.infinite(fcs)] <- max(fcs[!is.infinite(fcs)])+1
#logfcs <- log(fcs)
plot(fcs,
     -log(as.numeric(zreportz$T4.multi.Enriched.APP.Bonferroni)),
     pch=nodepch,col="darkblue",cex=0.7,
     main="T4.multi.Enriched.APP.Bonferroni",
     xlab="log2(FoldChange)", ylab="-log(Bonferroni-adjusted P-value)"
)










#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ccmap")

library(ccdata)
library(ccmap)

# The command for the query itself is 
# 'l1000_res<-query_drugs(query_sig,l1000_es)', 
# where query sig is a named vector containing the effect sizes for each gene below the cutoff for significance. 
# Note that you will need to convert to HGNC symbols.


zreport <- read.table("all_pvalues.20190916_reformatted_pm.txt", sep = "\t", stringsAsFactors = FALSE)
zreportz <- zreport[2:dim(zreport)[1],2:dim(zreport)[2]]
rownames(zreportz) <- zreport[2:dim(zreport)[1],1]
colnames(zreportz) <- as.character(unlist(zreport[1,2:dim(zreport)[2]])) 

#  -1 (complete opposite, reverses gene signature) and 1 (identical, mimics signature).

cutoff <- 1e-10

query            <- zreportz[!is.na(as.numeric(zreportz$T3.multi.Enriched.PS1andUNDEF.Bonferroni)),][ as.numeric(zreportz$T3.multi.Enriched.PS1andUNDEF.Bonferroni)[ !is.na(as.numeric(zreportz$T3.multi.Enriched.PS1andUNDEF.Bonferroni)) ] < cutoff , ]
query            <- query[ order(query$T3.multi.Enriched.PS1andUNDEF.Bonferroni) , ]
queryCmap        <- as.numeric(zreportz[rownames(query),"T3.FC_ps1_andUndefined"])
names(queryCmap) <- query$Symbols #rownames(query)

l1000_res_ps1 <-query_drugs(queryCmap,drug_info = c("cmap", "l1000"), sorted = TRUE)
head( l1000_res_ps1 )
tail( l1000_res_ps1 )

query            <- zreportz[!is.na(as.numeric(zreportz$T3.multi.Enriched.APPandUNDEF.Bonferroni)),][ as.numeric(zreportz$T3.multi.Enriched.APPandUNDEF.Bonferroni)[ !is.na(as.numeric(zreportz$T3.multi.Enriched.PS1andUNDEF.Bonferroni)) ] < cutoff , ]
query            <- query[ order(query$T3.multi.Enriched.APPandUNDEF.Bonferroni) , ]
queryCmap        <- as.numeric(zreportz[rownames(query),"T3.FC_app_andUndefined"])
names(queryCmap) <- query$Symbols #rownames(query)

l1000_res_app <-query_drugs(queryCmap,drug_info = c("cmap", "l1000"), sorted = TRUE)
head( l1000_res_app )
tail( l1000_res_app )

write.table(l1000_res_ps1, "l1000_res_ps1_1e-10.tsv", sep="\t")
write.table(l1000_res_app, "l1000_res_app_1e-10.tsv", sep="\t")

cutoff <- 1e-50

query            <- zreportz[!is.na(as.numeric(zreportz$T3.multi.Enriched.PS1andUNDEF.Bonferroni)),][ as.numeric(zreportz$T3.multi.Enriched.PS1andUNDEF.Bonferroni)[ !is.na(as.numeric(zreportz$T3.multi.Enriched.PS1andUNDEF.Bonferroni)) ] < cutoff , ]
query            <- query[ order(query$T3.multi.Enriched.PS1andUNDEF.Bonferroni) , ]
queryCmap        <- as.numeric(zreportz[rownames(query),"T3.FC_ps1_andUndefined"])
names(queryCmap) <- query$Symbols #rownames(query)

l1000_res_ps1 <-query_drugs(queryCmap,drug_info = c("cmap", "l1000"), sorted = TRUE)
head( l1000_res_ps1 )
tail( l1000_res_ps1 )

query            <- zreportz[!is.na(as.numeric(zreportz$T3.multi.Enriched.APPandUNDEF.Bonferroni)),][ as.numeric(zreportz$T3.multi.Enriched.APPandUNDEF.Bonferroni)[ !is.na(as.numeric(zreportz$T3.multi.Enriched.PS1andUNDEF.Bonferroni)) ] < cutoff , ]
query            <- query[ order(query$T3.multi.Enriched.APPandUNDEF.Bonferroni) , ]
queryCmap        <- as.numeric(zreportz[rownames(query),"T3.FC_app_andUndefined"])
names(queryCmap) <- query$Symbols #rownames(query)

l1000_res_app <-query_drugs(queryCmap,drug_info = c("cmap", "l1000"), sorted = TRUE)
head( l1000_res_app )
tail( l1000_res_app )

write.table(l1000_res_ps1, "l1000_res_ps1_1e-50.tsv", sep="\t")
write.table(l1000_res_app, "l1000_res_app_1e-50.tsv", sep="\t")
 

