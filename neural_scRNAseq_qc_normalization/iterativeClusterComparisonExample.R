

setwd("C:\\Users\\jamis\\Desktop")
library(vegan)
library(gplots)
 
##########
# COUNTS # 
##########
intron.abund <- read.table("BRAIN_in_ex/intron_intergenic.Gx_25998.Cx_50/INTRON_INTERGENIC.Gx_25998.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.txt", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
#intron.abund <- read.table("BRAIN_in_ex/intron_intergenic.Gx_1000.Cx_50/INTRON_INTERGENIC.Gx_1000.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.txt", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
colnames(intron.abund) <- intron.abund[1, ]
intron.abund           <- intron.abund[-1,]
rownames(intron.abund) <- intron.abund[, 1]
intron.abund           <- intron.abund[,-1]

##########
# COUNTS # 
##########
covariate.matrix <- read.table("BRAIN_in_ex/20161026_covariate_table.format.txt", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
colnames(covariate.matrix) <- covariate.matrix[1, ]
covariate.matrix           <- covariate.matrix[-1,]
rownames(covariate.matrix) <- covariate.matrix[, 1]
covariate.matrix           <- covariate.matrix[,-1]
# glia  inh  exc 
# 177   646  970


####################################################
# SIGNIF. OF SEPARATION BETWEEN 2 EXAMPLE CLUSTERS # 
####################################################
# pull all ids associated with cluster id 1
rn.a <- rownames(covariate.matrix)[ covariate.matrix$Cluster_ID_20161007 == 28 ]
rn.a <- rn.a[!is.na(rn.a)]
# and cluster 2
rn.b <- rownames(covariate.matrix)[ covariate.matrix$Cluster_ID_20161007 == 33 ]
rn.b <- rn.b[!is.na(rn.b)]

dim(intron.abund[rn.a,])
dim(intron.abund[rn.b,])
 
# metadata  and counts for cluster id 1
md.a <- covariate.matrix[rn.a,1:134] 
        #intron.abund[rn.a,1:49]
ct.a <- intron.abund[rn.a,50:dim(intron.abund)[1]]
ct.a <- sapply(ct.a, as.numeric)
rownames(ct.a) <- rn.a
ct.rel.a <- decostand(ct.a, method="total")
# and 2
md.b <- covariate.matrix[rn.b,1:134] 
ct.b <- intron.abund[rn.b,50:dim(intron.abund)[1]]
ct.b <- sapply(ct.b, as.numeric)
rownames(ct.b) <- rn.b
ct.rel.b <- decostand(ct.b, method="total")

compare_matrix <- as.data.frame(rbind(ct.rel.a, ct.rel.b))
compare_matrix[,"compare"] <- rep(NA, dim(compare_matrix)[1]) 
compare_matrix[rn.a,]$compare <- "Cluster28"
compare_matrix[rn.b,]$compare <- "Cluster33"

# significance of GENE association (lm)
L <- length(colnames(compare_matrix)[ colnames(compare_matrix) != "compare" ])
gene_pvalues <- matrix(nrow=L,ncol=6)
colnames(gene_pvalues) <- c("pvalue.Raw", "pvalue.FDR", "isEnrichedinA", "slope", "intercept", "fStat")
rownames(gene_pvalues) <- colnames(compare_matrix)[ colnames(compare_matrix) != "compare" ]
gi <- 1
for( g in colnames(compare_matrix)[ colnames(compare_matrix) != "compare" ] ){
  if(gi %% 20 == 0){ print(paste(g,(gi/L)*100,"%")) }
  reg <- lm( as.numeric( compare_matrix[,g] ) ~ as.factor(compare_matrix$compare) + 1 )
         #lm( as.factor(compare_matrix$compare) ~ as.numeric( compare_matrix[,g] ) )
  residual <- resid(reg)
  if (t(coef(summary(reg)))[2,1] > 0 & is.nan(t(coef(summary(reg)))[2,1]) == 0){
    cpus.lm2 <- as.matrix( anova(reg) )
    fstat    <- summary(reg)$fstatistic
    slope    <- unname(coef(reg)[2])
    int      <- unname(coef(reg)[1])
    posneg   <- as.numeric(slope > 0)
    raw_pval <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
  }else{
    raw_pval <- slope <- posneg <- int <- NA
  }
  gene_pvalues[gi,] <- c(raw_pval, NA, posneg, slope, int, unname(fstat[1]))
  gi <- gi +1
}
gene_pvalues[,"pvalue.FDR"] <- p.adjust(gene_pvalues[,"pvalue.Raw"], method="fdr")

# significance of METADATA association (lm)
metadata_pvalues <- matrix(nrow=134,ncol=6)
colnames(metadata_pvalues) <- c("pvalue.Raw", "pvalue.FDR", "isEnrichedinA", "slope", "intercept", "fStat")
rownames(metadata_pvalues) <- colnames(covariate.matrix)[1:134]
for(m in colnames(covariate.matrix)[1:22]){
  print(m)
  # factored metadata: varied enough to compare?
  cmpr.L <- length(unique(covariate.matrix[,m]))
  if(cmpr.L < 4 | cmpr.L > dim(covariate.matrix)[1]-3){
    pop.m <- covariate.matrix[ !is.na( covariate.matrix[,m] ) , ]
    pop.m.wAbund <- rownames(pop.m)[ rownames(pop.m) %in% rownames(compare_matrix) ]
    hm <- sapply(compare_matrix[pop.m.wAbund,],as.numeric)
    rownames(hm) <- pop.m.wAbund
    hm[is.na(hm)] <- 0
    for(g in colnames(hm)){
      reg <- lm( as.numeric( hm[,g] ) ~ as.factor(compare_matrix$compare) + 1 )
    }
    #cr <- rownames(hm) %in% rn.a
    #cr[cr == TRUE] <- "darkviolet"
    #cr[cr == FALSE] <- "darkorange"
    #rc <- as.numeric(as.factor(covariate.matrix[ pop.m.wAbund , m ]))
    #rc <- c("maroon","darkgreen")[rc]
    #heatmap.2( log( hm[,names(sort(gene_pvalues[,"pvalue.FDR"])[1:50])] +0.000000001), 
    #           trace = "none", col=redblue(20), mar=c(5,15), keysize = 0.75, 
    #           colRow = cr, cexRow = 0.5, cexCol=0.3, 
    #           RowSideColors = rc)
  }else{
    print("Skipped")
  }
}
for(m in colnames(covariate.matrix)[23:134]){
  # numerical metadata
  pop.m <- covariate.matrix[ !is.na( covariate.matrix[,m] ) , ]
  pop.m.wAbund <- rownames(pop.m)[ rownames(pop.m) %in% rownames(compare_matrix) ]
  hm <- sapply(compare_matrix[pop.m.wAbund,],as.numeric)
  for(g in colnames(hm)){
    reg <- lm( hm[,g] ~ as.numeric(covariate.matrix[pop.m.wAbund,m]) + 1 )
  }
}

sum_mat <- matrix(nrow=8,ncol=length(colnames(hm))-1)
# test case
x <- 1
for(m in c("ACTB_Ct",
           "% Exact Duplicates(Raw Seq.)",
           "cDNA_Pico_AdjConc",
           "% GC Content(Pretrimmed)",
           "percentTrimmedOverRawReads",
           "BP Count(Raw Seq.)",
           "[MitoCore] Average Total Count of (Non-Zero) Core Genes(Core Genes)",
           "Std Dev. Phred Score(Raw Seq.)"
)){
  pop.m <- covariate.matrix[ !is.na( covariate.matrix[,m] ) , ]
  pop.m.wAbund <- rownames(pop.m)[ rownames(pop.m) %in% rownames(compare_matrix) ]
  hm <- sapply(compare_matrix[pop.m.wAbund,],as.numeric)
  y <- 1
  for(g in colnames(hm)[colnames(hm) != "compare"]){
    reg <- lm( hm[,g] ~ as.factor(covariate.matrix[pop.m.wAbund,"Cluster_ID_20161007"]) + as.numeric(covariate.matrix[pop.m.wAbund,m]) + 1 )
    if (t(coef(summary(reg)))[2,1] > 0 & is.nan(t(coef(summary(reg)))[2,1]) == 0){
      fstat    <- summary(reg)$fstatistic
      raw_pval <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
    }else{
      raw_pval <- NA
    }
    sum_mat[x,y] <- raw_pval
    y <- y + 1
  }
  x <- x + 1
}

rownames(sum_mat) <- c("ACTB_Ct",
                      "% Exact Duplicates(Raw Seq.)",
                      "cDNA_Pico_AdjConc",
                      "% GC Content(Pretrimmed)",
                      "percentTrimmedOverRawReads",
                      "BP Count(Raw Seq.)",
                      "[MitoCore] Average Total Count of (Non-Zero) Core Genes(Core Genes)",
                      "Std Dev. Phred Score(Raw Seq.)")
colnames(sum_mat) <- colnames(hm)[colnames(hm) != "compare"]

mins <- list()
mi <- 1
for(g in 1:dim(sum_mat)[2]){
  mins[mi] <- min(sum_mat[,g][!is.na(sum_mat[,g])])
  mi <- mi + 1
}
sum_mat[ , order(unlist(mins)) ][,1:100]
 
heatmap.2( -log( sum_mat[ , order(unlist(mins)) ][,1:100] ), 
                      trace = "none", col=redblue(20), mar=c(5,15), keysize = 0.75, 
                      #colRow = cr, 
                      cexRow = 0.5, cexCol=0.3, 
                      #RowSideColors = rc
           )     
