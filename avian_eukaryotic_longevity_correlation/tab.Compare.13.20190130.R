
library(colorRamps)
library(gplots)
setwd ("~/Desktop/develop.Avian")
 
list1 <- c(paste0("RERUN.","alignment.","CAPER.","logwMLSresid_v_relAbs.report.tsv"),              #   "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv", 
           paste0("RERUN.","alignment.","LM.","logwMLSresid_v_relAbs.report.tsv"),                 #   "alignmentPlotsAndTables/alignment.10.LM.logwMLSresid_v_relAbsreport.tsv",
           paste0("RERUN.","alignment.","LM.","logwMLSresid_v_absAbs.report.tsv"),                 #   "alignmentPlotsAndTables/alignment.11.LM.logwMLSresid_v_absAbsreport.tsv",
           paste0("RERUN.","alignment.","LM.","logwMLSresid_v_binaryAbs.report.tsv"),              #   "alignmentPlotsAndTables/alignment.12.LM.logwMLSresid_v_binaryAbsreport.tsv",
           paste0("RERUN.","alignment.","LM.","wMLSresid_v_relAbs.report.tsv"),                    #   "alignmentPlotsAndTables/alignment.13.LM.wMLSresid_v_relAbsreport.tsv",
           paste0("RERUN.","alignment.","LM.","wMLSresid_v_absAbs.report.tsv"),                    #   "alignmentPlotsAndTables/alignment.14.LM.wMLSresid_v_absAbsreport.tsv",
           paste0("RERUN.","alignment.","LM.","wMLSresid_v_binaryAbs.report.tsv"),                 #   "alignmentPlotsAndTables/alignment.15.LM.wMLSresid_v_binaryAbsreport.tsv",
           paste0("RERUN.","alignment.","LM.","MLS_v_relAbs.report.tsv"),                          #   "alignmentPlotsAndTables/alignment.16.LM.MLS_v_relAbsreport.tsv",
           paste0("RERUN.","alignment.","LM.","MLS_v_absAbs.report.tsv"),                          #   "alignmentPlotsAndTables/alignment.17.LM.MLS_v_absAbsreport.tsv",
           paste0("RERUN.","alignment.","LM.","MLS_v_binaryAbs.report.tsv"),                       #   "alignmentPlotsAndTables/alignment.18.LM.MLS_v_binaryAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport.tsv",
           paste0("RERUN.","alignment.","CAPER.","logwMLSresid_v_absAbs.report.tsv"),              #   "alignmentPlotsAndTables/alignment.2.CAPER.logwMLSresid_v_absAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.20.extremeLongevity_sd2_hi.logwMLSresid_v_relAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.21.extremeLongevity_sd2_lo.logwMLSresid_v_relAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.22.extremeLongevity_sd1_all.logwMLSresid_v_relAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.23.extremeLongevity_sd1_hi.logwMLSresid_v_relAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.24.extremeLongevity_sd1_lo.logwMLSresid_v_relAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.26.extremeLongevity_sd1andsd2_hi.logwMLSresid_v_relAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.27.extremeLongevity_sd1andsd2_lo.logwMLSresid_v_relAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.28.extremeLongevity_sd2_all.logwMLSresid_v_binaryAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.29.extremeLongevity_sd2_hi.logwMLSresid_v_binaryAbsreport.tsv",
           paste0("RERUN.","alignment.","CAPER.","logwMLSresid_v_binaryAbs.report.tsv"),           #   "alignmentPlotsAndTables/alignment.3.CAPER.logwMLSresid_v_binaryAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.30.extremeLongevity_sd2_lo.logwMLSresid_v_binaryAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.31.extremeLongevity_sd1_all.logwMLSresid_v_binaryAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.32.extremeLongevity_sd1_hi.logwMLSresid_v_binaryAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.33.extremeLongevity_sd1_lo.logwMLSresid_v_binaryAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.34.extremeLongevity_sd1andsd2_all.logwMLSresid_v_binaryAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.35.extremeLongevity_sd1andsd2_hi.logwMLSresid_v_binaryAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.36.extremeLongevity_sd1andsd2_lo.logwMLSresid_v_binaryAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.37.Family_v_relAbsreport.tsv","alignmentPlotsAndTables/alignment.38.Family_v_absAbsreport.tsv",
           "alignmentPlotsAndTables/alignment.39.Family_v_binaryAbsreport.tsv",
           paste0("RERUN.","alignment.","CAPER.","wMLSresid_v_relAbs.report.tsv"),                #   "alignmentPlotsAndTables/alignment.4.CAPER.wMLSresid_v_relAbsreport.tsv",
           paste0("RERUN.","alignment.","CAPER.","wMLSresid_v_absAbs.report.tsv"),                #   "alignmentPlotsAndTables/alignment.5.CAPER.wMLSresid_v_absAbsreport.tsv",
           paste0("RERUN.","alignment.","CAPER.","wMLSresid_v_binaryAbs.report.tsv"),             #   "alignmentPlotsAndTables/alignment.6.CAPER.wMLSresid_v_binaryAbsreport.tsv",
           paste0("RERUN.","alignment.","CAPER.","MLS_v_relAbs.report.tsv"),                      #   "alignmentPlotsAndTables/alignment.7.CAPER.MLS_v_relAbsreport.tsv",
           paste0("RERUN.","alignment.","CAPER.","MLS_v_absAbs.report.tsv"),                      #   "alignmentPlotsAndTables/alignment.8.CAPER.MLS_v_absAbsreport.tsv",
           paste0("RERUN.","alignment.","CAPER.","MLS_v_binaryAbs.report.tsv"),                   #   "alignmentPlotsAndTables/alignment.9.CAPER.MLS_v_binaryAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","CAPER.","logwMLSresid_v_relAbs.report.tsv"),                #   "deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","LM.","logwMLSresid_v_relAbs.report.tsv"),                   #   "deNovoPlotsAndTables/denovo_highConfidence.10.LM.logwMLSresid_v_relAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","LM.","logwMLSresid_v_absAbs.report.tsv"),                   #   "deNovoPlotsAndTables/denovo_highConfidence.11.LM.logwMLSresid_v_absAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","LM.","logwMLSresid_v_binaryAbs.report.tsv"),                 #   "deNovoPlotsAndTables/denovo_highConfidence.12.LM.logwMLSresid_v_binaryAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","LM.","wMLSresid_v_relAbs.report.tsv"),                      #   "deNovoPlotsAndTables/denovo_highConfidence.13.LM.wMLSresid_v_relAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","LM.","wMLSresid_v_absAbs.report.tsv")                      #   "deNovoPlotsAndTables/denovo_highConfidence.14.LM.wMLSresid_v_absAbsreport.tsv")
          )
list2 <- c(paste0("RERUN.","denovo_highConfidence.","LM.","wMLSresid_v_binaryAbs.report.tsv"),                   #   "deNovoPlotsAndTables/denovo_highConfidence.15.LM.wMLSresid_v_binaryAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","LM.","MLS_v_relAbs.report.tsv"),                            #   "deNovoPlotsAndTables/denovo_highConfidence.16.LM.MLS_v_relAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","LM.","MLS_v_absAbs.report.tsv"),                            #   "deNovoPlotsAndTables/denovo_highConfidence.17.LM.MLS_v_absAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","LM.","MLS_v_binaryAbs.report.tsv"),                         #   "deNovoPlotsAndTables/denovo_highConfidence.18.LM.MLS_v_binaryAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","CAPER.","wMLSresid_v_absAbs.report.tsv")                    # "deNovoPlotsAndTables/denovo_highConfidence.2.CAPER.logwMLSresid_v_absAbsreport.tsv")
          )
list3 <- c("deNovoPlotsAndTables/denovo_highConfidence.20.extremeLongevity_sd2_hi.logwMLSresid_v_relAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.21.extremeLongevity_sd2_lo.logwMLSresid_v_relAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.22.extremeLongevity_sd1_all.logwMLSresid_v_relAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.23.extremeLongevity_sd1_hi.logwMLSresid_v_relAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.24.extremeLongevity_sd1_lo.logwMLSresid_v_relAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport.tsv")
list4 <- c("deNovoPlotsAndTables/denovo_highConfidence.26.extremeLongevity_sd1andsd2_hi.logwMLSresid_v_relAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.27.extremeLongevity_sd1andsd2_lo.logwMLSresid_v_relAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.28.extremeLongevity_sd2_all.logwMLSresid_v_binaryAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.29.extremeLongevity_sd2_hi.logwMLSresid_v_binaryAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","CAPER.","logwMLSresid_v_binaryAbs.report.tsv"),            #   "deNovoPlotsAndTables/denovo_highConfidence.3.CAPER.logwMLSresid_v_binaryAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.30.extremeLongevity_sd2_lo.logwMLSresid_v_binaryAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.31.extremeLongevity_sd1_all.logwMLSresid_v_binaryAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.32.extremeLongevity_sd1_hi.logwMLSresid_v_binaryAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.33.extremeLongevity_sd1_lo.logwMLSresid_v_binaryAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.34.extremeLongevity_sd1andsd2_all.logwMLSresid_v_binaryAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.35.extremeLongevity_sd1andsd2_hi.logwMLSresid_v_binaryAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.36.extremeLongevity_sd1andsd2_lo.logwMLSresid_v_binaryAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.37.Family_v_relAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.38.Family_v_absAbsreport.tsv",
           "deNovoPlotsAndTables/denovo_highConfidence.39.Family_v_binaryAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","CAPER.","wMLSresid_v_relAbs.report.tsv"),                 #   "deNovoPlotsAndTables/denovo_highConfidence.4.CAPER.wMLSresid_v_relAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","CAPER.","wMLSresid_v_absAbs.report.tsv"),                 #   "deNovoPlotsAndTables/denovo_highConfidence.5.CAPER.wMLSresid_v_absAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","CAPER.","wMLSresid_v_binaryAbs.report.tsv"),              #   "deNovoPlotsAndTables/denovo_highConfidence.6.CAPER.wMLSresid_v_binaryAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","CAPER.","MLS_v_relAbs.report.tsv"),                       #   "deNovoPlotsAndTables/denovo_highConfidence.7.CAPER.MLS_v_relAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","CAPER.","MLS_v_absAbs.report.tsv"),                       #   "deNovoPlotsAndTables/denovo_highConfidence.8.CAPER.MLS_v_absAbsreport.tsv",
           paste0("RERUN.","denovo_highConfidence.","CAPER.","MLS_v_binaryAbs.report.tsv")                     #   "deNovoPlotsAndTables/denovo_highConfidence.9.CAPER.MLS_v_binaryAbsreport.tsv") 
         )
           

# initial pass
model_initial_stats <- matrix(ncol=22,nrow=length(c(list1,list2,list3,list4)))
rownames(model_initial_stats) <- c(list1,list2,list3,list4)
colnames(model_initial_stats) <- c( 
                                    "count_nonNA",
                                    "count_nonNA_andZeroCutoff",
                                    "rawPvalue_signif_0.005__POS",
                                    "rawPvalue_signif_0.001__POS",
                                    "rawPvalue_signif_0.0005_POS",
                                    "fdrPvalue_signif_0.05___POS",
                                    "fdrPvalue_signif_0.01___POS",
                                    "rawPvalue_signif_0.005__NEG",
                                    "rawPvalue_signif_0.001__NEG",
                                    "rawPvalue_signif_0.0005_NEG",
                                    "fdrPvalue_signif_0.05___NEG",
                                    "fdrPvalue_signif_0.01___NEG",
                                    "rawPvalue_signif_0.005__ALL",
                                    "rawPvalue_signif_0.001__ALL",
                                    "rawPvalue_signif_0.0005_ALL",
                                    "fdrPvalue_signif_0.05___ALL",
                                    "fdrPvalue_signif_0.01___ALL",
                                    "rawPvalue_signif_0.005__ZeroCutoff",
                                    "rawPvalue_signif_0.001__ZeroCutoff",
                                    "rawPvalue_signif_0.0005_ZeroCutoff",
                                    "fdrPvalue_signif_0.05___ZeroCutoff",
                                    "fdrPvalue_signif_0.01___ZeroCutoff"
                                    )
significant_raw <- list()
all_groups_of_interest <- list()
for (i in c(list1,list2,list3,list4)) {
  
  print(i)
  #devel
  #i <- "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv"
  tab <- read.table(i, sep="\t", header=TRUE, fill=TRUE)
    
  model_initial_stats[i,"count_nonNA"]                  <- dim(tab)[1]
  model_initial_stats[i,"count_nonNA_andZeroCutoff"]    <- sum(tab$cutoff == 0)
  model_initial_stats[i,"rawPvalue_signif_0.005__POS"]  <- sum(tab$raw_pvalues <= 0.005 & tab$slopes > 0)
  model_initial_stats[i,"rawPvalue_signif_0.001__POS"]  <- sum(tab$raw_pvalues <= 0.001 & tab$slopes > 0)
  model_initial_stats[i,"rawPvalue_signif_0.0005_POS"]  <- sum(tab$raw_pvalues <= 0.0005 & tab$slopes > 0)
  model_initial_stats[i,"fdrPvalue_signif_0.05___POS"] <- sum(tab$fdr_pvalues <= 0.05 & tab$slopes > 0)
  model_initial_stats[i,"fdrPvalue_signif_0.01___POS"] <- sum(tab$fdr_pvalues <= 0.01 & tab$slopes > 0)
  model_initial_stats[i,"rawPvalue_signif_0.005__NEG"] <- sum(tab$raw_pvalues <= 0.005 & tab$slopes < 0)
  model_initial_stats[i,"rawPvalue_signif_0.001__NEG"] <- sum(tab$raw_pvalues <= 0.001 & tab$slopes < 0)
  model_initial_stats[i,"rawPvalue_signif_0.0005_NEG"] <- sum(tab$raw_pvalues <= 0.0005 & tab$slopes < 0)
  model_initial_stats[i,"fdrPvalue_signif_0.05___NEG"] <- sum(tab$fdr_pvalues <= 0.05 & tab$slopes < 0)
  model_initial_stats[i,"fdrPvalue_signif_0.01___NEG"] <- sum(tab$fdr_pvalues <= 0.01 & tab$slopes < 0)
  model_initial_stats[i,"rawPvalue_signif_0.005__ALL"]  <- sum(tab$raw_pvalues <= 0.005)
  model_initial_stats[i,"rawPvalue_signif_0.001__ALL"]  <- sum(tab$raw_pvalues <= 0.001)
  model_initial_stats[i,"rawPvalue_signif_0.0005_ALL"]  <- sum(tab$raw_pvalues <= 0.0005)
  model_initial_stats[i,"fdrPvalue_signif_0.05___ALL"] <- sum(tab$fdr_pvalues <= 0.05)
  model_initial_stats[i,"fdrPvalue_signif_0.01___ALL"] <- sum(tab$fdr_pvalues <= 0.01)
  model_initial_stats[i,"rawPvalue_signif_0.005__ZeroCutoff"] <- sum(tab$raw_pvalues <= 0.005 & tab$cutoff == 0)
  model_initial_stats[i,"rawPvalue_signif_0.001__ZeroCutoff"] <- sum(tab$raw_pvalues <= 0.001 & tab$cutoff == 0)
  model_initial_stats[i,"rawPvalue_signif_0.0005_ZeroCutoff"] <- sum(tab$raw_pvalues <= 0.0005 & tab$cutoff == 0)
  model_initial_stats[i,"fdrPvalue_signif_0.05___ZeroCutoff"] <- sum(tab$fdr_pvalues <= 0.05 & tab$cutoff == 0)
  model_initial_stats[i,"fdrPvalue_signif_0.01___ZeroCutoff"] <- sum(tab$fdr_pvalues <= 0.01 & tab$cutoff == 0)
  
  significant_raw <- append(significant_raw,        rownames(tab)[tab$raw_pvalues <= 0.001])
  all_groups_of_interest <- unlist(append(all_groups_of_interest, unlist(rownames(tab))))
  
}

write.table(model_initial_stats, file="model_initial_stats.tsv", sep = "\t")

all_convert <- list()
ni <- 1
for(n in unlist(all_groups_of_interest)){
  #print(n)
  #n <- 
  all_convert[ni] <- paste0("Cluster",as.numeric(strsplit( n , '_' )[[1]][2]))
  ni <- ni + 1
}
names(all_convert) <- all_groups_of_interest
  
# just curious - how many of these are interesting...  (50-120 per test)
png(filename = "model_initial_stats.png", height = 1000,width = 1000)
m <- model_initial_stats
#m <- as.data.frame(model_initial_stats)
#m <- m[order(rownames(m)),][1:78,]
colz <- rep("gray",length(rownames(m)))
colz[grepl("alignment", rownames(m))] <- "maroon"
colz[grepl("denovo", rownames(m))]    <- "darkblue"
colz[grepl("binary", rownames(m))]    <- "green"
colz[grepl("sd", rownames(m))]    <- "gray"
colz[grepl("LM", rownames(m))]    <- "white"
par(mar=c(4,42,4,4))
#barplot(m$count_nonNA, horiz = TRUE, names.arg = rownames(m), las=1, cex.names = 0.8, col = colz)
barplot(m[,"rawPvalue_signif_0.001__ALL"], horiz = TRUE, names.arg = rownames(m), las=1, cex.names = 0.8, col = colz)
dev.off()
  
summary_matrix <- matrix( nrow=length(unique(unlist(unname(all_convert)))) , ncol=(length(rownames(m))*6)+2 )
colnames(summary_matrix) <- c( "cutoff.denovo","cutoff.align",
                               paste0(rownames(m),".rawPval") , paste0(rownames(m),".rawPvalSortID"),
                               paste0(rownames(m),".fdrPval") , paste0(rownames(m),".fdrPvalSortID"),
                               paste0(rownames(m),".slope") ,   paste0(rownames(m),".intercept") 
                             )
rownames(summary_matrix) <- unique(unlist(unname(all_convert)))

for (i in rownames(m)) {
  
  print(i)
  #devel
  #i <- "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv"
  tab <- read.table(i, sep="\t", header=TRUE, fill=TRUE)
  rownames(tab) <- unlist(unname(all_convert[ rownames(tab) ]))
  
  if(grepl("alignment", i)==TRUE){
    summary_matrix[rownames(tab),"cutoff.align"] <- tab$cutoff
  }else{
    summary_matrix[rownames(tab),"cutoff.denovo"] <- tab$cutoff
  }
  summary_matrix[rownames(tab),paste0(i,".rawPval")]       <- tab$raw_pvalues
  #summary_matrix[rownames(tab),paste0(i,".rawPvalSortID")] <- order(tab$raw_pvalues)
  summary_matrix[rownames(tab),paste0(i,".fdrPval")]       <- tab$fdr_pvalues
  #summary_matrix[rownames(tab),paste0(i,".fdrPvalSortID")] <- order(tab$fdr_pvalues)
  summary_matrix[rownames(tab),paste0(i,".slope")]         <- tab$slopes
  summary_matrix[rownames(tab),paste0(i,".intercept")]     <- tab$intercepts
  
}
           
# how many ortholog groups got p-values in how many tests
rowSums( !is.na(summary_matrix[,paste0(rownames(m),".rawPval")]) )
 
summary_matrix.12.20190130 <- summary_matrix

##########

#load("20190118.Rdata")

#summary_matrix.12.20190130


# s <- summary_matrix[order(summary_matrix[,"RERUN.alignment.CAPER.logwMLSresid_v_relAbs.report.tsv.rawPvalSortID"]),
#                     grepl("rawPvalSortID",colnames(summary_matrix))][1:1000,]
# #s <- summary_matrix[order(summary_matrix[,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPvalSortID"]),
# #                    grepl("rawPvalSortID",colnames(summary_matrix))][1:1000,]
# #s <- summary_matrix[order(summary_matrix[,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPvalSortID"]),
# #                    grepl("rawPvalSortID",colnames(summary_matrix))][1:1000,]
# s <- s[,colSums(s, na.rm=TRUE) > 0]
# s[is.na(s)] <- 2*max(s, na.rm = TRUE)
# 
# colz <- rep("gray",length(colnames(s)))
# colz[grepl("alignment", colnames(s))] <- "maroon"
# colz[grepl("denovo", colnames(s))]    <- "darkblue"
# colz[grepl("binary", colnames(s))]    <- "green"
# colz[grepl("sd", colnames(s))]    <- "gray"
# colz[grepl("LM", colnames(s))]    <- "white"
# 
# png(filename = "dev3.png", height = 1000,width = 1000)
# heatmap.2(t(s), trace='none', col=c(rev(gray.colors(150)),matlab.like2(150)),
#           margins = c(10,40),
#           keysize=0.7, key.par = list(cex=0.5), RowSideColors = colz)
# dev.off()
 
#d <- summary_matrix[grepl('MCL_0', rownames(summary_matrix)), grepl('rawPvalSortID',colnames(summary_matrix))]
#d <- d[,colSums(d, na.rm=TRUE) > 0][1:100,]
#d[is.na(d)] <- 2*max(d, na.rm = TRUE)
#png(filename = "dev2.png", height = 1000,width = 1000)
#heatmap.2(t(d), trace='none', col=c(rev(gray.colors(150)),matlab.like2(150)),
#          margins = c(10,40),
#          keysize=0.7, key.par = list(cex=0.5))
#dev.off()

A <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Anas_platyrhynchos.id.txt", sep=" ", header=FALSE, fill=TRUE)
B <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Calypte_anna.id.txt", sep=" ", header=FALSE, fill=TRUE)
C <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Charadrius_vociferus.id.txt", sep=" ", header=FALSE, fill=TRUE)
D <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Columba_livia.id.txt", sep=" ", header=FALSE, fill=TRUE)
E <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Corvus_brachyrhynchos.id.txt", sep=" ", header=FALSE, fill=TRUE)
F <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Nipponia_nippon.id.txt", sep=" ", header=FALSE, fill=TRUE)
G <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Peregrine.id.txt", sep=" ", header=FALSE, fill=TRUE)
H <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Picoides_pubescens.id.txt", sep=" ", header=FALSE, fill=TRUE)
I <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Struthio_camelus.id.txt", sep=" ", header=FALSE, fill=TRUE)
X <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Taeniopygia_guttata.rough.id.txt", sep="\t", header=FALSE, fill=TRUE)
Y <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Gallus_gallus.rough.id.txt", sep="\t", header=FALSE, fill=TRUE)
Z <- read.table("/Users/apple/Desktop/DESKTOP_CLEANUP_2/Meleagris_gallopavo.rough.id.txt", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
gff_all <- rbind(A, B, C, D, E, F, G, H, I)
all_ids_rough <- c(as.character(X$V2), as.character(Y$V2), as.character(Z$V2))
all_ids_gff   <- c(as.character(A$V2), as.character(B$V2), as.character(C$V2), as.character(D$V2), as.character(E$V2), as.character(F$V2), as.character(G$V2), as.character(H$V2))

load(file="20190101_3pm.Rdata")
rpkm.nu2cmpr <- rpkm[ rpkm$V1 != "Gene" , ]

si <- 1
s <- summary_matrix[order(summary_matrix[,"RERUN.alignment.LM.logwMLSresid_v_relAbs.report.tsv.rawPval"]),]
s_withLabels <- as.data.frame(s)
s_withLabels$proteinIDs <- rep(NA,dim(s_withLabels)[1])
s_withLabels$transcriptIDs <- rep(NA,dim(s_withLabels)[1])
s_withLabels[,c(as.character(mls$AlignLabel))] <- rep(0,dim(s_withLabels)[1])
# for(ss in rownames(s)){
#   if(si %% 100 == 0)( print(si/dim(s_withLabels)[1]) )
#   clustID             <- strsplit(ss,'r')[[1]][2]
#   orthologClusterIDs  <- ref11_mcl[clustID,][ ref11_mcl[clustID,] != "" ]
#   transcriptIDs <- proteinIDs <- queryIDs <- symbols <- list()
#   mi <- 1
#   for(o in orthologClusterIDs){
#     mat <- matrix()
#     if(o %in% rpkm.nu2cmpr$V1){
#       mat <- rpkm.nu2[ rpkm.nu2cmpr$V1 == o , ]
#       mat[,5] <- rep(o,dim(mat)[1])
#       mat[,6] <- rep(NA,dim(mat)[1])
#       mat[,7] <- rep(clustID,dim(mat)[1])
#     }else if(o %in% transcript2prot$V2){
#       o2 <- transcript2prot[ transcript2prot$V2 == o , 1 ]
#       if( o2 %in% rpkm.nu2cmpr$V1 ){
#         mat <- rpkm.nu2[ rpkm.nu2cmpr$V1 == o2 , ]
#         mat[,5] <- rep(o,dim(mat)[1])
#         mat[,6] <- rep(o2,dim(mat)[1])
#         mat[,7] <- rep(clustID,dim(mat)[1])
#       }else{
#         mat <- matrix(nrow=1,ncol=7)
#         mat[,5] <- rep(o,dim(mat)[1])
#         mat[,6] <- rep(o2,dim(mat)[1])
#         mat[,7] <- rep(clustID,dim(mat)[1])
#       }
#     }
#     proteinIDs[mi]    <- paste0(unique(mat[,5]),collapse=",")
#     transcriptIDs[mi] <- paste0(unique(mat[,1]),collapse=",")
#     queryIDs[mi]  <- paste0(unique(mat[,3]),collapse=",")
#     mi <- mi + 1
#   }
#   keepers <- list()
#   for(p in queryIDs){ keepers <- append(keepers, unique(strsplit(p,","))[[1]]) }
#   keepers <- unlist(keepers)[unlist(keepers) != "NA"]
#   #as.numeric(c(as.character(mls$AlignLabel)) %in% keepers)
#   s_withLabels[ss,c(as.character(mls$AlignLabel))] <- as.numeric(c(as.character(mls$AlignLabel)) %in% keepers)
# 
#   keepers <- list()
#   for(p in proteinIDs){ keepers <- append(keepers, unique(strsplit(p,","))[[1]]) }
#   keepers <- unlist(keepers)[unlist(keepers) != "NA"]
#   s_withLabels[ss,"proteinIDs"]    <-  paste0(unique(keepers),collapse=",")
# 
#   #print(s_withLabels[ss,"proteinIDs"])
# 
#   keepers <- list()
#   for(p in transcriptIDs){ keepers <- append(keepers, unique(strsplit(p,","))[[1]]) }
#   keepers <- unlist(keepers)[unlist(keepers) != "NA"]
#   s_withLabels[ss,"transcriptIDs"] <- paste0(unique(keepers),collapse=",")
# 
#   #print(s_withLabels[ss,"transcriptIDs"])
# 
#   si <- si + 1
# }


# save.image(file="20190113.5.Rdata")
save.image(file="20190130.1.Rdata")

# load("20190113.5.Rdata")

s_withLabels$geneSymbols <- rep(NA,dim(s_withLabels)[1])

#symbols[mi]       <-
#for(n in transcripts_in_cluster){
# geneSymbols <- list()
# gi <- 1
# for(s in s_withLabels[,"transcriptIDs"]){
#   if(gi %% 250 == 0)( print(gi/dim(s_withLabels)[1]) )
#   geneSymbols[gi] <- ""
#   for(n in strsplit(s,",")[[1]]){
#     #gene symbol
#     if( n %in% all_ids_gff ){
#           geneSymbols[gi] <- paste( geneSymbols[gi] , as.character(unique(gff_all[gff_all$V2 %in% n,]$V3)), collapse="," )
#           #print(paste(gi,n,as.character(unique(gff_all[gff_all$V2 %in% n,]$V3))))
#     }else if( n %in% all_ids_rough){
#       if(n %in% as.character(X$V2)){
#           geneSymbols[gi] <- paste( geneSymbols[gi] , as.character( X[as.character(X$V2)==n,19] ), collapse="," )
#          # print(paste(gi,n,as.character( X[as.character(X$V2)==n,19] )))
#       }else if(n %in% as.character(Y$V2)){
#           geneSymbols[gi] <- paste( geneSymbols[gi] , as.character( X[as.character(Y$V2)==n,17] ), collapse="," )
#           #print(paste(gi,n,as.character( X[as.character(Y$V2)==n,17] )))
#       }else if (n %in% as.character(Z$V2)){
#           geneSymbols[gi] <- paste( geneSymbols[gi] , as.character( Z[as.character(Z$V2)==n,19] ), collapse="," )
#           #print(paste(gi,n, as.character(Z[as.character(Z$V2)==n,19] )))
#       }
#     }else{
#     }
#     #print(paste (gi,n,geneSymbols[gi]))
#   }
#   gi <- gi + 1
# }

s_withLabels$geneSymbols <- unlist(geneSymbols)

#write.table(s_withLabels, file="summary_matrix.20190121.tsv", sep = '\t', quote = FALSE)
write.table(s_withLabels, file="summary_matrix.20190130.tsv", sep = '\t', quote = FALSE)

#cmpr <- s_withLabels[ order(s_withLabels$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval`) ,  ]
#cmpr_signif_ClustIDs_0.01 <- cmpr[ , "deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.001
#order_cmpr1 <- order(as.numeric(tmpc$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval`))


png(filename = "denovo.fdrSortedPval100_compareAcrossModels.AlignCutoffLabel.20190130.png", height = 1000,width = 1400)
relevant_comparison_columns <- colnames(s_withLabels)[ grepl("fdrPval", colnames(s_withLabels)) ][ colSums(!is.na(s_withLabels[ , colnames(s_withLabels)[ grepl("fdrPval", colnames(s_withLabels)) ] ])) > 0 ][1:57]
raw_pvalues       <- s_withLabels[ order(s_withLabels$`RERUN.alignment.CAPER.logwMLSresid_v_relAbs.report.tsv.fdrPval`) ,
                                    relevant_comparison_columns]
raw_pvalue_orders <- raw_pvalues
for(n in colnames( raw_pvalue_orders )){
  raw_pvalue_orders[ , n ] <- order(raw_pvalue_orders[,n] , na.last = TRUE)
}
heatmap.2( t(as.matrix(raw_pvalue_orders[ 1:100 , ])) ,
           trace="none",
           col=c(matlab.like2(150)),
           margins = c(10,60),
           keysize=0.7, key.par = list(cex=0.5),
           Colv = FALSE,
           ColSideColors = greenred(50)[as.numeric( s_withLabels[ rownames(raw_pvalue_orders)[1:100] , "cutoff.align" ])+1]
           )
dev.off()

png(filename = "denovo.rawSortedPval100_compareAcrossModels.DeNovoCutoffLabel.20190130.png", height = 1000,width = 1400)
heatmap.2( t(as.matrix(raw_pvalue_orders[ 1:100 , ])) ,trace="none",col=c(matlab.like2(150)), margins = c(10,60),keysize=0.7, key.par = list(cex=0.5), Colv = FALSE,
           ColSideColors = greenred(50)[as.numeric( s_withLabels[ rownames(raw_pvalue_orders)[1:100] , "cutoff.denovo" ])+1]
)
dev.off()

relevant_comparison_columns <-  colnames(s_withLabels)[ grepl("fdrPval", colnames(s_withLabels)) ][ colSums(!is.na(s_withLabels[ , colnames(s_withLabels)[ grepl("fdrPval", colnames(s_withLabels)) ] ])) > 0 ][1:57]
raw_pvalues       <- s_withLabels[  order(s_withLabels$`alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval`) , relevant_comparison_columns]
raw_pvalue_orders <- raw_pvalues
for(n in colnames( raw_pvalue_orders )){
  raw_pvalue_orders[ , n ] <- order(raw_pvalue_orders[,n] , na.last = TRUE)
}
#png(filename = "align.rawSortedPval100_compareAcrossModels.AlignCutoffLabel.20190130.png", height = 1000,width = 1400)
#heatmap.2( t(as.matrix(raw_pvalue_orders[ 1:100 , ])) ,trace="none",col=c(matlab.like2(150)), margins = c(10,60),keysize=0.7, key.par = list(cex=0.5), Colv = FALSE,
#           ColSideColors = greenred(50)[as.numeric( s_withLabels[ rownames(raw_pvalue_orders)[1:100] , "cutoff.align" ])+1]
#)
#dev.off()
#png(filename = "align.rawSortedPval100_compareAcrossModels.DeNovoCutoffLabel.20190130.png", height = 1000,width = 1400)
#heatmap.2( t(as.matrix(raw_pvalue_orders[ 1:100 , ])) ,trace="none",col=c(matlab.like2(150)), margins = c(10,60),keysize=0.7, key.par = list(cex=0.5), Colv = FALSE,
#           ColSideColors = greenred(50)[as.numeric( s_withLabels[ rownames(raw_pvalue_orders)[1:100] , "cutoff.denovo" ])+1]
#)
#dev.off()

s.clean <- s_withLabels[ , colSums( !is.na(s_withLabels)) > 0 ]
s.clean <- s.clean[ , !grepl("SortID", colnames(s.clean)) ]
write.table( s.clean , file = "summary_matrix.clean.20190130.tsv", sep = "\t" )

# s_withLabels[ , grepl(".1.CAPER.", colnames(s_withLabels)) ]




#################################################################################################################################





#s.clean[ order(s.clean$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval`) , ]$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval`

# "alignmentPlotsAndTables/alignment.4.CAPER.wMLSresid_v_relAbsreport.tsv.rawPval"




# 
# 
# load("20190104_2pm.Rdata")
# # ^ reload the new environment
# 
# alignment.abundance_matrix.abs <- ag_count.hum#[ shared_abundances <= 46, ]
# alignment.abundance_matrix.rel <- ag_count.hum.rel#[ shared_abundances <= 46, ]
# alignment.shared_abundances    <- shared_abundances#[ shared_abundances <= 46 ]
# 
# denovo_highConfidence.abundance_matrix.abs <- og_counts.high_conf.abs
# denovo_highConfidence.abundance_matrix.rel <- og_counts.high_conf.rel
# denovo_highConfidence.shared_abundances    <- shared_abundances.high_conf
# 
# denovo_lowConfidence.abundance_matrix.abs <- og_counts.low_conf.abs
# denovo_lowConfidence.abundance_matrix.rel <- og_counts.low_conf.rel
# denovo_lowConfidence.shared_abundances    <- shared_abundances.low_conf
# 
# rownames( alignment.abundance_matrix.abs ) <- rownames( denovo_highConfidence.abundance_matrix.abs ) <- rownames( denovo_lowConfidence.abundance_matrix.abs ) <- paste0("Cluster",c(1:dim(alignment.abundance_matrix.abs)[1]))
# rownames( alignment.abundance_matrix.rel ) <- rownames( denovo_highConfidence.abundance_matrix.rel ) <- rownames( denovo_lowConfidence.abundance_matrix.rel ) <- paste0("Cluster",c(1:dim(alignment.abundance_matrix.abs)[1]))
# names( alignment.shared_abundances ) <- names( denovo_highConfidence.shared_abundances ) <- names( denovo_lowConfidence.shared_abundances ) <- paste0("Cluster",c(1:dim(alignment.abundance_matrix.abs)[1]))
# 
# dist.bray<-vegdist(t(emotion_abund_matrix), method="bray")
# 
# library("vegan")
# #plot.new()
# #par(mfrow=c(5,10))
# dists <- matrix(nrow=50,ncol=6)
# colnames(dists) <- c("alignment_V_deNovo_highConfidence",            "alignment_V_deNovo_highConfidence.comparableCount",
#                      "alignment_V_deNovolowConfidence",              "alignment_V_deNovolowConfidence.comparableCount",
#                      "deNovo_highConfidence_vs_deNovo_lowConfidence","deNovo_highConfidence_vs_deNovo_lowConfidence.comparableCount"
#                      ) # ,"alignPop","dnHiPop","dnLoPop")
# rownames(dists) <- colnames(alignment.abundance_matrix.rel)
# dists.bray <- dists
# #png(filename = "q.png", height = 5000,width = 5000)
# for(q in 1:50){
#   #if(q %% 200 ==TRUE) print(q/dim(alignment.abundance_matrix.rel)[1])
#   matQ       <- matrix(nrow=dim(alignment.abundance_matrix.rel)[1],ncol=3)
#   matQ[,1]   <- unname(unlist(alignment.abundance_matrix.rel[,q]))
#   matQ[,2]   <- unname(unlist(denovo_highConfidence.abundance_matrix.rel[,q]))
#   matQ[,3]   <- unname(unlist(denovo_lowConfidence.abundance_matrix.rel[,q]))
#   dist1      <- vegdist(t(matQ), method = "euclidean")
#   dist1.bray <- vegdist(t(matQ), method = "bray")
#   #dist1t      <- vegdist(matQ, method = "euclidean")
#   #heatmap.2(as.matrix(dist1.bray), #main=q,
#   #          trace='none') #, col=c(rev(gray.colors(150)),matlab.like2(150)),keysize=0.2) # margins = c(1,1)
#   dists[q,]      <- c( as.matrix(dist1)[1,2]      , sum(matQ[,1] > 0 & matQ[,2] > 0),
#                        as.matrix(dist1)[1,3]      , sum(matQ[,1] > 0 & matQ[,3] > 0),
#                        as.matrix(dist1)[2,3]      , sum(matQ[,2] > 0 & matQ[,3] > 0) )
#   dists.bray[q,] <- c( as.matrix(dist1.bray)[1,2]      , sum(matQ[,1] > 0 & matQ[,2] > 0),
#                        as.matrix(dist1.bray)[1,3]      , sum(matQ[,1] > 0 & matQ[,3] > 0),
#                        as.matrix(dist1.bray)[2,3]      , sum(matQ[,2] > 0 & matQ[,3] > 0) )
# }
# #dev.off()
# png(filename = "reference_distance_matrix.Euclidean.png", height = 1000,width = 1400)
# heatmap.2(t(dists[order(dists[,1]),c(1,3,5)]),
#           col=redgreen, trace='n',margins=c(25,30),
#           cexRow=1.4,
#           cexCol=2, keysize = 0.4)
# dev.off()
# write.table(t(dists[order(dists[,1]),c(1,3,5)]), file="reference_distance_matrix.Euclidean.png")
# 
# png(filename = "reference_distance_matrix.Bray.png", height = 1000,width = 1400)
# heatmap.2(t(dists[order(dists.bray[,1]),c(1,3,5)]),
#           col=redgreen, trace='n',margins=c(25,30),
#           cexRow=1.4,
#           cexCol=2, keysize = 0.4)
# dev.off()
# write.table(t(dists[order(dists.bray[,1]),c(1,3,5)]), file="reference_distance_matrix.Bray.png")
# 
# 
# 
# 
# 
# 
# 
# dists.trans <- matrix(nrow=dim(alignment.abundance_matrix.rel)[1],ncol=6)
# colnames(dists.trans) <- c("alignment_V_deNovo_highConfidence",            "alignment_V_deNovo_highConfidence.comparableCount",
#                      "alignment_V_deNovolowConfidence",              "alignment_V_deNovolowConfidence.comparableCount",
#                      "deNovo_highConfidence_vs_deNovo_lowConfidence","deNovo_highConfidence_vs_deNovo_lowConfidence.comparableCount"
# ) # ,"alignPop","dnHiPop","dnLoPop")
# rownames(dists.trans) <- rownames(alignment.abundance_matrix.rel)
# dists.trans.bray <- dists
# #png(filename = "q.png", height = 5000,width = 5000)
# for(q in 1:dim(alignment.abundance_matrix.rel)[1]){
#   if(q %% 200 ==TRUE) print(q/dim(alignment.abundance_matrix.rel)[1])
#   matQ       <- matrix(nrow=50,ncol=3)
#   matQ[,1]   <- unname(unlist(alignment.abundance_matrix.rel[q,]))
#   matQ[,2]   <- unname(unlist(denovo_highConfidence.abundance_matrix.rel[q,]))
#   matQ[,3]   <- unname(unlist(denovo_lowConfidence.abundance_matrix.rel[q,]))
#   dists1      <- vegdist(t(matQ), method = "euclidean")
#   dists1bray <- vegdist(t(matQ), method = "bray")
# 
#   #dist1t      <- vegdist(matQ, method = "euclidean")
#   #heatmap.2(as.matrix(dist1.bray), #main=q,
#   #          trace='none') #, col=c(rev(gray.colors(150)),matlab.like2(150)),keysize=0.2) # margins = c(1,1)
#   dists.trans[q,] <- c( as.matrix(dists1)[1,2]      , sum(matQ[,1] > 0 & matQ[,2] > 0),
#                        as.matrix(dists1)[1,3]      , sum(matQ[,1] > 0 & matQ[,3] > 0),
#                        as.matrix(dists1)[2,3]      , sum(matQ[,2] > 0 & matQ[,3] > 0 ))
# 
#   #if(sum(is.na(dists1bray)) == 0){
#   #  dists.trans.bray[q,] <- c( as.matrix(dists.trans.bray)[1,2]      , sum(matQ[,1] > 0 & matQ[,2] > 0),
#   #                       as.matrix(dists.trans.bray)[1,3]      , sum(matQ[,1] > 0 & matQ[,3] > 0),
#   #                       as.matrix(dists.trans.bray)[2,3]      , sum(matQ[,2] > 0 & matQ[,3] > 0) )
#   #}
# 
# }
# png(filename = "alignment_V_deNovo_highConfidence_SimilarAbundance_distance.Euclidean.png", height = 1000,width = 1400)
# heatmap.2( dists.trans[ rev(order(unname(dists.trans[ , "alignment_V_deNovo_highConfidence" ]))) , c(1,3,5) ][ 1:100 , ],
#            margins=c(30,10),
#            cexCol=1.4)
# dev.off()
# write.table(dists.trans[ rev(order(unname(dists.trans[ , "alignment_V_deNovo_highConfidence" ]))) , ],
#             file="transcript_distance_matrix.Euclidean.tsv")
# 
# s_withLabels2 <- as.data.frame(s_withLabels)
# s_withLabels2$alignment_vs_denovoHi.dist <- rep(NA,dim(s_withLabels2)[1])
# s_withLabels2$alignment_vs_denovoLo.dist <- rep(NA,dim(s_withLabels2)[1])
# s_withLabels2$denovoHi_vs_denovoLo.dist <- rep(NA,dim(s_withLabels2)[1])
# s_withLabels2[ rownames(dists.trans)[rownames(dists.trans) %in% rownames(s_withLabels2)] , "alignment_vs_denovoHi.dist"] <- dists.trans[ rownames(dists.trans)[rownames(dists.trans) %in% rownames(s_withLabels)] , 1 ]
# s_withLabels2[ rownames(dists.trans)[rownames(dists.trans) %in% rownames(s_withLabels2)] , "alignment_vs_denovoLo.dist"] <- dists.trans[ rownames(dists.trans)[rownames(dists.trans) %in% rownames(s_withLabels)] , 3 ]
# s_withLabels2[ rownames(dists.trans)[rownames(dists.trans) %in% rownames(s_withLabels2)] , "denovoHi_vs_denovoLo.dist"] <- dists.trans[ rownames(dists.trans)[rownames(dists.trans) %in% rownames(s_withLabels)] , 5 ]
# write.table(s_withLabels2, file="summary_matrix.clean.dist.20190122.tsv", sep="\t")
# 
#            #dev.off()
# 
# 
# 
# 
# # s_withLabels$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval`
# r1 <- rownames(s_withLabels[ order(s_withLabels2$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval`) , ])
# r2 <- rownames(s_withLabels[ order(s_withLabels2$`alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval`) , ])
# intop1k <- r1[1:1000][r1[1:1000] %in% r2[1:1000]]
# order(r1)[r1 %in% intop1k]
# 
# r1f <- rownames(s_withLabels[ order(s_withLabels2$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval`) , ])
# r2f <- rownames(s_withLabels[ order(s_withLabels2$`alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval`) , ])
# intop1kf <- r1f[1:1000][r1f[1:1000] %in% r2f[1:1000]]
# 
# inBoth <- intop1k[intop1k %in% intop1kf]
# order(r1)[r1 %in% inBoth]
# 
# rD <- s_withLabels2[inBoth,]$cutoff.denovo
# rA <- s_withLabels2[inBoth,]$cutoff.align
# names(rD) <- names(rA) <- rownames(s_withLabels2[inBoth,])
# inBothAndNonZero.d <- inBoth[inBoth %in% names(rD[rD > 0])]
# inBothAndNonZero.a <- inBoth[inBoth %in% names(rA[rA > 0])]
# 
# inBothAndNonZero <- inBothAndNonZero.d[inBothAndNonZero.d %in% inBothAndNonZero.a]
# order(r2)[r2 %in% inBothAndNonZero]
# 
# inBoth.a <- r2[1:1000][r2[1:1000] %in% r2f[1:1000]]
# inBothAndNonZero.a <- inBoth.a[inBoth.a %in% names(rA[rA > 0])]
# order(r1)[r1 %in% inBothAndNonZero.a]
# 
# 
# r5 <- s_withLabels2$`alignmentPlotsAndTables/alignment.3.CAPER.logwMLSresid_v_binaryAbsreport.tsv.rawPval`
# r6 <- s_withLabels2$`alignmentPlotsAndTables/alignment.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport.tsv.rawPval`
# inBothAndNonZero[ inBothAndNonZero %in% r5[1:1000] ]
# inBothAndNonZero[ inBothAndNonZero %in% r6[1:1000] ]
# 
# r3 <- s_withLabels2[inBoth,]$alignment_vs_denovoHi.dist


# save.image(file="20190118.Rdata")

load("20190118.Rdata")
summary_matrix.12.20190130

s_withLabels3 <- s_withLabels2

ri <- 1
for(rn in rownames(s_withLabels3)){
    
  if( ri %% 200 == 0 ){ print ( ri/length(rownames(s_withLabels3)) ) }
  for(w in c(".rawPval",".fdrPval",".slope",".intercept")){  
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","CAPER.","logwMLSresid_v_relAbs.report.tsv",w) ]              #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.10.LM.logwMLSresid_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","LM.","logwMLSresid_v_relAbs.report.tsv",w) ]   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.11.LM.logwMLSresid_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","LM.","logwMLSresid_v_absAbs.report.tsv",w) ]                 #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.12.LM.logwMLSresid_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","LM.","logwMLSresid_v_binaryAbs.report.tsv",w) ]              #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.13.LM.wMLSresid_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","LM.","wMLSresid_v_relAbs.report.tsv",w) ]                    #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.14.LM.wMLSresid_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","LM.","wMLSresid_v_absAbs.report.tsv",w) ]                    #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.15.LM.wMLSresid_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","LM.","wMLSresid_v_binaryAbs.report.tsv",w) ]                 #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.16.LM.MLS_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","LM.","MLS_v_relAbs.report.tsv",w) ]                          #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.17.LM.MLS_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","LM.","MLS_v_absAbs.report.tsv",w) ]                          #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.18.LM.MLS_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","LM.","MLS_v_binaryAbs.report.tsv",w) ]                       #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.2.CAPER.logwMLSresid_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","CAPER.","logwMLSresid_v_absAbs.report.tsv",w) ]              #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.3.CAPER.logwMLSresid_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","CAPER.","logwMLSresid_v_binaryAbs.report.tsv",w) ]           #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.4.CAPER.wMLSresid_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","CAPER.","wMLSresid_v_relAbs.report.tsv",w) ]                #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.5.CAPER.wMLSresid_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","CAPER.","wMLSresid_v_absAbs.report.tsv",w) ]                #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.6.CAPER.wMLSresid_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","CAPER.","wMLSresid_v_binaryAbs.report.tsv",w) ]             #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.7.CAPER.MLS_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","CAPER.","MLS_v_relAbs.report.tsv",w) ]                      #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.8.CAPER.MLS_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","CAPER.","MLS_v_absAbs.report.tsv",w) ]                      #   
    
    s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.9.CAPER.MLS_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","alignment.","CAPER.","MLS_v_binaryAbs.report.tsv",w) ]                   #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","CAPER.","logwMLSresid_v_relAbs.report.tsv",w) ]                #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.10.LM.logwMLSresid_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","LM.","logwMLSresid_v_relAbs.report.tsv",w) ]                   #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.11.LM.logwMLSresid_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","LM.","logwMLSresid_v_absAbs.report.tsv",w) ]                   #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.12.LM.logwMLSresid_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","LM.","logwMLSresid_v_binaryAbs.report.tsv",w) ]                 #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.13.LM.wMLSresid_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","LM.","wMLSresid_v_relAbs.report.tsv",w) ]                      #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.14.LM.wMLSresid_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","LM.","wMLSresid_v_absAbs.report.tsv",w) ]                      #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.15.LM.wMLSresid_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","LM.","wMLSresid_v_binaryAbs.report.tsv",w) ]                   #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.16.LM.MLS_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","LM.","MLS_v_relAbs.report.tsv",w) ]                            #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.17.LM.MLS_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","LM.","MLS_v_absAbs.report.tsv",w) ]                            #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.18.LM.MLS_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","LM.","MLS_v_binaryAbs.report.tsv",w) ]                         #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.2.CAPER.logwMLSresid_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","CAPER.","wMLSresid_v_absAbs.report.tsv",w) ]                    # 
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.3.CAPER.logwMLSresid_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","CAPER.","logwMLSresid_v_binaryAbs.report.tsv",w) ]            #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.4.CAPER.wMLSresid_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","CAPER.","wMLSresid_v_relAbs.report.tsv",w) ]                 #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.5.CAPER.wMLSresid_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","CAPER.","wMLSresid_v_absAbs.report.tsv",w) ]                 #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.6.CAPER.wMLSresid_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","CAPER.","wMLSresid_v_binaryAbs.report.tsv",w) ]              #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.7.CAPER.MLS_v_relAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","CAPER.","MLS_v_relAbs.report.tsv",w) ]                       #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.8.CAPER.MLS_v_absAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","CAPER.","MLS_v_absAbs.report.tsv",w) ]                       #   
    
    s_withLabels3[ rn , paste0("deNovoPlotsAndTables/denovo_highConfidence.9.CAPER.MLS_v_binaryAbsreport.tsv",w) ] <- summary_matrix.12.20190130[ rn , paste0("RERUN.","denovo_highConfidence.","CAPER.","MLS_v_binaryAbs.report.tsv",w) ]                     #   
     
    if(w == ".rawPval"){
      #print(paste(
      #            rn,
      #            s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv",w) ],
      #            s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.4.CAPER.wMLSresid_v_relAbsreport.tsv",w) ],
      #            s_withLabels3[ rn , paste0("alignmentPlotsAndTables/alignment.7.CAPER.MLS_v_relAbsreport.tsv",w) ]
      #           )
      #    )
    }
  }
  ri <- ri + 1 
}
  
write.table(s_withLabels3, file="20190201_AM.tsv", sep = "\t")

s_withLabels4 <- s_withLabels3[,!grepl("SortID", colnames(s_withLabels3))]

write.table(s_withLabels4, file="20190201_AM.2.tsv", sep = "\t")


# 
# 
# 
#  
# #sum_mat.pop[,"fdr_pvalues"] <- p.adjust( as.numeric(sum_mat.pop[,3] ) , method="fdr" )
# 
# # save.image("20190201.Rdata")
# 
# s_withLabels2.orig <- s_withLabels2
# s_withLabels2 <- s_withLabels3
# 
# ddd <- s_withLabels2[ !is.na(s_withLabels2$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval`) , ]
# 
# www                   <- s_withLabels2[ , 3:80 ]
# www[is.na(www)]       <- 1
# colSums(www < 0.005)
# 
# 
# 
# sum(www[,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.0001)
# sum(www[,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.0005)
# sum(www[,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.001)
# sum(www[,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.005)
# sum(s_withLabels2[ , "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.05, na.rm=TRUE)
# sum(s_withLabels2[ , "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.01, na.rm=TRUE)
# 
# 
# sum(s_withLabels2[s_withLabels2$cutoff.align == 0,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.0001, na.rm=TRUE)
# sum(s_withLabels2[s_withLabels2$cutoff.align == 0,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.0005, na.rm=TRUE)
# sum(s_withLabels2[s_withLabels2$cutoff.align == 0,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.001, na.rm=TRUE)
# sum(s_withLabels2[s_withLabels2$cutoff.align == 0,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.005, na.rm=TRUE)
# sum(s_withLabels2[s_withLabels2$cutoff.align == 0 , "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.05, na.rm=TRUE)
# sum(s_withLabels2[s_withLabels2$cutoff.align == 0 , "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.01, na.rm=TRUE)
# 
# 
# sum(s_withLabels2[s_withLabels2$cutoff.align <= 4,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.0001, na.rm=TRUE)
# sum(s_withLabels2[s_withLabels2$cutoff.align  <= 4,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.0005, na.rm=TRUE)
# sum(s_withLabels2[s_withLabels2$cutoff.align  <= 4,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.001, na.rm=TRUE)
# sum(s_withLabels2[s_withLabels2$cutoff.align  <= 4,"alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval"] < 0.005, na.rm=TRUE)
# sum(s_withLabels2[s_withLabels2$cutoff.align  <= 4 , "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.05, na.rm=TRUE)
# sum(s_withLabels2[s_withLabels2$cutoff.align  <= 4 , "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.01, na.rm=TRUE)
# 
# sum( s_withLabels2[,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.0001, na.rm=TRUE)
# sum( s_withLabels2[,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.0005, na.rm=TRUE)
# sum( s_withLabels2[,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.001, na.rm=TRUE)
# sum( s_withLabels2[,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.005, na.rm=TRUE)
# sum( s_withLabels2[,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.05, na.rm=TRUE)
# sum( s_withLabels2[,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.01, na.rm=TRUE)
# 
# 
# sum( s_withLabels2[s_withLabels2$cutoff.denovo == 0,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.0001, na.rm=TRUE)
# sum( s_withLabels2[s_withLabels2$cutoff.denovo == 0,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.0005, na.rm=TRUE)
# sum( s_withLabels2[s_withLabels2$cutoff.denovo == 0,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.001, na.rm=TRUE)
# sum( s_withLabels2[s_withLabels2$cutoff.denovo == 0,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.005, na.rm=TRUE)
# sum( s_withLabels2[s_withLabels2$cutoff.denovo == 0,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.05, na.rm=TRUE)
# sum( s_withLabels2[s_withLabels2$cutoff.denovo == 0,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.01, na.rm=TRUE)
# 
# sum( s_withLabels2[s_withLabels2$cutoff.denovo <= 12,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.0001, na.rm=TRUE)
# sum( s_withLabels2[s_withLabels2$cutoff.denovo <= 12,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.0005, na.rm=TRUE)
# sum( s_withLabels2[s_withLabels2$cutoff.denovo <= 12,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.001, na.rm=TRUE)
# sum( s_withLabels2[s_withLabels2$cutoff.denovo <= 12,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.005, na.rm=TRUE)
# sum( s_withLabels2[s_withLabels2$cutoff.denovo <= 12,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.05, na.rm=TRUE)
# sum( s_withLabels2[s_withLabels2$cutoff.denovo <= 12,"deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval" ] < 0.01, na.rm=TRUE)
# 
# #extremeLongevity_sd1_hi.logwMLSresid_v_relAbsreport.tsv.fdrPval
# #extremeLongevity_sd1_hi.logwMLSresid_v_relAbsreport.tsv.fdrPval
# 
# 
# 
# 
# colSums( s_withLabels2[, grepl("fdrPval",colnames(s_withLabels2)) ] < 0.05, na.rm=TRUE)[ !grepl( "LM", names(colSums( s_withLabels2[, grepl("fdrPval",colnames(s_withLabels2)) ] < 0.05, na.rm=TRUE) )) ][ colSums( s_withLabels2[, grepl("fdrPval",colnames(s_withLabels2)) ] < 0.05, na.rm=TRUE)[ !grepl( "LM", names(colSums( s_withLabels2[, grepl("fdrPval",colnames(s_withLabels2)) ] < 0.05, na.rm=TRUE) )) ] > 0  ]
# 
# 
# colSums( s_withLabels2[, grepl("fdrPval",colnames(s_withLabels2)) ] < 0.01, na.rm=TRUE)[ !grepl( "LM", names(colSums( s_withLabels2[, grepl("fdrPval",colnames(s_withLabels2)) ] < 0.01, na.rm=TRUE) )) ][ colSums( s_withLabels2[, grepl("fdrPval",colnames(s_withLabels2)) ] < 0.05, na.rm=TRUE)[ !grepl( "LM", names(colSums( s_withLabels2[, grepl("fdrPval",colnames(s_withLabels2)) ] < 0.05, na.rm=TRUE) )) ] > 0  ]
# 
#  
# s_withLabels
# 
# 
# sq <- s_withLabels2[ !is.na(s_withLabels2[ "deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ]) , ]
# sq <- sq[ sq[ "deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.005 , ]
# denovolist       <- unique(unlist(strsplit(sq$geneSymbols," "))[ unlist(strsplit(sq$geneSymbols," ")) != "" ])
# denovolist.huref <- denovolist[ !grepl("_", denovolist) ]
# denovolist.hurefClean <- list()
# for(b in strsplit(denovolist.huref,"-")){
#   denovolist.hurefClean <- append(denovolist.hurefClean, b[1])
# }
# denovolist.hurefClean <- unlist(denovolist.hurefClean)
# 
# sq <- s_withLabels2[ !is.na(s_withLabels2[ "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ]) , ]
# sq <- sq[ sq[ "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval" ] < 0.005 , ]
# alignlist       <- unique(unlist(strsplit(sq$geneSymbols," "))[ unlist(strsplit(sq$geneSymbols," ")) != "" ])
# alignlist.huref <- alignlist[ !grepl("_", alignlist) ]
# alignlist.hurefClean <- list()
# for(b in strsplit(alignlist.huref,"-")){
#   alignlist.hurefClean <- append(alignlist.hurefClean, b[1])
# }
# alignlist.hurefClean <- unlist(alignlist.hurefClean)
# 
# 
# sum( denovolist.hurefClean %in% alignlist.hurefClean )
# 
# write.table( as.matrix(denovolist.hurefClean) , file = "denovolist.hurefClean.txt" )
# write.table( as.matrix(alignlist.hurefClean) , file = "alignlist.hurefClean.txt" )
# write.table( as.matrix(denovolist.hurefClean[denovolist.hurefClean %in% alignlist.hurefClean]) , file = "alignANDdenovolist.hurefClean.txt" )
# 
# 
# 
# s_small <- s_withLabels2[ order(s_withLabels2$`alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval`), c( "deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval",
#                                                                                                                                        "deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval",
#                                                                                                                                        "deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.slope",
#                                                                                                                                        "deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.intercept",
#                                                                                                                                        "cutoff.denovo",
#                                                                                                                                        "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval",
#                                                                                                                                        "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval",
#                                                                                                                                        "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.slope",
#                                                                                                                                        "alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.intercept",
#                                                                                                                                        "cutoff.align",
#                                                                                                                                        "alignment_vs_denovoHi.dist",
#                                                                                                                                        "proteinIDs",                                                                                                               
#                                                                                                                                        "transcriptIDs",
#                                                                                                                                        "geneSymbols"
# )
# ]
# 
# write.table( s_small , file = "briefSum.20190130.tsv", sep = "\t" )
# write.table( s_withLabels2 , file = "longSum.20190130.tsv", sep = "\t" )
#  
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 


# START 1/27/2019
#library(colorRamps)
#library(gplots)
#setwd ("~/Desktop/develop.Avian")
#load(file="20190118.Rdata")

Ma.B.MLS     <- read.table("elife-19130-table1-data1-v2.B.MLS.txt", sep="\t", header=TRUE, fill=TRUE)
Ma.D.resMLS  <- read.table("elife-19130-table1-data1-v2.D.resMLS.txt", sep="\t", header=TRUE, fill=TRUE)

Avian.CAPER.DeNovo.logMLSresid <- rownames(s_withLabels4[ order(s_withLabels4$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval`) , ])
Avian.CAPER.Align.logMLSresid  <- rownames(s_withLabels4[ order(s_withLabels4$`alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval`) , ])

Avian.CAPER.DeNovo.logMLSresid.FDR  <- rownames(s_withLabels4[ order(s_withLabels4$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval`) , ])
Avian.CAPER.Align.logMLSresid.FDR   <- rownames(s_withLabels4[ order(s_withLabels4$`alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.fdrPval`) , ])

sum(toupper(as.character(Ma.B.MLS$Symbol)) %in% unlist(strsplit(s_withLabels4[ Avian.CAPER.DeNovo.logMLSresid , "geneSymbols" ], " "))[ unlist(strsplit(s_withLabels4[ Avian.CAPER.DeNovo.logMLSresid , "geneSymbols" ], " ")) != "" ])

shared_with_ma <- toupper(as.character(Ma.B.MLS$Symbol))[ toupper(as.character(Ma.B.MLS$Symbol)) %in% unlist(strsplit(s_withLabels4[ Avian.CAPER.DeNovo.logMLSresid , "geneSymbols" ], " "))[ unlist(strsplit(s_withLabels2[ Avian.CAPER.DeNovo.logMLSresid , "geneSymbols" ], " ")) != "" ] ]

Ma.B.MLS.sorted <- Ma.B.MLS[ toupper(as.character(Ma.B.MLS$Symbol)) %in% shared_with_ma , ][ order(Ma.B.MLS[ toupper(as.character(Ma.B.MLS$Symbol)) %in% shared_with_ma , "p.value.all" ]) , ]
Ma.D.resMLS.sorted <- Ma.B.MLS[ toupper(as.character(Ma.D.resMLS$Symbol)) %in% shared_with_ma , ][ order(Ma.D.resMLS[ toupper(as.character(Ma.D.resMLS$Symbol)) %in% shared_with_ma , "p.value.all" ]) , ]
rownames_shared_with_ma <- list()
rownames_shared_with_ma.concat <- list()
si <- 1
for(s in shared_with_ma){
  rownames_shared_with_ma <- append( x = rownames_shared_with_ma, values = rownames(s_withLabels4[ grepl( s , s_withLabels4$geneSymbols ) , ]) )
  rownames_shared_with_ma.concat[si] <- paste(rownames(s_withLabels4[ grepl( s , s_withLabels4$geneSymbols ) , ]))
  si <- si + 1 
}

rownames_shared_with_ma.concat <- unlist(rownames_shared_with_ma.concat)
names(rownames_shared_with_ma.concat) <- shared_with_ma

ma_sum_matrix <- matrix(nrow=length(rownames(s_withLabels4)),
                          #length(shared_with_ma),
                        ncol=23
                       )
#rownames(ma_sum_matrix) <- c(1:length(shared_with_ma))
rownames(ma_sum_matrix) <- c(rownames(s_withLabels4))
colnames(ma_sum_matrix) <- c("huref_gene_symbol",
                             "OrthoMCL_ClusterID",
                             "align_queryZero",
                             "denovo_queryZero",
                             "alignment_vs_denovoHi.dist",
                             "MLS.p_value.Ma",
                             "MLS.p_value.Align",
                             "MLS.p_value.DeNovo",
                             "resMLS.p_value.Ma",
                             "resMLS.p_value.Align",
                             "resMLS.p_value.DeNovo",
                             "reslogMLS.p_value.Align",
                             "reslogMLS.p_value.DeNovo",
                             "reslogMLS_binary.p_value.Align" , 
                             "reslogMLS_binary.p_value.DeNovo" , 
                             "reslogMLS_sd2.p_value.Align" , 
                             "reslogMLS_sd2.p_value.DeNovo" , 
                             "reslogMLS_sd1.p_value.Align" , 
                             "reslogMLS_sd1.p_value.DeNovo",
                             "proteinIDs",                                                                                                                   
                             "transcriptIDs",
                             "ma_rank.MLS",
                             "ma_rank.resMLS"
                            )
si <- 1
#for(s in shared_with_ma){
for(id in rownames(s_withLabels4)){
  
  if(si %% 50 == TRUE){ print(si/sum(s_withLabels4$geneSymbols != "")) }
  
  for(sym in unlist(strsplit(s_withLabels4[ id , "geneSymbols" ]," "))[unlist(strsplit(s_withLabels4[ id , "geneSymbols" ]," ")) != ""]){
    
    #print(paste(si, sym))
    
    #ma_sum_matrix[si,"huref_gene_symbol"]             <- s
    ma_sum_matrix[id,"huref_gene_symbol"]             <- s <- sym
    
    #id <- unname(rownames_shared_with_ma.concat[s])
    ma_sum_matrix[id,"OrthoMCL_ClusterID"]            <- id
    
    ma_sum_matrix[id,"ma_rank.MLS"]                   <- match(s, toupper(Ma.B.MLS.sorted$Symbol))
  
    ma_sum_matrix[id,"proteinIDs"]                    <- s_withLabels4$proteinIDs[ rownames(s_withLabels4) == id ]                                                                                           
    ma_sum_matrix[id,"transcriptIDs"]                 <- s_withLabels4$transcriptIDs[ rownames(s_withLabels4) == id ]
    
    ma_sum_matrix[id,"ma_rank.resMLS"]                <- match(s, toupper(Ma.D.resMLS.sorted$Symbol))
    
    ma_sum_matrix[id,"align_queryZero"]               <- s_withLabels4[rownames(s_withLabels4) == id,]$cutoff.align
    ma_sum_matrix[id,"denovo_queryZero"]              <- s_withLabels4[rownames(s_withLabels4) == id,]$cutoff.denovo
    ma_sum_matrix[id,"alignment_vs_denovoHi.dist"]    <- s_withLabels4[rownames(s_withLabels4) == id,]$alignment_vs_denovoHi.dist
    
    if(sum(toupper(Ma.B.MLS$Symbol) %in% s) > 0){
      ma_sum_matrix[id,"MLS.p_value.Ma"]              <- Ma.B.MLS[ toupper(Ma.B.MLS$Symbol) == s , "p.value.all" ]
    }
    if(sum(toupper(Ma.D.resMLS$Symbol) %in% s) > 0){
      ma_sum_matrix[id,"resMLS.p_value.Ma"]           <- Ma.D.resMLS[ toupper(Ma.D.resMLS$Symbol) == s , "p.value.all" ]
    }
    
    ma_sum_matrix[id,"MLS.p_value.Align"]               <- s_withLabels4$`alignmentPlotsAndTables/alignment.7.CAPER.MLS_v_relAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ]
    ma_sum_matrix[id,"MLS.p_value.DeNovo"]              <- s_withLabels4$`deNovoPlotsAndTables/denovo_highConfidence.7.CAPER.MLS_v_relAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ]
    ma_sum_matrix[id,"resMLS.p_value.Align"]            <- s_withLabels4$`alignmentPlotsAndTables/alignment.4.CAPER.wMLSresid_v_relAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ]
    ma_sum_matrix[id,"resMLS.p_value.DeNovo"]           <- s_withLabels4$`deNovoPlotsAndTables/denovo_highConfidence.4.CAPER.wMLSresid_v_relAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ]
    ma_sum_matrix[id,"reslogMLS.p_value.Align"]         <- s_withLabels4$`alignmentPlotsAndTables/alignment.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ]
    ma_sum_matrix[id,"reslogMLS.p_value.DeNovo"]        <- s_withLabels4$`deNovoPlotsAndTables/denovo_highConfidence.1.CAPER.logwMLSresid_v_relAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ] 
    
    ma_sum_matrix[id,"reslogMLS_binary.p_value.Align"]  <- s_withLabels4$`alignmentPlotsAndTables/alignment.3.CAPER.logwMLSresid_v_binaryAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ] 
    ma_sum_matrix[id,"reslogMLS_binary.p_value.DeNovo"] <- s_withLabels4$`deNovoPlotsAndTables/denovo_highConfidence.3.CAPER.logwMLSresid_v_binaryAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ] 
    
    ma_sum_matrix[id,"reslogMLS_sd2.p_value.Align"]  <- s_withLabels4$`alignmentPlotsAndTables/alignment.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ] 
    ma_sum_matrix[id,"reslogMLS_sd2.p_value.DeNovo"] <- s_withLabels4$`deNovoPlotsAndTables/denovo_highConfidence.19.extremeLongevity_sd2_all.logwMLSresid_v_relAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ] 
    
    ma_sum_matrix[id,"reslogMLS_sd1.p_value.Align"]  <- s_withLabels4$`alignmentPlotsAndTables/alignment.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ] 
    ma_sum_matrix[id,"reslogMLS_sd1.p_value.DeNovo"] <- s_withLabels4$`deNovoPlotsAndTables/denovo_highConfidence.25.extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbsreport.tsv.rawPval`[ rownames(s_withLabels4) == id ] 
    
    si <- si + 1
    
  }
  
}
 
write.table(ma_sum_matrix, "sum_matrix_AM_20190201.tsv", sep = "\t")


masummatreduced <- ma_sum_matrix[,colSums(!is.na(ma_sum_matrix))>0]
rownames(masummatreduced) <- masummatreduced[,"OrthoMCL_ClusterID"]
#masummatreduced <- masummatreduced[,!grepl("rank", colnames(masummatreduced))]
masummatreduced.ordered <- masummatreduced[ order(as.double(masummatreduced[ , "reslogMLS.p_value.Align"])) , ] 

masummatreduced.ordered.df                                            <- as.data.frame(masummatreduced.ordered)[, !grepl("rank",colnames(masummatreduced.ordered))]

masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "reslogMLS.p_value.Align" ]))) , ]
masummatreduced.ordered.df$reslogMLS.rank.Align <- rep(NA, dim(masummatreduced.ordered.df)[1])
z <- rownames( masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df$reslogMLS.p_value.Align) , ] )
masummatreduced.ordered.df[ z , "reslogMLS.rank.Align" ] <- order(as.double(as.character(masummatreduced.ordered.df[ z, "reslogMLS.p_value.Align" ])))

masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "MLS.p_value.DeNovo" ]))) , ]
masummatreduced.ordered.df$MLS.rank.DeNovo <- rep(NA, dim(masummatreduced.ordered.df)[1])
z <- rownames( masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df$MLS.p_value.DeNovo) , ] )
masummatreduced.ordered.df[ z , "MLS.rank.DeNovo" ] <- order(as.double(as.character(masummatreduced.ordered.df[ z, "MLS.p_value.DeNovo" ])))

#
masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "MLS.p_value.Ma" ]))) , ]
masummatreduced.ordered.df$MLS.rank.Ma <- rep(NA, dim(masummatreduced.ordered.df)[1])
z <- rownames( masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df$MLS.p_value.Ma) , ] )
masummatreduced.ordered.df[ z , "MLS.rank.Ma" ] <- order(as.double(as.character(masummatreduced.ordered.df[ z, "MLS.p_value.Ma" ])))

#
masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "resMLS.p_value.Ma" ]))) , ]
masummatreduced.ordered.df$resMLS.rank.Ma <- rep(NA, dim(masummatreduced.ordered.df)[1])
z <- rownames( masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df$resMLS.p_value.Ma) , ] )
masummatreduced.ordered.df[ z , "resMLS.rank.Ma" ] <- order(as.double(as.character(masummatreduced.ordered.df[ z, "resMLS.p_value.Ma" ])))

#
masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "resMLS.p_value.Align" ]))) , ]
masummatreduced.ordered.df$resMLS.rank.Align <- rep(NA, dim(masummatreduced.ordered.df)[1])
z <- rownames( masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df$resMLS.p_value.Align) , ] )
masummatreduced.ordered.df[ z , "resMLS.rank.Align" ] <- order(as.double(as.character(masummatreduced.ordered.df[ z, "resMLS.p_value.Align" ])))

#
masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "resMLS.p_value.DeNovo" ]))) , ]
masummatreduced.ordered.df$resMLS.rank.DeNovo <- rep(NA, dim(masummatreduced.ordered.df)[1])
z <- rownames( masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df$resMLS.p_value.DeNovo) , ] )
masummatreduced.ordered.df[ z , "resMLS.rank.DeNovo" ] <- order(as.double(as.character(masummatreduced.ordered.df[ z, "resMLS.p_value.DeNovo" ])))

#
masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "reslogMLS.p_value.Align" ]))) , ]
masummatreduced.ordered.df$reslogMLS.rank.Align <- rep(NA, dim(masummatreduced.ordered.df)[1])
z <- rownames( masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df$reslogMLS.p_value.Align) , ] )
masummatreduced.ordered.df[ z , "reslogMLS.rank.Align" ] <- order(as.double(as.character(masummatreduced.ordered.df[ z, "reslogMLS.p_value.Align" ])))

#
masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "reslogMLS.p_value.DeNovo" ]))) , ]
masummatreduced.ordered.df$reslogMLS.rank.DeNovo <- rep(NA, dim(masummatreduced.ordered.df)[1])
z <- rownames( masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df$reslogMLS.p_value.DeNovo) , ] )
masummatreduced.ordered.df[ z , "reslogMLS.rank.DeNovo" ] <- order(as.double(as.character(masummatreduced.ordered.df[ z, "reslogMLS.p_value.DeNovo" ])))

#
#masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "reslogMLS_binary.p_value.Align" ]))) , ]
#masummatreduced.ordered.df$reslogMLS_binary.rank.Align <- rep(NA, dim(masummatreduced.ordered.df)[1])
#z <- rownames( masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df$reslogMLS_binary.p_value.Align) , ] )
#masummatreduced.ordered.df[ z , "reslogMLS_binary.rank.Align" ] <- order(as.double(as.character(masummatreduced.ordered.df[ z, "reslogMLS_binary.p_value.Align" ])))

#
masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "reslogMLS_binary.p_value.DeNovo" ]))) , ]
masummatreduced.ordered.df$reslogMLS_binary.rank.DeNovo <- rep(NA, dim(masummatreduced.ordered.df)[1])
z <- rownames( masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df$reslogMLS_binary.p_value.DeNovo) , ] )
masummatreduced.ordered.df[ z , "reslogMLS_binary.rank.DeNovo" ] <- order(as.double(as.character(masummatreduced.ordered.df[ z, "reslogMLS_binary.p_value.DeNovo" ])))
 
masummatreduced.ordered.df                      <- masummatreduced.ordered.df[ order(as.double(as.character(masummatreduced.ordered.df[ , "reslogMLS.p_value.Align" ]))) , ]
write.table(masummatreduced.ordered.df, file="OGsignif.orderedOn.align_rank.reslogMLS.tsv", sep = "\t")

masummatreduced.ordered.ma <- masummatreduced.ordered.df[ !is.na(masummatreduced.ordered.df[ , "MLS.rank.Ma"]) | !is.na(masummatreduced.ordered.df[ , "resMLS.rank.Ma"]) , ]
write.table(masummatreduced.ordered.ma, file="OGsignif_sharedWithMa.ordered.orderedOn.align_rank.MLS.tsv", sep = "\t")

png(filename = "ranks.reslogMLS.p_value.Align.20190130.png", height = 1000,width = 1000)
masummatreduced.ordered.ma                      <- masummatreduced.ordered.ma[ order(as.double(as.character(masummatreduced.ordered.ma[ , "reslogMLS.p_value.Align" ]))) , ]
heatmap.2(log(as.matrix(masummatreduced.ordered.ma[, grepl("rank",colnames(masummatreduced.ordered.ma))])),
          margins = c(15,10),
          keysize=0.7, key.par = list(cex=0.5),
          cexCol = 0.8, cexRow=0.1,
          Rowv = NA, Colv=NA,
          tracecol = "black",
          scale = "column"
          )
dev.off()

png(filename = "ranks.resMLS.p_value.Ma.png", height = 1000,width = 1000)
masummatreduced.ordered.ma                      <- masummatreduced.ordered.ma[ order(as.double(as.character(masummatreduced.ordered.ma[ , "resMLS.p_value.Ma" ]))) , ]
heatmap.2(log(as.matrix(masummatreduced.ordered.ma[, grepl("rank",colnames(masummatreduced.ordered.ma))])),
          margins = c(15,10),
          keysize=0.7, key.par = list(cex=0.5),
          cexCol = 0.8, cexRow=0.1,
          Rowv = NA, Colv=NA,
          tracecol = "black",
          scale = "column"
)
dev.off()

png(filename = "ranks.MLS.p_value.Ma.png", height = 1000,width = 1000)
masummatreduced.ordered.ma                      <- masummatreduced.ordered.ma[ order(as.double(as.character(masummatreduced.ordered.ma[ , "MLS.p_value.Ma" ]))) , ]
heatmap.2(log(as.matrix(masummatreduced.ordered.ma[, grepl("rank",colnames(masummatreduced.ordered.ma))])),
          margins = c(15,10),
          keysize=0.7, key.par = list(cex=0.5),
          cexCol = 0.8, cexRow=0.1,
          Rowv = NA, Colv=NA,
          tracecol = "black",
          scale = "column"
)
dev.off()
 

png(filename = "ranks.reslogMLS.p_value.DeNovo.png", height = 1000,width = 1000)
masummatreduced.ordered.ma                      <- masummatreduced.ordered.ma[ order(as.double(as.character(masummatreduced.ordered.ma[ , "reslogMLS.p_value.DeNovo" ]))) , ]
heatmap.2(log(as.matrix(masummatreduced.ordered.ma[, grepl("rank",colnames(masummatreduced.ordered.ma))])),
          margins = c(15,10),
          keysize=0.7, key.par = list(cex=0.5),
          cexCol = 0.8, cexRow=0.1,
          Rowv = NA, Colv=NA,
          tracecol = "black",
          scale = "column"
)
dev.off()

png(filename = "ranks.reslogMLS.p_value.DeNovo.png", height = 1000,width = 1000)
masummatreduced.ordered.ma                      <- masummatreduced.ordered.ma[ order(as.double(as.character(masummatreduced.ordered.ma[ , "reslogMLS_binary.rank.DeNovo" ]))) , ]
heatmap.2(log(as.matrix(masummatreduced.ordered.ma[, grepl("rank",colnames(masummatreduced.ordered.ma))])),
          margins = c(15,10),
          keysize=0.7, key.par = list(cex=0.5),
          cexCol = 0.8, cexRow=0.1,
          Rowv = NA, Colv=NA,
          tracecol = "black",
          scale = "column"
)
dev.off()

# save.image(file="20190205.Rdata")
