
setwd("C:\\Users\\jamis\\Desktop\\birdAgain")
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

m1 <- read.table( "C:\\Users\\jamis\\Desktop\\AvianPostPHD\\model_performance.2020610b.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )
m1 <- m1[-1,-1]

m2  <- read.table("C:\\Users\\jamis\\Desktop\\old_desktop_20200521\\model.p.20200212.TSV",       sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)

colnames(m1)

colnames(m2)

#######################################################################################
#
# Table 1 :  Minimal Table Requests for Oct 4 copy
# 
# Model:          Compare LM and Caper                (Omit MDMR)
#                 Show "model 2" : Mass ~ MaxLifesPan (Omit "model 1" : Lifespan ~ Mass )
# Abundance:      Numeric only                        (Omit Rank-order)
#                 Compare MLS, MLSW, MLSLW
#                 Compare De Novo and Alignment       (omit shared)
# Directionality: Merge pos + neg
# Metadata:       "Validated" AnAge                   (Omit curated metadata)
# Counts (x3):   1) Signif < 0.05
#                2) Signif < 0.005
#                3) FDR Signif < 0.05
#
#######################################################################################


# signif_report_OGnamesBelowCutoff_directional.<rowname>.<cutoff>
# signif_report_GeneNamesBelowCutoff_directional.<rowname>.<cutoff>
# signif_report_TopGobelowCutoff_directional.<cutoff>.tsv

ma.mls  = read.table("C:\\Users\\jamis\\Downloads\\Ma.elife-19130-table1-data1-v2.B.MLS.txt", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )
ma.mls.cutoff.signif = ma.mls$Symbol[ ma.mls$p.value.all < 0.05 ]

ma.mlsw = read.table("C:\\Users\\jamis\\Downloads\\Ma.elife-19130-table1-data1-v2.D.resMLS.txt", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )
ma.mlsw.cutoff.signif = ma.mlsw$Symbol[ ma.mlsw$p.value.all < 0.05 ]

ma.cutoff.signif = toupper(unique(c(ma.mls.cutoff.signif,ma.mlsw.cutoff.signif)))
 
table1_reportMatrix = matrix(nrow=9,ncol=8)
table1_reportMatrix = matrix(nrow=999,ncol=8)
rownames(table1_reportMatrix)=c("og_list_0.05_raw","gene_list_0.05_raw","gene_ma_overlap_list_0.05_raw",
                                "og_list_0.005_raw","gene_list_0.005_raw","gene_ma_overlap_list_0.005_raw",
                                "og_list_0.05_fdr","gene_list_0.05_fdr","gene_ma_overlap_list_0.05_fdr")

colnames(table1_reportMatrix) = 
  c("9.AW_A~abund_ref.cts.rel.LM.slope", 
    "9.AW_A~abund_ref.cts.rel.CAPER.slope",
    "12.lnAW_A~abund_ref.cts.rel.LM.slope",
    "12.lnAW_A~abund_ref.cts.rel.CAPER.slope",
    "33.AW_A~abund_dno.cts.rel.LM.slope",
    "33.AW_A~abund_dno.cts.rel.CAPER.slope",
    "36.lnAW_A~abund_dno.cts.rel.LM.slope",
    "36.lnAW_A~abund_dno.cts.rel.CAPER.slope"
  )
#table1_reportMatrix[1,1] = sum( m1[,colz[[1]]] < 0.05 )

for(modelname in 
    c("9.AW_A~abund_ref.cts.rel.LM.slope", 
      "9.AW_A~abund_ref.cts.rel.CAPER.slope",
      "12.lnAW_A~abund_ref.cts.rel.LM.slope",
      "12.lnAW_A~abund_ref.cts.rel.CAPER.slope",
      "33.AW_A~abund_dno.cts.rel.LM.slope",
      "33.AW_A~abund_dno.cts.rel.CAPER.slope",
      "36.lnAW_A~abund_dno.cts.rel.LM.slope",
      "36.lnAW_A~abund_dno.cts.rel.CAPER.slope"
    )){

    cutoff = 5e-05
    tmp.og.names   = read.table(paste0("C:\\Users\\jamis\\Desktop\\AvianPostPHD\\signif_report\\signif_report_OGnamesBelowCutoff_directional.",modelname,".",cutoff,".tsv"), sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )
    tmp.gene.names = read.table(paste0("C:\\Users\\jamis\\Desktop\\AvianPostPHD\\signif_report\\signif_report_GeneNamesBelowCutoff_directional.",modelname,".",cutoff,".tsv"), sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )
    
    nu_tmp.og.names = unlist(strsplit( tmp.og.names[2,], " " ))[2:33]
    names(nu_tmp.og.names) = c("raw_pos_0.01",	"raw_neg_0.01",	
    "fdr_pos_0.01",	"fdr_neg_0.01",	"raw_pos_0.05",	"raw_neg_0.05",	
    "fdr_pos_0.05",	"fdr_neg_0.05",	"raw_pos_0.001",	"raw_neg_0.001",	
    "fdr_pos_0.001",	"fdr_neg_0.001",	"raw_pos_0.005",	"raw_neg_0.005",	
    "fdr_pos_0.005",	"fdr_neg_0.005",	"raw_pos_1e-04",	"raw_neg_1e-04",	
    "fdr_pos_1e-04",	"fdr_neg_1e-04",	"raw_pos_5e-04",	"raw_neg_5e-04",	
    "fdr_pos_5e-04",	"fdr_neg_5e-04",	"raw_pos_1e-05",	"raw_neg_1e-05",	
    "fdr_pos_1e-05",	"fdr_neg_1e-05",	"raw_pos_5e-05",	"raw_neg_5e-05",
    "fdr_pos_5e-05",	"fdr_neg_5e-05")
    
    og_list_0.05_raw = strsplit( paste(unname( nu_tmp.og.names[  c("raw_neg_0.05", "raw_pos0.05"  ) ] ), collapse=",") , "," )[[1]]
    og_list_0.05_raw = toupper(og_list_0.05_raw[og_list_0.05_raw != "NULL"])
    og_ma_overlap_list_0.05_raw = og_list_0.05_raw[og_list_0.05_raw %in% ma.cutoff.signif ]
    
    og_list_0.005_raw = strsplit( paste(unname( nu_tmp.og.names[ c("raw_neg_0.005","raw_pos_0.005") ] ), collapse=",") , "," )[[1]]
    og_list_0.005_raw = toupper(og_list_0.005_raw[og_list_0.005_raw != "NULL"])
    og_ma_overlap_list_0.005_raw = og_list_0.005_raw[og_list_0.005_raw %in% ma.cutoff.signif ]
    
    og_list_0.05_fdr = strsplit( paste(unname( nu_tmp.og.names[ c("fdr_neg_0.05","fdr_pos_0.05") ] ), collapse=",") , "," )[[1]]
    og_list_0.05_fdr = toupper(og_list_0.05_fdr[og_list_0.05_fdr != "NULL"])
    
    nu_tmp.gene.names = unlist(strsplit( tmp.gene.names[2,], " " ))[2:33]
    names(nu_tmp.gene.names) = c("raw_pos_0.01",	"raw_neg_0.01",	
                               "fdr_pos_0.01",	"fdr_neg_0.01",	"raw_pos_0.05",	"raw_neg_0.05",	
                               "fdr_pos_0.05",	"fdr_neg_0.05",	"raw_pos_0.001",	"raw_neg_0.001",	
                               "fdr_pos_0.001",	"fdr_neg_0.001",	"raw_pos_0.005",	"raw_neg_0.005",	
                               "fdr_pos_0.005",	"fdr_neg_0.005",	"raw_pos_1e-04",	"raw_neg_1e-04",	
                               "fdr_pos_1e-04",	"fdr_neg_1e-04",	"raw_pos_5e-04",	"raw_neg_5e-04",	
                               "fdr_pos_5e-04",	"fdr_neg_5e-04",	"raw_pos_1e-05",	"raw_neg_1e-05",	
                               "fdr_pos_1e-05",	"fdr_neg_1e-05",	"raw_pos_5e-05",	"raw_neg_5e-05",
                               "fdr_pos_5e-05",	"fdr_neg_5e-05")
    gene_list_0.05_raw = strsplit( paste(unname( nu_tmp.gene.names[  c("raw_neg_0.05", "raw_pos0.05"  ) ] ), collapse=",") , "," )[[1]]
    gene_list_0.05_raw = toupper(gene_list_0.05_raw[gene_list_0.05_raw != "NULL"])
    gene_ma_overlap_list_0.05_raw = gene_list_0.05_raw[gene_list_0.05_raw %in% ma.cutoff.signif ]
    
    gene_list_0.005_raw = strsplit( paste(unname( nu_tmp.gene.names[ c("raw_neg_0.005","raw_pos_0.005") ] ), collapse=",") , "," )[[1]]
    gene_list_0.005_raw = toupper(gene_list_0.005_raw[gene_list_0.005_raw != "NULL"])
    gene_ma_overlap_list_0.005_raw = gene_list_0.005_raw[gene_list_0.005_raw %in% ma.cutoff.signif ]
    
    gene_list_0.05_fdr = strsplit( paste(unname( nu_tmp.gene.names[ c("fdr_neg_0.05","fdr_pos_0.05") ] ), collapse=",") , "," )[[1]]
    gene_list_0.05_fdr = toupper(gene_list_0.05_fdr[gene_list_0.05_fdr != "NULL"])
    gene_ma_overlap_list_0.05_fdr = gene_list_0.05_fdr[gene_list_0.05_fdr %in% ma.cutoff.signif ]
    
    table1_reportMatrix[  , modelname     ] = c( length(og_list_0.05_raw) , length(gene_list_0.05_raw), length(gene_ma_overlap_list_0.05_raw), 
                                                length(og_list_0.005_raw) , length(gene_list_0.005_raw), length(gene_ma_overlap_list_0.005_raw), 
                                                 length(og_list_0.05_fdr) , length(gene_list_0.05_fdr), length(gene_ma_overlap_list_0.05_fdr)
                                               )
    
    
}

#table2_reportMatrix[ , modelname      ] = c()

# append og_correlations

sum( og_correlation$LOGImputedRelAbs.pval[!is.na(og_correlation$LOGImputedRelAbs.pval)] < 0.05 )
#[1] 13284
sum( og_correlation$LOGImputedRelAbs.pval[!is.na(og_correlation$LOGImputedRelAbs.pval)] < 0.005 )
#[1] 12164
sum( p.adjust(og_correlation$LOGImputedRelAbs.pval[!is.na(og_correlation$LOGImputedRelAbs.pval)],method="fdr") < 0.05 )
# 

#############################################################################################################################################################

# 
# colz=list()
# colz.fdr=list()
# #################################################
# #
# # m1: select 12 models for testing: 
# # 
# # row 2: Phenotype = MLSW
# # 
# #  set 1: reference alignment
# #   col 1: phylo_contrast
# #   col 2: LM
# colz[1] = colnames(m1)[grepl("AW_A.abund_ref.cts.rel.CAPER.p_val",colnames(m1))][1]
# colz[2] = colnames(m1)[grepl("AW_A.abund_ref.cts.rel.LM.p_val"   ,colnames(m1))][1]
# colz.fdr[1] = colnames(m1)[grepl("AW_A.abund_ref.cts.rel.CAPER.p_val.FDR",colnames(m1))][1]
# colz.fdr[2] = colnames(m1)[grepl("AW_A.abund_ref.cts.rel.LM.p_val.FDR"   ,colnames(m1))][1]
# #  set 2: de novo
# #   col 1: phylo_contrast
# #   col 2: LM
# colz[3] = colnames(m1)[grepl("AW_A.abund_dno.cts.rel.CAPER.p_val",colnames(m1))][1]
# colz[4] = colnames(m1)[grepl("AW_A.abund_dno.cts.rel.LM.p_val"   ,colnames(m1))][1]
# colz.fdr[3] = colnames(m1)[grepl("AW_A.abund_dno.cts.rel.CAPER.p_val.FDR",colnames(m1))][1]
# colz.fdr[4] = colnames(m1)[grepl("AW_A.abund_dno.cts.rel.LM.p_val.FDR"   ,colnames(m1))][1]
# # 
# #################################################
# #
# # row 3: Phenotype = MLSLW
# # 
# #  set 1: reference alignment
# #   col 1: phylo_contrast
# #   col 2: LM
# colz[5] = colnames(m1)[grepl("AW_A.abund_ref.cts.rel.CAPER.p_val",colnames(m1))][2]
# colz[6] = colnames(m1)[grepl("AW_A.abund_ref.cts.rel.LM.p_val"   ,colnames(m1))][2]
# colz.fdr[5] = colnames(m1)[grepl("AW_A.abund_ref.cts.rel.CAPER.p_val.FDR",colnames(m1))][2]
# colz.fdr[6] = colnames(m1)[grepl("AW_A.abund_ref.cts.rel.LM.p_val.FDR"   ,colnames(m1))][2]
# #  set 2: de novo
# #   col 1: phylo_contrast
# #   col 2: LM
# colz[7] = colnames(m1)[grepl("AW_A.abund_dno.cts.rel.CAPER.p_val",colnames(m1))][2]
# colz[8] = colnames(m1)[grepl("AW_A.abund_dno.cts.rel.LM.p_val"   ,colnames(m1))][2]
# colz.fdr[7] = colnames(m1)[grepl("AW_A.abund_dno.cts.rel.CAPER.p_val.FDR",colnames(m1))][2]
# colz.fdr[8] = colnames(m1)[grepl("AW_A.abund_dno.cts.rel.LM.p_val.FDR"   ,colnames(m1))][2]
# # 
# #################################################

#######################################################################################
#
# TABLE 2 #
# BASed on x02-compareSignifAndRenrTopGO*.R
#
#######################################################################################

setwd("C:\\Users\\jamis\\Desktop\\AvianPostPHD")
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

# pull in environment from x01_modelExe_20200612.R
# load("deveL-complete-202200613.Rdata")
load("named_ogs.20200612.Rdata")
# model_performance.orig = model_performance[-1,-1]

write.table(as.matrix(cmpq), file = "names_ogs.20200612.20201010.tsv", sep="\t")

report_table2 = as.data.frame(matrix())

query_list = c(
  "9.AW_A~abund_ref.cts.rel.CAPER",
  "9.AW_A~abund_ref.cts.rel.LM", 
  "12.lnAW_A~abund_ref.cts.rel.CAPER",
  "12.lnAW_A~abund_ref.cts.rel.LM",
  "33.AW_A~abund_dno.cts.rel.CAPER",
  "33.AW_A~abund_dno.cts.rel.LM",
  "36.lnAW_A~abund_dno.cts.rel.CAPER",
  "36.lnAW_A~abund_dno.cts.rel.LM"
)
for(query in query_list){

  roword.query  = rownames(model_performance)[order( model_performance[,paste0(query,'.p_val.FDR')] )]
  model_performance.query = model_performance[roword.query,]
    
  for(w in c(1:7)){

    og =  rownames(model_performance.query)[w]
    report_table2[paste0(query,"_Ref_Phy_Name"),         paste0("OG",w)] = og
    if(og %in% names(cmpq)){
      report_table2[paste0(query,"_Ref_Phy_HurefName"),    paste0("OG",w)] = unname(cmpq[og])[[1]]
    }
    report_table2[paste0(query,"_Ref_Phy_SlopeOfAssoc"), paste0("OG",w)] = model_performance.query[ , paste0(query,'.','slope') ][w]
    report_table2[paste0(query,"Ref_Phy_p-val"),        paste0("OG",w)]  = model_performance.query[ , paste0(query,'.','p_val') ][w]
    report_table2[paste0(query,"Ref_Phy_FDR"),          paste0("OG",w)]  = model_performance.query[ , paste0(query,'.','p_val.FDR') ][w]
    report_table2[paste0(query,"_Coeff"), paste0("OG",w)]                 = og_correlation[og,'CCC.estRho']

  }
  
}
report_table2 = report_table2[-1,-1]
write.table(report_table2, file="table2_fdrSort.tsv", sep="\t")
 
 



report_table2 = as.data.frame(matrix())

query_list = c(
  "9.AW_A~abund_ref.cts.rel.CAPER",
  "9.AW_A~abund_ref.cts.rel.LM", 
  "12.lnAW_A~abund_ref.cts.rel.CAPER",
  "12.lnAW_A~abund_ref.cts.rel.LM",
  "33.AW_A~abund_dno.cts.rel.CAPER",
  "33.AW_A~abund_dno.cts.rel.LM",
  "36.lnAW_A~abund_dno.cts.rel.CAPER",
  "36.lnAW_A~abund_dno.cts.rel.LM"
)
for(query in query_list){
  
  roword.query  = rownames(model_performance)[order( model_performance[,paste0(query,'.p_val')] )]
  model_performance.query = model_performance[roword.query,]
  
  for(w in c(1:7)){
    
    og =  rownames(model_performance.query)[w]
    report_table2[paste0(query,"_Ref_Phy_Name"),         paste0("OG",w)] = og
    if(og %in% names(cmpq)){
      report_table2[paste0(query,"_Ref_Phy_HurefName"),    paste0("OG",w)] = unname(cmpq[og])[[1]]
    }
    report_table2[paste0(query,"_Ref_Phy_SlopeOfAssoc"), paste0("OG",w)] = model_performance.query[ , paste0(query,'.','slope') ][w]
    report_table2[paste0(query,"Ref_Phy_p-val"),        paste0("OG",w)] = model_performance.query[ , paste0(query,'.','p_val') ][w]
    report_table2[paste0(query,"Ref_Phy_FDR"),          paste0("OG",w)] = model_performance.query[ , paste0(query,'.','p_val.FDR') ][w]
    report_table2[paste0(query,"_Coeff"), paste0("OG",w)]                 = og_correlation[og,'CCC.estRho']

  }
  
}
report_table2 = report_table2[-1,-1]
write.table(report_table2, file="table2_pSort.tsv", sep="\t")

#######################################################################################
#
# table 3 # 
# 
#######################################################################################

bp.pos <- read.table( "C:\\Users\\jamis\\Desktop\\AvianPostPHD\\TopGO_BP.raw_pos_12.lnAW_A~abund_ref.cts.rel.CAPER.slope.0.05.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )
bp.neg <- read.table( "C:\\Users\\jamis\\Desktop\\AvianPostPHD\\TopGO_BP.raw_neg_12.lnAW_A~abund_ref.cts.rel.CAPER.slope.0.05.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )

mf.pos <- read.table( "C:\\Users\\jamis\\Desktop\\AvianPostPHD\\TopGO_MF.raw_pos_12.lnAW_A~abund_ref.cts.rel.CAPER.slope.0.05.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )
mf.neg <- read.table( "C:\\Users\\jamis\\Desktop\\AvianPostPHD\\TopGO_MF.raw_neg_12.lnAW_A~abund_ref.cts.rel.CAPER.slope.0.05.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )

cc.pos <- read.table( "C:\\Users\\jamis\\Desktop\\AvianPostPHD\\TopGO_CC.raw_pos_12.lnAW_A~abund_ref.cts.rel.CAPER.slope.0.05.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )
cc.neg <- read.table( "C:\\Users\\jamis\\Desktop\\AvianPostPHD\\TopGO_CC.raw_neg_12.lnAW_A~abund_ref.cts.rel.CAPER.slope.0.05.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE )
 
bp.cut = sort(c(bp.pos$classicKS, bp.neg$classicKS))[5]
bp.pos.top = bp.pos[ bp.pos$classicKS <= bp.cut , ]
bp.neg.top = bp.neg[ bp.pos$classicKS <= bp.cut , ]
bp.pos.top[,'dir']=rep("pos",nrow(bp.pos.top))
bp.neg.top[,'dir']=rep("neg",nrow(bp.neg.top))
bp.top = rbind(bp.pos.top,bp.neg.top)
bp.top = bp.top[,c("Term","classicKS","dir")]
bp.top = bp.top[order(bp.top[,'classicKS']),]

mf.cut = sort(c(mf.pos$classicKS, mf.neg$classicKS))[5]
mf.pos.top = mf.pos[ mf.pos$classicKS <= mf.cut , ]
mf.neg.top = mf.neg[ mf.neg$classicKS <= mf.cut , ]
mf.pos.top[,'dir']=rep("pos",nrow(mf.pos.top))
mf.neg.top[,'dir']=rep("neg",nrow(mf.neg.top))
mf.top = rbind(mf.pos.top,mf.neg.top)
mf.top = mf.top[,c("Term","classicKS","dir")]
mf.top = mf.top[order(mf.top[,'classicKS']),]

cc.cut = sort(c(cc.pos$classicKS, cc.neg$classicKS))[5]
cc.pos.top = cc.pos[ cc.pos$classicKS <= cc.cut , ]
cc.neg.top = cc.neg[ cc.neg$classicKS <= cc.cut , ]
cc.pos.top[,'dir']=rep("pos",nrow(cc.pos.top))
cc.neg.top[,'dir']=rep("neg",nrow(cc.neg.top))
cc.top = rbind(cc.pos.top,cc.neg.top)
cc.top = cc.top[,c("Term","classicKS","dir")]
cc.top = cc.top[order(cc.top[,'classicKS']),]
 
View(rbind(bp.top,mf.top,cc.top))
                 
ip.pos <- read.table( "C:\\Users\\jamis\\Desktop\\positive_IPA_example_DaF.txt", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE )
ip.pos[,'dir']=rep("pos",nrow(ip.pos))
ip.neg <- read.table( "C:\\Users\\jamis\\Desktop\\negative_IPA_example_DaF.txt", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE )
ip.neg[,'dir']=rep("neg",nrow(ip.neg))
ip.top=rbind(ip.pos[3:8,],ip.neg[3:8,])
ip.top=ip.top[order(ip.top$V3),]                 
ip.top[1:5,]

#######################################################################################
#
# Figure 4
# 
#######################################################################################
#
table2_og_ids = c("EOG090F01WE","EOG090F0233","EOG090F07NB","EOG090F06CH","EOG090F02OY","EOG090F003H","EOG090F00UZ","EOG090F022R","EOG090F08Y7","EOG090F038X","EOG090F009Y","EOG090F01F5","EOG090F089O","EOG090F001K","EOG090F04K0","EOG090F0CIG","EOG090F0AUS","EOG090F02SH","EOG090F01JL","EOG090F0233","EOG090F028Q","EOG090F009Y","EOG090F022R","EOG090F02SH","EOG090F032U","EOG090F0433","EOG090F06N2","EOG090F086R")

# png("metadataInitialValidation_CAPER.MLS.logAW_A.20200528.png", width=600, height=1000)
# f <- as.formula("MLS02 ~ logAW")
# compdata <- comparative.data(phylotree, compare_me, phylotreeLabel)
# reg   <- crunch(f, data=compdata)
# reg1 <- reg
# fit1 <- fit
# residual_report[ names(reg$contrast.data$studentResid) , 1 ] <- unname(reg$contrast.data$studentResid)
# plot(phylotree, type = "phylogram", main="ln( AW (A) )",cex = 0.5,x.lim = c(0,140)
# )
# sn <-  c( "A","B","C","D","E","F","G","H","J","L","M","O","P","Q","R","S","T","U","V","W","X","Y","Z",
#           "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY" )
# centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
# nodelabels(text=phylotree$node.label,cex=0.9,
#            bg=redblue(centered_color_range)[ round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1 ],
#            frame="circle"
# )
# axis(side = 1, at = seq(0, 140, 20), labels = FALSE, lwd = 2)
# mtext(seq(0, 140, 20), side = 1, at = seq(0, 140, 20), line = 1, las = 2)
# dev.off()
# 
# png("metadataInitialValidation_CAPER.MLS.logAW_A.LEGEND.20200528.png", width=600, height=200)
# centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
# plot(1:91,rep(1,91),col=redblue(centered_color_range),pch=15, cex = 2, ylab="", xlab="Relative Phylogeny Contrast (%)")
# dev.off()
# 
# residual_comparison.lm[,"lnAW_A"] <- resi
# residual_comparison.caper[,"lnAW_A"] <- unname(reg$contrast.data$studentResid)
# residual_comparison.dir_and_sd_status[,"lnAW_A_sd1"] <- abs(resi) > sd1
# residual_comparison.dir_and_sd_status[,"lnAW_A_sd2"] <- abs(resi) > sd2
# residual_comparison.dir_and_sd_status[,"lnAW_A_direction"] <- resi >= 0

for(t in table2_og_ids){
  
  png(filename = paste0(t,".png"), width=600, height=600)
  
  fit <- lm( log( abund_dno.cts.rel[ t , ] ) ~ compare_me[,"lnAW_A"] )
    #lm( log( abund_dno.cts.rel[ t , ] ) ~ resid( fit.mlsw05 ) )
  xmin <- min(compare_me[,"lnAW_A"] ) #min(resid( fit.mlsw05 ))
  xmax <- max(compare_me[,"lnAW_A"] ) #max(resid( fit.mlsw05 ))
  ymin <- min(c(  log( abund_dno.cts.rel[ t , ] )  ,  log( abund_ref.cts.rel[ t , ] )  ))
  ymax <- max(c(  log( abund_dno.cts.rel[ t , ] )  ,  log( abund_ref.cts.rel[ t , ] )  ))
  plot( compare_me[,"lnAW_A"], # resid( fit.mlsw05 ),
        log( abund_dno.cts.rel[ t , ] ),
        pch=1,
        main=t,
        cex=2,
        xlab=" ", # "logMLS Residual (Anage AW)",
        ylab=" ", # "log(relative abundance)",
        xlim=c(xmin,xmax),
        ylim=c(ymin,ymax)
  )
  sd2 <- sd(abs(fit$residuals))*2
  sd1 <- sd(abs(fit$residuals))
  abline(fit$coefficients[1],fit$coefficients[2], lty=2)
  abline(fit$coefficients[1],fit$coefficients[2], lty=2)
  abline(fit$coefficients[1]+sd2,fit$coefficients[2], col="gray", lty=4)
  abline(fit$coefficients[1]-sd2,fit$coefficients[2], col="gray", lty=4)
  abline(fit$coefficients[1]+sd1,fit$coefficients[2], col="lightblue", lty=4)
  abline(fit$coefficients[1]-sd1,fit$coefficients[2], col="lightblue", lty=4) 
  
  par(new=TRUE)
  fit <- lm( log( abund_ref.cts.rel[ t , ] ) ~ compare_me[,"lnAW_A"] )
  #lm( log( abund_dno.cts.rel[ t , ] ) ~ resid( fit.mlsw05 ) )
  plot( compare_me[,"lnAW_A"], # resid( fit.mlsw05 ),
        log( abund_ref.cts.rel[ t , ] ),
        pch=16,
        main=t,
        cex=2,
        xlab="log(MLS) Residual (Anage AW)",
        ylab="ln(Relative Abundance)",
        xlim=c(xmin,xmax),
        ylim=c(ymin,ymax)
  )
  sd2 <- sd(abs(fit$residuals))*2
  sd1 <- sd(abs(fit$residuals))
  abline(fit$coefficients[1],fit$coefficients[2], lty=1)
  abline(fit$coefficients[1],fit$coefficients[2], lty=1)
  abline(fit$coefficients[1]+sd2,fit$coefficients[2], col="gray", lty=4)
  abline(fit$coefficients[1]-sd2,fit$coefficients[2], col="gray", lty=4)
  abline(fit$coefficients[1]+sd1,fit$coefficients[2], col="lightblue", lty=4)
  abline(fit$coefficients[1]-sd1,fit$coefficients[2], col="lightblue", lty=4) 
  
  dev.off()
}
 


#######################################################################################
#
# Figure 1
# 
#######################################################################################

png("figure1.png", width=1000,height=1294)

plot.new()
#m <- rbind(c(1:6))
#layout(m,c(0.60,0.08,0.08,0.08,0.08,0.08))
m <- rbind(c(1:4))
layout(m,c(0.64,0.12,0.12,0.12))
# > par()$mar
# [1] 5.1 4.1 4.1 2.1
par(mar=c(0,0,0,0))
plot(phylotree,
     type = "phylogram", 
     #main="Avian Phylogeny",
     #tip.color = unlist(cl),
     cex = 1.75, #cex = 1.4, 
     cex.main=2,
     x.lim = c(0,140),
     label.offset=0.5
     
)

f <- as.formula("AW ~ MLS04")
compdata <- comparative.data(phylotree, compare_me, phylotreeLabel)
reg   <- crunch(f, data=compdata)
nodelabels(text=phylotree$node.label,cex=1.2,bg="white",frame="circle",
           col = rev(redgreen(100)[1:75])[round(as.numeric(unlist(abs(reg$mod$residuals)/max(abs(reg$mod$residuals)) ))*75)],
             #rev(redgreen(100))[round(as.numeric(unlist(abs(reg_LM_t_LOGw$mod$residuals)/max(abs(reg_LM_t_LOGw$mod$residuals)) ))*100)],
)
namez = query_metadata[ match(phylotree$tip.label, query_metadata$PhylotreeLabel), 'MASHrefsimple' ]
idz = as.numeric(as.factor( query_metadata[ match(phylotree$tip.label, query_metadata$PhylotreeLabel), 'MASHrefsimple' ] ))
refz = primary.colors(26)[3:22][ idz ]
textcol=rep("white",length(refz))
textcol[idz %in% c(6,5,7,15,14,16)] = "black"
 
#  query_metadata$MASHrefsimple
bp = barplot( rep(1,length(phylotree$tip.label)), col=refz, horiz=TRUE, axes = FALSE, main = "\n\nBest Ref", cex.main=2 )
text(rep(0.5,length(phylotree$tip.label)),
     bp,
     idz,
     col=textcol,
     cex=1.75
     )
bp = barplot( rep(1,length(phylotree$tip.label)), col="white", border="white", horiz=TRUE, axes = FALSE, main = "\n\nMLS", cex.main=2 )
text(rep(0.5,length(phylotree$tip.label)),
     bp,
     compare_me[ match(phylotree$tip.label, query_metadata$PhylotreeLabel) , "MLS04" ],
     cex=1.75
)
bp = barplot( rep(1,length(phylotree$tip.label)), col="white", border="white", horiz=TRUE, axes = FALSE, main = "\n\nln(AW)", cex.main=2 )
text(rep(0.5,length(phylotree$tip.label)),
     bp,
     round( compare_me[ match(phylotree$tip.label, query_metadata$PhylotreeLabel) , "logAW" ], digits=3 ),
     cex=1.75
)

dev.off()
 





#######################################################################################
#
png("figure2.png", width=1000, height=1000)
plot.new()
par(cex.lab=2)
par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
par(oma=c(0,0,0,0))
fit <- lm(compare_me[,"MLS04"] ~ compare_me[,"logAW"])
sd2 <- sd(abs(fit$residuals))*2
sd1 <- sd(abs(fit$residuals))
resi <- resid(fit) 
text_col = rep("black", 49)
#text_col[ abs(resi) > sd1 & abs(resi) < sd2 ] <- "blue"
#text_col[ abs(resi) > sd2 ]                   <- "black"
#text_col[ abs(resi) < sd1 ]                   <- "gray"
plot(                   
  compare_me[,"logAW"] ,
  compare_me[,"MLS04"] ,
  pch=16, #xlim=c(0,85) , ylim=c(0,12), 
  col=text_col, 
  ylab = "Maximum Lifespan (Years)",
  xlab="ln( Adult Weight (Grams) ) ",
  #main="AW (A)", 
  cex=2, cex.main=2, cex.axis=2) 
abline(fit$coefficients[1],fit$coefficients[2])
abline(fit$coefficients[1]+sd2,fit$coefficients[2], col="black", lty=4, lwd=2)
abline(fit$coefficients[1]-sd2,fit$coefficients[2], col="black", lty=4, lwd=2)
abline(fit$coefficients[1]+sd1,fit$coefficients[2], col="black", lty=3, lwd=2)
abline(fit$coefficients[1]-sd1,fit$coefficients[2], col="black", lty=3, lwd=2) 

update_lab = compare_me$TreeFormat
update_lab[ abs(resi) < sd1 ] <- ""
pos_lab = rep(4,length(update_lab))
pos_lab[resi > 0] = 2
pos_lab[22] = 2
pos_lab[46] = 1
pos_lab[47] = 3
update_lab[update_lab != ""] = query_metadata[ match( query_metadata$PhylotreeLabel, compare_me$TreeFormat ) , "Common.Name" ][update_lab != ""]
update_lab[update_lab == "European House Sparrow"] = "European\nHouse Sparrow"
update_lab[update_lab == "Double Crested Cormorant"] = "Double Crested\nCormorant\n"
update_lab[update_lab == "Cedar waxwing"] = "                        Cedar Waxwing"
update_lab[update_lab == "Horned Lark"] = "Horned Lark            "
  
text( compare_me[,"logAW"] ,
      compare_me[,"MLS04"] ,
      update_lab,
      pos=pos_lab,
      cex=1.5
    )

text( compare_me[,"logAW"] ,
      compare_me[,"MLS04"] ,
      update_lab,
      pos=pos_lab,
      cex=1.5
)
dev.off()
 
 
# Figure 2
# 
#######################################################################################













