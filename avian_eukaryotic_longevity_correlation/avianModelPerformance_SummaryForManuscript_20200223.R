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

# 2/9/2020
# J. McCorrisnoyo
# Hello all,
# I've attached the newest p-value and correlation analysis matrices. A summary of these results is provided in slides 1-29 of the attached powerpoint presentation.
# As a highest priority before the U24 call (updates pending, slides 30+) I will be summarizing the count of significant genes in the most highly correlated metadata 
# contexts (O+A BMG, A AW) in the positive and negative direction. I will favor non-imputed abundances in these example slides.  When evaluating the positive and negative 
# top significant genes I will provide additional feedback on obvious human reference equivalents.

# figures_20200212.R

pvalue_report      <- read.table("model.p.20200212.tsv", sep="\t", stringsAsFactors = FALSE, fill=TRUE)

residual_report    <- read.table("resid.p.20200212.tsv", sep="\t", stringsAsFactors = FALSE, fill=TRUE)

correlation_report <- read.table("correlation_report.20200212.tsv", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
correlation_report <- correlation_report[-1,-1]
 
# https://www.orthodb.org/v9.1/download/odb9v1_OG_xrefs.tab.gz
odb9v1_OG_xrefs <- read.table("odb9v1_OG_xrefs.tab", sep="\t", stringsAsFactors = FALSE, fill=TRUE)
#odb9v1_OG_xrefs <- odb9v1_OG_xrefs[-1,-1]
 
# save.image("formatted_base_requirements.Rdata")

###############################################################################################################################################

# abundance comparisons

# import
query_metadata       <- read.table("FinalSpeciesList.20191220e.txt",       sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)[1:49,]
query_metadata.anage <- read.table("anage.data.build14.Oct.2017.aves2.txt", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)

# fix bad scientific names
# Yellow Warbler  ORIG="Dendroica petechia"   ANAGE="Setophaga petechia"
# House Finch     ORIG="Carpodacus mexicanus" ANAGE="Haemorhous mexicanus"
query_metadata.anage[ query_metadata.anage$Scientific.name == "Setophaga petechia"   , "Scientific.name" ] <- "Dendroica petechia"
query_metadata.anage[ query_metadata.anage$Scientific.name == "Haemorhous mexicanus" , "Scientific.name" ] <- "Carpodacus mexicanus"

phylotree <- read.tree(text="(Ostrich:113,((Ruby_throated_Hummingbird:79.6,(((Double_crested_Cormorant:71.5,Sandhill_Crane:71.5)AT:3,(Rock_Dove:26.3,Mourning_Dove:26.3)AU:48.2)AV:4,((Killdeer:67.9,((Woodcock:7.6,Spotted_Sandpiper:7.6)AP:55.8,(Caspian_Tern:22.5,(Ring_billed_Gull:1.58,Herring_Gull:1.58)AN:20.92)AO:40.9)AQ:4.5)AR:11.5,((Great_Horned_Owl:77.7,((Red_tailed_Hawk:24.2,Coopers_Hawk:24.2)AJ:51.7,(Downy_Woodpecker:14.3,(Red_bellied_Woodpecker:7.05,Hairy_Woodpecker:7.05)AH:7.15)AI:61.6)AK:1.8)AL:0.7,((American_Crow:37.8,Red_eyed_Vireo:37.8)AF:12.3,((Tufted_Titmouse:42.2,(Horned_Lark:34.9,(Tree_Swallow:18.8,Barn_Swallow:18.8)AB:16.1)AC:7.3)AD:3.4,((House_Sparrow:29.3,(House_Finch:26,((Northern_Cardinal:21.5,(Brown_headed_Cowbird:10.5,(Yellow_Warbler:5.25,(Yellow_throated_Warbler:4.09,Yellow_rumped_Warbler:4.09)N:1.06)O:5.25)P:11)Q:1.4,(Common_Grackle:12,(Chipping_Sparrow:9.97,(Song_Sparrow:8.46,American_Tree_Sparrow:8.46)K:1.51)L:2.03)M:10.9)R:3.1)S:3.3)T:15.8,(Cedar_waxwing:44.6,((White_breasted_Nuthatch:32.2,(Carolina_Wren:14.3,House_wren:14.3)U:17.9)V:6.6,(American_Robin:28.1,(Gray_Catbird:21.9,Starling:21.9)W:6.2)X:10.7)Y:5.8)Z:0.5)AA:0.5)AE:4.5)AG:28.3)AM:1)AS:0.1)AW:0.1)AX:18,((Pheasant:14.6,(Turkey:13.3,Ruffed_Grouse:13.3)H:1.3)I:57.8,((Canada_Goose:13.5,Mute_Swan:13.5)F:10.6,(Wood_Duck:12.3,(Northern_Shoveler:7.34,(Gadwall:5.77,(Mallard:4.24,(Green_winged_Teal:3.99,Pintail_Duck:3.99)A:1.25)B:1.53)C:1.6)D:4.93)E:11.8)G:48.3)J:25.2)AY:15.4)AZ;")
phylotree       <- drop.tip(phylotree, c("Ostrich","Yellow_throated_Warbler","Pheasant"))
phylotree       <- drop.tip(phylotree, c("American_Tree_Sparrow"))
ref_pallette <- c("darkred","red","darkorange","darkgoldenrod1","yellowgreen","darkgreen","aquamarine3","blue","darkblue","blueviolet","grey","black","orange","pink","green","yellow","lightblue","brown","violet","orange"
)

query_metadata$AnAge.BMG <- query_metadata$Body.Mass.Grams #rep(NA, dim(query_metadata)[1])
query_metadata$AnAge.MLS <- query_metadata$MLS.Years #rep(NA, dim(query_metadata)[1])
toggle <- rep(0,49)
tt <- 1
for(s in query_metadata$Scientific.Name){
  if(s %in% query_metadata.anage$Scientific.name){ # all
    if( !is.na(query_metadata.anage[query_metadata.anage$Scientific.name == s,"Body.mass.g"]) ){
      print(paste(s, "A TRUE"))
      query_metadata[query_metadata$Scientific.Name == s,"AnAge.BMG"]   <- query_metadata.anage[query_metadata.anage$Scientific.name == s,"Body.mass.g"]
      toggle[tt] <- 1
    }
    if(!is.na(query_metadata.anage[query_metadata.anage$Scientific.name == s,"Maximum.longevity.yrs"]) ){
      print(paste(s, "B TRUE"))
      query_metadata[query_metadata$Scientific.Name == s,"AnAge.MLS"]   <- query_metadata.anage[query_metadata.anage$Scientific.name == s,"Maximum.longevity.yrs"]
      if(toggle[tt] == 1){ toggle[tt] <- 3 }
      else{ toggle[tt] <- 2 }
    }else{
      print(s)
    }
  }
  tt <- tt + 1
}

query_metadata.anage_plus <- query_metadata.anage
query_metadata.anage_plus$plussed <- rep(NA,dim(query_metadata.anage_plus)[1])
for(q in which(toggle==2)){
  # update BMG only
  query_metadata.anage_plus[query_metadata.anage_plus$Scientific.name == query_metadata[q,"Scientific.Name"],"Body.mass.g"] <- query_metadata[q,"Body.Mass.Grams"]
  query_metadata.anage_plus[query_metadata.anage_plus$Scientific.name == query_metadata[q,"Scientific.Name"],"plussed"] <- 1
}
for(q in which(toggle==3)){
  # update neither
}
for(q in which(toggle==0)){
  # update both MLS and BMG
  # and add a new row
  query_metadata.anage_plus[query_metadata[q,"Scientific.Name"],"Body.mass.g"] <- query_metadata[q,"Body.Mass.Grams"]
  query_metadata.anage_plus[query_metadata[q,"Scientific.Name"],"Maximum.longevity.yrs"] <- query_metadata[q,"MLS.Years"]
  query_metadata.anage_plus[query_metadata[q,"Scientific.Name"],"Scientific.name"] <- query_metadata[q,"Scientific.Name"]
  query_metadata.anage_plus[query_metadata[q,"Scientific.Name"],"plussed"] <- 1
}

abund_ref <- read.table("Ref.OG.abundance2.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
abund_ref <- abund_ref[-dim(abund_ref)[1],]

# clean out single quotes
abund_dno <- read.table("Denovo.OG.abundance3.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
abund_dno <- abund_dno[-dim(abund_dno)[1],]


abund_ref.cts <- as.matrix(abund_ref[,query_metadata$AbundanceLabel])
rownames(abund_ref.cts) <- abund_ref$OrthoDB_OG_ID

abund_dno.cts <- as.matrix(abund_dno[,query_metadata$AbundanceLabel])
rownames(abund_dno.cts) <- abund_dno$OrthoDB_OG_ID

c <- rownames(abund_ref.cts)[ rownames(abund_ref.cts)  %in% rownames(abund_dno.cts) ]
a <- rownames(abund_ref.cts)[ !rownames(abund_ref.cts) %in% rownames(abund_dno.cts) ]
b <- rownames(abund_dno.cts)[ !rownames(abund_dno.cts)  %in% rownames(abund_ref.cts) ]
# a : the de novo set entirely contains the de novo OG set
# b : there are 10,767 unique OGs in the de novo OG set
# c : 7,983 OGs are shared and comparable betwen the sets and subject for correlation analysis

# normalize within sample
abund_ref.cts.rel <- t(decostand(t(abund_ref.cts), method="total"))
rownames(abund_ref.cts.rel) <- abund_ref$OrthoDB_OG_ID
abund_dno.cts.rel <- t(decostand(t(abund_dno.cts), method="total"))
rownames(abund_dno.cts.rel) <- abund_dno$OrthoDB_OG_ID


ovrlps <- query_metadata$Scientific.Name[ query_metadata$Scientific.Name %in% query_metadata.anage$Scientific.name ]
newL <- list()
newB <- list()
l <- 1
#for(o in ovrlps){
for(o in query_metadata$Scientific.Name){
  if(o %in% ovrlps){
    newL[l] <- query_metadata.anage[query_metadata.anage$Scientific.name == o,]$Maximum.longevity.yrs
    newB[l] <- query_metadata.anage[query_metadata.anage$Scientific.name == o,]$Body.mass.g
  }else{
    newL[l] <- NA
    newB[l] <- NA
  }
  l <- l + 1
}
newL <- unlist(newL)
newB <- unlist(newB)

newL2 <- list()
newB2 <- list()
l <- 1
for(o in query_metadata$Scientific.Name){
  newL2[l] <- query_metadata.anage_plus[query_metadata.anage_plus$Scientific.name == o,]$Maximum.longevity.yrs
  newB2[l] <- query_metadata.anage_plus[query_metadata.anage_plus$Scientific.name == o,]$Body.mass.g
  l <- l + 1
}
newL2 <- unlist(newL2)
newB2 <- unlist(newB2)

aw1 <- list()
l <- 1
for(o in query_metadata$Scientific.Name){
  aw1[l] <- query_metadata.anage_plus[query_metadata.anage_plus$Scientific.name == o,]$Adult.weight.g
  l <- l + 1
}
aw1 <- unlist(aw1)

# metadata
compare_me     <- as.data.frame(matrix(nrow=49))
# MLS
# orig
compare_me[,"MLS01"] <- as.numeric(query_metadata$MLS.Years)
# orig+anage
compare_me[,"MLS02"] <- as.numeric(query_metadata$AnAge.MLS)
# anage
#compare_me[,"MLS03"] <- newL
# anage+orig
compare_me[,"MLS04"] <- as.numeric(newL2)
# BMG
# orig
compare_me[,"BMG01"] <- as.numeric(query_metadata$Body.Mass.Grams)
# orig+anage
compare_me[,"BMG02"] <- as.numeric(query_metadata$AnAge.BMG)
# anage
#compare_me[,"BMG03"] <- newB
# anage+orig
compare_me[,"BMG04"] <- as.numeric(newB2)
# anage adult weight
compare_me[,"AW"] <- as.numeric(aw1)
# BMG
# orig
compare_me[,"logBMG01"] <- as.numeric(log(query_metadata$Body.Mass.Grams))
# orig+anage
compare_me[,"logBMG02"] <- as.numeric(log(query_metadata$AnAge.BMG))
# anage
#compare_me[,"BMG03"] <- newB
# anage+orig
compare_me[,"logBMG04"] <- as.numeric(log(newB2))
# anage adult weight
compare_me[,"logAW"] <- as.numeric(log(aw1))
# MLSW - leaf residuals
# orig
fit.mlsw01 <- lm(compare_me[,"BMG01"] ~ compare_me[,"MLS01"])
compare_me[,"MLSW01"] <- as.numeric( resid( fit.mlsw01 ) )
# orig+anage
fit.mlsw02 <- lm(compare_me[,"BMG02"] ~ compare_me[,"MLS02"])
compare_me[,"MLSW02"] <- as.numeric( resid( fit.mlsw02 ) )
# anage
#fit.mlsw03 <- lm(compare_me[,"BMG03"] ~ compare_me[,"MLS03"])
#compare_me[names(resid( fit.mlsw03 )),"MLSW03"] <- resid( fit.mlsw03 )
# anage+orig
fit.mlsw04 <- lm(compare_me[,"BMG04"] ~ compare_me[,"MLS04"])
compare_me[,"MLSW04"] <- as.numeric( resid( fit.mlsw04 ) )
# anage adult weight
fit.mlsw05 <- lm(compare_me[,"AW"] ~ compare_me[,"MLS04"])
compare_me[,"MLSW05"] <- as.numeric( resid( fit.mlsw05 ) )
# MLSLW
# orig
fit.mlslw01 <- lm(log(compare_me[,"BMG01"]) ~ compare_me[,"MLS01"])
compare_me[,"MLSLW01"] <- as.numeric( resid( fit.mlslw01 ) )
# orig+anage
fit.mlslw02 <- lm(log(compare_me[,"BMG02"]) ~ compare_me[,"MLS02"])
compare_me[,"MLSLW02"] <- as.numeric( resid( fit.mlslw02 ) )
# anage
#if we include this column in the compdata we lose tips in our comparison comp.data
#fit.mlslw03 <- lm(log10(compare_me[,"BMG03"]) ~ compare_me[,"MLS03"])
#compare_me[names(resid( fit.mlsw03 )),"MLSW03"] <- resid( fit.mlslw03 )
# anage+orig
fit.mlslw04 <- lm(log(compare_me[,"BMG04"]) ~ compare_me[,"MLS04"])
compare_me[,"MLSLW04"] <- as.numeric( resid( fit.mlslw04 ) )
# anage adult weight
fit.mlsw05 <- lm(compare_me[,"logAW"] ~ compare_me[,"MLS04"])
compare_me[,"MLSLW05"] <- as.numeric( resid( fit.mlsw05 ) )

compare_me <- compare_me[,-1]
compare_me$TreeFormat <- query_metadata$PhylotreeLabel

# save.image("formatted_base_requirements2.Rdata")

###############################################################################################################################################
d
model.p.dos <- pvalue_report
rsq <- function (x, y) cor(x, y) ^ 2
 
t_correlation <- as.data.frame(matrix())

mmmmm <- model.p.dos[ !is.na(model.p.dos$`X11.f_05.CAPER.p_val.FDR`), ]
# align
#[1] "EOG090F03QS" "EOG090F0C08" "EOG090F02SE"

png( "examples.strictly_corr2.20200227.png" , width=400, height=1200 )

plot.new()
par(mfrow=c(3,1))
tt <- 1
#for(t in rownames( mmmmm[ mmmmm$`X11.f_05.CAPER.p_val.FDR` < 0.005 , ] )){
for(t in rownames( mmmmm[ mmmmm$`X23.f_05.CAPER.p_val.FDR` < 0.005 , ] )){
  fit <- lm( log( abund_dno.cts.rel[ t , ] ) ~ resid( fit.mlsw05 ) )
  xmin <- min(resid( fit.mlsw05 ))
  xmax <- max(resid( fit.mlsw05 ))
  ymin <- min(c(  log( abund_dno.cts.rel[ t , ] )  ,  log( abund_ref.cts.rel[ t , ] )  ))
  ymax <- max(c(  log( abund_dno.cts.rel[ t , ] )  ,  log( abund_ref.cts.rel[ t , ] )  ))
  plot( resid( fit.mlsw05 ),
        log( abund_dno.cts.rel[ t , ] ),
        pch=1,
        main=t,
        cex=2,
        xlab="logMLS Residual (Anage AW)",
        ylab="log(relative abundance)",
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
  fit <- lm( log( abund_ref.cts.rel[ t , ] ) ~ resid( fit.mlsw05 ) )
  plot( resid( fit.mlsw05 ),
        log( abund_ref.cts.rel[ t , ] ),
        pch=16,
        main=t,
        cex=2,
        xlab="logMLS Residual (Anage AW)",
        ylab="log(relative abundance)",
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
  
  dontomit <- abund_dno.cts.rel[t,] != 0
  dontomit[ abund_ref.cts.rel[t,] == 0 ] <- FALSE
  r2 <- rsq(log(abund_dno.cts.rel[t,])[dontomit], log(abund_ref.cts.rel[t,])[dontomit])
  fstat    <- summary(fit)$fstatistic
  pval    <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
  #if(tt == 2){ legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd(abs(fit$residuals)),3)),paste0("i: ",round(fit$coefficients[1],3)), paste0("r2: ", round(r2,3)), paste0("log(p): ", round(log(pval),3)))) }
  #else{legend("topright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd(abs(fit$residuals)),3)),paste0("i: ",round(fit$coefficients[1],3)), paste0("r2: ", round(r2,3)), paste0("log(p): ", round(log(pval),3)))) }
  t_correlation[ t ,  "RelAbs.slope"  ]     <- fit$coefficients[2]
  t_correlation[ t ,  "RelAbs.intercept"  ] <- fit$coefficients[3]
  t_correlation[ t ,  "RelAbs.sd"  ]        <- sd(abs(fit$residuals))
  t_correlation[ t ,  "RelAbs.r2"  ]        <- r2
  t_correlation[ t ,  "RelAbs.pval"  ]      <- pval
  tt <- tt + 1
}

dev.off()
 
#t <- "EOG090F03QS"

corr_mat <- matrix(nrow=2,ncol=4)

corr_mat[1,1] <- correlation_report[t,"ccc11.01.lm_est.rho"]
corr_mat[2,1] <- correlation_report[t,"ccc11.01.caper_est.rho"]

corr_mat[1,2] <- correlation_report[t,"ccc11.02.lm_est.rho"]
corr_mat[2,2] <- correlation_report[t,"ccc11.02.caper_est.rho"]

corr_mat[1,3] <- correlation_report[t,"ccc11.04.lm_est.rho"]
corr_mat[2,3] <- correlation_report[t,"ccc11.04.caper_est.rho"]

corr_mat[1,4] <- correlation_report[t,"ccc11.05.lm_est.rho"]
corr_mat[2,4] <- correlation_report[t,"ccc11.05.caper_est.rho"]

rep_mat <- matrix(nrow=8,ncol=3)

rep_mat[1,1] <- pvalue_report[t,"X11.f_02.LM.slope"]
rep_mat[1,2] <- pvalue_report[t,"X11.f_02.LM.p_val"]
rep_mat[1,3] <- pvalue_report[t,"X11.f_02.LM.p_val.FDR"]

rep_mat[2,1] <- pvalue_report[t,"X23.f_02.LM.slope"]
rep_mat[2,2] <- pvalue_report[t,"X23.f_02.LM.p_val"]
rep_mat[2,3] <- pvalue_report[t,"X23.f_02.LM.p_val.FDR"]

rep_mat[3,1] <- pvalue_report[t,"X11.f_02.CAPER.slope"]
rep_mat[3,2] <- pvalue_report[t,"X11.f_02.CAPER.p_val"]
rep_mat[3,3] <- pvalue_report[t,"X11.f_02.CAPER.p_val.FDR"]

rep_mat[4,1] <- pvalue_report[t,"X23.f_02.CAPER.slope"]
rep_mat[4,2] <- pvalue_report[t,"X23.f_02.CAPER.p_val"]
rep_mat[4,3] <- pvalue_report[t,"X23.f_02.CAPER.p_val.FDR"]
 
rep_mat[5,2] <- pvalue_report[t,"X11.f_02.MDMR_mash.pv"]
rep_mat[5,3] <- pvalue_report[t,"X11.f_02.MDMR_mash.pv.FDR"]

rep_mat[6,2] <- pvalue_report[t,"X23.f_02.MDMR_mash.pv"]
rep_mat[6,3] <- pvalue_report[t,"X23.f_02.MDMR_mash.pv.FDR"]

rep_mat[7,2] <- pvalue_report[t,"X11.f_02.MDMR_phylo.pv"]
rep_mat[7,3] <- pvalue_report[t,"X11.f_02.MDMR_phylo.pv.FDR"]

rep_mat[8,2] <- pvalue_report[t,"X23.f_02.MDMR_phylo.pv"]
rep_mat[8,3] <- pvalue_report[t,"X23.f_02.MDMR_phylo.pv.FDR"]

View( t(rep_mat) )
 
rep_mat <- matrix(nrow=8,ncol=3)

rep_mat[1,1] <- pvalue_report[t,"X11.f_05.LM.slope"]
rep_mat[1,2] <- pvalue_report[t,"X11.f_05.LM.p_val"]
rep_mat[1,3] <- pvalue_report[t,"X11.f_05.LM.p_val.FDR"]

rep_mat[2,1] <- pvalue_report[t,"X23.f_05.LM.slope"]
rep_mat[2,2] <- pvalue_report[t,"X23.f_05.LM.p_val"]
rep_mat[2,3] <- pvalue_report[t,"X23.f_05.LM.p_val.FDR"]

rep_mat[3,1] <- pvalue_report[t,"X11.f_05.CAPER.slope"]
rep_mat[3,2] <- pvalue_report[t,"X11.f_05.CAPER.p_val"]
rep_mat[3,3] <- pvalue_report[t,"X11.f_05.CAPER.p_val.FDR"]

rep_mat[4,1] <- pvalue_report[t,"X23.f_05.CAPER.slope"]
rep_mat[4,2] <- pvalue_report[t,"X23.f_05.CAPER.p_val"]
rep_mat[4,3] <- pvalue_report[t,"X23.f_05.CAPER.p_val.FDR"]

rep_mat[5,2] <- pvalue_report[t,"X11.f_05.MDMR_mash.pv"]
rep_mat[5,3] <- pvalue_report[t,"X11.f_05.MDMR_mash.pv.FDR"]

rep_mat[6,2] <- pvalue_report[t,"X23.f_05.MDMR_mash.pv"]
rep_mat[6,3] <- pvalue_report[t,"X23.f_05.MDMR_mash.pv.FDR"]

rep_mat[7,2] <- pvalue_report[t,"X11.f_05.MDMR_phylo.pv"]
rep_mat[7,3] <- pvalue_report[t,"X11.f_05.MDMR_phylo.pv.FDR"]

rep_mat[8,2] <- pvalue_report[t,"X23.f_05.MDMR_phylo.pv"]
rep_mat[8,3] <- pvalue_report[t,"X23.f_05.MDMR_phylo.pv.FDR"]

View( t(rep_mat) )

odb9v1_OG_xrefs[odb9v1_OG_xrefs[,1]==t,]





mmmmm <- model.p.dos[ !is.na(model.p.dos$`X23.f_05.CAPER.p_val.FDR`), ]
# denovo
#    "EOG090F055R" "EOG090F0AQK" "EOG090F02SE"

plot.new()
par(mfrow=c(3,1))
for(t in rev(rownames( mmmmm[ mmmmm$`23.f_05.CAPER.p_val.FDR` < 0.005 , ] ))){
  fitttt <- lm( log( abund_dno.cts.rel[ t , ] ) ~ resid( fit.mlsw05 ) )
  plot( resid( fit.mlsw05 ),
        log( abund_dno.cts.rel[ t , ] ),
        pch=16,
        main=t
  )
  abline(fitttt$coefficients[1],fitttt$coefficients[2])
}
