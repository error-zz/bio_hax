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

# import
query_metadata       <- read.table("input_x01_modelExe/FinalSpeciesList.20191220e.txt",       sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)[1:49,]
query_metadata.anage <- read.table("input_x01_modelExe/anage.data.build14.Oct.2017.aves2.txt", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)

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

residual_report <- as.data.frame(matrix())

# DEBUG SIMPLE PLOTS + CAPTURE RESIDUALS + SD1 + 2 #

#########################################################################################################################################################
# METADATA VALIDATION ###################################################################################################################################
#########################################################################################################################################################

############################################
# BUILD COMP.DATA OF ALL METRICS FOR MOELS #
############################################

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
compare_me[,"MLS01"] <- as.numeric(query_metadata$MLS.Years)
compare_me[,"MLS02"] <- as.numeric(query_metadata$AnAge.MLS)
compare_me[,"MLS04"] <- as.numeric(newL2)
# BMG
compare_me[,"BMG01"] <- as.numeric(query_metadata$Body.Mass.Grams)
compare_me[,"BMG02"] <- as.numeric(query_metadata$AnAge.BMG)
compare_me[,"BMG04"] <- as.numeric(newB2)
compare_me[,"AW"] <- as.numeric(aw1)
# logBMG
compare_me[,"logBMG01"] <- as.numeric(log(query_metadata$Body.Mass.Grams))
compare_me[,"logBMG02"] <- as.numeric(log(query_metadata$AnAge.BMG))
compare_me[,"logBMG04"] <- as.numeric(log(newB2))
compare_me[,"logAW"] <- as.numeric(log(aw1))
# phylotree labels
compare_me[,"phylotreeLabel"] <- phylotree$tip.label

compare_me <- compare_me[,-1]
compare_me$TreeFormat <- query_metadata$PhylotreeLabel



residual_comparison.lm    <- as.data.frame(matrix(nrow=49,ncol=6))
residual_comparison.caper <- as.data.frame(matrix(nrow=48,ncol=6))
colnames(residual_comparison.lm) <- colnames(residual_comparison.caper) <- c("BMG_O","BMG_A","AW_A","lnBMG_O","lnBMG_A","lnAW_A")
rownames(residual_comparison.lm) <- compare_me$TreeFormat
rownames(residual_comparison.caper) <- sn <-  c( "A","B","C","D","E","F","G","H","J","L","M","O","P","Q","R","S","T","U","V","W","X","Y","Z",
                                                 "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY" )

residual_comparison.dir_and_sd_status <- as.data.frame(matrix(nrow=49,ncol=18))
rownames( residual_comparison.dir_and_sd_status ) <- compare_me$TreeFormat
colnames( residual_comparison.dir_and_sd_status ) <- c("BMG_O_sd1","BMG_O_sd2","BMG_O_direction",
                                                       "BMG_A_sd1","BMG_A_sd2","BMG_A_direction",
                                                       "AW_A_sd1","AW_A_sd2","AW_A_direction",
                                                       "lnBMG_O_sd1","lnBMG_O_sd2","lnBMG_O_direction",
                                                       "lnBMG_A_sd1","lnBMG_A_sd2","lnBMG_A_direction",
                                                       "lnAW_A_sd1","lnAW_A_sd2","lnAW_A_direction"
)

#########################################################################################################################################################
# MLS # 1 of 3 ##########################################################################################################################################
#########################################################################################################################################################

png("metadataInitialValidation_LM.MLS.BMG_O.20200528.png", width=600, height=600)
fit <- lm(compare_me[,"MLS01"] ~ compare_me[,"BMG01"])
sd2 <- sd(abs(fit$residuals))*2
sd1 <- sd(abs(fit$residuals))
resi <- resid(fit) 
text_col = rep("black", 49)
text_col[ abs(resi) > sd1 & abs(resi) < sd2 ] <- "blue"
text_col[ abs(resi) > sd2 ]                   <- "black"
text_col[ abs(resi) < sd1 ]                   <- "gray"
plot(                   
  compare_me[,"BMG01"] ,
  compare_me[,"MLS01"] ,
  pch=16, #xlim=c(0,85) , ylim=c(0,12), 
  col=text_col, ylab = "Maximum Lifespan (Years) [O]",xlab="Body Mass (Grams) [O]",
  main="BMG (O)", cex=2, cex.main=2, cex.axis=2) 
abline(fit$coefficients[1],fit$coefficients[2])
abline(fit$coefficients[1]+sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]-sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]+sd1,fit$coefficients[2], col="lightblue", lty=4)
abline(fit$coefficients[1]-sd1,fit$coefficients[2], col="lightblue", lty=4) 
legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd1,3)),paste0("i: ",round(fit$coefficients[1],3))))
#barplot(sort(resid(fit)), horiz=TRUE)
query_metadata[ resi > sd2,      "01_ORIG.SD2.neg"] <- 1
query_metadata[ resi < -sd2,      "01_ORIG.SD2.pos"] <- 1
query_metadata[ abs(resi) < sd2, "01_ORIG.SD2.all"] <- 1
query_metadata[ resi > sd1,      "01_ORIG.SD1.neg"] <- 1
query_metadata[ resi < -sd1,      "01_ORIG.SD1.pos"] <- 1
query_metadata[ abs(resi) > sd1, "01_ORIG.SD1.all"] <- 1
query_metadata[ , "01_ORIG.lm_residual"] <- resi
dev.off()

png("metadataInitialValidation_CAPER.MLS.BMG_O.20200528.png", width=600, height=1000)
f <- as.formula("MLS01 ~ BMG01")
compdata <- comparative.data(phylotree, compare_me, phylotreeLabel)
reg   <- crunch(f, data=compdata)
reg1 <- reg
fit1 <- fit
residual_report[ names(reg$contrast.data$studentResid) , 1 ] <- unname(reg$contrast.data$studentResid)
plot(phylotree, type = "phylogram", main="BMG (O)",cex = 0.5,x.lim = c(0,140)
)
sn <-  c( "A","B","C","D","E","F","G","H","J","L","M","O","P","Q","R","S","T","U","V","W","X","Y","Z",
          "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY" )
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
nodelabels(text=phylotree$node.label,cex=0.9,
           bg=redblue(centered_color_range)[ round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1 ],
           frame="circle"
)
axis(side = 1, at = seq(0, 140, 20), labels = FALSE, lwd = 2)
mtext(seq(0, 140, 20), side = 1, at = seq(0, 140, 20), line = 1, las = 2)
dev.off()

png("metadataInitialValidation_CAPER.MLS.BMG_O.LEGEND.20200528.png", width=600, height=200)
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
plot(1:91,rep(1,91),col=redblue(centered_color_range),pch=15, cex = 2, ylab="", xlab="Relative Phylogeny Contrast (%)")
dev.off()

residual_comparison.lm[,"BMG_O"] <- resi
residual_comparison.caper[,"BMG_O"] <- unname(reg$contrast.data$studentResid)
residual_comparison.dir_and_sd_status[,"BMG_O_sd1"] <- abs(resi) > sd1
residual_comparison.dir_and_sd_status[,"BMG_O_sd2"] <- abs(resi) > sd2
residual_comparison.dir_and_sd_status[,"BMG_O_direction"] <- resi >= 0

#########################################################################################################################################################
# MLS # 2 of 3 ##########################################################################################################################################
#########################################################################################################################################################

png("metadataInitialValidation_LM.MLS.BMG_A.20200528.png", width=600, height=600)
fit <- lm(compare_me[,"MLS02"] ~ compare_me[,"BMG02"])
sd2 <- sd(abs(fit$residuals))*2
sd1 <- sd(abs(fit$residuals))
resi <- resid(fit) 
text_col = rep("black", 49)
text_col[ abs(resi) > sd1 & abs(resi) < sd2 ] <- "blue"
text_col[ abs(resi) > sd2 ]                   <- "black"
text_col[ abs(resi) < sd1 ]                   <- "gray"
plot(                   
  compare_me[,"BMG02"] ,
  compare_me[,"MLS02"] ,
  pch=16, #xlim=c(0,85) , ylim=c(0,12), 
  col=text_col, ylab = "Maximum Lifespan (Years) [A]",xlab="Body Mass (Grams) [A]",
  main="BMG (A)", cex=2, cex.main=2, cex.axis=2) 
abline(fit$coefficients[1],fit$coefficients[2])
abline(fit$coefficients[1]+sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]-sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]+sd1,fit$coefficients[2], col="lightblue", lty=4)
abline(fit$coefficients[1]-sd1,fit$coefficients[2], col="lightblue", lty=4) 
legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd1,3)),paste0("i: ",round(fit$coefficients[1],3))))
#barplot(sort(resid(fit)), horiz=TRUE)
query_metadata[ resi > sd2,      "01_ORIG.SD2.neg"] <- 1
query_metadata[ resi < -sd2,      "01_ORIG.SD2.pos"] <- 1
query_metadata[ abs(resi) < sd2, "01_ORIG.SD2.all"] <- 1
query_metadata[ resi > sd1,      "01_ORIG.SD1.neg"] <- 1
query_metadata[ resi < -sd1,      "01_ORIG.SD1.pos"] <- 1
query_metadata[ abs(resi) > sd1, "01_ORIG.SD1.all"] <- 1
query_metadata[ , "01_ORIG.lm_residual"] <- resi
dev.off()

png("metadataInitialValidation_CAPER.MLS.BMG_A.20200528.png", width=600, height=1000)
f <- as.formula("MLS02 ~ BMG02")
compdata <- comparative.data(phylotree, compare_me, phylotreeLabel)
reg   <- crunch(f, data=compdata)
reg1 <- reg
fit1 <- fit
residual_report[ names(reg$contrast.data$studentResid) , 1 ] <- unname(reg$contrast.data$studentResid)
plot(phylotree, type = "phylogram", main="BMG (A)",cex = 0.5,x.lim = c(0,140)
)
sn <-  c( "A","B","C","D","E","F","G","H","J","L","M","O","P","Q","R","S","T","U","V","W","X","Y","Z",
          "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY" )
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
nodelabels(text=phylotree$node.label,cex=0.9,
           bg=redblue(centered_color_range)[ round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1 ],
           frame="circle"
)
axis(side = 1, at = seq(0, 140, 20), labels = FALSE, lwd = 2)
mtext(seq(0, 140, 20), side = 1, at = seq(0, 140, 20), line = 1, las = 2)
dev.off()

png("metadataInitialValidation_CAPER.MLS.BMG_A.LEGEND.20200528.png", width=600, height=200)
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
plot(1:91,rep(1,91),col=redblue(centered_color_range),pch=15, cex = 2, ylab="", xlab="Relative Phylogeny Contrast (%)")
dev.off()

residual_comparison.lm[,"BMG_A"] <- resi
residual_comparison.caper[,"BMG_A"] <- unname(reg$contrast.data$studentResid)
residual_comparison.dir_and_sd_status[,"BMG_A_sd1"] <- abs(resi) > sd1
residual_comparison.dir_and_sd_status[,"BMG_A_sd2"] <- abs(resi) > sd2
residual_comparison.dir_and_sd_status[,"BMG_A_direction"] <- resi >= 0

#########################################################################################################################################################
# MLS # 3 of 3 ##########################################################################################################################################
#########################################################################################################################################################

png("metadataInitialValidation_LM.MLS.AW_A.20200528.png", width=600, height=600)
fit <- lm(compare_me[,"MLS02"] ~ compare_me[,"AW"])
sd2 <- sd(abs(fit$residuals))*2
sd1 <- sd(abs(fit$residuals))
resi <- resid(fit) 
text_col = rep("black", 49)
text_col[ abs(resi) > sd1 & abs(resi) < sd2 ] <- "blue"
text_col[ abs(resi) > sd2 ]                   <- "black"
text_col[ abs(resi) < sd1 ]                   <- "gray"
plot(                   
  compare_me[,"AW"] ,
  compare_me[,"MLS02"] ,
  pch=16, #xlim=c(0,85) , ylim=c(0,12), 
  col=text_col, ylab = "Maximum Lifespan (Years) [A]",xlab="Adult Weight (Grams) [A]",
  main="AW (A)", cex=2, cex.main=2, cex.axis=2) 
abline(fit$coefficients[1],fit$coefficients[2])
abline(fit$coefficients[1]+sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]-sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]+sd1,fit$coefficients[2], col="lightblue", lty=4)
abline(fit$coefficients[1]-sd1,fit$coefficients[2], col="lightblue", lty=4) 
legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd1,3)),paste0("i: ",round(fit$coefficients[1],3))))
#barplot(sort(resid(fit)), horiz=TRUE)
query_metadata[ resi > sd2,      "01_ORIG.SD2.neg"] <- 1
query_metadata[ resi < -sd2,      "01_ORIG.SD2.pos"] <- 1
query_metadata[ abs(resi) < sd2, "01_ORIG.SD2.all"] <- 1
query_metadata[ resi > sd1,      "01_ORIG.SD1.neg"] <- 1
query_metadata[ resi < -sd1,      "01_ORIG.SD1.pos"] <- 1
query_metadata[ abs(resi) > sd1, "01_ORIG.SD1.all"] <- 1
query_metadata[ , "01_ORIG.lm_residual"] <- resi
dev.off()

png("metadataInitialValidation_CAPER.MLS.AW_A.20200528.png", width=600, height=1000)
f <- as.formula("MLS02 ~ BMG01")
compdata <- comparative.data(phylotree, compare_me, phylotreeLabel)
reg   <- crunch(f, data=compdata)
reg1 <- reg
fit1 <- fit
residual_report[ names(reg$contrast.data$studentResid) , 1 ] <- unname(reg$contrast.data$studentResid)
plot(phylotree, type = "phylogram", main="AW (A)",cex = 0.5,x.lim = c(0,140)
)
sn <-  c( "A","B","C","D","E","F","G","H","J","L","M","O","P","Q","R","S","T","U","V","W","X","Y","Z",
          "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY" )
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
nodelabels(text=phylotree$node.label,cex=0.9,
           bg=redblue(centered_color_range)[ round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1 ],
           frame="circle"
)
axis(side = 1, at = seq(0, 140, 20), labels = FALSE, lwd = 2)
mtext(seq(0, 140, 20), side = 1, at = seq(0, 140, 20), line = 1, las = 2)
dev.off()

png("metadataInitialValidation_CAPER.MLS.AW_A.LEGEND.20200528.png", width=600, height=200)
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
plot(1:91,rep(1,91),col=redblue(centered_color_range),pch=15, cex = 2, ylab="", xlab="Relative Phylogeny Contrast (%)")
dev.off()

residual_comparison.lm[,"AW_A"] <- resi
residual_comparison.caper[,"AW_A"] <- unname(reg$contrast.data$studentResid)
residual_comparison.dir_and_sd_status[,"AW_A_sd1"] <- abs(resi) > sd1
residual_comparison.dir_and_sd_status[,"AW_A_sd2"] <- abs(resi) > sd2
residual_comparison.dir_and_sd_status[,"AW_A_direction"] <- resi >= 0

#########################################################################################################################################################
# logMLS # 1 of 3 ##########################################################################################################################################
#########################################################################################################################################################

png("metadataInitialValidation_LM.MLS.logBMG_O.20200528.png", width=600, height=600)
fit <- lm(compare_me[,"MLS01"] ~ compare_me[,"logBMG01"])
sd2 <- sd(abs(fit$residuals))*2
sd1 <- sd(abs(fit$residuals))
resi <- resid(fit) 
text_col = rep("black", 49)
text_col[ abs(resi) > sd1 & abs(resi) < sd2 ] <- "blue"
text_col[ abs(resi) > sd2 ]                   <- "black"
text_col[ abs(resi) < sd1 ]                   <- "gray"
plot(                   
  compare_me[,"logBMG01"] ,
  compare_me[,"MLS01"] ,
  pch=16, #xlim=c(0,85) , ylim=c(0,12), 
  col=text_col, ylab = "Maximum Lifespan (Years) [O]",xlab="Natural Log of Body Mass (Grams) [O]",
  main="ln( BMG (O) )", cex=2, cex.main=2, cex.axis=2) 
abline(fit$coefficients[1],fit$coefficients[2])
abline(fit$coefficients[1]+sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]-sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]+sd1,fit$coefficients[2], col="lightblue", lty=4)
abline(fit$coefficients[1]-sd1,fit$coefficients[2], col="lightblue", lty=4) 
legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd1,3)),paste0("i: ",round(fit$coefficients[1],3))))
#barplot(sort(resid(fit)), horiz=TRUE)
query_metadata[ resi > sd2,      "01_ORIG.SD2.neg"] <- 1
query_metadata[ resi < -sd2,      "01_ORIG.SD2.pos"] <- 1
query_metadata[ abs(resi) < sd2, "01_ORIG.SD2.all"] <- 1
query_metadata[ resi > sd1,      "01_ORIG.SD1.neg"] <- 1
query_metadata[ resi < -sd1,      "01_ORIG.SD1.pos"] <- 1
query_metadata[ abs(resi) > sd1, "01_ORIG.SD1.all"] <- 1
query_metadata[ , "01_ORIG.lm_residual"] <- resi
dev.off()

png("metadataInitialValidation_CAPER.MLS.logBMG_O.20200528.png", width=600, height=1000)
f <- as.formula("MLS01 ~ logBMG01")
compdata <- comparative.data(phylotree, compare_me, phylotreeLabel)
reg   <- crunch(f, data=compdata)
reg1 <- reg
fit1 <- fit
residual_report[ names(reg$contrast.data$studentResid) , 1 ] <- unname(reg$contrast.data$studentResid)
plot(phylotree, type = "phylogram", main="ln( BMG (O) )",cex = 0.5,x.lim = c(0,140)
)
sn <-  c( "A","B","C","D","E","F","G","H","J","L","M","O","P","Q","R","S","T","U","V","W","X","Y","Z",
          "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY" )
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
nodelabels(text=phylotree$node.label,cex=0.9,
           bg=redblue(centered_color_range)[ round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1 ],
           frame="circle"
)
axis(side = 1, at = seq(0, 140, 20), labels = FALSE, lwd = 2)
mtext(seq(0, 140, 20), side = 1, at = seq(0, 140, 20), line = 1, las = 2)
dev.off()

png("metadataInitialValidation_CAPER.MLS.logBMG_O.LEGEND.20200528.png", width=600, height=200)
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
plot(1:91,rep(1,91),col=redblue(centered_color_range),pch=15, cex = 2, ylab="", xlab="Relative Phylogeny Contrast (%)")
dev.off()

residual_comparison.lm[,"lnBMG_O"] <- resi
residual_comparison.caper[,"lnBMG_O"] <- unname(reg$contrast.data$studentResid)
residual_comparison.dir_and_sd_status[,"lnBMG_O_sd1"] <- abs(resi) > sd1
residual_comparison.dir_and_sd_status[,"lnBMG_O_sd2"] <- abs(resi) > sd2
residual_comparison.dir_and_sd_status[,"lnBMG_O_direction"] <- resi >= 0

#########################################################################################################################################################
# logMLS # 2 of 3 ##########################################################################################################################################
#########################################################################################################################################################

png("metadataInitialValidation_LM.MLS.logBMG_A.20200528.png", width=600, height=600)
fit <- lm(compare_me[,"MLS02"] ~ compare_me[,"logBMG02"])
sd2 <- sd(abs(fit$residuals))*2
sd1 <- sd(abs(fit$residuals))
resi <- resid(fit) 
text_col = rep("black", 49)
text_col[ abs(resi) > sd1 & abs(resi) < sd2 ] <- "blue"
text_col[ abs(resi) > sd2 ]                   <- "black"
text_col[ abs(resi) < sd1 ]                   <- "gray"
plot(                   
  compare_me[,"logBMG02"] ,
  compare_me[,"MLS02"] ,
  pch=16, #xlim=c(0,85) , ylim=c(0,12), 
  col=text_col, ylab = "Maximum Lifespan (Years) [A]",xlab="Natural Log of Body Mass (Grams) [A]",
  main="ln( BMG (A) )", cex=2, cex.main=2, cex.axis=2) 
abline(fit$coefficients[1],fit$coefficients[2])
abline(fit$coefficients[1]+sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]-sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]+sd1,fit$coefficients[2], col="lightblue", lty=4)
abline(fit$coefficients[1]-sd1,fit$coefficients[2], col="lightblue", lty=4) 
legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd1,3)),paste0("i: ",round(fit$coefficients[1],3))))
#barplot(sort(resid(fit)), horiz=TRUE)
query_metadata[ resi > sd2,      "01_ORIG.SD2.neg"] <- 1
query_metadata[ resi < -sd2,      "01_ORIG.SD2.pos"] <- 1
query_metadata[ abs(resi) < sd2, "01_ORIG.SD2.all"] <- 1
query_metadata[ resi > sd1,      "01_ORIG.SD1.neg"] <- 1
query_metadata[ resi < -sd1,      "01_ORIG.SD1.pos"] <- 1
query_metadata[ abs(resi) > sd1, "01_ORIG.SD1.all"] <- 1
query_metadata[ , "01_ORIG.lm_residual"] <- resi
dev.off()

png("metadataInitialValidation_CAPER.MLS.logBMG_A.20200528.png", width=600, height=1000)
f <- as.formula("MLS02 ~ logBMG02")
compdata <- comparative.data(phylotree, compare_me, phylotreeLabel)
reg   <- crunch(f, data=compdata)
reg1 <- reg
fit1 <- fit
residual_report[ names(reg$contrast.data$studentResid) , 1 ] <- unname(reg$contrast.data$studentResid)
plot(phylotree, type = "phylogram", main="ln( BMG (A) )",cex = 0.5,x.lim = c(0,140)
)
sn <-  c( "A","B","C","D","E","F","G","H","J","L","M","O","P","Q","R","S","T","U","V","W","X","Y","Z",
          "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY" )
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
nodelabels(text=phylotree$node.label,cex=0.9,
           bg=redblue(centered_color_range)[ round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1 ],
           frame="circle"
)
axis(side = 1, at = seq(0, 140, 20), labels = FALSE, lwd = 2)
mtext(seq(0, 140, 20), side = 1, at = seq(0, 140, 20), line = 1, las = 2)
dev.off()

png("metadataInitialValidation_CAPER.MLS.logBMG_A.LEGEND.20200528.png", width=600, height=200)
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
plot(1:91,rep(1,91),col=redblue(centered_color_range),pch=15, cex = 2, ylab="", xlab="Relative Phylogeny Contrast (%)")
dev.off()

residual_comparison.lm[,"lnBMG_A"] <- resi
residual_comparison.caper[,"lnBMG_A"] <- unname(reg$contrast.data$studentResid)
residual_comparison.dir_and_sd_status[,"lnBMG_A_sd1"] <- abs(resi) > sd1
residual_comparison.dir_and_sd_status[,"lnBMG_A_sd2"] <- abs(resi) > sd2
residual_comparison.dir_and_sd_status[,"lnBMG_A_direction"] <- resi >= 0

#########################################################################################################################################################
# logMLS # 3 of 3 ##########################################################################################################################################
#########################################################################################################################################################

png("metadataInitialValidation_LM.MLS.logAW_A.20200528.png", width=600, height=600)
fit <- lm(compare_me[,"MLS02"] ~ compare_me[,"logAW"])
sd2 <- sd(abs(fit$residuals))*2
sd1 <- sd(abs(fit$residuals))
resi <- resid(fit) 
text_col = rep("black", 49)
text_col[ abs(resi) > sd1 & abs(resi) < sd2 ] <- "blue"
text_col[ abs(resi) > sd2 ]                   <- "black"
text_col[ abs(resi) < sd1 ]                   <- "gray"
plot(                   
  compare_me[,"logAW"] ,
  compare_me[,"MLS02"] ,
  pch=16, #xlim=c(0,85) , ylim=c(0,12), 
  col=text_col, ylab = "Maximum Lifespan (Years) [A]",xlab="Natural Log of Adult Weight (Grams) [A]",
  main="ln( AW (A) )", cex=2, cex.main=2, cex.axis=2) 
abline(fit$coefficients[1],fit$coefficients[2])
abline(fit$coefficients[1]+sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]-sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]+sd1,fit$coefficients[2], col="lightblue", lty=4)
abline(fit$coefficients[1]-sd1,fit$coefficients[2], col="lightblue", lty=4) 
legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd1,3)),paste0("i: ",round(fit$coefficients[1],3))))
#barplot(sort(resid(fit)), horiz=TRUE)
query_metadata[ resi > sd2,      "01_ORIG.SD2.neg"] <- 1
query_metadata[ resi < -sd2,      "01_ORIG.SD2.pos"] <- 1
query_metadata[ abs(resi) < sd2, "01_ORIG.SD2.all"] <- 1
query_metadata[ resi > sd1,      "01_ORIG.SD1.neg"] <- 1
query_metadata[ resi < -sd1,      "01_ORIG.SD1.pos"] <- 1
query_metadata[ abs(resi) > sd1, "01_ORIG.SD1.all"] <- 1
query_metadata[ , "01_ORIG.lm_residual"] <- resi
dev.off()

png("metadataInitialValidation_CAPER.MLS.logAW_A.20200528.png", width=600, height=1000)
f <- as.formula("MLS02 ~ logAW")
compdata <- comparative.data(phylotree, compare_me, phylotreeLabel)
reg   <- crunch(f, data=compdata)
reg1 <- reg
fit1 <- fit
residual_report[ names(reg$contrast.data$studentResid) , 1 ] <- unname(reg$contrast.data$studentResid)
plot(phylotree, type = "phylogram", main="ln( AW (A) )",cex = 0.5,x.lim = c(0,140)
)
sn <-  c( "A","B","C","D","E","F","G","H","J","L","M","O","P","Q","R","S","T","U","V","W","X","Y","Z",
          "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY" )
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
nodelabels(text=phylotree$node.label,cex=0.9,
           bg=redblue(centered_color_range)[ round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1 ],
           frame="circle"
)
axis(side = 1, at = seq(0, 140, 20), labels = FALSE, lwd = 2)
mtext(seq(0, 140, 20), side = 1, at = seq(0, 140, 20), line = 1, las = 2)
dev.off()

png("metadataInitialValidation_CAPER.MLS.logAW_A.LEGEND.20200528.png", width=600, height=200)
centered_color_range <- round(mean(as.numeric(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100) + abs(min(round((reg$contrast.data$studentResid[  phylotree$node.label  ] ) / max(abs( reg$contrast.data$studentResid[  phylotree$node.label  ] ))*100)))+1))*2)
plot(1:91,rep(1,91),col=redblue(centered_color_range),pch=15, cex = 2, ylab="", xlab="Relative Phylogeny Contrast (%)")
dev.off()

residual_comparison.lm[,"lnAW_A"] <- resi
residual_comparison.caper[,"lnAW_A"] <- unname(reg$contrast.data$studentResid)
residual_comparison.dir_and_sd_status[,"lnAW_A_sd1"] <- abs(resi) > sd1
residual_comparison.dir_and_sd_status[,"lnAW_A_sd2"] <- abs(resi) > sd2
residual_comparison.dir_and_sd_status[,"lnAW_A_direction"] <- resi >= 0





write.table(residual_comparison.lm, file="residual_comparison.lm.tsv", sep='\t')
write.table(residual_comparison.caper, file="residual_comparison.caper.tsv", sep='\t')
write.table(residual_comparison.dir_and_sd_status, file="residual_comparison.dir_and_sd_status.tsv", sep='\t')



#################################################################################################################################################################################################

#################################################################################################################################################################################################

#################################################################################################################################################################################################




png("phylotree_distance.20200530.png",width=1000,height=1000)
heatmap.2( as.matrix(distTips(phylotree, tips = "all")), trace='none', margins = c(10,10) )
dev.off()

# xs
mash.dit           <- read.table("input_x01_modelExe/mash.dist.tsv", sep="\t")
rownames(mash.dit) <- mash.dit[,1]
mash.dit           <- mash.dit[,-1]
colnames(mash.dit) <- rownames(mash.dit)

png("MASH_distance.20200211.png",width=1000,height=1000)
heatmap.2( as.matrix(mash.dit), trace='none', margins = c(10,10), main="MASH (all species)" )
dev.off()

qa  <- query_metadata$AbundanceLabel
qas <- list()
qc  <- 1
for(a in qa){
  qas[qc] <- substr(a, 2, nchar(a))
  qc      <- qc + 1
}
qas <- unlist(qas)

png("MASH_subset_distance.png",width=1000,height=1000)
heatmap.2( as.matrix(mash.dit)[qas,qas], trace='none', margins = c(10,10), main="MASH (all species)" )
dev.off()




# MDMR

# bmg 

cmpr.phyl <- as.matrix(distTips(phylotree, tips = "all"))[query_metadata$PhylotreeLabel,query_metadata$PhylotreeLabel]
cmpr.mash <- as.matrix(mash.dit)[qas,qas]

lm_mdmr_validation_summary <- matrix(ncol=28,nrow=3)
#rownames(lm_mdmr_validation_summary) <- compare_me$TreeFormat
colnames(lm_mdmr_validation_summary) <- c(   "meta_BMG_O.mash","meta_BMG_A.mash","meta_AW_A.mash",
                                             "meta_lnBMG_O.mash","meta_lnBMG_A.mash","meta_lnAW_A.mash",
                                             "meta_BMG_O.phylo","meta_BMG_A.phylo","meta_AW_A.phylo",
                                             "meta_lnBMG_O.phylo","meta_lnBMG_A.phylo","meta_lnAW_A.phylo",
                                             "BMG_O.mash","BMG_A.mash","AW_A.mash",
                                             "lnBMG_O.mash","lnBMG_A.mash","lnAW_A.mash",
                                             "BMG_O.phylo","BMG_A.phylo","AW_A.phylo",
                                             "lnBMG_O.phylo","lnBMG_A.phylo","lnAW_A.phylo",
                                             "meta_MLS_O.mash","meta_MLS_A.mash",
                                             "meta_MLS_O.phylo","meta_MLS_A.phylo")

lm_mdmr_validation_summary[1,"meta_BMG_O.mash"]   <- mdmr(X = compare_me$BMG01 , D = cmpr.mash, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_BMG_A.mash"]   <- mdmr(X = compare_me$BMG02 , D = cmpr.mash, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_AW_A.mash"]    <- mdmr(X = compare_me$AW ,    D = cmpr.mash, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_lnBMG_O.mash"] <- mdmr(X = compare_me$logBMG01 , D = cmpr.mash, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_lnBMG_A.mash"] <- mdmr(X = compare_me$logBMG02 , D = cmpr.mash, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_lnAW_A.mash"]  <- mdmr(X = compare_me$logAW ,    D = cmpr.mash, seed = 13579)$p.prec[1,1]

lm_mdmr_validation_summary[1,"meta_BMG_O.phylo"]   <- mdmr(X = compare_me$BMG01 , D = cmpr.phyl, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_BMG_A.phylo"]   <- mdmr(X = compare_me$BMG02 , D = cmpr.phyl, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_AW_A.phylo"]    <- mdmr(X = compare_me$AW ,    D = cmpr.phyl, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_lnBMG_O.phylo"] <- mdmr(X = compare_me$logBMG01 , D = cmpr.phyl, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_lnBMG_A.phylo"] <- mdmr(X = compare_me$logBMG02 , D = cmpr.phyl, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_lnAW_A.phylo"]  <- mdmr(X = compare_me$logAW ,    D = cmpr.phyl, seed = 13579)$p.prec[1,1]

lm_mdmr_validation_summary[1,"meta_MLS_O.mash"]   <- mdmr(X = compare_me$MLS01 , D = cmpr.mash, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_MLS_A.mash"]   <- mdmr(X = compare_me$MLS02 , D = cmpr.mash, seed = 13579)$p.prec[1,1]

lm_mdmr_validation_summary[1,"meta_MLS_O.phylo"]   <- mdmr(X = compare_me$MLS01 , D = cmpr.phyl, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"meta_MLS_A.phylo"]   <- mdmr(X = compare_me$MLS02 , D = cmpr.phyl, seed = 13579)$p.prec[1,1]

lm_mdmr_validation_summary[1,"BMG_O.mash"]      <- mdmr(X = residual_comparison.lm$BMG_O , D = cmpr.mash, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"BMG_A.mash"]      <- mdmr(X = residual_comparison.lm$BMG_A , D = cmpr.mash, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"AW_A.mash"]       <- mdmr(X = residual_comparison.lm$AW_A , D = cmpr.mash, seed = 13579)$p.prec[1,1]

lm_mdmr_validation_summary[1,"lnBMG_O.mash"]    <- mdmr(X = residual_comparison.lm$lnBMG_O , D = cmpr.mash, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"lnBMG_A.mash"]    <- mdmr(X = residual_comparison.lm$lnBMG_A , D = cmpr.mash, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"lnAW_A.mash"]     <- mdmr(X = residual_comparison.lm$lnAW_A , D = cmpr.mash, seed = 13579)$p.prec[1,1]

lm_mdmr_validation_summary[1,"BMG_O.phylo"]     <- mdmr(X = residual_comparison.lm$BMG_O , D = cmpr.phyl, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"BMG_A.phylo"]     <- mdmr(X = residual_comparison.lm$BMG_A , D = cmpr.phyl, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"AW_A.phylo"]      <- mdmr(X = residual_comparison.lm$AW_A , D = cmpr.phyl, seed = 13579)$p.prec[1,1]

lm_mdmr_validation_summary[1,"lnBMG_O.phylo"]   <- mdmr(X = residual_comparison.lm$lnBMG_O , D = cmpr.phyl, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"lnBMG_A.phylo"]   <- mdmr(X = residual_comparison.lm$lnBMG_A , D = cmpr.phyl, seed = 13579)$p.prec[1,1]
lm_mdmr_validation_summary[1,"lnAW_A.phylo"]    <- mdmr(X = residual_comparison.lm$lnAW_A , D = cmpr.phyl, seed = 13579)$p.prec[1,1]



lm_mdmr_validation_summary[2,"meta_BMG_O.mash"]   <- mdmr(X = compare_me$BMG01 , D = cmpr.mash, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_BMG_A.mash"]   <- mdmr(X = compare_me$BMG02 , D = cmpr.mash, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_AW_A.mash"]    <- mdmr(X = compare_me$AW ,    D = cmpr.mash, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_lnBMG_O.mash"] <- mdmr(X = compare_me$logBMG01 , D = cmpr.mash, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_lnBMG_A.mash"] <- mdmr(X = compare_me$logBMG02 , D = cmpr.mash, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_lnAW_A.mash"]  <- mdmr(X = compare_me$logAW ,    D = cmpr.mash, seed = 13579)$stat[1,1]

lm_mdmr_validation_summary[2,"meta_BMG_O.phylo"]   <- mdmr(X = compare_me$BMG01 , D = cmpr.phyl, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_BMG_A.phylo"]   <- mdmr(X = compare_me$BMG02 , D = cmpr.phyl, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_AW_A.phylo"]    <- mdmr(X = compare_me$AW ,    D = cmpr.phyl, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_lnBMG_O.phylo"] <- mdmr(X = compare_me$logBMG01 , D = cmpr.phyl, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_lnBMG_A.phylo"] <- mdmr(X = compare_me$logBMG02 , D = cmpr.phyl, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_lnAW_A.phylo"]  <- mdmr(X = compare_me$logAW ,    D = cmpr.phyl, seed = 13579)$stat[1,1]

lm_mdmr_validation_summary[2,"meta_MLS_O.mash"]   <- mdmr(X = compare_me$MLS01 , D = cmpr.mash, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_MLS_A.mash"]   <- mdmr(X = compare_me$MLS02 , D = cmpr.mash, seed = 13579)$stat[1,1]

lm_mdmr_validation_summary[2,"meta_MLS_O.phylo"]   <- mdmr(X = compare_me$MLS01 , D = cmpr.phyl, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"meta_MLS_A.phylo"]   <- mdmr(X = compare_me$MLS02 , D = cmpr.phyl, seed = 13579)$stat[1,1]


lm_mdmr_validation_summary[2,"BMG_O.mash"]      <- mdmr(X = residual_comparison.lm$BMG_O , D = cmpr.mash, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"BMG_A.mash"]      <- mdmr(X = residual_comparison.lm$BMG_A , D = cmpr.mash, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"AW_A.mash"]       <- mdmr(X = residual_comparison.lm$AW_A , D = cmpr.mash, seed = 13579)$stat[1,1]

lm_mdmr_validation_summary[2,"lnBMG_O.mash"]    <- mdmr(X = residual_comparison.lm$lnBMG_O , D = cmpr.mash, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"lnBMG_A.mash"]    <- mdmr(X = residual_comparison.lm$lnBMG_A , D = cmpr.mash, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"lnAW_A.mash"]     <- mdmr(X = residual_comparison.lm$lnAW_A , D = cmpr.mash, seed = 13579)$stat[1,1]

lm_mdmr_validation_summary[2,"BMG_O.phylo"]     <- mdmr(X = residual_comparison.lm$BMG_O , D = cmpr.phyl, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"BMG_A.phylo"]     <- mdmr(X = residual_comparison.lm$BMG_A , D = cmpr.phyl, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"AW_A.phylo"]      <- mdmr(X = residual_comparison.lm$AW_A , D = cmpr.phyl, seed = 13579)$stat[1,1]

lm_mdmr_validation_summary[2,"lnBMG_O.phylo"]   <- mdmr(X = residual_comparison.lm$lnBMG_O , D = cmpr.phyl, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"lnBMG_A.phylo"]   <- mdmr(X = residual_comparison.lm$lnBMG_A , D = cmpr.phyl, seed = 13579)$stat[1,1]
lm_mdmr_validation_summary[2,"lnAW_A.phylo"]    <- mdmr(X = residual_comparison.lm$lnAW_A , D = cmpr.phyl, seed = 13579)$stat[1,1]




lm_mdmr_validation_summary[3,"meta_BMG_O.mash"]   <- mdmr(X = compare_me$BMG01 , D = cmpr.mash, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_BMG_A.mash"]   <- mdmr(X = compare_me$BMG02 , D = cmpr.mash, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_AW_A.mash"]    <- mdmr(X = compare_me$AW ,    D = cmpr.mash, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_lnBMG_O.mash"] <- mdmr(X = compare_me$logBMG01 , D = cmpr.mash, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_lnBMG_A.mash"] <- mdmr(X = compare_me$logBMG02 , D = cmpr.mash, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_lnAW_A.mash"]  <- mdmr(X = compare_me$logAW ,    D = cmpr.mash, seed = 13579)$pr.sq[1,1]

lm_mdmr_validation_summary[3,"meta_BMG_O.phylo"]   <- mdmr(X = compare_me$BMG01 , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_BMG_A.phylo"]   <- mdmr(X = compare_me$BMG02 , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_AW_A.phylo"]    <- mdmr(X = compare_me$AW ,    D = cmpr.phyl, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_lnBMG_O.phylo"] <- mdmr(X = compare_me$logBMG01 , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_lnBMG_A.phylo"] <- mdmr(X = compare_me$logBMG02 , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_lnAW_A.phylo"]  <- mdmr(X = compare_me$logAW ,    D = cmpr.phyl, seed = 13579)$pr.sq[1,1]

lm_mdmr_validation_summary[3,"meta_MLS_O.mash"]   <- mdmr(X = compare_me$MLS01 , D = cmpr.mash, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_MLS_A.mash"]   <- mdmr(X = compare_me$MLS02 , D = cmpr.mash, seed = 13579)$pr.sq[1,1]

lm_mdmr_validation_summary[3,"meta_MLS_O.phylo"]   <- mdmr(X = compare_me$MLS01 , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"meta_MLS_A.phylo"]   <- mdmr(X = compare_me$MLS02 , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]


lm_mdmr_validation_summary[3,"BMG_O.mash"]      <- mdmr(X = residual_comparison.lm$BMG_O , D = cmpr.mash, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"BMG_A.mash"]      <- mdmr(X = residual_comparison.lm$BMG_A , D = cmpr.mash, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"AW_A.mash"]       <- mdmr(X = residual_comparison.lm$AW_A , D = cmpr.mash, seed = 13579)$pr.sq[1,1]

lm_mdmr_validation_summary[3,"lnBMG_O.mash"]    <- mdmr(X = residual_comparison.lm$lnBMG_O , D = cmpr.mash, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"lnBMG_A.mash"]    <- mdmr(X = residual_comparison.lm$lnBMG_A , D = cmpr.mash, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"lnAW_A.mash"]     <- mdmr(X = residual_comparison.lm$lnAW_A , D = cmpr.mash, seed = 13579)$pr.sq[1,1]

lm_mdmr_validation_summary[3,"BMG_O.phylo"]     <- mdmr(X = residual_comparison.lm$BMG_O , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"BMG_A.phylo"]     <- mdmr(X = residual_comparison.lm$BMG_A , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"AW_A.phylo"]      <- mdmr(X = residual_comparison.lm$AW_A , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]

lm_mdmr_validation_summary[3,"lnBMG_O.phylo"]   <- mdmr(X = residual_comparison.lm$lnBMG_O , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"lnBMG_A.phylo"]   <- mdmr(X = residual_comparison.lm$lnBMG_A , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]
lm_mdmr_validation_summary[3,"lnAW_A.phylo"]    <- mdmr(X = residual_comparison.lm$lnAW_A , D = cmpr.phyl, seed = 13579)$pr.sq[1,1]


write.table(lm_mdmr_validation_summary, file="lm_mdmr_validation_summary.tsv", sep='\t')


#################################################################################################################################################################################################

#################################################################################################################################################################################################

#################################################################################################################################################################################################


# abundance comparisons

abund_ref <- read.table("input_x01_modelExe/Ref.OG.abundance2.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
abund_ref <- abund_ref[-dim(abund_ref)[1],]

# clean out single quotes
abund_dno <- read.table("input_x01_modelExe/Denovo.OG.abundance3.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
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

# rank 1 imputation with svd
abund_ref.cts.rel.impute <- abund_ref.cts.rel
abund_ref.cts.rel.impute[ abund_ref.cts.rel.impute == 0 ] <- NA
abund_ref.cts.rel.impute <- impute.svd(abund_ref.cts.rel.impute, 1)$x
rownames(abund_ref.cts.rel.impute) <- abund_ref$OrthoDB_OG_ID
abund_dno.cts.rel.impute <- abund_dno.cts.rel
abund_dno.cts.rel.impute[ abund_dno.cts.rel.impute == 0 ] <- NA
abund_dno.cts.rel.impute <- impute.svd(abund_dno.cts.rel.impute, 1)$x
rownames(abund_dno.cts.rel.impute) <- abund_dno$OrthoDB_OG_ID

abund_ref.cts.rel.LOGimpute <- log(abund_ref.cts.rel)
abund_ref.cts.rel.LOGimpute[ abund_ref.cts.rel == 0 ] <- NA
abund_ref.cts.rel.LOGimpute <- impute.svd(abund_ref.cts.rel.LOGimpute, 1)$x
rownames(abund_ref.cts.rel.LOGimpute) <- abund_ref$OrthoDB_OG_ID
abund_dno.cts.rel.LOGimpute <- log(abund_dno.cts.rel)
abund_dno.cts.rel.LOGimpute[ abund_dno.cts.rel == 0 ] <- NA
abund_dno.cts.rel.LOGimpute <- impute.svd(abund_dno.cts.rel.LOGimpute, 1)$x
rownames(abund_dno.cts.rel.LOGimpute) <- abund_dno$OrthoDB_OG_ID

og_correlation <- as.data.frame(matrix())

rsq <- function (x, y) cor(x, y) ^ 2

oo <- 1


for(og in rownames(abund_ref.cts.rel)){
  
  if( oo %% 20 == 0 ){  print(paste(og, (oo/length(rownames(abund_ref.cts.rel)))*100, '%' ))  }
  
  og_correlation[ og ,  "denovo.count" ] <- sum(abund_dno.cts.rel[og,] > 0)
  og_correlation[ og ,  "alignm.count" ] <- sum(abund_ref.cts.rel[og,] > 0)
  
  dontomit <- abund_dno.cts.rel[og,] != 0
  dontomit[ abund_ref.cts.rel[og,] == 0 ] <- FALSE
  
  if(sum(dontomit) > 1){
    
    png(paste0("corr/correlation_of_OGabs.",og,".png"),height=800,width=600)
    plot.new()
    par(mfrow=c(4,2))
    plot(log(abund_ref.cts.rel[og,]), log(abund_dno.cts.rel[og,]),main=paste(og, "RelAbs, Zeroes Omitted"), ylab = "log(De Novo Relative Abundance) (Non-Imputed)", xlab="log(Reference-derived Abundance) (Non-Imputed)",pch=16)
    fit <- lm( log(abund_dno.cts.rel[og,])[dontomit] ~ log(abund_ref.cts.rel[og,])[dontomit]  )
    abline(fit$coefficients[1],fit$coefficients[2])
    r2 <- rsq(log(abund_dno.cts.rel[og,])[dontomit], log(abund_ref.cts.rel[og,])[dontomit])
    fstat    <- summary(fit)$fstatistic
    pval    <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
    legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd(abs(fit$residuals)),3)),paste0("i: ",round(fit$coefficients[1],3)), paste0("r2: ", round(r2,3)), paste0("log(p): ", round(log(pval),3))))
    og_correlation[ og ,  "RelAbs.slope"  ]     <- fit$coefficients[2]
    og_correlation[ og ,  "RelAbs.intercept"  ] <- fit$coefficients[3]
    og_correlation[ og ,  "RelAbs.sd"  ]        <- sd(abs(fit$residuals))
    og_correlation[ og ,  "RelAbs.r2"  ]        <- r2
    og_correlation[ og ,  "RelAbs.pval"  ]      <- pval
    
    plot( log(abund_ref.cts.rel.impute[og,]), log(abund_dno.cts.rel.impute[og,]), main=paste(og, "RelAbs (Imputed>Log)"), ylab = "log(De Novo-derived Relative Abundance) (Imputed)", xlab = "log(Reference-derived Relative Abundance) (Imputed)",pch=16)
    fit <- lm( log(abund_dno.cts.rel.impute[og,]) ~ log(abund_ref.cts.rel.impute[og,]) )
    abline(fit$coefficients[1],fit$coefficients[2])
    r2 <- rsq(log(abund_dno.cts.rel.impute[og,]), log(abund_ref.cts.rel.impute[og,]))
    fstat    <- summary(fit)$fstatistic
    pval    <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
    legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd(abs(fit$residuals)),3)),paste0("i: ",round(fit$coefficients[1],3)), paste0("r2: ", round(r2,3)), paste0("log(p): ", round(log(pval),3))))
    og_correlation[ og ,  "ImputedRelAbs.slope"  ]     <- fit$coefficients[2]
    og_correlation[ og ,  "ImputedRelAbs.intercept"  ] <- fit$coefficients[3]
    og_correlation[ og ,  "ImputedRelAbs.sd"  ]        <- sd(abs(fit$residuals))
    og_correlation[ og ,  "ImputedRelAbs.r2"  ]        <- r2
    og_correlation[ og ,  "ImputedRelAbs.pval"  ]      <- pval
    
    plot( abund_ref.cts.rel.LOGimpute[og,], abund_dno.cts.rel.LOGimpute[og,], main=paste(og, "RelAbs (Log>Imputed)"), ylab = "log(De Novo-derived Relative Abundance) (Imputed)", xlab = "log(Reference-derived Relative Abundance) (Imputed)",pch=16)
    fit <- lm( abund_dno.cts.rel.LOGimpute[og,] ~ abund_ref.cts.rel.LOGimpute[og,] )
    abline(fit$coefficients[1],fit$coefficients[2])
    r2 <- rsq(abund_dno.cts.rel.LOGimpute[og,], abund_ref.cts.rel.LOGimpute[og,])
    fstat    <- summary(fit)$fstatistic
    pval    <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
    legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd(abs(fit$residuals)),3)),paste0("i: ",round(fit$coefficients[1],3)), paste0("r2: ", round(r2,3)), paste0("log(p): ", round(log(pval),3))))
    og_correlation[ og ,  "LOGImputedRelAbs.slope"  ]     <- fit$coefficients[2]
    og_correlation[ og ,  "LOGImputedRelAbs.intercept"  ] <- fit$coefficients[3]
    og_correlation[ og ,  "LOGImputedRelAbs.sd"  ]        <- sd(abs(fit$residuals))
    og_correlation[ og ,  "LOGImputedRelAbs.r2"  ]        <- r2
    og_correlation[ og ,  "LOGImputedRelAbs.pval"  ]      <- pval
    
    plot(order(abund_ref.cts.rel[og,][dontomit]), order(abund_dno.cts.rel[og,][dontomit]), main=paste(og, "Rank Order"), ylab = "Order (De Novo Abundance)", xlab="Order (Reference Abundance)",pch=16)
    fit <- lm( order(abund_dno.cts.rel[og,][dontomit]) ~ order(abund_ref.cts.rel[og,][dontomit]) )
    abline(fit$coefficients[1],fit$coefficients[2])
    r2 <- rsq(log(abund_dno.cts.rel.impute[og,]), log(abund_ref.cts.rel.impute[og,]))
    fstat    <- summary(fit)$fstatistic
    pval    <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
    legend("bottomright", legend= c(paste0("s: ",round(fit$coefficients[2],3)),paste0("sd: ",round(sd(abs(fit$residuals)),3)),paste0("i: ",round(fit$coefficients[1],3))))
    og_correlation[ og ,  "Order.slope"  ]     <- fit$coefficients[2]
    og_correlation[ og ,  "Order.intercept"  ] <- fit$coefficients[3]
    og_correlation[ og ,  "Order.sd"  ]        <- sd(abs(fit$residuals))
    og_correlation[ og ,  "Order.ImputedRelAbs.r2"  ]        <- r2
    og_correlation[ og ,  "Order.pval"  ]      <- pval
    
    resid1.order <- resid(lm(order(log(abund_dno.cts.rel[og,])[dontomit]) ~ order(log(abund_ref.cts.rel[og,])[dontomit])))
    resid2.numbr <- resid(lm(log(abund_dno.cts.rel[og,])[dontomit] ~ log(abund_ref.cts.rel[og,])[dontomit]))
    ccc01 <- CCC( resid1.order , resid2.numbr )
    tmp.mean <- mean(ccc01$blalt$delta)
    tmp.sd <- sqrt(var(ccc01$blalt$delta))
    plot(ccc01$blalt$mean, main="CCC",pch=16)
    abline(h = tmp.mean, lty = 1, col = "gray")
    abline(h = tmp.mean - (2 * tmp.sd), lty = 2, col = "gray")
    abline(h = tmp.mean + (2 * tmp.sd), lty = 2, col = "gray")
    legend(x = "topleft", legend = c("Mean difference",                                  "Mean difference +/ 2SD"), lty = c(1,2), bty = "n")
    legend("bottomright", legend=paste(names(ccc01$rho.c),round(unname(ccc01$rho.c), 3)))
    og_correlation[ og ,  "CCC.mean"  ]        <- tmp.mean
    og_correlation[ og ,  "CCC.sd"  ]          <- tmp.sd
    og_correlation[ og ,  "CCC.estRho"  ]      <- unname(ccc01$rho.c)[1]
    og_correlation[ og ,  "CCC.lwr.ci"  ]      <- unname(ccc01$rho.c)[2]
    og_correlation[ og ,  "CCC.upr.ci"  ]      <- unname(ccc01$rho.c)[3]
    
    dev.off()
    
  }else{
    print(paste(og,"ZERO CASE"))
  }
  oo <- oo + 1
  
}

#save.image("develo20200531_4pm.Rdata")
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
load("develo20200531_4pm.Rdata")

# save.image("20200210.Rdata")

og_correlation <- og_correlation[-1,-1]

write.table(og_correlation, file = "og_correlation.20200531.tsv", sep="\t")

for(t in 1:dim(og_correlation)[1]){
  og_correlation[t,"min.count"] <- min( og_correlation[t,]$denovo.count, og_correlation[t,]$alignm.count )
}

#  END INPUT RE-PARSING ##########################################################################################################################################################################
##################################################################################################################################################################################################


formulaz <- c(  as.formula("MLS01   ~ abund_ref.cts.rel"), 
                as.formula("MLS02   ~ abund_ref.cts.rel"),
                as.formula("MLS01   ~ abund_ref.cts.rel.ord") , 
                as.formula("MLS02   ~ abund_ref.cts.rel.ord")  ,
                as.formula("MLS01   ~ abund_ref.cts.rel.LOGimpute"),  
                as.formula("MLS02   ~ abund_ref.cts.rel.LOGimpute") , 
                as.formula("BMG_O   ~ abund_ref.cts.rel")  ,
                as.formula("BMG_A   ~ abund_ref.cts.rel")   ,
                as.formula("AW_A    ~ abund_ref.cts.rel")   ,
                as.formula("lnBMG_O ~ abund_ref.cts.rel") ,
                as.formula("lnBMG_A ~ abund_ref.cts.rel") ,
                as.formula("lnAW_A  ~ abund_ref.cts.rel"),
                as.formula("BMG_O   ~ abund_ref.cts.rel.ord"),  
                as.formula("BMG_A   ~ abund_ref.cts.rel.ord") ,  
                as.formula("AW_A    ~ abund_ref.cts.rel.ord")  , 
                as.formula("lnBMG_O ~ abund_ref.cts.rel.ord") ,
                as.formula("lnBMG_A ~ abund_ref.cts.rel.ord") ,
                as.formula("lnAW_A  ~ abund_ref.cts.rel.ord"),
                as.formula("BMG_O   ~ abund_ref.cts.rel.LOGimpute"),  
                as.formula("BMG_A   ~ abund_ref.cts.rel.LOGimpute") ,  
                as.formula("AW_A    ~ abund_ref.cts.rel.LOGimpute")  , 
                as.formula("lnBMG_O ~ abund_ref.cts.rel.LOGimpute") ,
                as.formula("lnBMG_A ~ abund_ref.cts.rel.LOGimpute") ,
                as.formula("lnAW_A  ~ abund_ref.cts.rel.LOGimpute") ,
                as.formula("MLS01   ~ abund_dno.cts.rel"), 
                as.formula("MLS02   ~ abund_dno.cts.rel"),
                as.formula("MLS01   ~ abund_dno.cts.rel.ord") , 
                as.formula("MLS02   ~ abund_dno.cts.rel.ord")  ,
                as.formula("MLS01   ~ abund_dno.cts.rel.LOGimpute"),  
                as.formula("MLS02   ~ abund_dno.cts.rel.LOGimpute") , 
                as.formula("BMG_O   ~ abund_dno.cts.rel")  ,
                as.formula("BMG_A   ~ abund_dno.cts.rel")   ,
                as.formula("AW_A    ~ abund_dno.cts.rel")   ,
                as.formula("lnBMG_O ~ abund_dno.cts.rel") ,
                as.formula("lnBMG_A ~ abund_dno.cts.rel") ,
                as.formula("lnAW_A  ~ abund_dno.cts.rel"),
                as.formula("BMG_O   ~ abund_dno.cts.rel.ord"),  
                as.formula("BMG_A   ~ abund_dno.cts.rel.ord") ,  
                as.formula("AW_A    ~ abund_dno.cts.rel.ord")  , 
                as.formula("lnBMG_O ~ abund_dno.cts.rel.ord") ,
                as.formula("lnBMG_A ~ abund_dno.cts.rel.ord") ,
                as.formula("lnAW_A  ~ abund_dno.cts.rel.ord"),
                as.formula("BMG_O   ~ abund_dno.cts.rel.LOGimpute"),  
                as.formula("BMG_A   ~ abund_dno.cts.rel.LOGimpute") ,  
                as.formula("AW_A    ~ abund_dno.cts.rel.LOGimpute")  , 
                as.formula("lnBMG_O ~ abund_dno.cts.rel.LOGimpute") ,
                as.formula("lnBMG_A ~ abund_dno.cts.rel.LOGimpute") ,
                as.formula("lnAW_A  ~ abund_dno.cts.rel.LOGimpute")
)


for(d in colnames(residual_comparison.lm)){
  compare_me[,d] <- residual_comparison.lm[,d]
}

model_performance <- as.data.frame(matrix())
resid_performance <- as.data.frame(matrix())
mc <- 1
###for( m in 49:0 ){
m <- 49
oc <- 1
for( o in rownames(og_correlation[ og_correlation$min.count==m, ]) ){
  print(paste(oc, m, "(",round((oc/length(rownames(og_correlation[ og_correlation$min.count==m, ])))*100,2),"%)    ----    ",oc,o,"(",round((oc/length(rownames(og_correlation[ og_correlation$min.count==m, ])))*100,2),"%)"))
  
  compare_me$abund_ref.cts.rel <- as.numeric(abund_ref.cts.rel[o,query_metadata$AbundanceLabel])
  compare_me$abund_dno.cts.rel <- as.numeric(abund_dno.cts.rel[o,])
  compare_me$abund_ref.cts.rel.ord <- as.numeric(order(abund_ref.cts.rel[o,]))
  compare_me$abund_dno.cts.rel.ord <- as.numeric(order(abund_dno.cts.rel[o,]))
  compare_me$abund_ref.cts.rel.LOGimpute <- as.numeric(abund_ref.cts.rel.LOGimpute[o,])
  compare_me$abund_dno.cts.rel.LOGimpute <- as.numeric(abund_dno.cts.rel.LOGimpute[o,])
  #compare_me$abund_ref.cts.rel.LOGimpute.ord <- as.numeric(order(abund_ref.cts.rel.LOGimpute[o,]))
  #compare_me$abund_dno.cts.rel.LOGimpute.ord <- as.numeric(order(abund_dno.cts.rel.LOGimpute[o,]))
  compare_me$abund_ref.cts.binary <- as.factor( as.numeric(abund_ref.cts.rel[o,query_metadata$AbundanceLabel]) > 0 )
  compare_me$abund_dno.cts.binary <- as.factor( abund_dno.cts.rel[o,] > 0 )
  
  tryCatch( compdata <- comparative.data(phylotree, compare_me, TreeFormat) ,  error = function(e) message(paste("ROBUST FAIL")))
  
  ###################
  #     numeric     #
  ###################
  fc <- 1
  for(formz in formulaz){
    if(fc %in% c(1:24) ){ abund_type <- "ref." }
    if(fc %in% c(25:48)){ abund_type <- "dno." }
    
    formzPrint <- paste0( as.character(formz)[2] , "~" , as.character(formz)[3] )
    print(paste("!",o,formzPrint))
    
    # lm
    reg.lm.f_01   <- lm    (formz, data=compare_me)
    if(!is.na(reg.lm.f_01)){
      if(!is.na(coef(reg.lm.f_01))){
        if (t(coef(summary(reg.lm.f_01)))[2,1] > 0 & is.nan(t(coef(summary(reg.lm.f_01)))[2,1]) == 0){
          cpus.lm2 <- as.matrix( anova(reg.lm.f_01) )
          fstat    <- summary(reg.lm.f_01)$fstatistic
          model_performance[o,paste0(fc,".",formzPrint,".LM.p_val")] <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
          model_performance[o,paste0(fc,".",formzPrint,".LM.slope")] <- coef(pgls(formz,data=compdata))[2]
          model_performance[o,paste0(fc,".",formzPrint,".LM.itcpt")] <- coef(pgls(formz,data=compdata))[1]
          model_performance[o,paste0(fc,".",formzPrint,"LM.ispos")] <- as.numeric(coef(pgls(formz,data=compdata))[2] > 0)
          resid_performance[o,paste0(fc,".",formzPrint,"LM","-",query_metadata$Common.Name)] <- unname(resid(reg.lm.f_01))
        }else{ model_performance[o,paste0(fc,".",formzPrint,".LM.p_val")] <- model_performance[o,paste0(fc,".",formzPrint,".LM.slope")]    <- model_performance[o,paste0(fc,".",formzPrint,".LM.itcpt")]   <- model_performance[o,paste0(fc,".",formzPrint,".LM.ispos")]      <- NA }
      }else{   model_performance[o,paste0(fc,".",formzPrint,".LM.p_val")] <- model_performance[o,paste0(fc,".",formzPrint,".LM.slope")]    <- model_performance[o,paste0(fc,".",formzPrint,".LM.itcpt")]   <- model_performance[o,paste0(fc,".",formzPrint,".LM.ispos")]      <- NA }
    }else{     model_performance[o,paste0(fc,".",formzPrint,".LM.p_val")] <- model_performance[o,paste0(fc,".",formzPrint,".LM.slope")]    <- model_performance[o,paste0(fc,".",formzPrint,".LM.itcpt")]   <- model_performance[o,paste0(fc,".",formzPrint,".LM.ispos")]      <- NA }
    
    # phylogeny contrast fit
    reg.ph.f_01   <- crunch(formz, data=compdata) #brunch(f, data=compdata, robust = 4) ,  error = function(e) message(paste("FAIL",og)))
    if(!is.na(reg.ph.f_01)){
      if(!is.na(coef(reg.ph.f_01))){
        if (t(coef(summary(reg.ph.f_01)))[2,1] > 0 & is.nan(t(coef(summary(reg.ph.f_01)))[2,1]) == 0){
          cpus.lm2 <- as.matrix( anova(reg.ph.f_01) )
          fstat    <- summary(reg.ph.f_01)$fstatistic
          model_performance[o,paste0(fc,".",formzPrint,".CAPER.p_val")] <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
          model_performance[o,paste0(fc,".",formzPrint,".CAPER.slope")] <- coef(pgls(formz,data=compdata))[2]
          model_performance[o,paste0(fc,".",formzPrint,".CAPER.itcpt")] <- coef(pgls(formz,data=compdata))[1]
          model_performance[o,paste0(fc,".",formzPrint,".CAPER.ispos")] <- as.numeric(coef(pgls(formz,data=compdata))[2] > 0)
          resid_performance[o,paste0(fc,".",formzPrint,".CAPER","-",names(resid(reg.ph.f_01)))] <- unname(resid(reg.ph.f_01))
        }else{ model_performance[o,paste0(fc,".",formzPrint,".CAPER.p_val")] <- model_performance[o,paste0(fc,".",formzPrint,".CAPER.slope")]    <- model_performance[o,paste0(fc,".",formzPrint,".CAPER.itcpt")]   <- model_performance[o,paste0(fc,".",formzPrint,".CAPER.ispos")]      <- NA }
      }else{ model_performance[o,paste0(fc,".",formzPrint,".CAPER.p_val")] <- model_performance[o,paste0(fc,".",formzPrint,".CAPER.slope")]    <- model_performance[o,paste0(fc,".",formzPrint,".CAPER.itcpt")]   <- model_performance[o,paste0(fc,".",formzPrint,"..CAPER.ispos")]      <- NA }
    }else{ model_performance[o,paste0(fc,".",formzPrint,".CAPER.p_val")] <- model_performance[o,paste0(fc,".",formzPrint,".CAPER.slope")]    <- model_performance[o,paste0(fc,".",formzPrint,".CAPER.itcpt")]   <- model_performance[o,paste0(fc,".",formzPrint,".CAPER.ispos")]      <- NA }
    
    # run mdmr comparisons - lm 
    resi.f_01 <- resid(reg.lm.f_01) 
    names(resi.f_01) <- query_metadata$PhylotreeLabel
    yy.f_01 <- list()
    y <- 1
    for(t in rownames(as.matrix(distTips(phylotree, tips = "all")))){
      yy.f_01[y] <- resi.f_01[t]
      y <- y + 1
    }
    yy.f_01 <- unlist(yy.f_01)
    mdmr.res.f_01.phyl <- mdmr(X = yy.f_01, D = as.matrix(cmpr.phyl), seed = 13579)
    mdmr.res.f_01.mash <- mdmr(X = yy.f_01, D = as.matrix(cmpr.mash), seed = 13579)
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_phylo.stat")]    <- mdmr.res.f_01.phyl$stat[1,1]
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_phylo.pr.sql")]  <- mdmr.res.f_01.phyl$pr.sq[1,1]
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_phylo.p_val")]      <- mdmr.res.f_01.phyl$pv[1,1]
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_phylo.p_val.FDR")]      <- p.adjust( model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_phylo.p_val")] , method="fdr" )
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_phylo.p.prec")]  <- mdmr.res.f_01.phyl$p.prec[1,1]
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_phylo.df")]      <- mdmr.res.f_01.phyl$df[1,1]
    
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_mash.stat")]    <- mdmr.res.f_01.mash$stat[1,1]
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_mash.pr.sql")]  <- mdmr.res.f_01.mash$pr.sq[1,1]
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_mash.p_val")]      <- mdmr.res.f_01.mash$pv[1,1]
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_mash.p.prec")]  <- mdmr.res.f_01.mash$p.prec[1,1]
    model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_mash.df")]      <- mdmr.res.f_01.mash$df[1,1]
    
    
    # run alignment/denovo correlation analysess
    if(fc > 24){
      
      # lm context
      formzPrint_A <- paste0( as.character(formulaz[fc-24][[1]])[2] , "~" , as.character(formulaz[fc-24][[1]])[3] )
      formzPrint_B <- paste0( as.character(formz)[2] , "~" , as.character(formz)[3] )
      
      cq.lm <- CCC(  as.numeric(unname( resid_performance[o,paste0(fc-24,".",formzPrint_A,"LM","-",query_metadata$Common.Name)] ) ) ,
                     as.numeric(unname( resid_performance[o,paste0(fc,".",formzPrint_B,"LM","-",query_metadata$Common.Name)] ) ) 
      )
      model_performance[ o , paste0(fc,".",formzPrint_B,".lm_est.rho") ] <- cq.lm$rho.c[1]
      model_performance[ o , paste0(fc,".",formzPrint_B,".lm_lwr.ci") ]  <- cq.lm$rho.c[2]
      model_performance[ o , paste0(fc,".",formzPrint_B,".lm_upr.ci") ]  <- cq.lm$rho.c[3]
      model_performance[ o , paste0(fc,".",formzPrint_B,".lm_mean") ]    <- mean(cq.lm$blalt$delta)
      model_performance[ o , paste0(fc,".",formzPrint_B,".lm_stdev") ]   <- sqrt(var(cq.lm$blalt$delta))
      
      #caper context
      cq.caper <- CCC(  as.numeric(unname( resid_performance[o,paste0(fc-24,".",formzPrint_A,".CAPER","-",names(resid(reg.ph.f_01)))] ) ),
                        as.numeric(unname( resid_performance[o,paste0(fc,".",formzPrint_B,".CAPER","-",names(resid(reg.ph.f_01)))] ) )
      )
      model_performance[ o , paste0(fc,".",formzPrint_B,".caper_est.rho") ] <- cq.lm$rho.c[1]
      model_performance[ o , paste0(fc,".",formzPrint_B,".caper_lwr.ci") ]  <- cq.lm$rho.c[2]
      model_performance[ o , paste0(fc,".",formzPrint_B,".caper_upr.ci") ]  <- cq.lm$rho.c[3]
      model_performance[ o , paste0(fc,".",formzPrint_B,".caper_mean") ]    <- mean(cq.lm$blalt$delta)
      model_performance[ o , paste0(fc,".",formzPrint_B,".caper_stdev") ]   <- sqrt(var(cq.lm$blalt$delta))
      
    }
    
    
    
    fc <- fc + 1
  } # END OF FOR LOOP
  oc <- oc + 1
}







# save.image("develo20200609.Rdata")
# # save.image("develo20200601_8am.Rdata")

write.table(model_performance, file="model_performance.2020610.tsv", sep="\t")
write.table(resid_performance, file="resid_performance.2020610.tsv", sep="\t")

# FDR adjust p-values
pp <-  colnames(model_performance)[ grepl('p_val', colnames(model_performance)) ]
for(p in pp){
  model_performance[ ,  paste0(p,".FDR") ] <- p.adjust( model_performance[ , p ] , method="fdr" )
}

#model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_mash.p_val.FDR")]      <- p.adjust( model_performance[o,paste0(fc,".",formzPrint,".lm_MDMR_mash.p_val")] , method="fdr" )







#################################################################################################################################################################################################

#################################################################################################################################################################################################

#################################################################################################################################################################################################


# fix mislabelled correlation terms

mc <- 1
###for( m in 49:0 ){
m <- 49
oc <- 1
for( o in rownames(og_correlation[ og_correlation$min.count==m, ]) ){
  
  print(paste(oc, m, "(",round((oc/length(rownames(og_correlation[ og_correlation$min.count==m, ])))*100,2),"%)    ----    ",oc,o,"(",round((oc/length(rownames(og_correlation[ og_correlation$min.count==m, ])))*100,2),"%)"))
  
  fc <- 1
  for(formz in formulaz){
    if(fc %in% c(1:24) ){ abund_type <- "ref." }
    if(fc %in% c(25:48)){ abund_type <- "dno." }
    
    # run alignment/denovo correlation analysess
    if(fc > 24){
      #formzPrint <- paste0( as.character(formz)[2] , "~" , as.character(formz)[3] )
      #print(paste("correlation re-run!",o,formzPrint))
      
      # lm context
      formzPrint_A <- paste0( as.character(formulaz[fc-24][[1]])[2] , "~" , as.character(formulaz[fc-24][[1]])[3] )
      formzPrint_B <- paste0( as.character(formz)[2] , "~" , as.character(formz)[3] )
      
      cq.lm <- CCC(  as.numeric(unname( resid_performance[o,paste0(fc-24,".",formzPrint_A,"LM","-",query_metadata$Common.Name)] ) ) ,
                     as.numeric(unname( resid_performance[o,paste0(fc,".",formzPrint_B,"LM","-",query_metadata$Common.Name)] ) ) 
      )
      model_performance[ o , paste0(fc,".",formzPrint_B,".lm_est.rho") ] <- cq.lm$rho.c[1]
      model_performance[ o , paste0(fc,".",formzPrint_B,".lm_lwr.ci") ]  <- cq.lm$rho.c[2]
      model_performance[ o , paste0(fc,".",formzPrint_B,".lm_upr.ci") ]  <- cq.lm$rho.c[3]
      model_performance[ o , paste0(fc,".",formzPrint_B,".lm_mean") ]    <- mean(cq.lm$blalt$delta)
      model_performance[ o , paste0(fc,".",formzPrint_B,".lm_stdev") ]   <- sqrt(var(cq.lm$blalt$delta))
      
      #caper context
      cq.caper <- CCC(  as.numeric(unname( resid_performance[o,paste0(fc-24,".",formzPrint_A,".CAPER","-",names(resid(reg.ph.f_01)))] ) ),
                        as.numeric(unname( resid_performance[o,paste0(fc,".",formzPrint_B,".CAPER","-",names(resid(reg.ph.f_01)))] ) )
      )
      model_performance[ o , paste0(fc,".",formzPrint_B,".caper_est.rho") ] <- cq.caper$rho.c[1]
      model_performance[ o , paste0(fc,".",formzPrint_B,".caper_lwr.ci") ]  <- cq.caper$rho.c[2]
      model_performance[ o , paste0(fc,".",formzPrint_B,".caper_upr.ci") ]  <- cq.caper$rho.c[3]
      model_performance[ o , paste0(fc,".",formzPrint_B,".caper_mean") ]    <- mean(cq.caper$blalt$delta)
      model_performance[ o , paste0(fc,".",formzPrint_B,".caper_stdev") ]   <- sqrt(var(cq.caper$blalt$delta))
      
    }
    fc <- fc + 1
  } 
  oc <- oc + 1
}

write.table(model_performance, file="model_performance.2020610b.tsv", sep="\t")
write.table(resid_performance, file="resid_performance.2020610b.tsv", sep="\t")


#################################################################################################################################################################################################

#################################################################################################################################################################################################

#################################################################################################################################################################################################



#par(mfrow=c())
par(mfrow=c(1,1),
    mar=c(2,2,2,2))
# pp <-  colnames(model_performance)[ grepl('p_val', colnames(model_performance)) ]
cc <- colnames(model_performance)[ grepl('slope', colnames(model_performance)) ]

belowCutoff <- as.data.frame(matrix())
belowCutoff.fdr <- as.data.frame(matrix())

for(c in cc){
  
  ty <- grepl( 'FDR' , c )
  if( ty == TRUE){
    prefix <- substr(c,1,nchar(c)-10)
    yylab <- '-log( FDR-adjusted(p) )'
  }else{
    prefix <- substr(c,1,nchar(c)-6)
    yylab <- '-log( p )'
  }
  pValue   <- paste0( prefix , '.p_val' )
  fdrValue <- paste0( pValue , ".FDR" )
  xmax <- max(model_performance[ , c][ !is.na(model_performance[ , c])  ])
  
  estrhocols <- colnames(model_performance)[ grepl('est.rho',colnames(model_performance)) ]
  iter <- strsplit(c,'\\.')[[1]][1]
  
  #if( sum( grepl(iter,estrhocols) ) > 0 ){
  if(iter <=24){  
    lmid <- estrhocols[ grepl(as.character(as.numeric(iter)+24),estrhocols) ][ grepl('lm',estrhocols[ grepl( as.character(as.numeric(iter)+24) ,estrhocols) ]) ]
    cpid <- estrhocols[ grepl(as.character(as.numeric(iter)+24),estrhocols) ][ grepl('caper',estrhocols[ grepl(as.character(as.numeric(iter)+24) ,estrhocols) ]) ]
  }else{
    lmid <- estrhocols[ grepl(iter,estrhocols) ][ grepl('lm',   estrhocols[ grepl(iter,estrhocols) ]) ]
    cpid <- estrhocols[ grepl(iter,estrhocols) ][ grepl('caper',estrhocols[ grepl(iter,estrhocols) ]) ]
  }
  
  png(paste0("volcano/rawP.",prefix,".png"), width=600, height=600)
  
  if(grepl('.LM',prefix)==TRUE){    
    colmin <- min(round(model_performance[,lmid]*100)[!is.na(round(model_performance[,lmid]*100))])
    colmax <- max(round(model_performance[,lmid]*100)[!is.na(round(model_performance[,lmid]*100))])
    if(colmin < 0){
      colcol <- round(model_performance[,lmid]*100)[!is.na(round(model_performance[,lmid]*100))]+abs(colmin)
    }else{
      colcol <- round(model_performance[,lmid]*100)[!is.na(round(model_performance[,lmid]*100))]
    }
    colcol <- redgreen(max(colcol))[colcol]
    #green2red(min(round(model_performance[,lmid]*100)[!is.na(round(model_performance[,lmid]*100))]))[ round(model_performance[,lmid]*100) - min(round(model_performance[,lmid]*100)[!is.na(round(model_performance[,lmid]*100))]) ]
  }else{
    colmin <- min(round(model_performance[,cpid]*100)[!is.na(round(model_performance[,cpid]*100))])
    colmax <- max(round(model_performance[,cpid]*100)[!is.na(round(model_performance[,cpid]*100))])
    if(colmin < 0){
      colcol <- round(model_performance[,cpid]*100)[!is.na(round(model_performance[,cpid]*100))]+abs(colmin)
    }else{
      colcol <- round(model_performance[,lmid]*100)[!is.na(round(model_performance[,lmid]*100))]
    }
    colcol <- redgreen(max(colcol))[colcol]
    #green2red(min(round(model_performance[,cpid]*100)[!is.na(round(model_performance[,cpid]*100))]))[ round(model_performance[,cpid]*100) - min(round(model_performance[,cpid]*100)[!is.na(round(model_performance[,cpid]*100))]) ]
  }
  
  plot(      model_performance[ , c] ,
             -log( model_performance[ , pValue] ),
             main= paste( prefix , ": Raw P-Values" ),
             xlab="Slope of Association",
             ylab=yylab,
             xlim=c(-xmax,xmax),
             col=colcol,
             pch=16
  )
  
  for(l in c(0.05,0.005,0.0005,0.00005)){
    abline(h=-log(l),col="gray")
    text(x=-(xmax/4)*3,y=-log(l)+0.25,cex=0.5,col="darkgray",labels = l)  
  }
  
  legend("topleft",
         legend=c(colmin/100,((colmin+colmax)/2)/100,colmax/100),
         pch=16,
         col=c("red","black","green"),
         title = "CCC Rho.Q")
  
  dev.off()
  
  png(paste0("volcano/fdrP.",prefix,".png"), width=600, height=600)
  
  plot(      model_performance[ , c] ,
             -log( model_performance[ , fdrValue] ),
             main= paste( prefix , ": FDR-adjusted P-Values" ),
             xlab="Slope of Association",
             ylab=yylab,
             xlim=c(-xmax,xmax),
             col=colcol,
             pch=16
  )
  for(l in c(0.05,0.005,0.0005,0.00005)){
    abline(h=-log(l),col="gray")
    text(x=-(xmax/4)*3,y=-log(l)+0.25,cex=0.5,col="gray",labels = l)  
  }
  
  legend("topleft",
         legend=c(colmin/100,((colmin+colmax)/2)/100,colmax/100),
         pch=16,
         col=c("red","black","green"),
         title = "CCC Rho.Q")
  
  dev.off()
  
  belowCutoff[  pValue   ,     '0.01' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , pValue  ] < 0.01 )
  belowCutoff.fdr[  fdrValue , '0.01' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , fdrValue] < 0.01 )
  
  belowCutoff[  pValue   ,     '0.05' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , pValue  ] < 0.05 )
  belowCutoff.fdr[  fdrValue , '0.05' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , fdrValue] < 0.05 )
  
  belowCutoff[  pValue   ,     '0.001' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , pValue  ] < 0.001 )
  belowCutoff.fdr[  fdrValue , '0.001' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , fdrValue] < 0.001 )
  
  belowCutoff[  pValue   ,     '0.005' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , pValue  ] < 0.005 )
  belowCutoff.fdr[  fdrValue , '0.005' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , fdrValue] < 0.005 )
  
  belowCutoff[  pValue   ,     '0.0001' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , pValue  ] < 0.0001 )
  belowCutoff.fdr[  fdrValue , '0.0001' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , fdrValue] < 0.0001 )
  
  belowCutoff[  pValue   ,     '0.0005' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , pValue  ] < 0.0005 )
  belowCutoff.fdr[  fdrValue , '0.0005' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , fdrValue] < 0.0005 )
  
  belowCutoff[  pValue   ,     '0.00001' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , pValue  ] < 0.00001 )
  belowCutoff.fdr[  fdrValue , '0.00001' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , fdrValue] < 0.00001 )
  
  belowCutoff[  pValue   ,     '0.00005' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , pValue  ] < 0.00005 )
  belowCutoff.fdr[  fdrValue , '0.00005' ] <- sum( model_performance[ !is.na(model_performance[,pValue]) , fdrValue] < 0.00005 )
  #}
}

write.table(belowCutoff, file="belowCutoff.tsv", sep='\t')
write.table(belowCutoff.fdr, file="belowCutoff_FDR.tsv", sep='\t')



#################################################################################################################################################################################################

#################################################################################################################################################################################################

#################################################################################################################################################################################################


v9_v10_OGs_map     <-   read.table("../old_desktop_20200521/v9_v10_OGs_map.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
odb10v1_OG_xrefs     <- read.table("../old_desktop_20200521/odb10v1_OG_xrefs.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
odb10v1_OG2genes     <- read.table("../old_desktop_20200521/odb10v1_OG2genes.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)

# pull in match of indices to defined IDs
odb10v1_genes.2col    <- read.table("../AvianLongevity_ManuscriptData/odb10v1_genes.2col.tsv", sep=" ", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
# 
odb10v1_genes.2col <- odb10v1_genes.2col[ grepl(":",odb10v1_genes.2col[,1]) , ]
odb10v1_genes  <- odb10v1_genes.2col[,2]
names(odb10v1_genes) <- odb10v1_genes.2col[,1]


all_convert <- list()

for( c in cc[grepl("CAPER",cc)][ !grepl('.ord.',cc[grepl("CAPER",cc)]) ] ){
  #for(l in c(0.05,0.005,0.0005,0.00005)){ 
  #0.01,0.001,0.0001,0.00001 
  
  ty <- grepl( 'FDR' , c )
  if( ty == TRUE){
    prefix <- substr(c,1,nchar(c)-10)
    yylab <- '-log( FDR-adjusted(p) )'
  }else{
    prefix <- substr(c,1,nchar(c)-6)
    yylab <- '-log( p )'
  }
  pValue   <- paste0( prefix , '.p_val' )
  fdrValue <- paste0( pValue , ".FDR" )
  
  searchme <- rownames(model_performance)[ model_performance[,fdrValue] < 0.05 ]
  searchme <- searchme[ searchme %in% v9_v10_OGs_map[,2] ]
  
  print(paste(c, length(searchme), length(all_convert)))
  all_convert <- c(all_convert,searchme[!searchme %in% all_convert])
}



searchme <- all_convert
print( paste(c, length(searchme)) )

cmpq <- list()
ri <- 1
for(r in searchme){
  
  if(ri %% 10 == 0){
    print(  paste(c, ri, length(searchme), (ri/length(searchme)*100) )  )
  }
  
  if(ri > 1){
    # cmpq[ri,1] <-  v9_v10_OGs_map[ v9_v10_OGs_map[,2] == r , 3  ]
    tmpzzz <-  unique( odb10v1_OG2genes[ odb10v1_OG2genes[,1] %in%  v9_v10_OGs_map[ v9_v10_OGs_map[,2] == r , 3  ] , 2 ] )
    tmpzzz <-  tmpzzz[ tmpzzz %in% odb10v1_genes.2col[,1] ]
    tmpzzz <-  odb10v1_genes.2col[ odb10v1_genes.2col[,1] %in% tmpzzz , ]
    tmpzzz <-  table( tmpzzz[grepl("LOC",tmpzzz[,2])==FALSE,][,2] )
    cmpq[ri] <- paste(names(tmpzzz),collapse=',')
  }
  ri <- ri + 1
  
}
names(cmpq) <- searchme

# save.image("named_ogs.20200612.Rdata")

#}
#}

load( "named_ogs.20200612.Rdata" )
model_performance <- model_performance[-1,-1]
resid_performance <- resid_performance[-1,-1]

#################################################################################################################################################################################################

#################################################################################################################################################################################################

#################################################################################################################################################################################################


topDiffGenes <- function(allScore) {  return(allScore > 1)   }
#bm <- useMart("ensembl")
#bm <- useDataset("hsapiens_gene_ensembl", mart=bm)
bm <- useEnsembl(biomart = "ensembl", 
                 dataset = "hsapiens_gene_ensembl", 
                 mirror = "useast")


belowCutoff_directional <- as.data.frame(matrix())


signif_report_c <- as.data.frame(matrix())
signif_report_og <- as.data.frame(matrix())
signif_report_gene <- as.data.frame(matrix())
signif_report_topgo <- as.data.frame(matrix())
signif_report_topgo_strings <- as.data.frame(matrix())

for(c in cc[37:96]){
  
  ty <- grepl( 'FDR' , c )
  if( ty == TRUE){
    prefix <- substr(c,1,nchar(c)-10)
    yylab <- '-log( FDR-adjusted(p) )'
  }else{
    prefix <- substr(c,1,nchar(c)-6)
    yylab <- '-log( p )'
  }
  pValue   <- paste0( prefix , '.p_val' )
  fdrValue <- paste0( pValue , ".FDR" )
  xmax <- max(model_performance[ , c][ !is.na(model_performance[ , c])  ])
  
  estrhocols <- colnames(model_performance)[ grepl('est.rho',colnames(model_performance)) ]
  iter <- strsplit(c,'\\.')[[1]][1]
  
  
  model_performance.pos <- rownames(model_performance)[ model_performance[,c] > 0 ]
  model_performance.neg <- rownames(model_performance)[ model_performance[,c] < 0 ]
  
  cutoffs <- c(0.01,0.05,0.001,0.005,0.0001,0.0005,0.00001,0.00005)
  print(paste(nrow(model_performance)," total"))
  for(cutoff in cutoffs){
    model_performance.cut  <- model_performance[ model_performance[,pValue] < cutoff ,  ]
    model_performance.cut_fdr  <- model_performance[ model_performance[,fdrValue] < cutoff ,  ]
    
    model_performance.cut.pos <- model_performance.cut[ model_performance.cut[ , c ] > cutoff ,  ]
    model_performance.cut.neg <- model_performance.cut[ model_performance.cut[ , c ] < cutoff ,  ]
    model_performance.cut_fdr.pos <- model_performance.cut_fdr[ model_performance.cut_fdr[ , c ] > cutoff ,  ]
    model_performance.cut_fdr.neg <- model_performance.cut_fdr[ model_performance.cut_fdr[ , c ] < cutoff ,  ]
    
    signif_report_c[ c , paste0('raw_pos_',cutoff) ] <-  nrow(model_performance.cut.pos)
    signif_report_c[ c , paste0('raw_neg_',cutoff) ] <-  nrow(model_performance.cut.neg)
    signif_report_c[ c , paste0('fdr_pos_',cutoff) ] <-  nrow(model_performance.cut_fdr.pos)
    signif_report_c[ c , paste0('fdr_neg_',cutoff) ] <-  nrow(model_performance.cut_fdr.neg)
    
    signif_report_og[ c , paste0('raw_pos_',cutoff) ] <-  paste( names(cmpq[ rownames(model_performance.cut.pos) ]) , collapse="," )
    signif_report_og[ c , paste0('raw_neg_',cutoff) ] <-  paste( names(cmpq[ rownames(model_performance.cut.neg) ]) , collapse="," )
    signif_report_og[ c , paste0('fdr_pos_',cutoff) ] <-  paste( names(cmpq[ rownames(model_performance.cut_fdr.pos) ]) , collapse="," )
    signif_report_og[ c , paste0('fdr_neg_',cutoff) ] <-  paste( names(cmpq[ rownames(model_performance.cut_fdr.neg) ]) , collapse="," )
    
    u1 <- unique( strsplit( paste( unname(cmpq[ rownames(model_performance.cut.pos) ]) , collapse="," ) , ',' ) )[[1]]
    u1 <- u1[u1 != ""][ u1[u1 != ""] != "NULL" ]
    signif_report_gene[ c , paste0('raw_pos_',cutoff) ] <-  paste( u1 , collapse="," )
    u2 <- unique( strsplit( paste( unname(cmpq[ rownames(model_performance.cut.neg) ]) , collapse="," ) , ',' ) )[[1]]
    u2 <- u2[u2 != ""][ u2[u2 != ""] != "NULL" ]
    signif_report_gene[ c , paste0('raw_neg_',cutoff) ] <-  paste( u2 , collapse="," )
    u3 <- unique( strsplit( paste( unname(cmpq[ rownames(model_performance.cut_fdr.pos) ]) , collapse="," ) , ',' ) )[[1]]
    u3 <- u3[u3 != ""][ u3[u3 != ""] != "NULL" ]
    signif_report_gene[ c , paste0('fdr_pos_',cutoff) ] <-  paste( u3 , collapse="," )
    u4 <- unique( strsplit( paste( unname(cmpq[ rownames(model_performance.cut_fdr.neg) ]) , collapse="," ) , ',' ) )[[1]]
    u4 <- u4[u4 != ""][ u4[u4 != ""] != "NULL" ]
    signif_report_gene[ c , paste0('fdr_neg_',cutoff) ] <-  paste( u4 , collapse="," )
    
    print(  paste(  cutoff , "|raw|" , nrow(model_performance.cut) , 
                    "total" ,  nrow(model_performance.cut.pos) , 
                    "gr" ,  nrow(model_performance.cut.neg) , 
                    "lt     |fdr|" , nrow(model_performance.cut_fdr) , 
                    "total" ,  nrow(model_performance.cut_fdr.pos) , 
                    "gr" ,  nrow(model_performance.cut_fdr.neg) , "lt" ) )
    
    u1 <- u1[!is.na(u1)]
    u2 <- u2[!is.na(u2)]
    u3 <- u3[!is.na(u3)]
    u4 <- u4[!is.na(u4)]
    if(length(u1) > 1){
      results <- NA
      results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'go_id'), 
                       filter = "hgnc_symbol",
                       values = unlist(strsplit(u1,",")),
                       mart = bm)
      
      geneListX <- c(unlist(strsplit(u1,",")),
                     "FAKE")
      geneList <- factor(c(as.integer(rep(1,length(geneListX)-1)),0))
      names(geneList) <- geneListX
      
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
      
      queryFisher.classic.fisher.bp <- queryFisher.classic.ks.bp <- queryFisher.elim.ks.bp <- NA
      queryFisher.classic.fisher.mf <- queryFisher.classic.ks.mf <- queryFisher.elim.ks.mf <- NA
      queryFisher.classic.fisher.cc <- queryFisher.classic.ks.cc <- queryFisher.elim.ks.cc <- NA
      
      tryCatch( queryFisher.classic.fisher.bp <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.bp     <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.bp     <- runTest(queryGOdata.bp, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      tryCatch( queryFisher.classic.fisher.mf <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.mf     <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.mf     <- runTest(queryGOdata.mf, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      tryCatch( queryFisher.classic.fisher.cc <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.cc     <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.cc     <- runTest(queryGOdata.cc, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      allRes.bp <- allRes.mf <- allRes.cc <- NA
      tryCatch( allRes.bp <- GenTable( queryGOdata.bp, classicFisher = queryFisher.classic.fisher.bp, classicKS = queryFisher.classic.ks.bp, elimKS = queryFisher.elim.ks.bp, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.bp)) )   , error = function(e) message(paste("ROBUST FAIL")))        
      tryCatch( allRes.mf <- GenTable( queryGOdata.mf, classicFisher = queryFisher.classic.fisher.mf, classicKS = queryFisher.classic.ks.mf, elimKS = queryFisher.elim.ks.mf, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.mf)) )    , error = function(e) message(paste("ROBUST FAIL")))       
      tryCatch( allRes.cc <- GenTable( queryGOdata.cc, classicFisher = queryFisher.classic.fisher.cc, classicKS = queryFisher.classic.ks.cc, elimKS = queryFisher.elim.ks.cc, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.cc)) )    , error = function(e) message(paste("ROBUST FAIL")))       
      
      if(!is.na( allRes.bp)){
        write.table( allRes.bp , file = paste0("TopGO_BP.raw_pos_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'raw_pos_topGO_bp_ks0.5' ] <-  sum( as.numeric(allRes.bp$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'raw_pos_topGO_bp_ks0.5' ] <-  paste( paste( allRes.bp[ allRes.bp$classicKS < 0.05 , ][,1], allRes.bp[ allRes.bp$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
      
      if(!is.na( allRes.mf)){
        write.table( allRes.mf , file = paste0("TopGO_MF.raw_pos_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'raw_pos_topGO_mf_ks0.5' ] <-  sum( as.numeric(allRes.mf$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'raw_pos_topGO_mf_ks0.5' ] <-  paste( paste( allRes.mf[ allRes.mf$classicKS < 0.05 , ][,1], allRes.mf[ allRes.mf$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
      
      if(!is.na( allRes.cc)){
        write.table( allRes.cc , file = paste0("TopGO_CC.raw_pos_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'raw_pos_topGO_cc_ks0.5' ] <-  sum( as.numeric(allRes.cc$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'raw_pos_topGO_cc_ks0.5' ] <-  paste( paste( allRes.cc[ allRes.cc$classicKS < 0.05 , ][,1], allRes.cc[ allRes.cc$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
    }
    if(length(u2) > 1){
      results <- NA
      results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'go_id'), 
                       filter = "hgnc_symbol",
                       values = unlist(strsplit(u2,",")),
                       mart = bm)
      
      geneListX <- c(unlist(strsplit(u2,",")),
                     "FAKE")
      geneList <- factor(c(as.integer(rep(1,length(geneListX)-1)),0))
      names(geneList) <- geneListX
      
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
      
      queryFisher.classic.fisher.bp <- queryFisher.classic.ks.bp <- queryFisher.elim.ks.bp <- NA
      queryFisher.classic.fisher.mf <- queryFisher.classic.ks.mf <- queryFisher.elim.ks.mf <- NA
      queryFisher.classic.fisher.cc <- queryFisher.classic.ks.cc <- queryFisher.elim.ks.cc <- NA
      
      tryCatch( queryFisher.classic.fisher.bp <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.bp     <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.bp     <- runTest(queryGOdata.bp, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      tryCatch( queryFisher.classic.fisher.mf <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.mf     <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.mf     <- runTest(queryGOdata.mf, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      tryCatch( queryFisher.classic.fisher.cc <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.cc     <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.cc     <- runTest(queryGOdata.cc, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      allRes.bp <- allRes.mf <- allRes.cc <- NA
      tryCatch( allRes.bp <- GenTable( queryGOdata.bp, classicFisher = queryFisher.classic.fisher.bp, classicKS = queryFisher.classic.ks.bp, elimKS = queryFisher.elim.ks.bp, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.bp)) )   , error = function(e) message(paste("ROBUST FAIL")))        
      tryCatch( allRes.mf <- GenTable( queryGOdata.mf, classicFisher = queryFisher.classic.fisher.mf, classicKS = queryFisher.classic.ks.mf, elimKS = queryFisher.elim.ks.mf, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.mf)) )    , error = function(e) message(paste("ROBUST FAIL")))       
      tryCatch( allRes.cc <- GenTable( queryGOdata.cc, classicFisher = queryFisher.classic.fisher.cc, classicKS = queryFisher.classic.ks.cc, elimKS = queryFisher.elim.ks.cc, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.cc)) )    , error = function(e) message(paste("ROBUST FAIL")))       
      
      if(!is.na( allRes.bp)){
        write.table( allRes.bp , file = paste0("TopGO_BP.raw_neg_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'raw_neg_topGO_bp_ks0.5' ] <-  sum( as.numeric(allRes.bp$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'raw_neg_topGO_bp_ks0.5' ] <-  paste( paste( allRes.bp[ allRes.bp$classicKS < 0.05 , ][,1], allRes.bp[ allRes.bp$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
      
      if(!is.na( allRes.mf)){
        write.table( allRes.mf , file = paste0("TopGO_MF.raw_neg_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'raw_neg_topGO_mf_ks0.5' ] <-  sum( as.numeric(allRes.mf$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'raw_neg_topGO_mf_ks0.5' ] <-  paste( paste( allRes.mf[ allRes.mf$classicKS < 0.05 , ][,1], allRes.mf[ allRes.mf$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
      
      if(!is.na( allRes.cc)){
        write.table( allRes.cc , file = paste0("TopGO_CC.raw_neg_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'raw_neg_topGO_cc_ks0.5' ] <-  sum( as.numeric(allRes.cc$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'raw_neg_topGO_cc_ks0.5' ] <-  paste( paste( allRes.cc[ allRes.cc$classicKS < 0.05 , ][,1], allRes.cc[ allRes.cc$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
    }
    if(length(u3) > 1){
      results <- NA
      results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'go_id'), 
                       filter = "hgnc_symbol",
                       values = unlist(strsplit(u3,",")),
                       mart = bm)
      
      geneListX <- c(unlist(strsplit(u3,",")),
                     "FAKE")
      geneList <- factor(c(as.integer(rep(1,length(geneListX)-1)),0))
      names(geneList) <- geneListX
      
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
      
      queryFisher.classic.fisher.bp <- queryFisher.classic.ks.bp <- queryFisher.elim.ks.bp <- NA
      queryFisher.classic.fisher.mf <- queryFisher.classic.ks.mf <- queryFisher.elim.ks.mf <- NA
      queryFisher.classic.fisher.cc <- queryFisher.classic.ks.cc <- queryFisher.elim.ks.cc <- NA
      
      tryCatch( queryFisher.classic.fisher.bp <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.bp     <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.bp     <- runTest(queryGOdata.bp, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      tryCatch( queryFisher.classic.fisher.mf <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.mf     <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.mf     <- runTest(queryGOdata.mf, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      tryCatch( queryFisher.classic.fisher.cc <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.cc     <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.cc     <- runTest(queryGOdata.cc, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      allRes.bp <- allRes.mf <- allRes.cc <- NA
      tryCatch( allRes.bp <- GenTable( queryGOdata.bp, classicFisher = queryFisher.classic.fisher.bp, classicKS = queryFisher.classic.ks.bp, elimKS = queryFisher.elim.ks.bp, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.bp)) )   , error = function(e) message(paste("ROBUST FAIL")))        
      tryCatch( allRes.mf <- GenTable( queryGOdata.mf, classicFisher = queryFisher.classic.fisher.mf, classicKS = queryFisher.classic.ks.mf, elimKS = queryFisher.elim.ks.mf, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.mf)) )    , error = function(e) message(paste("ROBUST FAIL")))       
      tryCatch( allRes.cc <- GenTable( queryGOdata.cc, classicFisher = queryFisher.classic.fisher.cc, classicKS = queryFisher.classic.ks.cc, elimKS = queryFisher.elim.ks.cc, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.cc)) )    , error = function(e) message(paste("ROBUST FAIL")))       
      
      if(!is.na( allRes.bp)){
        write.table( allRes.bp , file = paste0("TopGO_BP.fdr_pos_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'fdr_pos_topGO_bp_ks0.5' ] <-  sum( as.numeric(allRes.bp$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'fdr_pos_topGO_bp_ks0.5' ] <-  paste( paste( allRes.bp[ allRes.bp$classicKS < 0.05 , ][,1], allRes.bp[ allRes.bp$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
      
      if(!is.na( allRes.mf)){
        write.table( allRes.mf , file = paste0("TopGO_MF.fdr_pos_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'fdr_pos_topGO_mf_ks0.5' ] <-  sum( as.numeric(allRes.mf$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'fdr_pos_topGO_mf_ks0.5' ] <-  paste( paste( allRes.mf[ allRes.mf$classicKS < 0.05 , ][,1], allRes.mf[ allRes.mf$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
      
      if(!is.na( allRes.cc)){
        write.table( allRes.cc , file = paste0("TopGO_CC.fdr_pos_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'fdr_pos_topGO_cc_ks0.5' ] <-  sum( as.numeric(allRes.cc$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'fdr_pos_topGO_cc_ks0.5' ] <-  paste( paste( allRes.cc[ allRes.cc$classicKS < 0.05 , ][,1], allRes.cc[ allRes.cc$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
    }
    if(length(u4) > 1){
      results <- NA
      results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'go_id'), 
                       filter = "hgnc_symbol",
                       values = unlist(strsplit(u4,",")),
                       mart = bm)
      
      geneListX <- c(unlist(strsplit(u4,",")),
                     "FAKE")
      geneList <- factor(c(as.integer(rep(1,length(geneListX)-1)),0))
      names(geneList) <- geneListX
      
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
      
      queryFisher.classic.fisher.bp <- queryFisher.classic.ks.bp <- queryFisher.elim.ks.bp <- NA
      queryFisher.classic.fisher.mf <- queryFisher.classic.ks.mf <- queryFisher.elim.ks.mf <- NA
      queryFisher.classic.fisher.cc <- queryFisher.classic.ks.cc <- queryFisher.elim.ks.cc <- NA
      
      tryCatch( queryFisher.classic.fisher.bp <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.bp     <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.bp     <- runTest(queryGOdata.bp, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      tryCatch( queryFisher.classic.fisher.mf <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.mf     <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.mf     <- runTest(queryGOdata.mf, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      tryCatch( queryFisher.classic.fisher.cc <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "fisher") , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.classic.ks.cc     <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      tryCatch( queryFisher.elim.ks.cc     <- runTest(queryGOdata.cc, algorithm = "elim",    statistic = "ks"    ) , error = function(e) message(paste("ROBUST FAIL")))
      
      allRes.bp <- allRes.mf <- allRes.cc <- NA
      tryCatch( allRes.bp <- GenTable( queryGOdata.bp, classicFisher = queryFisher.classic.fisher.bp, classicKS = queryFisher.classic.ks.bp, elimKS = queryFisher.elim.ks.bp, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.bp)) )   , error = function(e) message(paste("ROBUST FAIL")))        
      tryCatch( allRes.mf <- GenTable( queryGOdata.mf, classicFisher = queryFisher.classic.fisher.mf, classicKS = queryFisher.classic.ks.mf, elimKS = queryFisher.elim.ks.mf, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.mf)) )    , error = function(e) message(paste("ROBUST FAIL")))       
      tryCatch( allRes.cc <- GenTable( queryGOdata.cc, classicFisher = queryFisher.classic.fisher.cc, classicKS = queryFisher.classic.ks.cc, elimKS = queryFisher.elim.ks.cc, orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.cc)) )    , error = function(e) message(paste("ROBUST FAIL")))       
      
      if(!is.na( allRes.bp)){
        write.table( allRes.bp , file = paste0("TopGO_BP.fdr_neg_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'fdr_neg_topGO_bp_ks0.5' ] <-  sum( as.numeric(allRes.bp$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'fdr_neg_topGO_bp_ks0.5' ] <-  paste( paste( allRes.bp[ allRes.bp$classicKS < 0.05 , ][,1], allRes.bp[ allRes.bp$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
      
      if(!is.na( allRes.mf)){
        write.table( allRes.mf , file = paste0("TopGO_MF.raw_pos_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'raw_pos_topGO_mf_ks0.5' ] <-  sum( as.numeric(allRes.mf$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'raw_pos_topGO_mf_ks0.5' ] <-  paste( paste( allRes.mf[ allRes.mf$classicKS < 0.05 , ][,1], allRes.mf[ allRes.mf$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
      
      if(!is.na( allRes.cc)){
        write.table( allRes.cc , file = paste0("TopGO_CC.raw_pos_",c,".",cutoff,".tsv") , sep = '\t', row.names = FALSE )
        signif_report_topgo[ c , 'raw_pos_topGO_cc_ks0.5' ] <-  sum( as.numeric(allRes.cc$classicKS) < 0.05 )
        signif_report_topgo_strings[ c , 'raw_pos_topGO_cc_ks0.5' ] <-  paste( paste( allRes.cc[ allRes.cc$classicKS < 0.05 , ][,1], allRes.cc[ allRes.cc$classicKS < 0.05 , ][,2] ) , collapse = "," )
      }
    }
    
    
  }
  
}

#belowCutoff_directional.fdr <- as.data.frame(matrix())

signif_report_c <- as.data.frame(matrix())
signif_report_og <- as.data.frame(matrix())
signif_report_gene <- as.data.frame(matrix())
signif_report_topgo <- as.data.frame(matrix())
signif_report_topgo_strings <- as.data.frame(matrix())

write.table(  signif_report_c , file=paste0("signif_report/signif_report_belowCutoff_directional.",cutoff,".tsv")  )
write.table(  signif_report_og , file=paste0("signif_report/signif_report_OGnamesBelowCutoff_directional.",cutoff,".tsv")  )
write.table(  signif_report_gene , file=paste0("signif_report/signif_report_GeneNamesBelowCutoff_directional.",cutoff,".tsv")  )

write.table(  signif_report_topgo , file=paste0("signif_report/signif_report_TopGOBelowCutoff_directional.",cutoff,".tsv")  )
write.table(  signif_report_topgo_strings , file=paste0("signif_report/signif_report_TopGOStringsBelowCutoff_directional.",cutoff,".tsv")  )



# save.image("deveL-complete-202200613.Rdata")

#################################################################################################################################################################################################

#################################################################################################################################################################################################

#################################################################################################################################################################################################

load("deveL-complete-202200613.Rdata")



