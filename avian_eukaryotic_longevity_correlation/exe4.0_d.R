#setwd ("/Users/jmccorri/Desktop/exe.Avian")

library(ggplot2)
library(MASS)
library(vegan)
library(ggplot2)
library(Hmisc)
library(gplots)
library(ape)
library(caper)
setwd ("~/Desktop/develop.Avian")

######################################################################################################################################################################################################################################################
# og tables - high confidence
og_counts.high_conf           <- as.data.frame(read.table("backup20181230/OG.high_conf.abundance.tsv",     sep="\t", header=TRUE))
og_counts.high_conf           <- og_counts.high_conf[1:dim(og_counts.high_conf)[1]-1,]
rownames(og_counts.high_conf) <- as.character(og_counts.high_conf[,1])
og_counts.high_conf           <- og_counts.high_conf[,-1]
og_counts.high_conf           <- og_counts.high_conf[,-1]
og_counts.high_conf[is.na(og_counts.high_conf)] <- 0
og_counts.high_conf.abs       <- og_counts.high_conf[ , as.character(mls$AbundanceOrder) ]
remove(og_counts.high_conf)
og_counts.high_conf.rel <- t(decostand( t(og_counts.high_conf.abs), method="total" ))

shared_abundances.high_conf <- rowSums( og_counts.high_conf.abs == 0 )
png("denovo_og_shared_abundance_histogram.high_conf2.png", height=500,width=500)
barplot( unlist(table(shared_abundances.high_conf)) )
dev.off()  
write.table( table(shared_abundances.high_conf) , file = "denovo_og_shared_abundance_histogram.high_conf.tsv", sep = "\t", row.names = FALSE )

######################################################################################################################################################################################################################################################
# og tables - low confidence
og_counts.low_conf  <- as.data.frame(read.table("backup20181230/OG.abundance.tsv",     sep="\t", header=TRUE))
og_counts.low_conf           <- og_counts.low_conf[1:dim(og_counts.low_conf)[1]-1,]
rownames(og_counts.low_conf) <- as.character(og_counts.low_conf[,1])
og_counts.low_conf           <- og_counts.low_conf[,-1]
og_counts.low_conf           <- og_counts.low_conf[,-1]
og_counts.low_conf[is.na(og_counts.low_conf)] <- 0
og_counts.low_conf.abs       <- og_counts.low_conf[ , as.character(mls$AbundanceOrder) ]
remove(og_counts.low_conf)
og_counts.low_conf.rel <- t(decostand( t(og_counts.low_conf.abs), method="total" ))

shared_abundances.low_conf <- rowSums( og_counts.low_conf.abs == 0 )
png("denovo_og_shared_abundance_histogram.low_conf.png", height=500,width=500)
barplot( unlist(table(shared_abundances.low_conf)) )
dev.off()  
write.table( table(shared_abundances.low_conf) , file = "denovo_og_shared_abundance_histogram.low_conf.tsv", sep = "\t", row.names = FALSE )

######################################################################################################################################################################################################################################################
# import og ids
og_ids     <- read.table("backup20181230/OG.public_ids.txt", sep="\t", header=TRUE)
rownames(og_ids) <- og_ids[,1]

######################################################################################################################################################################################################################################################
# tree
phylotree <- read.tree(text="(Ostrich:113,((Ruby_throated_Hummingbird:79.6,(((Double_crested_Cormorant:71.5,Sandhill_Crane:71.5)AT:3,(Rock_Dove:26.3,Mourning_Dove:26.3)AU:48.2)AV:4,((Killdeer:67.9,((Woodcock:7.6,Spotted_Sandpiper:7.6)AP:55.8,(Caspian_Tern:22.5,(Ring_billed_Gull:1.58,Herring_Gull:1.58)AN:20.92)AO:40.9)AQ:4.5)AR:11.5,((Great_Horned_Owl:77.7,((Red_tailed_Hawk:24.2,Coopers_Hawk:24.2)AJ:51.7,(Downy_Woodpecker:14.3,(Red_bellied_Woodpecker:7.05,Hairy_Woodpecker:7.05)AH:7.15)AI:61.6)AK:1.8)AL:0.7,((American_Crow:37.8,Red_eyed_Vireo:37.8)AF:12.3,((Tufted_Titmouse:42.2,(Horned_Lark:34.9,(Tree_Swallow:18.8,Barn_Swallow:18.8)AB:16.1)AC:7.3)AD:3.4,((House_Sparrow:29.3,(House_Finch:26,((Northern_Cardinal:21.5,(Brown_headed_Cowbird:10.5,(Yellow_Warbler:5.25,(Yellow_throated_Warbler:4.09,Yellow_rumped_Warbler:4.09)N:1.06)O:5.25)P:11)Q:1.4,(Common_Grackle:12,(Chipping_Sparrow:9.97,(Song_Sparrow:8.46,American_Tree_Sparrow:8.46)K:1.51)L:2.03)M:10.9)R:3.1)S:3.3)T:15.8,(Cedar_waxwing:44.6,((White_breasted_Nuthatch:32.2,(Carolina_Wren:14.3,House_wren:14.3)U:17.9)V:6.6,(American_Robin:28.1,(Gray_Catbird:21.9,Starling:21.9)W:6.2)X:10.7)Y:5.8)Z:0.5)AA:0.5)AE:4.5)AG:28.3)AM:1)AS:0.1)AW:0.1)AX:18,((Pheasant:14.6,(Turkey:13.3,Ruffed_Grouse:13.3)H:1.3)I:57.8,((Canada_Goose:13.5,Mute_Swan:13.5)F:10.6,(Wood_Duck:12.3,(Northern_Shoveler:7.34,(Gadwall:5.77,(Mallard:4.24,(Green_winged_Teal:3.99,Pintail_Duck:3.99)A:1.25)B:1.53)C:1.6)D:4.93)E:11.8)G:48.3)J:25.2)AY:15.4)AZ;")
ref.ordered <- c("Taeniopygia_guttata","Taeniopygia_guttata","Taeniopygia_guttata","Taeniopygia_guttata","Taeniopygia_guttata","Calypte_anna","Taeniopygia_guttata","Taeniopygia_guttata","Corvus_brachyrhynchos","Taeniopygia_guttata","Charadrius_vociferus","Taeniopygia_guttata","Gallus_gallus","Taeniopygia_guttata","Taeniopygia_guttata","Taeniopygia_guttata","Picoides_pubescens","Charadrius_vociferus","Taeniopygia_guttata","Meleagris_gallopavo","Taeniopygia_guttata","Picoides_pubescens","Taeniopygia_guttata","Taeniopygia_guttata","Taeniopygia_guttata","Taeniopygia_guttata","Corvus_brachyrhynchos","Falcons_peregrine","Anas_platyrhynchos","Picoides_pubescens","Charadrius_vociferus","Anas_platyrhynchos","Anas_platyrhynchos","Nipponia_nippon","Taeniopygia_guttata","Taeniopygia_guttata","Taeniopygia_guttata","Gallus_gallus","Anas_platyrhynchos","Anas_platyrhynchos","Taeniopygia_guttata","Picoides_pubescens","Anas_platyrhynchos","Charadrius_vociferus","Falcons_peregrine","Charadrius_vociferus","Charadrius_vociferus","Columba_livia","Charadrius_vociferus","Anas_platyrhynchos","Charadrius_vociferus","Struthio_camelus","Anas_platyrhynchos")
ref.query_order <- c("Yellow_throated_Warbler",
                     "Yellow_rumped_Warbler",
                     "Horned_Lark","Cedar_waxwing",
                     "House_wren",
                     "Ruby_throated_Hummingbird","Carolina_Wren","White_breasted_Nuthatch",
                     "Red_eyed_Vireo","American_Tree_Sparrow","Killdeer","Yellow_Warbler","Ruffed_Grouse","Song_Sparrow","House_Finch","Chipping_Sparrow","Downy_Woodpecker","Spotted_Sandpiper",
                     "Tree_Swallow","Turkey","Tufted_Titmouse","Hairy_Woodpecker","Barn_Swallow","Brown_headed_Cowbird","American_Robin","Gray_Catbird","American_Crow","Coopers_Hawk","Northern_Shoveler",
                     "Red_bellied_Woodpecker","Woodcock","Gadwall","Wood_Duck","Double_crested_Cormorant","Common_Grackle",
                     "Starling",
                     "House_Sparrow",
                     "Pheasant",
                     "Green_winged_Teal",
                     "Pintail_Duck",
                     "Northern_Cardinal","Great_Horned_Owl","Mallard","Caspian_Tern","Red_tailed_Hawk","Mourning_Dove","Ring_billed_Gull","Rock_Dove","Sandhill_Crane","Canada_Goose","Herring_Gull",
                     "Ostrich","Mute_Swan") 

cl <- list()
cli <- 1
spectrum <- c("darkred","red","darkorange","darkgoldenrod1","yellowgreen","darkgreen","aquamarine3","blue","darkblue","blueviolet","violet","black")
names(spectrum) <- unique(ref.ordered)
for( w in phylotree$tip.label ){
  ro <- ref.ordered[ref.query_order == w]
  cl[cli] <- spectrum[ro]
  cli <- cli + 1
}
plot(phylotree,
     type = "phylogram", main="Query Bird Sample Phylogeny (A. Pickering)",
     tip.color = unlist(cl),
     cex = 1.4, x.lim = c(0,140)
)
# clean all metadata to match 50
phylotree       <- drop.tip(phylotree, c("Ostrich","Yellow_throated_Warbler","Pheasant"))
ref.ordered     <- ref.ordered    [ ! ref.query_order %in% c("Ostrich","Pheasant","Yellow_throated_Warbler") ]
ref.query_order <- ref.query_order[ ! ref.query_order %in% c("Ostrich","Pheasant","Yellow_throated_Warbler") ]
cl <- list()
cli <- 1
spectrum <- c("darkred","red","darkorange","darkgoldenrod1","yellowgreen","darkgreen","aquamarine3","blue","darkblue","blueviolet","violet","black")
names(spectrum) <- unique(ref.ordered)
for( w in phylotree$tip.label ){
  ro <- ref.ordered[ref.query_order == w]
  cl[cli] <- spectrum[ro]
  cli <- cli + 1
}
cl <- unlist(cl)
cl[5] <- "violet"
plot(phylotree,
     type = "phylogram", main="Query Bird Sample Phylogeny (A. Pickering)",
     tip.color = unlist(cl),
     cex = 0.7, x.lim = c(0,140)
)

# calculate residuals
f <- as.formula("A ~ C")
compare_me     <- as.data.frame(matrix(ncol=4,nrow=50))
compare_me[,1] <- mls$MLS.Years
compare_me[,2] <- rep(0,50)
compare_me[,3] <- log(mls$Body.Mass.Grams)
compare_me[,4] <-  mls$TreeFormat
colnames(compare_me) <- c("A","B","C","D")
compdata <- comparative.data(phylotree, compare_me, D)
reg_LM_t_LOGw <- crunch(f, data=compdata)
pvals <- list()
pp <- 1
if (t(coef(summary(reg_LM_t_LOGw)))[2,1] > 0 & is.nan(t(coef(summary(reg_LM_t_LOGw)))[2,1]) == 0){
  cpus.lm2 <- as.matrix( anova(reg_LM_t_LOGw) )
  fstat <- summary(reg_LM_t_LOGw)$fstatistic
  pvals <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
}else{
  pvals[pp] <- NA
}

q1 <- compare_me$A
q2 <- compare_me$C
plot(q1,q2,xlim = c(-5,45))
abline(lm(q2 ~ q1),col="red")
text(q1,q2,labels = compare_me$D,pos=2,cex=0.5)
lmlm  <- lm(q2 ~ q1)
lmres <- lmlm$residuals
names(lmres) <- compare_me$D 

png("avian50.MLS_residual.png", width = 600, height = 800)
par(mar=c(8,15,2,2))
barplot( sort(lmres), horiz=TRUE, las=2, cex.names = 0.7 )
dev.off()

png("avian50.longevity_lm.png", width = 600, height = 600)
q1 <- compare_me$A
q2 <- compare_me$C
fit <- lm(q1 ~ q2)
sd2 <- sd(abs(fit$residuals))*2
sd1 <- sd(abs(fit$residuals))
update_lab <- as.character(compare_me$D)
resi <- resid(fit) 
update_lab[ abs(resi) > sd1 & abs(resi) < sd2 ] <- compare_me$D[ abs(resi) > sd1 & abs(resi) < sd2 ]
update_lab[ abs(resi) > sd2 ]                   <- compare_me$D[ abs(resi) > sd2 ]
update_lab[ abs(resi) < sd1 ]                   <- ""
text_col = rep("black", 50)
text_col[ abs(resi) > sd1 & abs(resi) < sd2 ] <- "blue"
text_col[ abs(resi) > sd2 ]                   <- "black"
text_col[ abs(resi) < sd1 ]                   <- "gray"
plot(q2,q1,xlim = c(0,10), col=text_col,pch=16,xlab = "log(Mass) (Grams)",ylab="Maximum Lifespan (Years)" )
abline(fit$coefficients[1],fit$coefficients[2])
abline(fit$coefficients[1]+sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]-sd2,fit$coefficients[2], col="gray", lty=4)
abline(fit$coefficients[1]+sd1,fit$coefficients[2], col="lightblue", lty=4)
abline(fit$coefficients[1]-sd1,fit$coefficients[2], col="lightblue", lty=4)
update_lab <- as.character(compare_me$D)
pos_lab <- rep(4, 50)
pos_lab[  update_lab %in% c("Herring_Gull","Sandhill_Crane","Mute_Swan","House_Sparrow","Canada_Goose","Gray_Catbird","Ring_billed_Gull","Red_tailed_Hawk","Caspian_Tern") ] <- 2
pos_lab[  update_lab %in% c("Ruby_throated_Hummingbird","Woodcock","Common_Grackle") ] <- 3
pos_lab[  update_lab %in% c("Starling","Killdeer") ] <- 1
update_lab[update_lab == "Ruby_throated_Hummingbird"] <- "Ruby             \nthroated        \nHummingbird"
update_lab[ abs(resi) < sd1  ] <- ""
update_lab[ compare_me$D == "Woodcock" ] <- "                Woodcock"
update_lab <- gsub("_", " ", update_lab)
text(q2,q1,labels = update_lab,pos=pos_lab,cex=0.7, col=text_col)
dev.off() 

extremeLongevity_sd2_all <- compare_me$D[ abs(resi) > sd2  ]
extremeLongevity_sd2_hi  <- compare_me$D[ resi > sd2  ]
extremeLongevity_sd2_lo  <- compare_me$D[ resi < -sd2  ]

extremeLongevity_sd1_all <- compare_me$D[ abs(resi) > sd1 & abs(resi) < sd2 ]
extremeLongevity_sd1_hi  <- compare_me$D[ resi > sd1 & resi < sd2 ]
extremeLongevity_sd1_lo  <- compare_me$D[ resi < -sd1 & resi > -sd2 ]

extremeLongevity_sd1andsd2_all <- compare_me$D[ abs(resi) > sd1  ]
extremeLongevity_sd1andsd2_hi  <- compare_me$D[ resi > sd1  ]
extremeLongevity_sd1andsd2_lo  <- compare_me$D[ resi < -sd1  ]


######################################################################################################################################################################################################################################################
# DEVELOP PLOT # 
# REPLACE WITH VERSION FROM OTHER LAPTOP  ASAP!

# png(filename="avian50.MLS_vs_logWeight_Phylo.png", width = 1292, height = 1130 )
png(filename="avian50.MLS_vs_logWeight_Phylo.resize.png", width = 1682, height = 1300 )

plot.new()
m <- rbind(c(1:6))
layout(m,c(0.60,0.08,0.08,0.08,0.08,0.08))
plot(phylotree,
     type = "phylogram", main="Avian Phylogeny",
     tip.color = unlist(cl),
     cex = 2, #cex = 1.4, 
     cex.main=2,
     x.lim = c(0,120)
)
nodelabels(text=phylotree$node.label,cex=1.2,bg="white",frame="circle",
           col = rev(redgreen(100))[round(as.numeric(unlist(abs(reg_LM_t_LOGw$mod$residuals)/max(abs(reg_LM_t_LOGw$mod$residuals)) ))*100)],
)
spectrum <- spectrum[1:10]
axis(side = 1, at = seq(0, 140, 20), labels = FALSE, lwd = 2)
mtext(seq(0, 140, 20), side = 1, at = seq(0, 140, 20), line = 1, las = 2)
clade.matrix(phylotree)

# names(spectrum)
labelSpect <- c("Zebra Finch","Anna's Hummingbird","American Crow","Killdeer","Chicken","Downy Woodpecker","Turkey","Peregrine Falcon","Mallard","Crested Ibis"  )

legend("topleft", labelSpect, pch=22,pt.bg=c(spectrum,y.intersp=1.3),col="white", cex = 1.7, pt.cex=3.4, title = "Best Reference", bty = "n",lwd=2)
legend("topleft", labelSpect, pch=22,pt.bg=c(spectrum,y.intersp=1.3),col="white", cex = 1.7, pt.cex=3.4, title = "Best Reference", bty = "n",lwd=2)

par(mar=c(5.1,0.2,4.1,0.2))
barplot( mls[ match(phylotree$tip.label, mls$TreeFormat) , "MLS.Years" ], horiz = TRUE, col=unlist(cl), main="\nMLS\nYears", cex.main=2 )

q1 <- compare_me$A
q2 <- mls$Body.Mass.Grams
nonloglmres  <- resid(lm(q2 ~ q1))
barplot( mls[ match(phylotree$tip.label, mls$TreeFormat) , "Body.Mass.Grams" ], horiz = TRUE, col=unlist(cl), main="\nBody\nMass\nGrams", cex.main=2 )
barplot( nonloglmres, horiz=TRUE, las=2, main="\nMLS vs.\nBMG\nResidual", col=unlist(cl), names.arg = FALSE, cex.main=2 )

barplot( log( mls[ match(phylotree$tip.label, mls$TreeFormat) , "Body.Mass.Grams" ] ), horiz = TRUE, col=unlist(cl), main="log(BMG)", cex.main=2 )
barplot( lmres, horiz=TRUE, las=2, main="\nMLS vs.\nlog(BMG)\nResidual", col=unlist(cl), names.arg = FALSE, cex.main=2 )

dev.off()


# png(filename="avian50.MLS_vs_logWeight_Phylo.png", width = 1292, height = 1130 )
png(filename="avian50.MLS_vs_logWeight_Phylo.4print.png", width = 1682, height = 2176 )

plot.new()
m <- rbind(c(1:6))
layout(m,c(0.675,0.065,0.065,0.065,0.065,0.065))
plot(phylotree,
     type = "phylogram", main="Avian Phylogeny",
     tip.color = unlist(cl),
     cex = 3, #cex = 1.4, 
     cex.main=3, 
     x.lim = c(0,140)
)
nodelabels(text=phylotree$node.label,cex=2,bg="white",frame="circle",
           col = rev(redgreen(100))[round(as.numeric(unlist(abs(reg_LM_t_LOGw$mod$residuals)/max(abs(reg_LM_t_LOGw$mod$residuals)) ))*100)],
)
spectrum <- spectrum[1:10]
axis(side = 1, at = seq(0, 140, 20), labels = FALSE, lwd = 2)
mtext(seq(0, 140, 20), side = 1, at = seq(0, 140, 20), line = 1, las = 2)
clade.matrix(phylotree)

# names(spectrum)
labelSpect <- c("Zebra Finch","Anna's Hummingbird","American Crow","Killdeer","Chicken","Downy Woodpecker","Turkey","Peregrine Falcon","Mallard","Crested Ibis"  )

legend("topleft", labelSpect, pch=22,pt.bg=c(spectrum,y.intersp=1.3),col="white", cex = 2.5, pt.cex=5, title = "Best Reference", bty = "n",lwd=2)
legend("topleft", labelSpect, pch=22,pt.bg=c(spectrum,y.intersp=1.3),col="white", cex = 2.5, pt.cex=5, title = "Best Reference", bty = "n",lwd=2)

par(mar=c(5.1,0.2,4.1,0.2))
barplot( mls[ match(phylotree$tip.label, mls$TreeFormat) , "MLS.Years" ], horiz = TRUE, col=unlist(cl), main="\n\nMLS\nYears", cex.main=2.5 )

q1 <- compare_me$A
q2 <- mls$Body.Mass.Grams
nonloglmres  <- resid(lm(q2 ~ q1))
barplot( mls[ match(phylotree$tip.label, mls$TreeFormat) , "Body.Mass.Grams" ], horiz = TRUE, col=unlist(cl), main="\n\nBody\nMass\nGrams", cex.main=2.5 )
barplot( nonloglmres, horiz=TRUE, las=2, main="\n\nMLS vs.\nBMG\nResidual", col=unlist(cl), names.arg = FALSE, cex.main=2.5 )

barplot( log( mls[ match(phylotree$tip.label, mls$TreeFormat) , "Body.Mass.Grams" ] ), horiz = TRUE, col=unlist(cl), main="\nlog(BMG)", cex.main=2.5 )
barplot( lmres, horiz=TRUE, las=2, main="\n\nMLS vs.\nlog(BMG)\nResidual", col=unlist(cl), names.arg = FALSE, cex.main=2.5 )

dev.off()

# load("20190102/tmp20190101.Rdata")
# ^ to gather all above

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

# fix ag_count using RubyThroatedHummingbird

# load("20181210_10am.Rdata")

column_shared <- list()
c <- 1
for(a in colnames(ag_count)){
  column_shared[c] <- sum(ag_count[,a] > 0) #sum(ag_count.pop[,a] > 0)
  c <- c + 1
}

png(filename = "dev.png", height = 800,width = 800)
barplot(table(unlist(column_shared)))
dev.off()
 
rpkm.hum <- read.table("backup20181230/allWITHhummingbird.rpkm", sep=" ", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
rpkm.hum <- rpkm.hum[ rpkm.hum$V3 == "RubyThroatedHummingbird", ]

clusters_with_hummingbird <- list()
ci <- 1
ri <- 1
for(r in 1:dim(ref11_mcl)[1]){
  if(ci %% 200){ print(ci/dim(ref11_mcl)[1]) }
  if( sum(grepl("Aan", as.character(ref11_mcl[ 1, ]))) > 0 ){
    clusters_with_hummingbird[ri] <- r
    ri <- ri + 1
  }
  ci <- ci + 1
}
 
ag_count.hum <- as.data.frame(ag_count)
ag_count.hum[ "Ruby_throated_Hummingbird" , ] <- rep(0, dim(ag_count.hum)[2] )

ci <- 1
for( r in unlist(clusters_with_hummingbird) ){
     #ref11_mcl[ unlist(clusters_with_hummingbird) , ] ){
  if(ci %% 200){ print(ci/length(unlist(clusters_with_hummingbird))) }
  
  cont     <- ref11_mcl[ r , ][ref11_mcl[ r , ] != ""]
  cont.hum <- cont[grepl("Aan_", cont)]
  ag_count.hum[ "Ruby_throated_Hummingbird" , r ] <- sum( as.numeric(rpkm.hum[ rpkm.hum$V1 %in% cont.hum , ]$V2) )
  ci <- ci + 1
  
}


ag_count.hum <- ag_count.hum[ rownames(ag_count.hum) != "YellowThroatedWarbler" , ]

ag_count.hum.rel <- decostand( ag_count.hum, method="total" )

ag_count.hum     <- t(ag_count.hum)
ag_count.hum.rel <- t(ag_count.hum.rel)

ag_count.hum     <- ag_count.hum    [ , mls$AbundanceOrder]
ag_count.hum.rel <- ag_count.hum.rel[ , mls$AbundanceOrder]

column_shared.hum <- list()
c <- 1
for(a in rownames(ag_count.hum)){
  column_shared.hum[c] <- 50-sum(ag_count.hum[a,] > 0) #sum(ag_count.pop[,a] > 0)
  c <- c + 1
}

png(filename = "dev.png", height = 800,width = 800)
barplot(table(unlist(column_shared.hum)))
dev.off()

# save.image(file=paste0("20190104_2pm.Rdata"))

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

load("20190102/tmp20190101.Rdata")
# ^ reload the new environment

for(runtype in c("denovo_highConfidence","alignment","denovo_lowConfidence")){ 
  
  if(runtype=="alignment"){             shared_abundances    <- unlist(column_shared.hum)
                                        abundance_matrix.abs <- ag_count.hum[ shared_abundances <= 46, ]
                                        abundance_matrix.rel <- ag_count.hum.rel[ shared_abundances <= 46, ]
                                        shared_abundances    <- shared_abundances[ shared_abundances <= 46 ]
                                        rownames( abundance_matrix.abs ) <- rownames( abundance_matrix.rel ) <- paste0( "orthoMCL_",rownames(abundance_matrix.abs) )
                                         } 
  if(runtype=="denovo_highConfidence"){ abundance_matrix.abs <- og_counts.high_conf.abs 
                                        abundance_matrix.rel <- og_counts.high_conf.rel
                                        shared_abundances    <- shared_abundances.high_conf } 
  if(runtype=="denovo_lowConfidence"){  abundance_matrix.abs <- og_counts.low_conf.abs 
                                        abundance_matrix.rel <- og_counts.low_conf.rel 
                                        shared_abundances    <- shared_abundances.low_conf }
  colnames(abundance_matrix.abs) <- colnames(abundance_matrix.rel) <- mls$AbundanceOrder
  
  
    #for(modelx in c(1:39)){
    #ignore factor to factor comparison for now 
    
    for(modelx in c(1:27)){ #,37:38)){
      
      # Caper
      if(modelx == 1){  modelid <- "CAPER.logwMLSresid_v_relAbs" }
      if(modelx == 2){  modelid <- "CAPER.logwMLSresid_v_absAbs" }
      if(modelx == 3){  modelid <- "CAPER.logwMLSresid_v_binaryAbs" }
      if(modelx == 4){  modelid <- "CAPER.wMLSresid_v_relAbs" }
      if(modelx == 5){  modelid <- "CAPER.wMLSresid_v_absAbs" }
      if(modelx == 6){  modelid <- "CAPER.wMLSresid_v_binaryAbs" }
      if(modelx == 7){  modelid <- "CAPER.MLS_v_relAbs" }
      if(modelx == 8){  modelid <- "CAPER.MLS_v_absAbs" }
      if(modelx == 9){  modelid <- "CAPER.MLS_v_binaryAbs" }
      # LM
      if(modelx == 10){ modelid <- "LM.logwMLSresid_v_relAbs" }
      if(modelx == 11){ modelid <- "LM.logwMLSresid_v_absAbs" }
      if(modelx == 12){ modelid <- "LM.logwMLSresid_v_binaryAbs" }
      if(modelx == 13){ modelid <- "LM.wMLSresid_v_relAbs" }
      if(modelx == 14){ modelid <- "LM.wMLSresid_v_absAbs" }
      if(modelx == 15){ modelid <- "LM.wMLSresid_v_binaryAbs" }
      if(modelx == 16){ modelid <- "LM.MLS_v_relAbs" }
      if(modelx == 17){ modelid <- "LM.MLS_v_absAbs" }
      if(modelx == 18){ modelid <- "LM.MLS_v_binaryAbs" }
      # ExtremeLongevity - relative
      if(modelx == 19){ modelid <- "extremeLongevity_sd2_all.logwMLSresid_v_relAbs" }
      if(modelx == 20){ modelid <- "extremeLongevity_sd2_hi.logwMLSresid_v_relAbs" }
      if(modelx == 21){ modelid <- "extremeLongevity_sd2_lo.logwMLSresid_v_relAbs" }
      if(modelx == 22){ modelid <- "extremeLongevity_sd1_all.logwMLSresid_v_relAbs" }
      if(modelx == 23){ modelid <- "extremeLongevity_sd1_hi.logwMLSresid_v_relAbs" }
      if(modelx == 24){ modelid <- "extremeLongevity_sd1_lo.logwMLSresid_v_relAbs" }
      if(modelx == 25){ modelid <- "extremeLongevity_sd1andsd2_all.logwMLSresid_v_relAbs" }
      if(modelx == 26){ modelid <- "extremeLongevity_sd1andsd2_hi.logwMLSresid_v_relAbs" }
      if(modelx == 27){ modelid <- "extremeLongevity_sd1andsd2_lo.logwMLSresid_v_relAbs" }
      # ExtremeLongevity - binary
      if(modelx == 28){ modelid <- "extremeLongevity_sd2_all.logwMLSresid_v_binaryAbs" }
      if(modelx == 29){ modelid <- "extremeLongevity_sd2_hi.logwMLSresid_v_binaryAbs" }
      if(modelx == 30){ modelid <- "extremeLongevity_sd2_lo.logwMLSresid_v_binaryAbs" }
      if(modelx == 31){ modelid <- "extremeLongevity_sd1_all.logwMLSresid_v_binaryAbs" }
      if(modelx == 32){ modelid <- "extremeLongevity_sd1_hi.logwMLSresid_v_binaryAbs" }
      if(modelx == 33){ modelid <- "extremeLongevity_sd1_lo.logwMLSresid_v_binaryAbs" }
      if(modelx == 34){ modelid <- "extremeLongevity_sd1andsd2_all.logwMLSresid_v_binaryAbs" }
      if(modelx == 35){ modelid <- "extremeLongevity_sd1andsd2_hi.logwMLSresid_v_binaryAbs" }
      if(modelx == 36){ modelid <- "extremeLongevity_sd1andsd2_lo.logwMLSresid_v_binaryAbs" }
      if(modelx == 37){ modelid <- "Family_v_relAbs" }
      if(modelx == 38){ modelid <- "Family_v_absAbs" }
      if(modelx == 39){ modelid <- "Family_v_binaryAbs" }

      sum_mat <- matrix(ncol=7,nrow=length(rownames(abundance_matrix.abs)))
      colnames(sum_mat) <- c("cutoff","model","raw_pvalues","fdr_pvalues","posneg","slopes","intercepts")
      rownames(sum_mat) <- rownames(abundance_matrix.abs)
      
        for(cutoff in c(0:46)){

          runlist <- rownames(abundance_matrix.abs)[shared_abundances == cutoff]
            
          ogi <- 1
          for(og in runlist){

            if(ogi == 1 | ogi %% 20 == 0){  
              print(paste(runtype, ": (", modelx, ")", modelid, paste0(round((modelx/39)*100, digits=2),"%"), ": cutoff",cutoff, paste0( round((ogi/length( runlist ))*100, digits=2),"%")))
            }
              
            if(modelx %in% c(1,4,7,10,13,16,19:27,37)){ 
              # relative abundance
              compare_me     <- as.data.frame(matrix(ncol=22,nrow=sum(abundance_matrix.abs[og,] > 0)))
              compare_me[,1] <- as.numeric(lmres[ abundance_matrix.rel[og,] > 0 ])
              compare_me[,2] <- as.numeric(abundance_matrix.rel[og,][abundance_matrix.rel[og,] > 0])
            }else if(modelx %in% c(2,5,8,11,14,17,38)){
              # absolute abundance
              compare_me     <- as.data.frame(matrix(ncol=22,nrow=sum(abundance_matrix.abs[og,] > 0)))
              compare_me[,1] <- as.numeric(lmres[ abundance_matrix.abs[og,] > 0 ])
              compare_me[,2] <- as.numeric(abundance_matrix.abs[og,][abundance_matrix.abs[og,] > 0])
            }else{
              # binary abundance
              compare_me     <- as.data.frame(matrix(ncol=22,nrow=length(abundance_matrix.abs)))
              compare_me[,1] <- as.numeric(lmres)
              compare_me[,2] <- as.factor(abundance_matrix.abs[og,] > 0)
            }
            
            if(modelx %in% c(1,2,19:27)){
              #  log( mass )
              compare_me[,3] <- log(as.numeric(mls$Body.Mass.Grams[abundance_matrix.rel[og,] > 0]))
            }else if(modelx %in% c(4,5,7,8,13,14)){
              #  mass
              compare_me[,3] <- as.numeric(mls$Body.Mass.Grams[abundance_matrix.rel[og,] > 0])
            }else if(modelx %in% c(3, 12, 9, 28:36,39)){
              compare_me[,3] <- log(as.numeric(mls$Body.Mass.Grams))
            }else if(modelx %in% c(6, 15, 18)){
              compare_me[,3] <- as.numeric(mls$Body.Mass.Grams)
            }else{
              compare_me[,3] <- log(as.numeric(mls$Body.Mass.Grams[abundance_matrix.rel[og,] > 0]))
            }
            
            if(modelx %in% c(1,2,10,11,4,5,7,8,13,14,16,17,19:27,37:38)){
              compare_me[,4] <- as.numeric(mls$MLS.Years[abundance_matrix.rel[og,] > 0])
            }else{
              compare_me[,4] <- as.numeric(mls$MLS.Years)
            }
            
            #if(modelx %in% c(1,2,3,10,11,12,5,6,7,13,14,15,19:27)){
            if(modelx %in% c(1,2,10,11,4,5,7,8,13,14,16,17,19:27,37:38)){
                
              compare_me[,5 ] <- as.factor(mls$Family[abundance_matrix.rel[og,] > 0])
              compare_me[,6 ] <- as.character(mls$Source[abundance_matrix.rel[og,] > 0])
              compare_me[,7 ] <- as.character(mls$Query.is.Best.Ref[abundance_matrix.rel[og,] > 0])
              compare_me[,8 ] <- as.character(mls$TreeFormat[abundance_matrix.rel[og,] > 0])
              compare_me[,9 ] <- as.character(mls$AbundanceOrder[abundance_matrix.rel[og,] > 0])
              compare_me[,10] <- as.character(mls$Best.Reference.Common.Name[abundance_matrix.rel[og,] > 0])
              compare_me[,11] <- as.character(mls$Scientific.Name[abundance_matrix.rel[og,] > 0])
              compare_me[,12] <- as.character(mls$Common.Name[abundance_matrix.rel[og,] > 0])
              compare_me[,13] <- as.character(mls$AlignLabel[abundance_matrix.rel[og,] > 0])
              
              compare_me[,14] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd2_all )[abundance_matrix.rel[og,] > 0]
              compare_me[,15] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd2_hi )[abundance_matrix.rel[og,] > 0]
              compare_me[,16] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd2_lo )[abundance_matrix.rel[og,] > 0]
              
              compare_me[,17] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1_all )[abundance_matrix.rel[og,] > 0]
              compare_me[,18] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1_hi )[abundance_matrix.rel[og,] > 0]
              compare_me[,19] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1_lo )[abundance_matrix.rel[og,] > 0]
              
              compare_me[,20] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1andsd2_all )[abundance_matrix.rel[og,] > 0]
              compare_me[,21] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1andsd2_hi )[abundance_matrix.rel[og,] > 0]
              compare_me[,22] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1andsd2_lo )[abundance_matrix.rel[og,] > 0]
              
            }else{
              compare_me[,5 ] <- as.factor(mls$Family)
              compare_me[,6 ] <- as.character(mls$Source)
              compare_me[,7 ] <- as.character(mls$Query.is.Best.Ref)
              compare_me[,8 ] <- as.character(mls$TreeFormat)
              compare_me[,9 ] <- as.character(mls$AbundanceOrder)
              compare_me[,10] <- as.character(mls$Best.Reference.Common.Name)
              compare_me[,11] <- as.character(mls$Scientific.Name)
              compare_me[,12] <- as.character(mls$Common.Name)
              compare_me[,13] <- as.character(mls$AlignLabel)
              
              compare_me[,14] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd2_all )
              compare_me[,15] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd2_hi )
              compare_me[,16] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd2_lo )
              
              compare_me[,17] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1_all )
              compare_me[,18] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1_hi )
              compare_me[,19] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1_lo )
              
              compare_me[,20] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1andsd2_all )
              compare_me[,21] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1andsd2_hi )
              compare_me[,22] <- as.factor( mls$TreeFormat %in% extremeLongevity_sd1andsd2_lo )
            }
            colnames(compare_me) <- c("LMres.logW.v.MLS",
                                      "OG.Abundance",
                                      "Body.Mass.Grams",
                                      "MLS.Years",
                                      "Family",
                                      "Source",
                                      "Query.is.Best.Ref",
                                      "TreeFormat",
                                      "AbundanceOrder",
                                      "Best.Reference.Common.Name",
                                      "Scientific.Name",
                                      "Common.Name",
                                      "AlignLabel",
                                      "extremeLongevity_sd2_all","extremeLongevity_sd2_hi","extremeLongevity_sd2_lo",
                                      "extremeLongevity_sd1_all","extremeLongevity_sd1_hi","extremeLongevity_sd1_lo",
                                      "extremeLongevity_sd1andsd2_all","extremeLongevity_sd1andsd2_hi","extremeLongevity_sd1andsd2_lo")
            
            if(modelx %in% 1:18){     f <- as.formula("LMres.logW.v.MLS ~ OG.Abundance")  }
            if(modelx %in% c(19,28)){ f <- as.formula("extremeLongevity_sd2_all ~ OG.Abundance")  }
            if(modelx %in% c(20,29)){ f <- as.formula("extremeLongevity_sd2_hi ~ OG.Abundance")  }
            if(modelx %in% c(21,30)){ f <- as.formula("extremeLongevity_sd2_lo ~ OG.Abundance")  }
            if(modelx %in% c(22,31)){ f <- as.formula("extremeLongevity_sd1_all ~ OG.Abundance")  }
            if(modelx %in% c(23,32)){ f <- as.formula("extremeLongevity_sd1_hi ~ OG.Abundance")  }
            if(modelx %in% c(24,33)){ f <- as.formula("extremeLongevity_sd1_lo ~ OG.Abundance")  }
            if(modelx %in% c(25,34)){ f <- as.formula("extremeLongevity_sd1andsd2_all ~ OG.Abundance")  }
            if(modelx %in% c(26,35)){ f <- as.formula("extremeLongevity_sd1andsd2_hi ~ OG.Abundance")  }
            if(modelx %in% c(27,36)){ f <- as.formula("extremeLongevity_sd1andsd2_lo ~ OG.Abundance")  }
            if(modelx %in% c(37:39)){ f <- as.formula("Family ~ OG.Abundance")  }
            
            compdata <- comparative.data(phylotree, compare_me, TreeFormat)
            if(modelx %in% c(10:18)){
              reg   <- lm(f, data=compare_me)
            }else if (modelx %in% c(3,6,9,12,15,18,19:39)) {
              # (Burt, 1989)
              reg   <- brunch(f, data=compdata)
            }else{
              reg   <- crunch(f, data=compdata)
            }
            
            if (t(coef(summary(reg)))[2,1] > 0 & is.nan(t(coef(summary(reg)))[2,1]) == 0){
              cpus.lm2 <- as.matrix( anova(reg) )
              fstat    <- summary(reg)$fstatistic
              posneg   <- as.numeric(slope > 0)
              raw_pval <- unname(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
              if (modelx %in% c(19:39)) {
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
            sum_mat[ og , ] <- c(cutoff, modelid, raw_pval, NA, posneg, slope, int)
            
            #print(paste(cutoff, modelid, raw_pval, posneg, slope, int))
            
            ogi <- ogi + 1
            
          } # end for og
          
        } # end for cutoff
      
      #write.table(sum_mat, file=paste0(modelx,".",modelid,".tsv"), sep = "\t")
      
      sum_mat.pop <- sum_mat[ !is.na( sum_mat[,"raw_pvalues"] ) , ]
      sum_mat.pop[,"fdr_pvalues"] <- p.adjust( as.numeric(sum_mat.pop[,3] ) , method="fdr" )
      sum_mat.pop <- sum_mat.pop[ order(sum_mat.pop[,"raw_pvalues"]) , ]
      write.table(sum_mat.pop, file=paste0(runtype,".",modelx,".",modelid,"report.tsv"), sep = "\t")
      
      #if(modelx %in% c(1,2,10,11,4,5,7,8,13,14,16,17,19:27,37:38)){
        # colnum <- 47
      #}else{
        # colnum <- 50
      #}
        
      png(filename = paste0(runtype,".",modelx,".",modelid,".rawP.png"), height = 800,width = 800)
      
        maxX <- max(abs( c(min(as.numeric(sum_mat.pop[1:100,"slopes"])),max(as.numeric(sum_mat.pop[1:100,"slopes"]))) ))
        plot( as.numeric(sum_mat.pop[1:100,"slopes"]) , -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ),
              xlim=c(-maxX,maxX),
              pch=16, cex=2,
              col = redgreen(47)[ as.numeric(sum_mat.pop[,"cutoff"])[1:100]+1 ],
              main = paste(runtype,modelid),
              ylab = "-log( Raw P-value )",
              xlab = "Slope of Association"
              )
        
        abline(h=-log(0.05), col="purple", lty=3)
        text(maxX, -log(0.05), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.05) ),"significant @ 0.05\n"), pos=2, cex = 0.7)
        abline(h=-log(0.01), col="purple", lty=3)
        text(maxX, -log(0.01), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.01) ),"significant @ 0.01\n"), pos=2, cex = 0.7)
        abline(h=-log(0.005), col="purple", lty=3)
        text(maxX, -log(0.005), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.005) ),"significant @ 0.005\n"), pos=2, cex = 0.7)
        abline(h=-log(0.001), col="purple", lty=3)
        text(maxX, -log(0.001), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.001) ),"significant @ 0.001\n"), pos=2, cex = 0.7)
        abline(h=-log(0.0005), col="purple", lty=3)
        text(maxX, -log(0.0005), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.0005) ),"significant @ 0.0005\n"), pos=2, cex = 0.7)
        abline(h=-log(0.0001), col="purple", lty=3)
        text(maxX, -log(0.0001), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.0001) ),"significant @ 0.0001\n"), pos=2, cex = 0.7)
        
        legend("topleft", col = c("red","black","green"), title = "Queries Unrepresented by Cluster", legend = c(0,23,46), pch=16)
        
      dev.off()
       
      sum_mat.pop.fdrsort <- sum_mat.pop[ order(sum_mat.pop[,"fdr_pvalues"]) , ]
      png(filename = paste0(runtype,".",modelx,".",modelid,".fdrP.png"), height = 800,width = 800)
        maxX <- max(abs( c(min(as.numeric(sum_mat.pop.fdrsort[1:100,"slopes"])),max(as.numeric(sum_mat.pop.fdrsort[1:100,"slopes"]))) ))
        plot( as.numeric(sum_mat.pop.fdrsort[1:100,"slopes"]) , -log( as.numeric(sum_mat.pop.fdrsort[1:100,"fdr_pvalues"]) ) ,
              xlim=c(-maxX,maxX),
              pch=16, cex=2,
              col = redgreen(47)[ as.numeric(sum_mat.pop[,"cutoff"])[1:100]+1 ],
              main = paste(runtype,modelid),
              ylab = "-log( Raw P-value )",
              xlab = "Slope of Association"
              )
        abline(h=-log(0.05), col="purple", lty=3)
        text(maxX, -log(0.05), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.05) ),"significant @ 0.05\n"), pos=2, cex = 0.7)
        abline(h=-log(0.01), col="purple", lty=3)
        text(maxX, -log(0.01), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.01) ),"significant @ 0.01\n"), pos=2, cex = 0.7)
        abline(h=-log(0.005), col="purple", lty=3)
        text(maxX, -log(0.005), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.005) ),"significant @ 0.005\n"), pos=2, cex = 0.7)
        abline(h=-log(0.001), col="purple", lty=3)
        text(maxX, -log(0.001), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.001) ),"significant @ 0.001\n"), pos=2, cex = 0.7)
        abline(h=-log(0.0005), col="purple", lty=3)
        text(maxX, -log(0.0005), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.0005) ),"significant @ 0.0005\n"), pos=2, cex = 0.7)
        abline(h=-log(0.0001), col="purple", lty=3)
        text(maxX, -log(0.0001), col="purple", paste(sum( -log( as.numeric(sum_mat.pop[1:100,"raw_pvalues"]) ) > -log(0.0001) ),"significant @ 0.0001\n"), pos=2, cex = 0.7)
        
        legend("topleft", col = c("red","black","green"), title = "Queries Unrepresented by Cluster", legend = c(0,23,46), pch=16)
      dev.off()
      
    } # end for modelx
  
}# end for runtype












