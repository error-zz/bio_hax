
library(gplots)
library(colorRamps)
library(vegan)
library(ape)

#setwd ("/Users/apple/Desktop")
#clust_sorted.I14 <- read.table("clust_sorted.I14", sep="\t", header=FALSE)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
clust_sorted.I14 <- read.table( args[1], sep="\t", header=FALSE )

cols <- c("NumUniqRef","NumUniqQuery")
clust_info           <- matrix(0,ncol=2,nrow=length(unique(clust_sorted.I14$V4)))
clust_abund_ref      <- matrix(0,ncol=12,nrow=length(unique(clust_sorted.I14$V4)))
clust_abund_species  <- matrix(0,ncol=53,nrow=length(unique(clust_sorted.I14$V4)))
rownames(clust_info) <- rownames(clust_abund_ref) <- rownames(clust_abund_species) <- as.numeric(unique(clust_sorted.I14$V4))
colnames(clust_info) <- cols
colnames(clust_abund_ref) <- names(unlist(table(clust_sorted.I14[clust_sorted.I14$V4==2,]$V5)))
colnames(clust_abund_species) <- names(unlist(table(clust_sorted.I14[clust_sorted.I14$V4==2,]$V6)))

ci <- 1
for(c in rownames(clust_info)){
  if(ci %% 100 == 0){
    print(paste0( (ci/length(unique(clust_sorted.I14$V4)))*100, "%"))
  }
  subset <- clust_sorted.I14[clust_sorted.I14$V4==c,]
  clust_info[c,1] <- length(unique(clust_sorted.I14[clust_sorted.I14$V4==c,]$V6)) # table
  clust_info[c,2]   <- length(unique(clust_sorted.I14[clust_sorted.I14$V4==c,]$V5)) # table
  clust_abund_ref[c,] <- unname(unlist(table(clust_sorted.I14[clust_sorted.I14$V4==c,]$V5)))
  clust_abund_species[c,] <- unname(unlist(table(clust_sorted.I14[clust_sorted.I14$V4==c,]$V6)))
  # print(paste(c,clust_info[c,1],clust_info[as.factor(c),2],length(subset)))
  ci <- ci + 1
}

# save.image(file="clust_sorted.I14.Rdata")

png(filename=paste0("PerClustUniqQuery.png"), width = 1000, height = 1000 )
plot.new()
hist(clust_info[,1], xlab = "", breaks=52, 
     col ="darkgray",
     main="Number of Unique Query Species Represented per Cluster ID\n")
dev.off()

png(filename=paste0("PerClustUniqRef.png"), width = 1000, height = 1000 )
plot.new()
hist(clust_info[,2], xlab = "", breaks=11 ,
     col ="darkgray",
     main="Number of Unique Reference Species Represented per Cluster ID\n")
dev.off()


png(filename=paste0("ClustWithGrEq53queries.png"), width = 1000, height = 1000 )
plot.new()
plot_abund_subset <- clust_abund_ref[ rownames(clust_info[clust_info[,1] == 53,]) , ]
heatmap.2( 
  log(plot_abund_subset+0.1),
  col=c(rev(gray.colors(150)),matlab.like2(150)),
  trace="none",
  margins = c(15,15),
  keysize=0.7, key.par = list(cex=0.5),
  main="Ortholog Clusters representing all 53 query species\nWhere are most of the reference alignments?\n(To the most common best references.)"
)
dev.off()

png(filename=paste0("ClustWithGrEq1queries.png"), width = 1000, height = 1000 )
plot.new()
plot_abund_subset <- clust_abund_ref[ rownames(clust_info[clust_info[,1] == 1,]) , ]
heatmap.2( 
  log(plot_abund_subset+0.1),
  col=c(rev(gray.colors(150)),matlab.like2(150)),
  trace="none",
  margins = c(15,15),
  keysize=0.7, key.par = list(cex=0.5),
  main="Ortholog Clusters representing only 1 query species\nWhere are most of the reference alignments?"
)
dev.off()

png(filename=paste0("AllOrth_RefAlign.png"), width = 1000, height = 1000 )
plot.new()
plot_abund_subset <- clust_abund_ref
heatmap.2( 
  log(plot_abund_subset+0.1),
  col=c(rev(gray.colors(150)),matlab.like2(150)),
  trace="none",
  margins = c(15,15),
  keysize=0.7, key.par = list(cex=0.5),
  main="All Ortholog Clusters\nWhere are most of the reference alignments?"
)
dev.off()

png(filename=paste0("ClustWithGrEq53queries.png"), width = 1000, height = 1000 )
plot.new()
plot_abund_subset <- clust_abund_species[ rownames(clust_info[clust_info[,1] == 53,]) , ]
heatmap.2( 
  log(plot_abund_subset+0.1),
  col=c(rev(gray.colors(150)),matlab.like2(150)),
  trace="none",
  margins = c(15,15),
  keysize=0.7, key.par = list(cex=0.5),
  main="Ortholog Clusters representing all 53 query species\nWhich queries are most represented?"
)
dev.off()

png(filename=paste0("ClustWithGrEq53queries_gr212.png"), width = 1000, height = 1000 )
plot.new()
plot_abund_subset <- clust_abund_species[ rownames(clust_info[clust_info[,1] == 53,]) , ][ rowSums(clust_abund_species[ rownames(clust_info[clust_info[,1] == 53,]) , ]) > 212,  ]
heatmap.2( 
  log(plot_abund_subset+0.1),
  col=c(rev(gray.colors(150)),matlab.like2(150)),
  trace="none",
  margins = c(15,15),
  keysize=0.7, key.par = list(cex=0.5),
  main="Ortholog Clusters representing all 53 query species\nWhich queries are most represented?\nTotal Counts for Cluster ID > 212 (> 4 avg per query)"
)
dev.off()

png(filename=paste0("allClust.png"), width = 1000, height = 1000 )
plot.new()
plot_abund_subset <- clust_abund_species[ rownames(clust_info[clust_info[,1] == 1,]) , ]
heatmap.2( 
  log(plot_abund_subset+0.1),
  col=c(rev(gray.colors(150)),matlab.like2(150)),
  trace="none",
  margins = c(15,15),
  keysize=0.7, key.par = list(cex=0.5),
  main="Ortholog Clusters representing only 1 query species\nWhich queries are most represented?"
)
dev.off()

png(filename=paste0("allClust_gr212.png"), width = 1000, height = 1000 )
plot.new()
plot_abund_subset <- clust_abund_species[ rowSums(clust_abund_species) > 212 , ]
heatmap.2( 
  log(plot_abund_subset+0.1),
  col=c(rev(gray.colors(150)),matlab.like2(150)),
  trace="none",
  margins = c(15,15),
  keysize=0.7, key.par = list(cex=0.5),
  main="All Ortholog Clusters\nWhich queries are most represented??\nTotal Counts for Cluster ID > 212"
)
dev.off()

write.table(table(clust_info[,1]), "query_histo.tsv", sep = '\t')

write.table(table(clust_info[,2]), "ref_histo.tsv", sep = '\t')
# table(clust_info[,2])




expected_reference_ids <- c("ChippingSparrow","SongSparrow","AmericanTreeSparrow","Common_Grackle","YellowThroatedWarbler","YellowRumpedWarbler","YellowWarbler",
                            "BrownHeadedCowbird","NorthernCardinal","HouseFinch","HouseSparrow","CarolinaWren","Wren","WhiteBreastedNuthatch","AmericanRobin",
                            "GrayCatbird","Starling","CedarWaxwing","BarnSwallow","TreeSwallow","HornedLark","TuftedTitmouse")

# 53 #
#n <- 53
for(n in c(52,51,50,49,48,47,46)){
  #transcript_abund <- matrix(0, ncol=53,nrow=length(sort(unique(transcript_subset_for_cluster_n$V1))))
  transcript_abund <- matrix(0, ncol=53,nrow=length(sort(unique(clust_sorted.I14$V1))))
  colnames(transcript_abund) <- unique(clust_sorted.I14$V6)
  rownames(transcript_abund) <- sort(unique(clust_sorted.I14$V1))
  
  #clust_abund <- matrix(0, ncol=53,nrow=length(sort(unique(transcript_subset_for_cluster_n$V4))))
  clust_abund <- matrix(0, ncol=53,nrow=max(sort(unique(clust_sorted.I14$V4))))
  colnames(clust_abund) <- unique(clust_sorted.I14$V6)
  rownames(clust_abund) <- c(1:max(sort(unique(clust_sorted.I14$V4))))
    # as.factor(sort(unique(clust_sorted.I14$V4)))
   
  clust_sorted.I14.n_subset <- clust_sorted.I14[ clust_sorted.I14$V4 %in% rownames(clust_info[clust_info[,1] >= n , ]) , ]
  ci <- 1
  #for(c in 1:dim(clust_sorted.I14.n_subset)[1]){
  for(c in rev(unique(clust_sorted.I14[ clust_sorted.I14$V4 %in% rownames(clust_info[clust_info[,1] >= n , ]) , ]$V4))){
    # for each cluster id
    #if(ci %% 100 == 0){ 
      print(paste0( ci," : ",c," : ",(ci/dim(clust_sorted.I14.n_subset)[1])*100,"%" )) 
    #}
    ti <- 1
    for(t in unique(clust_sorted.I14.n_subset[ clust_sorted.I14.n_subset$V4 == as.factor(c) , ]$V1)){
      # print(paste0("  ",t))
      # for each transcript ID in each cluster ID
      for(query in unique(clust_sorted.I14.n_subset[ clust_sorted.I14.n_subset$V1 == t, ]$V6)){
        # print(paste0("     ",query))
        query_specific_transcript_sum <- sum(clust_sorted.I14.n_subset[ clust_sorted.I14.n_subset$V1 == t, ][ clust_sorted.I14.n_subset[ clust_sorted.I14.n_subset$V1 == t, ]$V6 == query,  ]$V3)
        transcript_abund[t,query] <- transcript_abund[t,query] + query_specific_transcript_sum
        clust_abund[c,query]      <- clust_abund[c,query]      + query_specific_transcript_sum
      }
      if(ti %% 100 == 0){ 
        print(paste0( "     ",ti," :",t)) 
      }
      ti <- ti + 1
    }
    ci <- ci + 1
  }
     
  write.table(transcript_abund, paste0("transcript_abund.", n, ".tsv"), sep = '\t')
  write.table(clust_abund, paste0("clust_abund.", n, ".tsv"), sep = '\t')
  save.image(file="clustsorted.Rdata")
  
  ca <- clust_abund[rowSums(clust_abund)>0,]
  ta <- transcript_abund[rowSums(transcript_abund)>0,]
  
  print(paste("DIM CA : ",n,":",dim(ca)))
  print(paste("DIM TA : ",n,":",dim(ta)))
  
  ca.rel <- t(decostand(t(ca), method='total'))
  ta.rel <- t(decostand(t(ta), method='total'))
  
  png(filename=paste0(n,".Clust_log_relAbs_all_in_gr20species.png"), width = 1400, height = 1400 )
  heatmap.2( 
    log(ca+0.001),
    col=c(rev(gray.colors(150)),matlab.like2(150)),
    trace="none",
    margins = c(15,15),
    keysize=0.7, key.par = list(cex=0.5),
    main="Summed Ortholog Cluster (Absolute Abundance)\nOrtholog Clusters Represented by > 20 Query Species"
  )
  dev.off() 
  
  png(filename=paste0("Clust_log_relAbs_gr0.075_in_gr20species.png"), width = 1400, height = 1400 )
  heatmap.2( 
    -log(ca.rel[rowSums(ca.rel)>0.075,]+0.01),
    col=c(rev(gray.colors(150)),matlab.like2(150)),
    trace="none",
    margins = c(15,15),
    keysize=0.7, key.par = list(cex=0.5),
    main="Summed Ortholog Cluster\n(-log(relativeAbundance) > 0.075)\nOrtholog Clusters Represented by > 20 Query Species"
  )
  dev.off() 
  
  png(filename=paste0(n,".Trans_RelAbs_gr24sharedQuery.png"), width = 1400, height = 1400 )
  heatmap.2( 
    #log(ta+0.001),
    -log( ta.rel[ rowSums(ta.rel==0) < 29 , ] + 0.005),
    col=c(rev(gray.colors(150)),matlab.like2(150)),
    trace="none",
    margins = c(15,15),
    keysize=0.7, key.par = list(cex=0.5),
    main="Summed Ortholog Transcript -log(relAbundance)\nTranscripts Represented by > 24 Query Species"
  )
  dev.off()
  
  png(filename=paste0("Trans_RelAbs_gr24sharedQuery_maxSum50.png"), width = 1400, height = 1400 )
  heatmap.2( 
    #log(ta+0.001),
    -log( ta.rel[ names(sort(rowSums(ta.rel[ rowSums(ta.rel==0) < 40 , ]))[1:50]) , ] + 0.01),
    col=c(rev(gray.colors(150)),matlab.like2(150)),
    trace="none",
    margins = c(15,15),
    keysize=0.7, key.par = list(cex=0.5),
    main="Summed Ortholog Transcript (Absolute Abundance)\nOrtholog Clusters Represented by > 24 Query Species"
  )
  dev.off()
  
  png(filename=paste0(n,".Trans_RelAbs_gr24sharedQuery_maxSum100.png"), width = 1400, height = 1400 )
  heatmap.2( 
    #log(ta+0.001),
    -log( ta.rel[ names(sort(rowSums(ta.rel[ rowSums(ta.rel==0) < 40 , ]))[1:100]) , ] + 0.01),
    col=c(rev(gray.colors(150)),matlab.like2(150)),
    trace="none",
    margins = c(15,15),
    keysize=0.7, key.par = list(cex=0.5),
    main="Summed Ortholog Transcript (Absolute Abundance)\nOrtholog Clusters Represented by > 24 Query Species"
  )
  dev.off()
  
  png(filename=paste0(n,".Trans_AbsAbs_gr24sharedQuery_maxSum100.png"), width = 1400, height = 1400 )
  heatmap.2( 
    #log(ta+0.001),
    log( ta[ names(sort(rowSums(ta.rel[ rowSums(ta.rel==0) < 40 , ]))[1:100]) , ] + 0.01),
    col=c(rev(gray.colors(150)),matlab.like2(150)),
    trace="none",
    margins = c(15,15),
    keysize=0.7, key.par = list(cex=0.5),
    main="Summed Ortholog Transcript (Absolute Abundance)\nOrtholog Clusters Represented by > 24 Query Species"
  )
  dev.off()
  
}

 
#rtmp <- rownames(clust_info[clust_info[,1] >= 45 , ])
#write.table(rtmp, file="c45.tsv")
#rtmp <- rownames(clust_info[clust_info[,1] >= 39 , ])
#write.table(rtmp, file="c39.tsv")
#rtmp <- rownames(clust_info[clust_info[,1] >= 30 , ])
#write.table(rtmp, file="c30.tsv")
#rtmp <- rownames(clust_info[clust_info[,1] >= 20 , ])
#write.table(rtmp, file="c20.tsv")

# clust_abund_species[ rownames(clust_info[clust_info[,1] >= 45 , ]) , ]
# clust_abund_species[ rownames(clust_info[clust_info[,1] >= 39 , ]) , ]
# clust_abund_species[ rownames(clust_info[clust_info[,1] >= 30 , ]) , ]
# clust_abund_species[ rownames(clust_info[clust_info[,1] >= 20 , ]) , ]



