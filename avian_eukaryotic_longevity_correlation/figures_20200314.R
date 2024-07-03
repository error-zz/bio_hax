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

load("postcorr_20200212.Rdata")

pp <-  colnames(model.p.dos)[ grepl('p_val$', colnames(model.p.dos)) ]
fd <-  colnames(model.p.dos)[ grepl('p_val.FDR', colnames(model.p.dos)) ] 

##########################################################################
#mmmmm <- model.p.dos[ !is.na(model.p.dos$`11.f_05.CAPER.p_val.FDR`), ]

all <- model.p.dos[ rownames( model.p.dos[ model.p.dos$`11.f_05.CAPER.p_val.FDR` < 0.05 , ] ) , ]
 
pos <- all[ all$`11.f_05.CAPER.ispos`== 1, ] 

neg <- all[ all$`11.f_05.CAPER.ispos`== 0, ] 
 


all2 <- model.p.dos[ rownames( model.p.dos[ model.p.dos$`23.f_05.CAPER.p_val.FDR` < 0.05 , ] ) , ]

pos2 <- all2[ all2$`23.f_05.CAPER.ispos`== 1, ] 

neg2 <- all2[ all2$`23.f_05.CAPER.ispos`== 0, ] 



all3 <- model.p.dos[ rownames( model.p.dos[ model.p.dos$`11.f_02.CAPER.p_val.FDR` < 0.05 , ] ) , ]

pos3 <- all3[ all3$`11.f_02.CAPER.ispos`== 1, ] 

neg3 <- all3[ all3$`11.f_02.CAPER.ispos`== 0, ] 

all4 <- model.p.dos[ rownames( model.p.dos[ model.p.dos$`23.f_02.CAPER.p_val.FDR` < 0.05 , ] ) , ]

pos4 <- all4[ all4$`23.f_02.CAPER.ispos`== 1, ] 

neg4 <- all4[ all4$`23.f_02.CAPER.ispos`== 0, ] 

 

sum( rownames(pos) %in% rownames(pos2) )-1
sum( rownames(neg) %in% rownames(neg2) )-1
sum( rownames(pos) %in% rownames(neg2) )-1


sum( rownames(pos) %in% rownames(pos3) )-1
sum( rownames(neg) %in% rownames(neg3) )-1
sum( rownames(pos) %in% rownames(neg3) )-1
 

sum( rownames(pos2) %in% rownames(pos3) )-1
sum( rownames(neg2) %in% rownames(neg3) )-1
sum( rownames(pos2) %in% rownames(neg3) )-1



sum( rownames(pos3) %in% rownames(pos4) )-1
sum( rownames(neg3) %in% rownames(neg4) )-1
sum( rownames(pos3) %in% rownames(neg4) )-1
 

sum( rownames(neg) %in% rownames(neg2) )
sum( rownames(neg2) %in% rownames(neg3) )
sum( rownames(neg3) %in% rownames(neg4) )
 
sum( rownames(pos3) %in% rownames(pos4) )


r <- rownames(pos)
r <- r[ r %in% rownames(pos2) ]
r <- r[ r %in% rownames(pos3) ]
allpos <- r[ r %in% rownames(pos4) ]


r <- rownames(neg)
r <- r[ r %in% rownames(neg2) ]
r <- r[ r %in% rownames(neg3) ]
allneg <- r[ r %in% rownames(neg4) ]

r <- rownames(all)
r <- r[ r %in% rownames(all2) ]
r <- r[ r %in% rownames(all3) ]
allall <- r[ r %in% rownames(all4) ]
  

allall[!is.na(allall)]
 

##########################################################################

# https://www.orthodb.org/?page=filelist

# v9_v10_OGs_map.tab.gz
# 1.	level tax id where both OGs are built
# 2.	previous OG id
# 3.	current OG id
# 4.	distance between the two OGs (0. - identical, 1.0 - totally unrelated)
# 5.	flag if the groups are "the best match" for each other

v9_v10_OGs_map     <- read.table("v9_v10_OGs_map.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)

pos.cur_og_id      <- v9_v10_OGs_map[ v9_v10_OGs_map[,2] %in% allpos , 3 ]
neg.cur_og_id      <- v9_v10_OGs_map[ v9_v10_OGs_map[,2] %in% allneg , 3 ]
all.cur_og_id      <- v9_v10_OGs_map[ v9_v10_OGs_map[,2] %in% allall , 3 ]
 


# get pathway info
# odb10v1_OG_xrefs.tab
odb10v1_OG_xrefs     <- read.table("odb10v1_OG_xrefs.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
pos.cur_og_path      <- odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] %in% pos.cur_og_id , c(2,3) ]
neg.cur_og_path      <- odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] %in% neg.cur_og_id , c(2,3) ]
all.cur_og_path      <- odb10v1_OG_xrefs[ odb10v1_OG_xrefs[,1] %in% all.cur_og_id , c(2,3) ]
 
#save.image("oops.Rdata")

#rownames(neg)[ !is.na(rownames(neg)) ]
#rownames(pos)[ !is.na(rownames(pos)) ]

# get gene identifier (numbers)
# odb10v1_OG2genes.tab
odb10v1_OG2genes     <- read.table("odb10v1_OG2genes.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
pos.cur_og_col      <- odb10v1_OG2genes[ odb10v1_OG2genes[,1] %in% pos.cur_og_id , 2 ]
neg.cur_og_col      <- odb10v1_OG2genes[ odb10v1_OG2genes[,1] %in% neg.cur_og_id , 2 ]
all.cur_og_col      <- odb10v1_OG2genes[ odb10v1_OG2genes[,1] %in% all.cur_og_id , 2 ]
 


#save.image("oops.Rdata")
load("oops.Rdata")
  
# get gene common name, interpro id, gene id, protein
# odb10v1_gene_xrefs.tab
# odb10v1_gene_xrefs     <- read.table("odb10v1_gene_xrefs.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
# ^ too big

# get ensembl identifier (numbers)
# odb10v1_genes.tab
#odb10v1_genes     <- read.table("odb10v1_genes.tab", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)



allall.cur_og_col      <- odb10v1_OG2genes[ odb10v1_OG2genes[,1] %in% all.cur_og_id , 2 ]




############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################


 
library("topGO")
library(biomaRt)
library(ALL)
library(colorRamps)
library(Rgraphviz)
topDiffGenes <- function(allScore) {
  return(allScore > 1)
}


bm <- useMart("ensembl")
bm <- useDataset("hsapiens_gene_ensembl", mart=bm)

geneListX <- c("AIG1","DYNLT1","EZR","GNS","IKBKE","LEMD3","MSRB3","NMBR","PEX3","RASSF3","RPL18A","RXYLT1","TBK1","WIF1","XPOT","FAKE")
geneList <- factor(c(as.integer(rep(1,length(geneListX)-1)),0))
names(geneList) <- geneListX

results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'go_id'), 
                                  filter = "hgnc_symbol",
                                  values = c("AIG1","DYNLT1","EZR","GNS","IKBKE","LEMD3","MSRB3","NMBR","PEX3","RASSF3","RPL18A","RXYLT1","TBK1","WIF1","XPOT"),
                                  mart = bm)

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

queryFisher.classic.fisher.bp <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "fisher")
queryFisher.classic.ks.bp     <- runTest(queryGOdata.bp, algorithm = "classic", statistic = "ks"    )
queryFisher.elim.ks.bp     <- runTest(queryGOdata.bp, algorithm = "elim",    statistic = "ks"    )

queryFisher.classic.fisher.mf <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "fisher")
queryFisher.classic.ks.mf     <- runTest(queryGOdata.mf, algorithm = "classic", statistic = "ks"    )
queryFisher.elim.ks.mf     <- runTest(queryGOdata.mf, algorithm = "elim",    statistic = "ks"    )

queryFisher.classic.fisher.cc <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "fisher")
queryFisher.classic.ks.cc     <- runTest(queryGOdata.cc, algorithm = "classic", statistic = "ks"    )
queryFisher.elim.ks.cc     <- runTest(queryGOdata.cc, algorithm = "elim",    statistic = "ks"    )

allRes.bp <- GenTable( queryGOdata.bp, classicFisher = queryFisher.classic.fisher, classicKS = queryFisher.classic.ks, elimKS = queryFisher.elim.ks, 
                    orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.bp)) )          

allRes.mf <- GenTable( queryGOdata.bp, classicFisher = queryFisher.classic.fisher, classicKS = queryFisher.classic.ks, elimKS = queryFisher.elim.ks, 
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.mf)) )          

allRes.cc <- GenTable( queryGOdata.bp, classicFisher = queryFisher.classic.fisher, classicKS = queryFisher.classic.ks, elimKS = queryFisher.elim.ks, 
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes=numNodes(graph(queryGOdata.cc)) )          
 
write.table(allRes.bp, file = "GO_table.BP.tsv", sep = '\t')
write.table(allRes.mf, file = "GO_table.MF.tsv", sep = '\t')
write.table(allRes.cc, file = "GO_table.CC.tsv", sep = '\t')
 
pValue.classic.bp <- score(queryFisher.classic.ks.bp)
pValue.elim.bp <- score(queryFisher.elim.ks.bp)[names(score(queryFisher.classic.ks.bp))]
gstat.bp <- termStat(queryGOdata.bp, names(pValue.classic.bp))
gSize.bp <- gstat.bp$Annotated / max(gstat.bp$Annotated) * 4

png( "pValue.classic.BP.png" , height=1000,width=1000 )
plot(pValue.classic.bp, pValue.elim.bp, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.bp, 
     col = redgreen(max(gstat.bp$Annotated))[ gstat.bp$Annotated ] )#gCol)
dev.off()

png( "LOGpValue.classic.BP.png" , height=1000,width=1000 )
plot(log(pValue.classic.bp), log(pValue.elim.bp), xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize.bp, 
     col = redgreen(max(gstat.bp$Annotated))[ gstat.bp$Annotated ] )#gCol)

sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
if(length(sel.go) > 0){
  write.table( cbind(termStat(queryGOdata, sel.go), elim = pValue.elim[sel.go], classic = pValue.classic[sel.go]) ,
               file ="highlightedKSoutlier.BP.tsv", sep = '\t')
}else{
  print("sel.go empty.")
}

png( "GO_tree.BP.png" , height=1000,width=1000 )
if(sum(allRes.bp$classicKS < 0.01) > 5){
  showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks), firstSigNodes = sum(allRes.bp$classicKS < 0.01), useInfo = 'all')
  title(paste0("Biological Processes\nCutoff p=0.01"))
}else if(sum(allRes$classicKS < 0.05) > 5 & sum(allRes$classicKS < 0.05) < 35){
  showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks), firstSigNodes = sum(allRes.bp$classicKS < 0.05), useInfo = 'all')
  title(paste0("Biological Processes\nCutoff p=0.05"))
}else{
  showSigOfNodes(queryGOdata.bp, score(queryFisher.elim.ks), firstSigNodes = 10, useInfo = 'all')
  title(paste0("Biological Processes\nCutoff count = 10 (no GO terms with p < 0.05)"))
}
dev.off() 


















queryGOdata.mf <- new( "topGOdata", 
                       description = "Simple session",
                       ontology = "MF",
                       allGenes = queryTest,
                       geneSel = topDiffGenes,
                       nodeSize = 2, # <- limit to annotations held over > 10 genes
                       annot=annFUN.gene2GO,
                       gene2GO=geneID2GO
)
queryGOdata.cc <- new( "topGOdata", 
                       description = "Simple session",
                       ontology = "CC",
                       allGenes = queryTest,
                       geneSel = topDiffGenes,
                       nodeSize = 2, # <- limit to annotations held over > 10 genes
                       annot=annFUN.gene2GO,
                       gene2GO=geneID2GO
)
                    
#
# results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'go_id'), 
#                  filter = "hgnc_symbol",
#                  values = names(queryTest), 
#                  mart = bm)


neg.go <- new("topGOdata", 
              desription="11.f_05.CAPER.p_val.FDR___neg",
              allGenes = c("AIG1","DYNLT1","EZR","GNS","IKBKE","LEMD3","MSRB3","NMBR","PEX3","RASSF3","RPL18A","RXYLT1","TBK1","WIF1","XPOT"),
              ontology="BP",
              ID="symbol"
             )
 