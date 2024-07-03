
#### MUTATION INPUT ###################

setwd ("/Users/apple/Desktop/BENG283.omics/FINAL_PROJECT/KIPAN")
mut_tab <- read.table ("mutation_matrix.top179genes.tsv", sep = "\t", header = T, comment.char = "", quote = "\"", stringsAsFactors = F)
mut_tab <- mut_tab[,-84]
rownames(mut_tab) <- mut_tab[,1]
mut_tab <- mut_tab[,-1]
mut_tab <- mut_tab[,-1]

meth_tab <- read.table ("methyl_cleanup.top179genes.tsv", sep = "\t", header = T, comment.char = "", quote = "\"", stringsAsFactors = F)
rownames(meth_tab) <- meth_tab[,1]
meth_tab <- meth_tab[,-1]
meth_tab <- meth_tab[,-1]
meth_tab <- meth_tab[,-1]
meth_tab <- t(meth_tab)
 
# MUTSIG - KMEANS ######################

# curate clustering
# kc2 <- kmeans(mut_tab, 2)
# kc6 <- kmeans(mut_tab, 6)

# pretty plot
# library(cluster)
# data(mut_tab)
# dissE <- daisy(mut_tab) 
# dE2   <- dissE^2
# sk2   <- silhouette(kc6$cl, dE2)
# plot(sk2, col = c("red", "orange", "yellow", "light blue", "blue", "green"), main = "MUTATION KMEANS k=6 (KIPAN)")
# sk2   <- silhouette(kc2$cl, dE2)
# plot(sk2, col = c("red", "blue"), main = "MUTATION KMEANS k=2 (KIPAN)")

# MUTSIG - HIERARCHICAL ##################

dist.mut_tab <- dist(mut_tab, method = "binary")
dist.mut_tab.single.link <- hclust(dist.mut_tab, method='ward.D')

# pretty plot
plot(dist.mut_tab.single.link, main = "Mutation Hierarchical Clustering (KIPAN, k=2)")
groups <- cutree(dist.mut_tab.single.link, k = 2)
rect.hclust(dist.mut_tab.single.link, k = 2, border = "blue")

plot(dist.mut_tab.single.link, main = "Mutation Hierarchical Clustering (KIPAN, k=6)")
groups <- cutree(dist.mut_tab.single.link, k = 6)
rect.hclust(dist.mut_tab.single.link, k = 6, border = "purple")

source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
A2Rplot(dist.mut_tab.single.link, k =6, boxes = FALSE, col.up = "gray50", col.down = c("red", "orange", "purple", "blue", "dark green", "#556270"), main = "Mutation Hierarchical Clustering (KIPAN, k=6)")

# METHYLATION - HIERARCHICAL ##################

dist.meth_tab <- dist(meth_tab, method = "euclidean")
dist.meth_tab.single.link <- hclust(dist.meth_tab, method='ward.D')

# pretty plot
plot(dist.meth_tab.single.link, main = "Methylation Hierarchical Clustering (KIPAN, k=2)")
groups <- cutree(dist.meth_tab.single.link, k = 2)
rect.hclust(dist.meth_tab.single.link, k = 2, border = "red")

plot(dist.meth_tab.single.link, main = "Methylation Hierarchical Clustering (KIPAN, k=6)")
groups <- cutree(dist.meth_tab.single.link, k = 6)
rect.hclust(dist.meth_tab.single.link, k = 6, border = "orange")

source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
A2Rplot(dist.meth_tab.single.link, k =6, boxes = FALSE, col.up = "gray50", col.down = c("red", "orange", "purple", "blue", "dark green", "#556270"), main = "Methylation Hierarchical Clustering (KIPAN, k=6)")

# NMF - METHYLATION #######################

library(NMF)
object <- nmfObject(meth_tab)
res <- nmf(object, 2)
coefmap(res, main = "Methylation NMF Clustering (KIPAN, k=4)")


setwd ("/Users/apple/Desktop/BENG283.omics/FINAL_PROJECT/KIPAN")
# cnv_tab <- read.table ("Gistic_Raw.txt", sep = "\t", header = T, comment.char = "", quote = "\"", stringsAsFactors = F)
cnv_tab <- read.table ("gistic_cleanup.top179genes.tsv", sep = "\t", header = T, comment.char = "", quote = "\"", stringsAsFactors = F)
cnv_tab <- cnv_tab[,-2]
cnv_tab <- cnv_tab[,-2]
rownames(cnv_tab) <- cnv_tab$Gene.Symbol
cnv_tab <- t(cnv_tab)
colnames(cnv_tab) <- cnv_tab[1,]
cnv_tab <- cnv_tab[-1,]
# hierarchical

dist.cnv_tab <- dist(cnv_tab, method = "euclidean")
dist.cnv_tab.single.link <- hclust(dist.cnv_tab, method='ward.D')

# pretty plot
plot(dist.cnv_tab.single.link, main = "Gistic Hierarchical Clustering (KIRP, k=4)")
groups <- cutree(dist.cnv_tab.single.link, k = 4)
rect.hclust(dist.cnv_tab.single.link, k = 4, border = "red")

source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
A2Rplot(dist.cnv_tab.single.link, k =4, boxes = FALSE, col.up = "gray50", col.down = c("red", "orange", "purple", "blue", "dark green", "#556270"), main = "Gistic Hierarchical Clustering (KIPAN, k=4)")

# nmf

library(NMF)
object <- nmfObject(cnv_tab)
res <- nmf(object, 2)
coefmap(res, main = "Gistic NMF Clustering (KIPAN, k=4)")

dev.off()
