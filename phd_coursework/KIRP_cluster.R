
#### MUTATION INPUT ###################

setwd ("/Users/apple/Desktop/BENG283.omics/FINAL_PROJECT/KIRP")
mut_tab <- read.table ("mutation_matrix.top29genes.tsv", sep = "\t", header = T, comment.char = "", quote = "\"", stringsAsFactors = F)
mut_tab <- mut_tab[,-84]
rownames(mut_tab) <- mut_tab[,1]
mut_tab <- mut_tab[,-1]
meth_tab <- read.table ("KIRP-TP.normalized.top29genes.txt", sep = "\t", header = T, comment.char = "", quote = "\"", stringsAsFactors = F)
rownames(meth_tab) <- meth_tab[,1]
meth_tab <- meth_tab[,-1]
meth_tab <- meth_tab[,-1]
# meth_tab <- t(meth_tab)

# MUTSIG - KMEANS ######################

# curate clustering
kc2 <- kmeans(mut_tab, 2)
kc6 <- kmeans(mut_tab, 6)

# pretty plot
library(cluster)
data(mut_tab)
dissE <- daisy(mut_tab) 
dE2   <- dissE^2
sk2   <- silhouette(kc6$cl, dE2)
plot(sk2, col = c("red", "orange", "yellow", "light blue", "blue", "green"))
sk2   <- silhouette(kc2$cl, dE2)
plot(sk2, col = c("red", "blue"))

# MUTSIG - HIERARCHICAL ##################

dist.mut_tab <- dist(mut_tab, method = "binary")
dist.mut_tab.single.link <- hclust(dist.mut_tab, method='ward.D')

# pretty plot
plot(dist.mut_tab.single.link)
groups <- cutree(dist.mut_tab.single.link, k = 2)
rect.hclust(dist.mut_tab.single.link, k = 2, border = "blue")

plot(dist.mut_tab.single.link)
groups <- cutree(dist.mut_tab.single.link, k = 6)
rect.hclust(dist.mut_tab.single.link, k = 6, border = "purple")

source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
A2Rplot(dist.mut_tab.single.link, k =6, boxes = FALSE, col.up = "gray50", col.down = c("red", "orange", "purple", "blue", "dark green", "#556270"))

# MUTSIG - HIERARCHICAL ##################

dist.meth_tab <- dist(meth_tab, method = "euclidean")
dist.meth_tab.single.link <- hclust(dist.meth_tab, method='ward.D')

# pretty plot
plot(dist.meth_tab.single.link)
groups <- cutree(dist.meth_tab.single.link, k = 2)
rect.hclust(dist.meth_tab.single.link, k = 2, border = "red")

plot(dist.meth_tab.single.link)
groups <- cutree(dist.meth_tab.single.link, k = 6)
rect.hclust(dist.meth_tab.single.link, k = 6, border = "orange")

source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
A2Rplot(dist.meth_tab.single.link, k =6, boxes = FALSE, col.up = "gray50", col.down = c("red", "orange", "purple", "blue", "dark green", "#556270"))

# NMF - METHYLATION #######################

library(NMF)
object <- nmfObject(meth_tab)
res <- nmf(object, 2)
coefmap(res)


