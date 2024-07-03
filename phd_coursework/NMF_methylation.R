
setwd ("/Users/apple/Desktop/BENG283.omics/FINAL_PROJECT/")

# meth_tab <- read.table ("mutation_matrix.top29genes.tsv", sep = "\t", header = T, comment.char = "", quote = "\"", stringsAsFactors = F)
meth_tab <- read.table ("KICH-TP.normalized.top29genes.tsv", sep = "\t", header = T, comment.char = "", quote = "\"", stringsAsFactors = F)
rownames(meth_tab) <- meth_tab[,1]

# delete empty columns
meth_tab <- meth_tab[,-1]
# meth_tab <- meth_tab[,-1]
# meth_tab <- meth_tab[,-37]

# install.packages("NMF")
# Load
library(NMF)

# res <- nmf(mat, 2:10, nrun = 200, seed = 123456)

object <- nmfObject(meth_tab)
res <- nmf(object, 2)
coefmap(res)

#plot(res)
#groups <- cutree(res, k = 2)
#rect.hclust(res, k = 2, border = "blue")

