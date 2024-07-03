
setwd ("/Users/apple/Desktop/BENG283.omics/FINAL_PROJECT/")

mut_tab <- read.table ("mutation_matrix.top29genes.tsv", sep = "\t", header = T, comment.char = "", quote = "\"", stringsAsFactors = F)
# mut_tab <- read.table ("KICH-TP.normalized.top29genes.tsv", sep = "\t", header = T, comment.char = "", quote = "\"", stringsAsFactors = F)
mut_tab <- t(mut_tab)
colnames(mut_tab) <- mut_tab[1,]

# delete empty columns
mut_tab <- mut_tab[-1,]
# mut_tab <- mut_tab[,-1]
# mut_tab <- mut_tab[,-1]
mut_tab <- mut_tab[-37,]

# install.packages("NMF")
# Load
library(NMF)

# remove zero rows
ind <- rowSums(mut_tab == 0) != ncol(mut_tab)
mut_tab <- mut_tab[ind, ]

# res <- nmf(mat, 2:10, nrun = 200, seed = 123456)
object <- nmfObject(mut_tab)
res <- nmf(object, 2)
coefmap(res)

#plot(res)
#groups <- cutree(res, k = 2)
#rect.hclust(res, k = 2, border = "blue")

