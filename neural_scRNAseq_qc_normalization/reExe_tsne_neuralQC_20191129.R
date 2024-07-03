
setwd("C:\\Users\\jamis\\Desktop\\neuralQC_20191108_weekend")
library(vegan)
library(umap)
library(Rtsne)
library(colorRamps)
 
# The basic data is contained in E_limp_ and I_limp_ for the exon and intron data, respectively.
# This data has been "quality-controlled" (by me), removing some samples with low transcript abundances and some genes with many 0-measurements.
# The resulting abundances have been log-transformed, and then the NANs (i.e., the log(0) entries) have been imputed using a rank-6 pca-based method.
# This process makes many of the same underlying assumptions used in the "MAST" model: namely, if transcripts are 'missing' (via a 0-measurment), we 
# fill them in assuming a maximum-likelihood gaussian model.

### import count tables
# exons
E_limp <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_limp_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
E_licl <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licl_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
E_licr <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licr_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
# introns
I_limp <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_limp_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
I_licl <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licl_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
I_licr <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licr_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
# remove all zero columns (just 1 each)
E_limp <- E_limp[ , colSums(is.na(E_limp)) != dim(E_limp)[1] ]
E_licl <- E_licl[ , colSums(is.na(E_licl)) != dim(E_licl)[1] ]
E_licr <- E_licr[ , colSums(is.na(E_licr)) != dim(E_licr)[1] ]
I_limp <- I_limp[ , colSums(is.na(I_limp)) != dim(I_limp)[1] ]
I_licl <- I_licl[ , colSums(is.na(I_licl)) != dim(I_licl)[1] ]
I_licr <- I_licr[ , colSums(is.na(I_licr)) != dim(I_licr)[1] ]
 
#save.image("base_import.R")
 
# I would be very interested to see what t-sne does to the I_corrected_ data! Also, you can stack the two matrices [ E_corrected , I_corrected_ ] next 
# to one another and try that too.

### 12 covariates
# exons
zeta_limp_E_un_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_limp_c012_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_limp_E_vn_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_limp_c012_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_E_un_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licl_c012_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_E_vn_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licl_c012_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_E_un_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licr_c012_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_E_vn_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licr_c012_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
# introns
zeta_limp_I_un_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_limp_c012_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_limp_I_vn_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_limp_c012_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_I_un_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licl_c012_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_I_vn_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licl_c012_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_I_un_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licr_c012_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_I_vn_c012   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licr_c012_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)

### 42 covariates
# exons
zeta_limp_E_un_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_limp_c042_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_limp_E_vn_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_limp_c042_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_E_un_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licl_c042_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_E_vn_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licl_c042_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_E_un_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licr_c042_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_E_vn_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licr_c042_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
# introns
zeta_limp_I_un_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_limp_c042_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_limp_I_vn_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_limp_c042_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_I_un_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licl_c042_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_I_vn_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licl_c042_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_I_un_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licr_c042_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_I_vn_c042   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licr_c042_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)

### 127 covariates
# exons
zeta_limp_E_un_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_limp_c127_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_limp_E_vn_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_limp_c127_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_E_un_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licl_c127_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_E_vn_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licl_c127_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_E_un_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licr_c127_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_E_vn_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/E_licr_c127_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
# introns
zeta_limp_I_un_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_limp_c127_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_limp_I_vn_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_limp_c127_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_I_un_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licl_c127_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licl_I_vn_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licl_c127_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_I_un_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licr_c127_zeta_un_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
zeta_licr_I_vn_c127   <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/I_licr_c127_zeta_vn_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)

# Covarite Ranks
C_rank_c127 <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/C_rank_c127_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
C_rank_c042 <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/C_rank_c042_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
C_rank_c012 <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/C_rank_c012_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)

# Sample and Cluster IDs
u_ID_ <-         as.list(read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/u_ID_sub_.nsv",       sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE))
u_CLabel_sub_ <- as.list(read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/str_CLabel_sub_.nsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE))

### run tsne and plot
colorbycluster <- primary.colors(46)[ as.numeric(unlist(u_CLabel_sub_)) ]
# exons
E_limp.tsne <- Rtsne(E_limp, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne <- Rtsne(E_licl, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne <- Rtsne(E_licr, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
png("20191109.E_limp.tsne.png", width=800,height=800)
plot(E_limp.tsne$Y, main="TSNE: Exon", col=colorbycluster,pch=16)
dev.off()
png("20191109.E_licl.tsne.png", width=800,height=800)
plot(E_licl.tsne$Y, main="TSNE: Exon", col=colorbycluster,pch=16)
dev.off()
png("20191109.E_licr.tsne.png", width=800,height=800)
plot(E_licr.tsne$Y, main="TSNE: Exon", col=colorbycluster,pch=16)
dev.off()
# introns
I_limp.tsne <- Rtsne(I_limp, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne <- Rtsne(I_licl, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne <- Rtsne(I_licr, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
png("20191109.I_limp.tsne.png", width=800,height=800)
plot(I_limp.tsne$Y, main="TSNE: Intron", col=colorbycluster,pch=16)
dev.off()
png("20191109.I_licl.tsne.png", width=800,height=800)
plot(I_licl.tsne$Y, main="TSNE: Intron", col=colorbycluster,pch=16)
dev.off()
png("20191109.I_licr.tsne.png", width=800,height=800)
plot(I_licr.tsne$Y, main="TSNE: Intron", col=colorbycluster,pch=16)
dev.off()

save.image("20191109.beforeCorrection.Rdata")

###########################################################################################################################################################

# Anyway, the linear model linking the covariates C_rank_ to the data E_limp_ is given by:
# E_lm_limp_ = [ (ones(n_u,1)*zeta_limp_E_un_(1,:) + C_rank_*zeta_limp_E_un_(2:end,:))*transpose(zeta_limp_E_vn_) ] ;

### 12 covariates
zz <- as.data.frame(matrix(ncol=1,nrow=dim(C_rank_c012)[1]))
zz[1:dim(C_rank_c012)[1],1] <- 1
# 1 column
E_lm_limp_c012.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c012[1,1]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_E_un_c012[2:dim(zeta_limp_E_un_c012)[1],1] ) )) %*% t(as.matrix(zeta_limp_E_vn_c012[,1]))
E_lm_licl_c012.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c012[1,1]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licl_E_un_c012)[1],1] ) )) %*% t(as.matrix(zeta_licl_E_vn_c012[,1]))
E_lm_licr_c012.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c012[1,1]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licr_E_un_c012)[1],1] ) )) %*% t(as.matrix(zeta_licr_E_vn_c012[,1]))
# 2 column
E_lm_limp_c012.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c012[1,1:2]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_E_un_c012[2:dim(zeta_limp_E_un_c012)[1],1:2] ) )) %*% t(as.matrix(zeta_limp_E_vn_c012[,1:2]))
E_lm_licl_c012.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c012[1,1:2]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licl_E_un_c012)[1],1:2] ) )) %*% t(as.matrix(zeta_licl_E_vn_c012[,1:2]))
E_lm_licr_c012.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c012[1,1:2]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licr_E_un_c012)[1],1:2] ) )) %*% t(as.matrix(zeta_licr_E_vn_c012[,1:2]))
# 3 column
E_lm_limp_c012.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c012[1,1:3]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_E_un_c012[2:dim(zeta_limp_E_un_c012)[1],1:3] ) )) %*% t(as.matrix(zeta_limp_E_vn_c012[,1:3]))
E_lm_licl_c012.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c012[1,1:3]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licl_E_un_c012)[1],1:3] ) )) %*% t(as.matrix(zeta_licl_E_vn_c012[,1:3]))
E_lm_licr_c012.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c012[1,1:3]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licr_E_un_c012)[1],1:3] ) )) %*% t(as.matrix(zeta_licr_E_vn_c012[,1:3]))
# 4 column
E_lm_limp_c012.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c012[1,1:4]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_E_un_c012[2:dim(zeta_limp_E_un_c012)[1],1:4] ) )) %*% t(as.matrix(zeta_limp_E_vn_c012[,1:4]))
E_lm_licl_c012.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c012[1,1:4]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licl_E_un_c012)[1],1:4] ) )) %*% t(as.matrix(zeta_licl_E_vn_c012[,1:4]))
E_lm_licr_c012.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c012[1,1:4]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licr_E_un_c012)[1],1:4] ) )) %*% t(as.matrix(zeta_licr_E_vn_c012[,1:4]))
# 5 column
E_lm_limp_c012.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c012[1,1:5]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_E_un_c012[2:dim(zeta_limp_E_un_c012)[1],1:5] ) )) %*% t(as.matrix(zeta_limp_E_vn_c012[,1:5]))
E_lm_licl_c012.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c012[1,1:5]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licl_E_un_c012)[1],1:5] ) )) %*% t(as.matrix(zeta_licl_E_vn_c012[,1:5]))
E_lm_licr_c012.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c012[1,1:5]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licr_E_un_c012)[1],1:5] ) )) %*% t(as.matrix(zeta_licr_E_vn_c012[,1:5]))
# 6 column
E_lm_limp_c012.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c012[1,1:6]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_E_un_c012[2:dim(zeta_limp_E_un_c012)[1],1:6] ) )) %*% t(as.matrix(zeta_limp_E_vn_c012[,1:6]))
E_lm_licl_c012.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c012[1,1:6]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licl_E_un_c012)[1],1:6] ) )) %*% t(as.matrix(zeta_licl_E_vn_c012[,1:6]))
E_lm_licr_c012.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c012[1,1:6]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_E_un_c012[2:dim(zeta_licr_E_un_c012)[1],1:6] ) )) %*% t(as.matrix(zeta_licr_E_vn_c012[,1:6]))


zz <- as.data.frame(matrix(ncol=1,nrow=dim(C_rank_c042)[1]))
zz[1:dim(C_rank_c042)[1],1] <- 1
# 1 column
E_lm_limp_c042.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c042[1,1]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_E_un_c042[2:dim(zeta_limp_E_un_c042)[1],1] ) )) %*% t(as.matrix(zeta_limp_E_vn_c042[,1]))
E_lm_licl_c042.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c042[1,1]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licl_E_un_c042)[1],1] ) )) %*% t(as.matrix(zeta_licl_E_vn_c042[,1]))
E_lm_licr_c042.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c042[1,1]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licr_E_un_c042)[1],1] ) )) %*% t(as.matrix(zeta_licr_E_vn_c042[,1]))
# 2 column
E_lm_limp_c042.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c042[1,1:2]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_E_un_c042[2:dim(zeta_limp_E_un_c042)[1],1:2] ) )) %*% t(as.matrix(zeta_limp_E_vn_c042[,1:2]))
E_lm_licl_c042.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c042[1,1:2]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licl_E_un_c042)[1],1:2] ) )) %*% t(as.matrix(zeta_licl_E_vn_c042[,1:2]))
E_lm_licr_c042.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c042[1,1:2]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licr_E_un_c042)[1],1:2] ) )) %*% t(as.matrix(zeta_licr_E_vn_c042[,1:2]))
# 3 column
E_lm_limp_c042.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c042[1,1:3]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_E_un_c042[2:dim(zeta_limp_E_un_c042)[1],1:3] ) )) %*% t(as.matrix(zeta_limp_E_vn_c042[,1:3]))
E_lm_licl_c042.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c042[1,1:3]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licl_E_un_c042)[1],1:3] ) )) %*% t(as.matrix(zeta_licl_E_vn_c042[,1:3]))
E_lm_licr_c042.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c042[1,1:3]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licr_E_un_c042)[1],1:3] ) )) %*% t(as.matrix(zeta_licr_E_vn_c042[,1:3]))
# 4 column
E_lm_limp_c042.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c042[1,1:4]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_E_un_c042[2:dim(zeta_limp_E_un_c042)[1],1:4] ) )) %*% t(as.matrix(zeta_limp_E_vn_c042[,1:4]))
E_lm_licl_c042.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c042[1,1:4]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licl_E_un_c042)[1],1:4] ) )) %*% t(as.matrix(zeta_licl_E_vn_c042[,1:4]))
E_lm_licr_c042.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c042[1,1:4]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licr_E_un_c042)[1],1:4] ) )) %*% t(as.matrix(zeta_licr_E_vn_c042[,1:4]))
# 5 column
E_lm_limp_c042.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c042[1,1:5]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_E_un_c042[2:dim(zeta_limp_E_un_c042)[1],1:5] ) )) %*% t(as.matrix(zeta_limp_E_vn_c042[,1:5]))
E_lm_licl_c042.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c042[1,1:5]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licl_E_un_c042)[1],1:5] ) )) %*% t(as.matrix(zeta_licl_E_vn_c042[,1:5]))
E_lm_licr_c042.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c042[1,1:5]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licr_E_un_c042)[1],1:5] ) )) %*% t(as.matrix(zeta_licr_E_vn_c042[,1:5]))
# 6 column
E_lm_limp_c042.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c042[1,1:6]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_E_un_c042[2:dim(zeta_limp_E_un_c042)[1],1:6] ) )) %*% t(as.matrix(zeta_limp_E_vn_c042[,1:6]))
E_lm_licl_c042.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c042[1,1:6]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licl_E_un_c042)[1],1:6] ) )) %*% t(as.matrix(zeta_licl_E_vn_c042[,1:6]))
E_lm_licr_c042.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c042[1,1:6]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_E_un_c042[2:dim(zeta_licr_E_un_c042)[1],1:6] ) )) %*% t(as.matrix(zeta_licr_E_vn_c042[,1:6]))


zz <- as.data.frame(matrix(ncol=1,nrow=dim(C_rank_c127)[1]))
zz[1:dim(C_rank_c127)[1],1] <- 1
# 1 column
E_lm_limp_c127.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c127[1,1]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_E_un_c127[2:dim(zeta_limp_E_un_c127)[1],1] ) )) %*% t(as.matrix(zeta_limp_E_vn_c127[,1]))
E_lm_licl_c127.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c127[1,1]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licl_E_un_c127)[1],1] ) )) %*% t(as.matrix(zeta_licl_E_vn_c127[,1]))
E_lm_licr_c127.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c127[1,1]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licr_E_un_c127)[1],1] ) )) %*% t(as.matrix(zeta_licr_E_vn_c127[,1]))
# 2 column
E_lm_limp_c127.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c127[1,1:2]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_E_un_c127[2:dim(zeta_limp_E_un_c127)[1],1:2] ) )) %*% t(as.matrix(zeta_limp_E_vn_c127[,1:2]))
E_lm_licl_c127.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c127[1,1:2]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licl_E_un_c127)[1],1:2] ) )) %*% t(as.matrix(zeta_licl_E_vn_c127[,1:2]))
E_lm_licr_c127.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c127[1,1:2]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licr_E_un_c127)[1],1:2] ) )) %*% t(as.matrix(zeta_licr_E_vn_c127[,1:2]))
# 3 column
E_lm_limp_c127.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c127[1,1:3]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_E_un_c127[2:dim(zeta_limp_E_un_c127)[1],1:3] ) )) %*% t(as.matrix(zeta_limp_E_vn_c127[,1:3]))
E_lm_licl_c127.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c127[1,1:3]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licl_E_un_c127)[1],1:3] ) )) %*% t(as.matrix(zeta_licl_E_vn_c127[,1:3]))
E_lm_licr_c127.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c127[1,1:3]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licr_E_un_c127)[1],1:3] ) )) %*% t(as.matrix(zeta_licr_E_vn_c127[,1:3]))
# 4 column
E_lm_limp_c127.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c127[1,1:4]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_E_un_c127[2:dim(zeta_limp_E_un_c127)[1],1:4] ) )) %*% t(as.matrix(zeta_limp_E_vn_c127[,1:4]))
E_lm_licl_c127.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c127[1,1:4]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licl_E_un_c127)[1],1:4] ) )) %*% t(as.matrix(zeta_licl_E_vn_c127[,1:4]))
E_lm_licr_c127.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c127[1,1:4]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licr_E_un_c127)[1],1:4] ) )) %*% t(as.matrix(zeta_licr_E_vn_c127[,1:4]))
# 5 column
E_lm_limp_c127.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c127[1,1:5]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_E_un_c127[2:dim(zeta_limp_E_un_c127)[1],1:5] ) )) %*% t(as.matrix(zeta_limp_E_vn_c127[,1:5]))
E_lm_licl_c127.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c127[1,1:5]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licl_E_un_c127)[1],1:5] ) )) %*% t(as.matrix(zeta_licl_E_vn_c127[,1:5]))
E_lm_licr_c127.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c127[1,1:5]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licr_E_un_c127)[1],1:5] ) )) %*% t(as.matrix(zeta_licr_E_vn_c127[,1:5]))
# 6 column
E_lm_limp_c127.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_E_un_c127[1,1:6]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_E_un_c127[2:dim(zeta_limp_E_un_c127)[1],1:6] ) )) %*% t(as.matrix(zeta_limp_E_vn_c127[,1:6]))
E_lm_licl_c127.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_E_un_c127[1,1:6]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licl_E_un_c127)[1],1:6] ) )) %*% t(as.matrix(zeta_licl_E_vn_c127[,1:6]))
E_lm_licr_c127.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_E_un_c127[1,1:6]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_E_un_c127[2:dim(zeta_licr_E_un_c127)[1],1:6] ) )) %*% t(as.matrix(zeta_licr_E_vn_c127[,1:6]))
 
# and the difference:
#   E_corrected_ = E_limp_ - E_lm_limp_
# would represent 'covariate-corrected' data, which could then be put into t-sne.
# 1 column
E_limp.corrected_c012.col1 <- E_limp - E_lm_limp_c012.col1
E_licl.corrected_c012.col1 <- E_licl - E_lm_licl_c012.col1
E_licr.corrected_c012.col1 <- E_licr - E_lm_licr_c012.col1
# 2 column
E_limp.corrected_c012.col2 <- E_limp - E_lm_limp_c012.col2
E_licl.corrected_c012.col2 <- E_licl - E_lm_licl_c012.col2
E_licr.corrected_c012.col2 <- E_licr - E_lm_licr_c012.col2
# 3 column
E_limp.corrected_c012.col3 <- E_limp - E_lm_limp_c012.col3
E_licl.corrected_c012.col3 <- E_licl - E_lm_licl_c012.col3
E_licr.corrected_c012.col3 <- E_licr - E_lm_licr_c012.col3
# 3 column
E_limp.corrected_c012.col4 <- E_limp - E_lm_limp_c012.col4
E_licl.corrected_c012.col4 <- E_licl - E_lm_licl_c012.col4
E_licr.corrected_c012.col4 <- E_licr - E_lm_licr_c012.col4
# 3 column
E_limp.corrected_c012.col5 <- E_limp - E_lm_limp_c012.col5
E_licl.corrected_c012.col5 <- E_licl - E_lm_licl_c012.col5
E_licr.corrected_c012.col5 <- E_licr - E_lm_licr_c012.col5
# 3 column
E_limp.corrected_c012.col6 <- E_limp - E_lm_limp_c012.col6
E_licl.corrected_c012.col6 <- E_licl - E_lm_licl_c012.col6
E_licr.corrected_c012.col6 <- E_licr - E_lm_licr_c012.col6

# 1 column
E_limp.corrected_c042.col1 <- E_limp - E_lm_limp_c042.col1
E_licl.corrected_c042.col1 <- E_licl - E_lm_licl_c042.col1
E_licr.corrected_c042.col1 <- E_licr - E_lm_licr_c042.col1
# 2 column
E_limp.corrected_c042.col2 <- E_limp - E_lm_limp_c042.col2
E_licl.corrected_c042.col2 <- E_licl - E_lm_licl_c042.col2
E_licr.corrected_c042.col2 <- E_licr - E_lm_licr_c042.col2
# 3 column
E_limp.corrected_c042.col3 <- E_limp - E_lm_limp_c042.col3
E_licl.corrected_c042.col3 <- E_licl - E_lm_licl_c042.col3
E_licr.corrected_c042.col3 <- E_licr - E_lm_licr_c042.col3
# 3 column
E_limp.corrected_c042.col4 <- E_limp - E_lm_limp_c042.col4
E_licl.corrected_c042.col4 <- E_licl - E_lm_licl_c042.col4
E_licr.corrected_c042.col4 <- E_licr - E_lm_licr_c042.col4
# 3 column
E_limp.corrected_c042.col5 <- E_limp - E_lm_limp_c042.col5
E_licl.corrected_c042.col5 <- E_licl - E_lm_licl_c042.col5
E_licr.corrected_c042.col5 <- E_licr - E_lm_licr_c042.col5
# 3 column
E_limp.corrected_c042.col6 <- E_limp - E_lm_limp_c042.col6
E_licl.corrected_c042.col6 <- E_licl - E_lm_licl_c042.col6
E_licr.corrected_c042.col6 <- E_licr - E_lm_licr_c042.col6

# 1 column
E_limp.corrected_c042.col1 <- E_limp - E_lm_limp_c042.col1
E_licl.corrected_c042.col1 <- E_licl - E_lm_licl_c042.col1
E_licr.corrected_c042.col1 <- E_licr - E_lm_licr_c042.col1
# 2 column
E_limp.corrected_c042.col2 <- E_limp - E_lm_limp_c042.col2
E_licl.corrected_c042.col2 <- E_licl - E_lm_licl_c042.col2
E_licr.corrected_c042.col2 <- E_licr - E_lm_licr_c042.col2
# 3 column
E_limp.corrected_c042.col3 <- E_limp - E_lm_limp_c042.col3
E_licl.corrected_c042.col3 <- E_licl - E_lm_licl_c042.col3
E_licr.corrected_c042.col3 <- E_licr - E_lm_licr_c042.col3
# 3 column
E_limp.corrected_c042.col4 <- E_limp - E_lm_limp_c042.col4
E_licl.corrected_c042.col4 <- E_licl - E_lm_licl_c042.col4
E_licr.corrected_c042.col4 <- E_licr - E_lm_licr_c042.col4
# 3 column
E_limp.corrected_c042.col5 <- E_limp - E_lm_limp_c042.col5
E_licl.corrected_c042.col5 <- E_licl - E_lm_licl_c042.col5
E_licr.corrected_c042.col5 <- E_licr - E_lm_licr_c042.col5
# 3 column
E_limp.corrected_c042.col6 <- E_limp - E_lm_limp_c042.col6
E_licl.corrected_c042.col6 <- E_licl - E_lm_licl_c042.col6
E_licr.corrected_c042.col6 <- E_licr - E_lm_licr_c042.col6

# 1 column
E_limp.corrected_c127.col1 <- E_limp - E_lm_limp_c127.col1
E_licl.corrected_c127.col1 <- E_licl - E_lm_licl_c127.col1
E_licr.corrected_c127.col1 <- E_licr - E_lm_licr_c127.col1
# 2 column
E_limp.corrected_c127.col2 <- E_limp - E_lm_limp_c127.col2
E_licl.corrected_c127.col2 <- E_licl - E_lm_licl_c127.col2
E_licr.corrected_c127.col2 <- E_licr - E_lm_licr_c127.col2
# 3 column
E_limp.corrected_c127.col3 <- E_limp - E_lm_limp_c127.col3
E_licl.corrected_c127.col3 <- E_licl - E_lm_licl_c127.col3
E_licr.corrected_c127.col3 <- E_licr - E_lm_licr_c127.col3
# 3 column
E_limp.corrected_c127.col4 <- E_limp - E_lm_limp_c127.col4
E_licl.corrected_c127.col4 <- E_licl - E_lm_licl_c127.col4
E_licr.corrected_c127.col4 <- E_licr - E_lm_licr_c127.col4
# 3 column
E_limp.corrected_c127.col5 <- E_limp - E_lm_limp_c127.col5
E_licl.corrected_c127.col5 <- E_licl - E_lm_licl_c127.col5
E_licr.corrected_c127.col5 <- E_licr - E_lm_licr_c127.col5
# 3 column
E_limp.corrected_c127.col6 <- E_limp - E_lm_limp_c127.col6
E_licl.corrected_c127.col6 <- E_licl - E_lm_licl_c127.col6
E_licr.corrected_c127.col6 <- E_licr - E_lm_licr_c127.col6
  
E_limp.tsne.corrected_c012.col1 <- Rtsne(E_limp.corrected_c012.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c012.col1 <- Rtsne(E_licl.corrected_c012.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c012.col1 <- Rtsne(E_licr.corrected_c012.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c012.col2 <- Rtsne(E_limp.corrected_c012.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c012.col2 <- Rtsne(E_licl.corrected_c012.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c012.col2 <- Rtsne(E_licr.corrected_c012.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c012.col3 <- Rtsne(E_limp.corrected_c012.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c012.col3 <- Rtsne(E_licl.corrected_c012.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c012.col3 <- Rtsne(E_licr.corrected_c012.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c012.col4 <- Rtsne(E_limp.corrected_c012.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c012.col4 <- Rtsne(E_licl.corrected_c012.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c012.col4 <- Rtsne(E_licr.corrected_c012.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c012.col5 <- Rtsne(E_limp.corrected_c012.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c012.col5 <- Rtsne(E_licl.corrected_c012.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c012.col5 <- Rtsne(E_licr.corrected_c012.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c012.col6 <- Rtsne(E_limp.corrected_c012.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c012.col6 <- Rtsne(E_licl.corrected_c012.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c012.col6 <- Rtsne(E_licr.corrected_c012.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c042.col1 <- Rtsne(E_limp.corrected_c042.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c042.col1 <- Rtsne(E_licl.corrected_c042.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c042.col1 <- Rtsne(E_licr.corrected_c042.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c042.col2 <- Rtsne(E_limp.corrected_c042.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c042.col2 <- Rtsne(E_licl.corrected_c042.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c042.col2 <- Rtsne(E_licr.corrected_c042.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c042.col3 <- Rtsne(E_limp.corrected_c042.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c042.col3 <- Rtsne(E_licl.corrected_c042.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c042.col3 <- Rtsne(E_licr.corrected_c042.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c042.col4 <- Rtsne(E_limp.corrected_c042.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c042.col4 <- Rtsne(E_licl.corrected_c042.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c042.col4 <- Rtsne(E_licr.corrected_c042.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c042.col5 <- Rtsne(E_limp.corrected_c042.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c042.col5 <- Rtsne(E_licl.corrected_c042.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c042.col5 <- Rtsne(E_licr.corrected_c042.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c042.col6 <- Rtsne(E_limp.corrected_c042.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c042.col6 <- Rtsne(E_licl.corrected_c042.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c042.col6 <- Rtsne(E_licr.corrected_c042.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c127.col1 <- Rtsne(E_limp.corrected_c127.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c127.col1 <- Rtsne(E_licl.corrected_c127.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c127.col1 <- Rtsne(E_licr.corrected_c127.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c127.col2 <- Rtsne(E_limp.corrected_c127.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c127.col2 <- Rtsne(E_licl.corrected_c127.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c127.col2 <- Rtsne(E_licr.corrected_c127.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c127.col3 <- Rtsne(E_limp.corrected_c127.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c127.col3 <- Rtsne(E_licl.corrected_c127.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c127.col3 <- Rtsne(E_licr.corrected_c127.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c127.col4 <- Rtsne(E_limp.corrected_c127.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c127.col4 <- Rtsne(E_licl.corrected_c127.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c127.col4 <- Rtsne(E_licr.corrected_c127.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c127.col5 <- Rtsne(E_limp.corrected_c127.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c127.col5 <- Rtsne(E_licl.corrected_c127.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c127.col5 <- Rtsne(E_licr.corrected_c127.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

E_limp.tsne.corrected_c127.col6 <- Rtsne(E_limp.corrected_c127.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licl.tsne.corrected_c127.col6 <- Rtsne(E_licl.corrected_c127.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
E_licr.tsne.corrected_c127.col6 <- Rtsne(E_licr.corrected_c127.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)


png("E_limp.corrected_c012.col1.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c012.col1$Y, main="TSNE: Exon (c012.col1)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c012.col1.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c012.col1$Y, main="TSNE: Exon (c012.col1)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c012.col1.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c012.col1$Y, main="TSNE: Exon (c012.col1)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c012.col2.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c012.col2$Y, main="TSNE: Exon (c012.col2)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c012.col2.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c012.col2$Y, main="TSNE: Exon (c012.col2)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c012.col2.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c012.col2$Y, main="TSNE: Exon (c012.col2)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c012.col3.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c012.col3$Y, main="TSNE: Exon (c012.col3)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c012.col3.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c012.col3$Y, main="TSNE: Exon (c012.col3)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c012.col3.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c012.col3$Y, main="TSNE: Exon (c012.col3)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c012.col4.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c012.col4$Y, main="TSNE: Exon (c012.col4)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c012.col4.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c012.col4$Y, main="TSNE: Exon (c012.col4)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c012.col4.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c012.col4$Y, main="TSNE: Exon (c012.col4)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c012.col5.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c012.col5$Y, main="TSNE: Exon (c012.col5)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c012.col5.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c012.col5$Y, main="TSNE: Exon (c012.col5)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c012.col5.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c012.col5$Y, main="TSNE: Exon (c012.col5)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c012.col6.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c012.col6$Y, main="TSNE: Exon (c012.col6)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c012.col6.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c012.col6$Y, main="TSNE: Exon (c012.col6)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c012.col6.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c012.col6$Y, main="TSNE: Exon (c012.col6)", col=colorbycluster,pch=16)
dev.off()

png("E_limp.corrected_c012.col_compare.tsne.png", width=1200,height=2400)
par(mfrow=c(6,3))
plot(E_limp.tsne.corrected_c012.col1$Y, main="TSNE: LIMP Exon (c012.col1)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c012.col1$Y, main="TSNE: LICL Exon (c012.col1)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c012.col1$Y, main="TSNE: LICR Exon (c012.col1)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c012.col2$Y, main="TSNE: LIMP Exon (c012.col2)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c012.col2$Y, main="TSNE: LICL Exon (c012.col2)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c012.col2$Y, main="TSNE: LICR Exon (c012.col2)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c012.col3$Y, main="TSNE: LIMP Exon (c012.col3)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c012.col3$Y, main="TSNE: LICL Exon (c012.col3)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c012.col3$Y, main="TSNE: LICR Exon (c012.col3)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c012.col4$Y, main="TSNE: LIMP Exon (c012.col4)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c012.col4$Y, main="TSNE: LICL Exon (c012.col4)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c012.col4$Y, main="TSNE: LICR Exon (c012.col4)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c012.col5$Y, main="TSNE: LIMP Exon (c012.col5)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c012.col5$Y, main="TSNE: LICL Exon (c012.col5)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c012.col5$Y, main="TSNE: LICR Exon (c012.col5)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c012.col6$Y, main="TSNE: LIMP Exon (c012.col6)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c012.col6$Y, main="TSNE: LICL Exon (c012.col6)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c012.col6$Y, main="TSNE: LICR Exon (c012.col6)", col=colorbycluster,pch=16)
dev.off()
par(mfrow=c(1,1))

png("E_limp.corrected_c042.col1.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c042.col1$Y, main="TSNE: Exon (c042.col1)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c042.col1.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c042.col1$Y, main="TSNE: Exon (c042.col1)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c042.col1.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c042.col1$Y, main="TSNE: Exon (c042.col1)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c042.col2.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c042.col2$Y, main="TSNE: Exon (c042.col2)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c042.col2.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c042.col2$Y, main="TSNE: Exon (c042.col2)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c042.col2.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c042.col2$Y, main="TSNE: Exon (c042.col2)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c042.col3.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c042.col3$Y, main="TSNE: Exon (c042.col3)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c042.col3.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c042.col3$Y, main="TSNE: Exon (c042.col3)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c042.col3.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c042.col3$Y, main="TSNE: Exon (c042.col3)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c042.col4.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c042.col4$Y, main="TSNE: Exon (c042.col4)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c042.col4.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c042.col4$Y, main="TSNE: Exon (c042.col4)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c042.col4.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c042.col4$Y, main="TSNE: Exon (c042.col4)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c042.col5.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c042.col5$Y, main="TSNE: Exon (c042.col5)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c042.col5.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c042.col5$Y, main="TSNE: Exon (c042.col5)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c042.col5.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c042.col5$Y, main="TSNE: Exon (c042.col5)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c042.col6.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c042.col6$Y, main="TSNE: Exon (c042.col6)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c042.col6.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c042.col6$Y, main="TSNE: Exon (c042.col6)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c042.col6.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c042.col6$Y, main="TSNE: Exon (c042.col6)", col=colorbycluster,pch=16)
dev.off()

png("E_limp.corrected_c042.col_compare.tsne.png", width=1200,height=2400)
par(mfrow=c(6,3))
plot(E_limp.tsne.corrected_c042.col1$Y, main="TSNE: LIMP Exon (c042.col1)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c042.col1$Y, main="TSNE: LICL Exon (c042.col1)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c042.col1$Y, main="TSNE: LICR Exon (c042.col1)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c042.col2$Y, main="TSNE: LIMP Exon (c042.col2)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c042.col2$Y, main="TSNE: LICL Exon (c042.col2)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c042.col2$Y, main="TSNE: LICR Exon (c042.col2)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c042.col3$Y, main="TSNE: LIMP Exon (c042.col3)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c042.col3$Y, main="TSNE: LICL Exon (c042.col3)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c042.col3$Y, main="TSNE: LICR Exon (c042.col3)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c042.col4$Y, main="TSNE: LIMP Exon (c042.col4)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c042.col4$Y, main="TSNE: LICL Exon (c042.col4)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c042.col4$Y, main="TSNE: LICR Exon (c042.col4)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c042.col5$Y, main="TSNE: LIMP Exon (c042.col5)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c042.col5$Y, main="TSNE: LICL Exon (c042.col5)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c042.col5$Y, main="TSNE: LICR Exon (c042.col5)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c042.col6$Y, main="TSNE: LIMP Exon (c042.col6)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c042.col6$Y, main="TSNE: LICL Exon (c042.col6)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c042.col6$Y, main="TSNE: LICR Exon (c042.col6)", col=colorbycluster,pch=16)
dev.off()
par(mfrow=c(1,1))

png("E_limp.corrected_c127.col1.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c127.col1$Y, main="TSNE: Exon (c127.col1)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c127.col1.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c127.col1$Y, main="TSNE: Exon (c127.col1)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c127.col1.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c127.col1$Y, main="TSNE: Exon (c127.col1)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c127.col2.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c127.col2$Y, main="TSNE: Exon (c127.col2)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c127.col2.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c127.col2$Y, main="TSNE: Exon (c127.col2)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c127.col2.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c127.col2$Y, main="TSNE: Exon (c127.col2)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c127.col3.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c127.col3$Y, main="TSNE: Exon (c127.col3)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c127.col3.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c127.col3$Y, main="TSNE: Exon (c127.col3)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c127.col3.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c127.col3$Y, main="TSNE: Exon (c127.col3)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c127.col4.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c127.col4$Y, main="TSNE: Exon (c127.col4)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c127.col4.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c127.col4$Y, main="TSNE: Exon (c127.col4)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c127.col4.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c127.col4$Y, main="TSNE: Exon (c127.col4)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c127.col5.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c127.col5$Y, main="TSNE: Exon (c127.col5)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c127.col5.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c127.col5$Y, main="TSNE: Exon (c127.col5)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c127.col5.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c127.col5$Y, main="TSNE: Exon (c127.col5)", col=colorbycluster,pch=16)
dev.off()
png("E_limp.corrected_c127.col6.tsne.png", width=800,height=800)
plot(E_limp.tsne.corrected_c127.col6$Y, main="TSNE: Exon (c127.col6)", col=colorbycluster,pch=16)
dev.off()
png("E_licl.corrected_c127.col6.tsne.png", width=800,height=800)
plot(E_licl.tsne.corrected_c127.col6$Y, main="TSNE: Exon (c127.col6)", col=colorbycluster,pch=16)
dev.off()
png("E_licr.corrected_c127.col6.tsne.png", width=800,height=800)
plot(E_licr.tsne.corrected_c127.col6$Y, main="TSNE: Exon (c127.col6)", col=colorbycluster,pch=16)
dev.off()

png("E_limp.corrected_c127.col_compare.tsne.png", width=1200,height=2400)
par(mfrow=c(6,3))
plot(E_limp.tsne.corrected_c127.col1$Y, main="TSNE: LIMP Exon (c127.col1)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c127.col1$Y, main="TSNE: LICL Exon (c127.col1)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c127.col1$Y, main="TSNE: LICR Exon (c127.col1)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c127.col2$Y, main="TSNE: LIMP Exon (c127.col2)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c127.col2$Y, main="TSNE: LICL Exon (c127.col2)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c127.col2$Y, main="TSNE: LICR Exon (c127.col2)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c127.col3$Y, main="TSNE: LIMP Exon (c127.col3)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c127.col3$Y, main="TSNE: LICL Exon (c127.col3)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c127.col3$Y, main="TSNE: LICR Exon (c127.col3)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c127.col4$Y, main="TSNE: LIMP Exon (c127.col4)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c127.col4$Y, main="TSNE: LICL Exon (c127.col4)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c127.col4$Y, main="TSNE: LICR Exon (c127.col4)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c127.col5$Y, main="TSNE: LIMP Exon (c127.col5)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c127.col5$Y, main="TSNE: LICL Exon (c127.col5)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c127.col5$Y, main="TSNE: LICR Exon (c127.col5)", col=colorbycluster,pch=16)
plot(E_limp.tsne.corrected_c127.col6$Y, main="TSNE: LIMP Exon (c127.col6)", col=colorbycluster,pch=16)
plot(E_licl.tsne.corrected_c127.col6$Y, main="TSNE: LICL Exon (c127.col6)", col=colorbycluster,pch=16)
plot(E_licr.tsne.corrected_c127.col6$Y, main="TSNE: LICR Exon (c127.col6)", col=colorbycluster,pch=16)
dev.off()
par(mfrow=c(1,1))

save.image("20191109.afterExonCorrection.Rdata")

###########################################################################################################################################################################################

# Anyway, the linear model linking the covariates C_rank_ to the data I_limp_ is given by:
# I_lm_limp_ = [ (ones(n_u,1)*zeta_limp_I_un_(1,:) + C_rank_*zeta_limp_I_un_(2:end,:))*transpose(zeta_limp_I_vn_) ] ;

### 12 covariates
zz <- as.data.frame(matrix(ncol=1,nrow=dim(C_rank_c012)[1]))
zz[1:dim(C_rank_c012)[1],1] <- 1
# 1 column
I_lm_limp_c012.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c012[1,1]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_I_un_c012[2:dim(zeta_limp_I_un_c012)[1],1] ) )) %*% t(as.matrix(zeta_limp_I_vn_c012[,1]))
I_lm_licl_c012.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c012[1,1]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licl_I_un_c012)[1],1] ) )) %*% t(as.matrix(zeta_licl_I_vn_c012[,1]))
I_lm_licr_c012.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c012[1,1]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licr_I_un_c012)[1],1] ) )) %*% t(as.matrix(zeta_licr_I_vn_c012[,1]))
# 2 column
I_lm_limp_c012.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c012[1,1:2]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_I_un_c012[2:dim(zeta_limp_I_un_c012)[1],1:2] ) )) %*% t(as.matrix(zeta_limp_I_vn_c012[,1:2]))
I_lm_licl_c012.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c012[1,1:2]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licl_I_un_c012)[1],1:2] ) )) %*% t(as.matrix(zeta_licl_I_vn_c012[,1:2]))
I_lm_licr_c012.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c012[1,1:2]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licr_I_un_c012)[1],1:2] ) )) %*% t(as.matrix(zeta_licr_I_vn_c012[,1:2]))
# 3 column
I_lm_limp_c012.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c012[1,1:3]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_I_un_c012[2:dim(zeta_limp_I_un_c012)[1],1:3] ) )) %*% t(as.matrix(zeta_limp_I_vn_c012[,1:3]))
I_lm_licl_c012.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c012[1,1:3]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licl_I_un_c012)[1],1:3] ) )) %*% t(as.matrix(zeta_licl_I_vn_c012[,1:3]))
I_lm_licr_c012.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c012[1,1:3]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licr_I_un_c012)[1],1:3] ) )) %*% t(as.matrix(zeta_licr_I_vn_c012[,1:3]))
# 4 column
I_lm_limp_c012.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c012[1,1:4]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_I_un_c012[2:dim(zeta_limp_I_un_c012)[1],1:4] ) )) %*% t(as.matrix(zeta_limp_I_vn_c012[,1:4]))
I_lm_licl_c012.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c012[1,1:4]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licl_I_un_c012)[1],1:4] ) )) %*% t(as.matrix(zeta_licl_I_vn_c012[,1:4]))
I_lm_licr_c012.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c012[1,1:4]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licr_I_un_c012)[1],1:4] ) )) %*% t(as.matrix(zeta_licr_I_vn_c012[,1:4]))
# 5 column
I_lm_limp_c012.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c012[1,1:5]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_I_un_c012[2:dim(zeta_limp_I_un_c012)[1],1:5] ) )) %*% t(as.matrix(zeta_limp_I_vn_c012[,1:5]))
I_lm_licl_c012.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c012[1,1:5]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licl_I_un_c012)[1],1:5] ) )) %*% t(as.matrix(zeta_licl_I_vn_c012[,1:5]))
I_lm_licr_c012.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c012[1,1:5]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licr_I_un_c012)[1],1:5] ) )) %*% t(as.matrix(zeta_licr_I_vn_c012[,1:5]))
# 6 column
I_lm_limp_c012.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c012[1,1:6]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_limp_I_un_c012[2:dim(zeta_limp_I_un_c012)[1],1:6] ) )) %*% t(as.matrix(zeta_limp_I_vn_c012[,1:6]))
I_lm_licl_c012.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c012[1,1:6]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licl_I_un_c012)[1],1:6] ) )) %*% t(as.matrix(zeta_licl_I_vn_c012[,1:6]))
I_lm_licr_c012.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c012[1,1:6]))) )) + as.matrix( C_rank_c012[,1:12] ) %*% as.matrix( zeta_licr_I_un_c012[2:dim(zeta_licr_I_un_c012)[1],1:6] ) )) %*% t(as.matrix(zeta_licr_I_vn_c012[,1:6]))


zz <- as.data.frame(matrix(ncol=1,nrow=dim(C_rank_c042)[1]))
zz[1:dim(C_rank_c042)[1],1] <- 1
# 1 column
I_lm_limp_c042.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c042[1,1]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_I_un_c042[2:dim(zeta_limp_I_un_c042)[1],1] ) )) %*% t(as.matrix(zeta_limp_I_vn_c042[,1]))
I_lm_licl_c042.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c042[1,1]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licl_I_un_c042)[1],1] ) )) %*% t(as.matrix(zeta_licl_I_vn_c042[,1]))
I_lm_licr_c042.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c042[1,1]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licr_I_un_c042)[1],1] ) )) %*% t(as.matrix(zeta_licr_I_vn_c042[,1]))
# 2 column
I_lm_limp_c042.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c042[1,1:2]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_I_un_c042[2:dim(zeta_limp_I_un_c042)[1],1:2] ) )) %*% t(as.matrix(zeta_limp_I_vn_c042[,1:2]))
I_lm_licl_c042.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c042[1,1:2]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licl_I_un_c042)[1],1:2] ) )) %*% t(as.matrix(zeta_licl_I_vn_c042[,1:2]))
I_lm_licr_c042.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c042[1,1:2]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licr_I_un_c042)[1],1:2] ) )) %*% t(as.matrix(zeta_licr_I_vn_c042[,1:2]))
# 3 column
I_lm_limp_c042.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c042[1,1:3]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_I_un_c042[2:dim(zeta_limp_I_un_c042)[1],1:3] ) )) %*% t(as.matrix(zeta_limp_I_vn_c042[,1:3]))
I_lm_licl_c042.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c042[1,1:3]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licl_I_un_c042)[1],1:3] ) )) %*% t(as.matrix(zeta_licl_I_vn_c042[,1:3]))
I_lm_licr_c042.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c042[1,1:3]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licr_I_un_c042)[1],1:3] ) )) %*% t(as.matrix(zeta_licr_I_vn_c042[,1:3]))
# 4 column
I_lm_limp_c042.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c042[1,1:4]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_I_un_c042[2:dim(zeta_limp_I_un_c042)[1],1:4] ) )) %*% t(as.matrix(zeta_limp_I_vn_c042[,1:4]))
I_lm_licl_c042.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c042[1,1:4]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licl_I_un_c042)[1],1:4] ) )) %*% t(as.matrix(zeta_licl_I_vn_c042[,1:4]))
I_lm_licr_c042.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c042[1,1:4]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licr_I_un_c042)[1],1:4] ) )) %*% t(as.matrix(zeta_licr_I_vn_c042[,1:4]))
# 5 column
I_lm_limp_c042.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c042[1,1:5]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_I_un_c042[2:dim(zeta_limp_I_un_c042)[1],1:5] ) )) %*% t(as.matrix(zeta_limp_I_vn_c042[,1:5]))
I_lm_licl_c042.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c042[1,1:5]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licl_I_un_c042)[1],1:5] ) )) %*% t(as.matrix(zeta_licl_I_vn_c042[,1:5]))
I_lm_licr_c042.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c042[1,1:5]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licr_I_un_c042)[1],1:5] ) )) %*% t(as.matrix(zeta_licr_I_vn_c042[,1:5]))
# 6 column
I_lm_limp_c042.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c042[1,1:6]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_limp_I_un_c042[2:dim(zeta_limp_I_un_c042)[1],1:6] ) )) %*% t(as.matrix(zeta_limp_I_vn_c042[,1:6]))
I_lm_licl_c042.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c042[1,1:6]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licl_I_un_c042)[1],1:6] ) )) %*% t(as.matrix(zeta_licl_I_vn_c042[,1:6]))
I_lm_licr_c042.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c042[1,1:6]))) )) + as.matrix( C_rank_c042[,1:42] ) %*% as.matrix( zeta_licr_I_un_c042[2:dim(zeta_licr_I_un_c042)[1],1:6] ) )) %*% t(as.matrix(zeta_licr_I_vn_c042[,1:6]))


zz <- as.data.frame(matrix(ncol=1,nrow=dim(C_rank_c127)[1]))
zz[1:dim(C_rank_c127)[1],1] <- 1
# 1 column
I_lm_limp_c127.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c127[1,1]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_I_un_c127[2:dim(zeta_limp_I_un_c127)[1],1] ) )) %*% t(as.matrix(zeta_limp_I_vn_c127[,1]))
I_lm_licl_c127.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c127[1,1]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licl_I_un_c127)[1],1] ) )) %*% t(as.matrix(zeta_licl_I_vn_c127[,1]))
I_lm_licr_c127.col1 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c127[1,1]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licr_I_un_c127)[1],1] ) )) %*% t(as.matrix(zeta_licr_I_vn_c127[,1]))
# 2 column
I_lm_limp_c127.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c127[1,1:2]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_I_un_c127[2:dim(zeta_limp_I_un_c127)[1],1:2] ) )) %*% t(as.matrix(zeta_limp_I_vn_c127[,1:2]))
I_lm_licl_c127.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c127[1,1:2]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licl_I_un_c127)[1],1:2] ) )) %*% t(as.matrix(zeta_licl_I_vn_c127[,1:2]))
I_lm_licr_c127.col2 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c127[1,1:2]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licr_I_un_c127)[1],1:2] ) )) %*% t(as.matrix(zeta_licr_I_vn_c127[,1:2]))
# 3 column
I_lm_limp_c127.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c127[1,1:3]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_I_un_c127[2:dim(zeta_limp_I_un_c127)[1],1:3] ) )) %*% t(as.matrix(zeta_limp_I_vn_c127[,1:3]))
I_lm_licl_c127.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c127[1,1:3]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licl_I_un_c127)[1],1:3] ) )) %*% t(as.matrix(zeta_licl_I_vn_c127[,1:3]))
I_lm_licr_c127.col3 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c127[1,1:3]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licr_I_un_c127)[1],1:3] ) )) %*% t(as.matrix(zeta_licr_I_vn_c127[,1:3]))
# 4 column
I_lm_limp_c127.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c127[1,1:4]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_I_un_c127[2:dim(zeta_limp_I_un_c127)[1],1:4] ) )) %*% t(as.matrix(zeta_limp_I_vn_c127[,1:4]))
I_lm_licl_c127.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c127[1,1:4]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licl_I_un_c127)[1],1:4] ) )) %*% t(as.matrix(zeta_licl_I_vn_c127[,1:4]))
I_lm_licr_c127.col4 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c127[1,1:4]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licr_I_un_c127)[1],1:4] ) )) %*% t(as.matrix(zeta_licr_I_vn_c127[,1:4]))
# 5 column
I_lm_limp_c127.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c127[1,1:5]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_I_un_c127[2:dim(zeta_limp_I_un_c127)[1],1:5] ) )) %*% t(as.matrix(zeta_limp_I_vn_c127[,1:5]))
I_lm_licl_c127.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c127[1,1:5]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licl_I_un_c127)[1],1:5] ) )) %*% t(as.matrix(zeta_licl_I_vn_c127[,1:5]))
I_lm_licr_c127.col5 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c127[1,1:5]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licr_I_un_c127)[1],1:5] ) )) %*% t(as.matrix(zeta_licr_I_vn_c127[,1:5]))
# 6 column
I_lm_limp_c127.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_limp_I_un_c127[1,1:6]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_limp_I_un_c127[2:dim(zeta_limp_I_un_c127)[1],1:6] ) )) %*% t(as.matrix(zeta_limp_I_vn_c127[,1:6]))
I_lm_licl_c127.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licl_I_un_c127[1,1:6]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licl_I_un_c127)[1],1:6] ) )) %*% t(as.matrix(zeta_licl_I_vn_c127[,1:6]))
I_lm_licr_c127.col6 <- ((  as.matrix( zz ) %*% t(as.matrix( unlist(unname(as.list(zeta_licr_I_un_c127[1,1:6]))) )) + as.matrix( C_rank_c127[,1:127] ) %*% as.matrix( zeta_licr_I_un_c127[2:dim(zeta_licr_I_un_c127)[1],1:6] ) )) %*% t(as.matrix(zeta_licr_I_vn_c127[,1:6]))

# and the difference:
#   I_corrected_ = I_limp_ - I_lm_limp_
# would represent 'covariate-corrected' data, which could then be put into t-sne.
# 1 column
I_limp.corrected_c012.col1 <- I_limp - I_lm_limp_c012.col1
I_licl.corrected_c012.col1 <- I_licl - I_lm_licl_c012.col1
I_licr.corrected_c012.col1 <- I_licr - I_lm_licr_c012.col1
# 2 column
I_limp.corrected_c012.col2 <- I_limp - I_lm_limp_c012.col2
I_licl.corrected_c012.col2 <- I_licl - I_lm_licl_c012.col2
I_licr.corrected_c012.col2 <- I_licr - I_lm_licr_c012.col2
# 3 column
I_limp.corrected_c012.col3 <- I_limp - I_lm_limp_c012.col3
I_licl.corrected_c012.col3 <- I_licl - I_lm_licl_c012.col3
I_licr.corrected_c012.col3 <- I_licr - I_lm_licr_c012.col3
# 3 column
I_limp.corrected_c012.col4 <- I_limp - I_lm_limp_c012.col4
I_licl.corrected_c012.col4 <- I_licl - I_lm_licl_c012.col4
I_licr.corrected_c012.col4 <- I_licr - I_lm_licr_c012.col4
# 3 column
I_limp.corrected_c012.col5 <- I_limp - I_lm_limp_c012.col5
I_licl.corrected_c012.col5 <- I_licl - I_lm_licl_c012.col5
I_licr.corrected_c012.col5 <- I_licr - I_lm_licr_c012.col5
# 3 column
I_limp.corrected_c012.col6 <- I_limp - I_lm_limp_c012.col6
I_licl.corrected_c012.col6 <- I_licl - I_lm_licl_c012.col6
I_licr.corrected_c012.col6 <- I_licr - I_lm_licr_c012.col6

# 1 column
I_limp.corrected_c042.col1 <- I_limp - I_lm_limp_c042.col1
I_licl.corrected_c042.col1 <- I_licl - I_lm_licl_c042.col1
I_licr.corrected_c042.col1 <- I_licr - I_lm_licr_c042.col1
# 2 column
I_limp.corrected_c042.col2 <- I_limp - I_lm_limp_c042.col2
I_licl.corrected_c042.col2 <- I_licl - I_lm_licl_c042.col2
I_licr.corrected_c042.col2 <- I_licr - I_lm_licr_c042.col2
# 3 column
I_limp.corrected_c042.col3 <- I_limp - I_lm_limp_c042.col3
I_licl.corrected_c042.col3 <- I_licl - I_lm_licl_c042.col3
I_licr.corrected_c042.col3 <- I_licr - I_lm_licr_c042.col3
# 3 column
I_limp.corrected_c042.col4 <- I_limp - I_lm_limp_c042.col4
I_licl.corrected_c042.col4 <- I_licl - I_lm_licl_c042.col4
I_licr.corrected_c042.col4 <- I_licr - I_lm_licr_c042.col4
# 3 column
I_limp.corrected_c042.col5 <- I_limp - I_lm_limp_c042.col5
I_licl.corrected_c042.col5 <- I_licl - I_lm_licl_c042.col5
I_licr.corrected_c042.col5 <- I_licr - I_lm_licr_c042.col5
# 3 column
I_limp.corrected_c042.col6 <- I_limp - I_lm_limp_c042.col6
I_licl.corrected_c042.col6 <- I_licl - I_lm_licl_c042.col6
I_licr.corrected_c042.col6 <- I_licr - I_lm_licr_c042.col6

# 1 column
I_limp.corrected_c042.col1 <- I_limp - I_lm_limp_c042.col1
I_licl.corrected_c042.col1 <- I_licl - I_lm_licl_c042.col1
I_licr.corrected_c042.col1 <- I_licr - I_lm_licr_c042.col1
# 2 column
I_limp.corrected_c042.col2 <- I_limp - I_lm_limp_c042.col2
I_licl.corrected_c042.col2 <- I_licl - I_lm_licl_c042.col2
I_licr.corrected_c042.col2 <- I_licr - I_lm_licr_c042.col2
# 3 column
I_limp.corrected_c042.col3 <- I_limp - I_lm_limp_c042.col3
I_licl.corrected_c042.col3 <- I_licl - I_lm_licl_c042.col3
I_licr.corrected_c042.col3 <- I_licr - I_lm_licr_c042.col3
# 3 column
I_limp.corrected_c042.col4 <- I_limp - I_lm_limp_c042.col4
I_licl.corrected_c042.col4 <- I_licl - I_lm_licl_c042.col4
I_licr.corrected_c042.col4 <- I_licr - I_lm_licr_c042.col4
# 3 column
I_limp.corrected_c042.col5 <- I_limp - I_lm_limp_c042.col5
I_licl.corrected_c042.col5 <- I_licl - I_lm_licl_c042.col5
I_licr.corrected_c042.col5 <- I_licr - I_lm_licr_c042.col5
# 3 column
I_limp.corrected_c042.col6 <- I_limp - I_lm_limp_c042.col6
I_licl.corrected_c042.col6 <- I_licl - I_lm_licl_c042.col6
I_licr.corrected_c042.col6 <- I_licr - I_lm_licr_c042.col6

# 1 column
I_limp.corrected_c127.col1 <- I_limp - I_lm_limp_c127.col1
I_licl.corrected_c127.col1 <- I_licl - I_lm_licl_c127.col1
I_licr.corrected_c127.col1 <- I_licr - I_lm_licr_c127.col1
# 2 column
I_limp.corrected_c127.col2 <- I_limp - I_lm_limp_c127.col2
I_licl.corrected_c127.col2 <- I_licl - I_lm_licl_c127.col2
I_licr.corrected_c127.col2 <- I_licr - I_lm_licr_c127.col2
# 3 column
I_limp.corrected_c127.col3 <- I_limp - I_lm_limp_c127.col3
I_licl.corrected_c127.col3 <- I_licl - I_lm_licl_c127.col3
I_licr.corrected_c127.col3 <- I_licr - I_lm_licr_c127.col3
# 3 column
I_limp.corrected_c127.col4 <- I_limp - I_lm_limp_c127.col4
I_licl.corrected_c127.col4 <- I_licl - I_lm_licl_c127.col4
I_licr.corrected_c127.col4 <- I_licr - I_lm_licr_c127.col4
# 3 column
I_limp.corrected_c127.col5 <- I_limp - I_lm_limp_c127.col5
I_licl.corrected_c127.col5 <- I_licl - I_lm_licl_c127.col5
I_licr.corrected_c127.col5 <- I_licr - I_lm_licr_c127.col5
# 3 column
I_limp.corrected_c127.col6 <- I_limp - I_lm_limp_c127.col6
I_licl.corrected_c127.col6 <- I_licl - I_lm_licl_c127.col6
I_licr.corrected_c127.col6 <- I_licr - I_lm_licr_c127.col6

I_limp.tsne.corrected_c012.col1 <- Rtsne(I_limp.corrected_c012.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c012.col1 <- Rtsne(I_licl.corrected_c012.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c012.col1 <- Rtsne(I_licr.corrected_c012.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c012.col2 <- Rtsne(I_limp.corrected_c012.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c012.col2 <- Rtsne(I_licl.corrected_c012.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c012.col2 <- Rtsne(I_licr.corrected_c012.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c012.col3 <- Rtsne(I_limp.corrected_c012.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c012.col3 <- Rtsne(I_licl.corrected_c012.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c012.col3 <- Rtsne(I_licr.corrected_c012.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c012.col4 <- Rtsne(I_limp.corrected_c012.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c012.col4 <- Rtsne(I_licl.corrected_c012.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c012.col4 <- Rtsne(I_licr.corrected_c012.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c012.col5 <- Rtsne(I_limp.corrected_c012.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c012.col5 <- Rtsne(I_licl.corrected_c012.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c012.col5 <- Rtsne(I_licr.corrected_c012.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c012.col6 <- Rtsne(I_limp.corrected_c012.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c012.col6 <- Rtsne(I_licl.corrected_c012.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c012.col6 <- Rtsne(I_licr.corrected_c012.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c042.col1 <- Rtsne(I_limp.corrected_c042.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c042.col1 <- Rtsne(I_licl.corrected_c042.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c042.col1 <- Rtsne(I_licr.corrected_c042.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c042.col2 <- Rtsne(I_limp.corrected_c042.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c042.col2 <- Rtsne(I_licl.corrected_c042.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c042.col2 <- Rtsne(I_licr.corrected_c042.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c042.col3 <- Rtsne(I_limp.corrected_c042.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c042.col3 <- Rtsne(I_licl.corrected_c042.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c042.col3 <- Rtsne(I_licr.corrected_c042.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c042.col4 <- Rtsne(I_limp.corrected_c042.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c042.col4 <- Rtsne(I_licl.corrected_c042.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c042.col4 <- Rtsne(I_licr.corrected_c042.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c042.col5 <- Rtsne(I_limp.corrected_c042.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c042.col5 <- Rtsne(I_licl.corrected_c042.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c042.col5 <- Rtsne(I_licr.corrected_c042.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c042.col6 <- Rtsne(I_limp.corrected_c042.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c042.col6 <- Rtsne(I_licl.corrected_c042.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c042.col6 <- Rtsne(I_licr.corrected_c042.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c127.col1 <- Rtsne(I_limp.corrected_c127.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c127.col1 <- Rtsne(I_licl.corrected_c127.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c127.col1 <- Rtsne(I_licr.corrected_c127.col1, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c127.col2 <- Rtsne(I_limp.corrected_c127.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c127.col2 <- Rtsne(I_licl.corrected_c127.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c127.col2 <- Rtsne(I_licr.corrected_c127.col2, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c127.col3 <- Rtsne(I_limp.corrected_c127.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c127.col3 <- Rtsne(I_licl.corrected_c127.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c127.col3 <- Rtsne(I_licr.corrected_c127.col3, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c127.col4 <- Rtsne(I_limp.corrected_c127.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c127.col4 <- Rtsne(I_licl.corrected_c127.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c127.col4 <- Rtsne(I_licr.corrected_c127.col4, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c127.col5 <- Rtsne(I_limp.corrected_c127.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c127.col5 <- Rtsne(I_licl.corrected_c127.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c127.col5 <- Rtsne(I_licr.corrected_c127.col5, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)

I_limp.tsne.corrected_c127.col6 <- Rtsne(I_limp.corrected_c127.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licl.tsne.corrected_c127.col6 <- Rtsne(I_licl.corrected_c127.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
I_licr.tsne.corrected_c127.col6 <- Rtsne(I_licr.corrected_c127.col6, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)


png("I_limp.corrected_c012.col1.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c012.col1$Y, main="TSNE: Intron (c012.col1)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c012.col1.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c012.col1$Y, main="TSNE: Intron (c012.col1)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c012.col1.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c012.col1$Y, main="TSNE: Intron (c012.col1)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c012.col2.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c012.col2$Y, main="TSNE: Intron (c012.col2)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c012.col2.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c012.col2$Y, main="TSNE: Intron (c012.col2)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c012.col2.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c012.col2$Y, main="TSNE: Intron (c012.col2)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c012.col3.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c012.col3$Y, main="TSNE: Intron (c012.col3)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c012.col3.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c012.col3$Y, main="TSNE: Intron (c012.col3)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c012.col3.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c012.col3$Y, main="TSNE: Intron (c012.col3)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c012.col4.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c012.col4$Y, main="TSNE: Intron (c012.col4)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c012.col4.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c012.col4$Y, main="TSNE: Intron (c012.col4)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c012.col4.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c012.col4$Y, main="TSNE: Intron (c012.col4)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c012.col5.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c012.col5$Y, main="TSNE: Intron (c012.col5)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c012.col5.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c012.col5$Y, main="TSNE: Intron (c012.col5)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c012.col5.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c012.col5$Y, main="TSNE: Intron (c012.col5)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c012.col6.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c012.col6$Y, main="TSNE: Intron (c012.col6)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c012.col6.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c012.col6$Y, main="TSNE: Intron (c012.col6)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c012.col6.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c012.col6$Y, main="TSNE: Intron (c012.col6)", col=colorbycluster,pch=16)
dev.off()

png("I_limp.corrected_c012.col_compare.tsne.png", width=1200,height=2400)
par(mfrow=c(6,3))
plot(I_limp.tsne.corrected_c012.col1$Y, main="TSNE: LIMP Intron (c012.col1)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c012.col1$Y, main="TSNE: LICL Intron (c012.col1)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c012.col1$Y, main="TSNE: LICR Intron (c012.col1)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c012.col2$Y, main="TSNE: LIMP Intron (c012.col2)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c012.col2$Y, main="TSNE: LICL Intron (c012.col2)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c012.col2$Y, main="TSNE: LICR Intron (c012.col2)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c012.col3$Y, main="TSNE: LIMP Intron (c012.col3)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c012.col3$Y, main="TSNE: LICL Intron (c012.col3)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c012.col3$Y, main="TSNE: LICR Intron (c012.col3)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c012.col4$Y, main="TSNE: LIMP Intron (c012.col4)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c012.col4$Y, main="TSNE: LICL Intron (c012.col4)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c012.col4$Y, main="TSNE: LICR Intron (c012.col4)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c012.col5$Y, main="TSNE: LIMP Intron (c012.col5)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c012.col5$Y, main="TSNE: LICL Intron (c012.col5)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c012.col5$Y, main="TSNE: LICR Intron (c012.col5)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c012.col6$Y, main="TSNE: LIMP Intron (c012.col6)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c012.col6$Y, main="TSNE: LICL Intron (c012.col6)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c012.col6$Y, main="TSNE: LICR Intron (c012.col6)", col=colorbycluster,pch=16)
dev.off()
par(mfrow=c(1,1))

png("I_limp.corrected_c042.col1.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c042.col1$Y, main="TSNE: Intron (c042.col1)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c042.col1.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c042.col1$Y, main="TSNE: Intron (c042.col1)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c042.col1.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c042.col1$Y, main="TSNE: Intron (c042.col1)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c042.col2.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c042.col2$Y, main="TSNE: Intron (c042.col2)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c042.col2.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c042.col2$Y, main="TSNE: Intron (c042.col2)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c042.col2.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c042.col2$Y, main="TSNE: Intron (c042.col2)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c042.col3.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c042.col3$Y, main="TSNE: Intron (c042.col3)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c042.col3.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c042.col3$Y, main="TSNE: Intron (c042.col3)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c042.col3.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c042.col3$Y, main="TSNE: Intron (c042.col3)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c042.col4.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c042.col4$Y, main="TSNE: Intron (c042.col4)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c042.col4.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c042.col4$Y, main="TSNE: Intron (c042.col4)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c042.col4.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c042.col4$Y, main="TSNE: Intron (c042.col4)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c042.col5.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c042.col5$Y, main="TSNE: Intron (c042.col5)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c042.col5.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c042.col5$Y, main="TSNE: Intron (c042.col5)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c042.col5.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c042.col5$Y, main="TSNE: Intron (c042.col5)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c042.col6.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c042.col6$Y, main="TSNE: Intron (c042.col6)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c042.col6.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c042.col6$Y, main="TSNE: Intron (c042.col6)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c042.col6.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c042.col6$Y, main="TSNE: Intron (c042.col6)", col=colorbycluster,pch=16)
dev.off()

png("I_limp.corrected_c042.col_compare.tsne.png", width=1200,height=2400)
par(mfrow=c(6,3))
plot(I_limp.tsne.corrected_c042.col1$Y, main="TSNE: LIMP Intron (c042.col1)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c042.col1$Y, main="TSNE: LICL Intron (c042.col1)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c042.col1$Y, main="TSNE: LICR Intron (c042.col1)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c042.col2$Y, main="TSNE: LIMP Intron (c042.col2)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c042.col2$Y, main="TSNE: LICL Intron (c042.col2)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c042.col2$Y, main="TSNE: LICR Intron (c042.col2)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c042.col3$Y, main="TSNE: LIMP Intron (c042.col3)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c042.col3$Y, main="TSNE: LICL Intron (c042.col3)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c042.col3$Y, main="TSNE: LICR Intron (c042.col3)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c042.col4$Y, main="TSNE: LIMP Intron (c042.col4)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c042.col4$Y, main="TSNE: LICL Intron (c042.col4)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c042.col4$Y, main="TSNE: LICR Intron (c042.col4)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c042.col5$Y, main="TSNE: LIMP Intron (c042.col5)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c042.col5$Y, main="TSNE: LICL Intron (c042.col5)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c042.col5$Y, main="TSNE: LICR Intron (c042.col5)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c042.col6$Y, main="TSNE: LIMP Intron (c042.col6)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c042.col6$Y, main="TSNE: LICL Intron (c042.col6)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c042.col6$Y, main="TSNE: LICR Intron (c042.col6)", col=colorbycluster,pch=16)
dev.off()
par(mfrow=c(1,1))

png("I_limp.corrected_c127.col1.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c127.col1$Y, main="TSNE: Intron (c127.col1)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c127.col1.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c127.col1$Y, main="TSNE: Intron (c127.col1)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c127.col1.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c127.col1$Y, main="TSNE: Intron (c127.col1)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c127.col2.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c127.col2$Y, main="TSNE: Intron (c127.col2)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c127.col2.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c127.col2$Y, main="TSNE: Intron (c127.col2)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c127.col2.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c127.col2$Y, main="TSNE: Intron (c127.col2)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c127.col3.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c127.col3$Y, main="TSNE: Intron (c127.col3)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c127.col3.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c127.col3$Y, main="TSNE: Intron (c127.col3)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c127.col3.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c127.col3$Y, main="TSNE: Intron (c127.col3)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c127.col4.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c127.col4$Y, main="TSNE: Intron (c127.col4)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c127.col4.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c127.col4$Y, main="TSNE: Intron (c127.col4)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c127.col4.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c127.col4$Y, main="TSNE: Intron (c127.col4)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c127.col5.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c127.col5$Y, main="TSNE: Intron (c127.col5)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c127.col5.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c127.col5$Y, main="TSNE: Intron (c127.col5)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c127.col5.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c127.col5$Y, main="TSNE: Intron (c127.col5)", col=colorbycluster,pch=16)
dev.off()
png("I_limp.corrected_c127.col6.tsne.png", width=800,height=800)
plot(I_limp.tsne.corrected_c127.col6$Y, main="TSNE: Intron (c127.col6)", col=colorbycluster,pch=16)
dev.off()
png("I_licl.corrected_c127.col6.tsne.png", width=800,height=800)
plot(I_licl.tsne.corrected_c127.col6$Y, main="TSNE: Intron (c127.col6)", col=colorbycluster,pch=16)
dev.off()
png("I_licr.corrected_c127.col6.tsne.png", width=800,height=800)
plot(I_licr.tsne.corrected_c127.col6$Y, main="TSNE: Intron (c127.col6)", col=colorbycluster,pch=16)
dev.off()

png("I_limp.corrected_c127.col_compare.tsne.png", width=1200,height=2400)
par(mfrow=c(6,3))
plot(I_limp.tsne.corrected_c127.col1$Y, main="TSNE: LIMP Intron (c127.col1)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c127.col1$Y, main="TSNE: LICL Intron (c127.col1)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c127.col1$Y, main="TSNE: LICR Intron (c127.col1)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c127.col2$Y, main="TSNE: LIMP Intron (c127.col2)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c127.col2$Y, main="TSNE: LICL Intron (c127.col2)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c127.col2$Y, main="TSNE: LICR Intron (c127.col2)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c127.col3$Y, main="TSNE: LIMP Intron (c127.col3)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c127.col3$Y, main="TSNE: LICL Intron (c127.col3)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c127.col3$Y, main="TSNE: LICR Intron (c127.col3)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c127.col4$Y, main="TSNE: LIMP Intron (c127.col4)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c127.col4$Y, main="TSNE: LICL Intron (c127.col4)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c127.col4$Y, main="TSNE: LICR Intron (c127.col4)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c127.col5$Y, main="TSNE: LIMP Intron (c127.col5)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c127.col5$Y, main="TSNE: LICL Intron (c127.col5)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c127.col5$Y, main="TSNE: LICR Intron (c127.col5)", col=colorbycluster,pch=16)
plot(I_limp.tsne.corrected_c127.col6$Y, main="TSNE: LIMP Intron (c127.col6)", col=colorbycluster,pch=16)
plot(I_licl.tsne.corrected_c127.col6$Y, main="TSNE: LICL Intron (c127.col6)", col=colorbycluster,pch=16)
plot(I_licr.tsne.corrected_c127.col6$Y, main="TSNE: LICR Intron (c127.col6)", col=colorbycluster,pch=16)
dev.off()
par(mfrow=c(1,1))

# save.image("20191109.afterExonAndIntronCorrection.Rdata")


# exe 1-7
# jump to here
load("20191109.afterExonAndIntronCorrection.Rdata")
























































###################################################################################




clust_par_absZ.cov.ex_c127 <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/AB_E_c127_C_absZ_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
clust_par_absZ.cov.ex_c127 <- clust_par_absZ.cov.ex_c127[ , !colSums(is.na(clust_par_absZ.cov.ex_c127)) > 0 ]

clust_par_absZ.cov.in_c127 <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/AB_I_c127_C_absZ_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
clust_par_absZ.cov.in_c127 <- clust_par_absZ.cov.in_c127[ , !colSums(is.na(clust_par_absZ.cov.in_c127)) > 0 ]

clust_par_relZ.cov.ex_c127 <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/AB_E_c127_C_rawZ_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
clust_par_relZ.cov.ex_c127 <- clust_par_relZ.cov.ex_c127[ , !colSums(is.na(clust_par_relZ.cov.ex_c127)) > 0 ]

clust_par_relZ.cov.in_c127 <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/AB_I_c127_C_rawZ_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
clust_par_relZ.cov.in_c127 <- clust_par_relZ.cov.in_c127[ , !colSums(is.na(clust_par_relZ.cov.in_c127)) > 0 ]

I_CLabel_pair_c127 <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/AB_I_c127_CLabel_pair_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
I_CLabel_pair_c127 <- unlist(as.list(I_CLabel_pair_c127[ , !colSums(is.na(I_CLabel_pair_c127)) > 0 ]))

E_CLabel_pair_c127 <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/AB_E_c127_CLabel_pair_.tsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
E_CLabel_pair_c127 <- unlist(as.list(E_CLabel_pair_c127[ , !colSums(is.na(E_CLabel_pair_c127)) > 0 ]))



#C_VariableName_ <- read.table("simple_rank_based_analysis_20191018/C_VariableName_.nsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
C_VariableName_ <- read.table("../../Google\ Drive/NeuralQC/simple_rank_based_analysis/dir_data/dir_mat/C_VariableName_.nsv", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
C_VariableName_ <- unlist(as.list(C_VariableName_[ , !colSums(is.na(C_VariableName_)) > 0 ]))

colnames(clust_par_relZ.cov.ex_c127) <- colnames(clust_par_relZ.cov.in_c127) <- C_VariableName_
rownames(clust_par_relZ.cov.ex_c127) <- E_CLabel_pair_c127
rownames(clust_par_relZ.cov.in_c127) <- I_CLabel_pair_c127
colnames(clust_par_absZ.cov.ex_c127) <- colnames(clust_par_absZ.cov.in_c127) <- C_VariableName_
rownames(clust_par_absZ.cov.ex_c127) <- E_CLabel_pair_c127
rownames(clust_par_absZ.cov.in_c127) <- I_CLabel_pair_c127

# I_CLabel_pair_.nsv
covariate.matrix <- read.table("../old_desktop_20191107/BRAIN_in_ex/20161026_covariate_table.format.txt", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
colnames(covariate.matrix) <- covariate.matrix[1, ]
covariate.matrix           <- covariate.matrix[-1,]
rownames(covariate.matrix) <- covariate.matrix[, 1]
covariate.matrix           <- covariate.matrix[,-1]

library(vegan)
library(gplots)
library(cluster)
library(pvclust)

# E_CLabel_pair_.nsv
# I_CLabel_pair_.nsv

# 1 : covariates (exons)

clust_par_relZ.cov.ex.dit_cov_c127 <- vegdist(t(clust_par_relZ.cov.ex_c127), method = "euclidean")
png("clust_par_relZ.cov.ex.dit_cov.20191111.png",width=1400,height=1400)
h <- heatmap.2( as.matrix(clust_par_relZ.cov.ex.dit_cov_c127), trace = "none", keysize=0.7, margins = c(20,20), col=rev(redgreen(150))  )
dev.off()

clust_par_absZ.cov.ex.dit_cov_c127 <- vegdist(t(clust_par_absZ.cov.ex_c127), method = "euclidean")
png("clust_par_absZ.cov.ex.dit_cov.20191111.png",width=1400,height=1400)
h <- heatmap.2( as.matrix(clust_par_absZ.cov.ex.dit_cov_c127), trace = "none", keysize=0.7, margins = c(20,20), col=rev(redgreen(150))  )
dev.off()

pvc.rel <- pvclust( as.matrix(clust_par_relZ.cov.ex_c127), method.dist = "cor", method.hclust = "average" )
png("clust_par_relZ.cov.ex.pvclust.20191111.png",width=1200,height=800)
plot(pvc.rel, cex=0.5)
pvrect(pvc.rel, alpha=0.97)
dev.off()

pvc.abs <- pvclust( as.matrix(clust_par_absZ.cov.ex_c127), method.dist = "cor", method.hclust = "average" )
png("clust_par_absZ.cov.ex.pvclust.20191111.png",width=1200,height=800)
plot(pvc.abs, cex=0.5)
pvrect(pvc.abs, alpha=0.97)
dev.off()

seplot(pvc.rel)
seplot(pvc.abs)

pvcr.abs <- pvpick(pvc.abs, alpha=0.99)
pvcr.rel <- pvpick(pvc.rel, alpha=0.99)
 

 


#save.image("onlyClusterPair_analysis.20191111.Rdata")
load("onlyClusterPair_analysis.20191111.Rdata")

orderedCovars <- unlist(pvcr.abs$clusters)
orderedCovars <- c(orderedCovars, unlist(pvcr.rel$clusters)[ !unlist(pvcr.rel$clusters) %in%  unlist(pvcr.abs$clusters) ] )
orderedCovars <- c(orderedCovars, unlist(C_VariableName_)[ !unlist(C_VariableName_) %in%  unlist(pvcr.abs$clusters) ] )
clust.pop <- matrix(nrow=2,ncol=length(orderedCovars))
rownames(clust.pop) <- c("abs","raw")
colnames(clust.pop) <- orderedCovars
clust.pop.col <- clust.pop
clust.pop.col[!is.na(clust.pop.col)] <- "white"
pi <- 1
for(p in pvcr.abs$clusters){
  for(q in p){
    clust.pop[1,q] <- pi
    clust.pop.col[1,q] <- primary.colors(27)[pi]
  }
  pi <- pi + 1
}
pi <- 1
for(p in pvcr.rel$clusters){
  for(q in p){
    clust.pop[2,q] <- pi
    clust.pop.col[2,q] <- primary.colors(27)[pi]
  }
  pi <- pi + 1
}
clust.pop[is.na(clust.pop)] <- 27
library(gplots)
png("clustered_contents.20191115.0.99.png", width=1400,height=600)
heatmap.2(clust.pop, col=unlist(c(primary.colors(26),"white")), 
          trace='none', Rowv = FALSE, Colv = TRUE, 
          margins = c(20,10), keysize=0.7)
dev.off()
write.table(clust.pop,"clustered_contents.20191115.0.99.tsv",sep="\t")
    
# Take the 21st cluster for example. AU p-value is estimated as 0.795, with SE (standard error) 0.733. By an analogue of standard normal 
# theory, the "true" AU p-value is roughly estimated to exist in between (AU - 2 * SE) and (AU + 2 * SE). From this rough inference, the 
# true AU p-value seems to be in between -0.671 and 2.261. However AU is defined to be between 0 and 1, so this inference is of no meaning 
# (Note: As shown in this example, the normal approximation fails for such large SE values. There is a better approximation for constructing 
# confidence intervals based on z-values, though not yet implemented in our package. This feature will be included in a future release). 

print(pvc.abs, which=pvcr.abs$edges, digits=10)
print(pvc.rel, which=pvcr.rel$edges, digits=10)

# save.image("onlyClustPairEnv.Rdata")

clust_of_interest <- 20
pvcr.abs$clusters[[20]]

clust_of_interest.covars <- pvcr.abs$clusters[[clust_of_interest]]
clust_par_relZ.cov.ex.dit_cov_c127 <- vegdist(clust_par_relZ.cov.ex_c127[,clust_of_interest.covars], method = "euclidean")

outliers <- sort(as.numeric(sort(unique(covariate.matrix[ covariate.matrix$IsOutlier_20161007 == "yes", "Cluster_ID_20160909" ]))))
exc <- sort(as.numeric(sort(unique(covariate.matrix[ covariate.matrix$Cluster_Grp_20161007 == "exc", "Cluster_ID_20160909" ]))))
inh <- sort(as.numeric(sort(unique(covariate.matrix[ covariate.matrix$Cluster_Grp_20161007 == "inh", "Cluster_ID_20160909" ]))))
glia <- sort(as.numeric(sort(unique(covariate.matrix[ covariate.matrix$Cluster_Grp_20161007 == "glia", "Cluster_ID_20160909" ]))))

cold <- list()
cildx <- list()
cildy <- list()
cd <- 1
for(cq in strsplit(rownames(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)),",")){
  if(as.numeric(cq[1]) %in% outliers || as.numeric(cq[2]) %in% outliers){
    cold[cd] <- "white"
  }else{
    cold[cd] <- "purple"
  }
  
  if(as.numeric(cq[1]) %in% exc){
    cildx[cd] <- "blue"
  }else if(as.numeric(cq[1]) %in% inh){
    cildx[cd] <- "orange"
  }else if(as.numeric(cq[1]) %in% glia){
    cildx[cd] <- "grey"
  }else{
    cildx[cd] <- "white"
  }
  if(as.numeric(cq[2]) %in% exc){
    cildy[cd] <- "blue"
  }else if(as.numeric(cq[2]) %in% inh){
    cildy[cd] <- "orange"
  }else if(as.numeric(cq[2]) %in% glia){
    cildy[cd] <- "grey"
  }else{
    cildy[cd] <- "white"
  }
  cd <- cd + 1
}

png("20191110.clust23.outlier.png", width=800,height=800)
heatmap.2( log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1), trace = "none", keysize=0.7, margins = c(10,10), col=rev(redgreen(150)),
           RowSideColors=unlist(cold)
           )
dev.off()

png("20191110.clust23.celltype.png", width=800,height=800)
heatmap.2( log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1), trace = "none", keysize=0.7, margins = c(10,10), col=rev(redgreen(150)),
           RowSideColors=unlist(cildx),
           ColSideColors=unlist(cildy)
)
dev.off()

#hh <- hclust(log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1))
pvc.clust_par_relZ.50      <- pvclust( log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1), method.dist = "cor", method.hclust = "average", nboot=50 )

#pvc.clust_par_relZ.250   <- pvclust( log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1), method.dist = "cor", method.hclust = "average", nboot=250 )
 
pvc.clust_par_relZ.50.pick <- pvpick(pvc.clust_par_relZ.50, alpha=0.99)
length(pvc.clust_par_relZ.50.pick$clusters)
#plot(sort(unlist(lapply(as.list(pvc.clust_par_relZ.50.pick$clusters), FUN=length))))
cmap <- matrix(nrow=46,ncol=46)
cmap[is.na(cmap)] <- 0
for(r in c(1:151)){
  for(s in strsplit( pvc.clust_par_relZ.50.pick$clusters[[r]] , "," )){
    cmap[as.numeric(s[1]),as.numeric(s[2])] <- cmap[as.numeric(s[1]),as.numeric(s[2])] + 1
    cmap[as.numeric(s[2]),as.numeric(s[1])] <- cmap[as.numeric(s[2]),as.numeric(s[1])] + 1
  }
  heatmap.2(cmap, col=redgreen(50), trace='none',
            ColSideColors = redgreen(max(colSums(cmap))+ 1)[ colSums(cmap)+1 ])
} 

#pvc.clust_par_relZ.pick <- pvpick(pvc.clust_par_relZ, alpha=0.9999)
#length(pvc.clust_par_relZ.pick$clusters)
#pvc.clust_par_relZ.pick <- pvpick(pvc.clust_par_relZ, alpha=0.9999999)
#length(pvc.clust_par_relZ.pick$clusters)
#pvc.clust_par_relZ.pick <- pvpick(pvc.clust_par_relZ, alpha=0.9999999999)
#length(pvc.clust_par_relZ.pick$clusters)

plot(pvc.clust_par_relZ.50, cex=3, )
pvrect(pvc.clust_par_relZ.50, alpha=0.70)
seplot(pvc.rel)
 
#pvc.clust_par_relZ.nL <- pvclust( as.matrix(clust_par_relZ.cov.ex.dit_cov_c127), method.dist = "cor", method.hclust = "average" )


#p <- prcomp(as.matrix(clust_par_relZ.cov.ex.dit_sam), scale=TRUE)


#plot(p$x)
 
clust_of_interest <- 24
clust_of_interest.covars <- names(p.pam$clustering[p.pam$clustering == clust_of_interest])
clust_par_relZ.cov.ex.dit_cov_c127 <- vegdist(clust_par_relZ.cov.ex_c127[,clust_of_interest.covars], method = "euclidean")
hy <- heatmap.2( log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1), trace = "none", keysize=0.7, margins = c(10,10), col=rev(redgreen(150))  )
  
par(mfrow=c(1,1))
ww_e_c127 <- hclust(as.dist(clust_par_relZ.cov.ex.dit_cov_c127), method="complete")

cutHere <- 10
table(unname(cutree(ww_e_c127, h=cutHere)))
hclust.clustcols <- primary.colors(max(cutree(ww_e_c127, h=cutHere)))[unname(cutree(ww_e_c127, h=cutHere))]
png("clust_par_relZ.cov.ex.dit_cov_c127.hierarchical.clust24.png", width=800,height=800)
heatmap.2( log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1 ), trace = "none", keysize=0.7, margins = c(5,5), col=rev(redgreen(150)),
           RowSideColors = hclust.clustcols,
           ColSideColors = hclust.clustcols
) 
dev.off()


cluster_assoc_density_matrix <- matrix(ncol=max(cutree(ww_e_c127, h=cutHere)), 
                                       nrow=46)
colnames(cluster_assoc_density_matrix) <- c(1:max(cutree(ww_e_c127, h=cutHere)))
rownames(cluster_assoc_density_matrix) <- c(1:46)
cluster_assoc_density_matrix[is.na(cluster_assoc_density_matrix)] <- 0
for( run_this_clust in 1:max(cutree(ww_e_c127, h=cutHere)) ){
  for( q in names(cutree(ww_e_c127, h=cutHere)[cutree(ww_e_c127, h=cutHere)==run_this_clust]) ){
    if(strsplit(q, ',')[[1]][1] != "NA" && strsplit(q, ',')[[1]][2] != "NA"){
      cluster_assoc_density_matrix[ strsplit(q, ',')[[1]][1] , run_this_clust ] <- cluster_assoc_density_matrix[ strsplit(q, ',')[[1]][1] , run_this_clust ] + 1
      cluster_assoc_density_matrix[ strsplit(q, ',')[[1]][2] , run_this_clust ] <- cluster_assoc_density_matrix[ strsplit(q, ',')[[1]][2] , run_this_clust ] + 1
    }
  }
}

#save.image("20191104devel.Rdata")
 
load("../old_desktop_20191107/20191104devel.Rdata")

sample_assoc_density_matrix <- matrix(ncol=max(cutree(ww_e_c127, h=cutHere)), 
                                       nrow=dim(E_licl)[1])

#u_CLabel_sub_.q <- unlist(u_CLabel_sub_)
#names(u_CLabel_sub_.q) <- u_ID_

colnames(sample_assoc_density_matrix) <- c(1:max(cutree(ww_e_c127, h=cutHere)))
rownames(sample_assoc_density_matrix) <- unlist(u_ID_)

sample_assoc_density_matrix[is.na(sample_assoc_density_matrix)] <- 0

for(x in rownames(sample_assoc_density_matrix)){
  for(y in colnames(sample_assoc_density_matrix)){
    sample_assoc_density_matrix[ x , y ] <- cluster_assoc_density_matrix[ unlist(u_CLabel_sub_)[ unlist(u_ID_sub_) == x ] , y ]
  }
}

#colscale = redgreen(  max(sample_assoc_density_matrix[,1][!is.na(sample_assoc_density_matrix[,1])])*2  )[38:74]
png("20191102.E_licl.tsne.COLORbyQCdensity.clust24.png", width=2000,height=500)
par(mfrow=c(1,6))
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,1]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,1]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,2]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,2]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,3]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,3]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,4]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,4]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,5]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,5]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,6]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,6]] )
dev.off()




png("20191102.E_licl.corrected_c127.tsne.COLORbyQCdensity.clust24.png", width=2000,height=500)
par(mfrow=c(1,6))
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,1]))))  )
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,1]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,2])))))
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,2]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,3])))))
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,3]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,4])))))
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,4]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,5])))))
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,5]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,6])))))
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,6]] )
dev.off()











clust_of_interest <- 13
clust_of_interest.covars <- names(p.pam$clustering[p.pam$clustering == clust_of_interest])
clust_par_relZ.cov.ex.dit_cov_c127 <- vegdist(clust_par_relZ.cov.ex_c127[,clust_of_interest.covars], method = "euclidean")
hy <- heatmap.2( log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1), trace = "none", keysize=0.7, margins = c(10,10), col=rev(redgreen(150))  )

par(mfrow=c(1,1))
ww_e_c127 <- hclust(as.dist(clust_par_relZ.cov.ex.dit_cov_c127), method="complete")

cutHere <- 10
table(unname(cutree(ww_e_c127, h=cutHere)))
hclust.clustcols <- rainbow(max(cutree(ww_e_c127, h=cutHere)))[unname(cutree(ww_e_c127, h=cutHere))]
png("clust_par_relZ.cov.ex.dit_cov_c127.hierarchical.clust13.png", width=800,height=800)
heatmap.2( log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1 ), trace = "none", keysize=0.7, margins = c(5,5), col=rev(redgreen(150)),
           RowSideColors = hclust.clustcols,
           ColSideColors = hclust.clustcols
) 
dev.off()


cluster_assoc_density_matrix <- matrix(ncol=max(cutree(ww_e_c127, h=cutHere)), 
                                       nrow=46)
colnames(cluster_assoc_density_matrix) <- c(1:max(cutree(ww_e_c127, h=cutHere)))
rownames(cluster_assoc_density_matrix) <- c(1:46)
cluster_assoc_density_matrix[is.na(cluster_assoc_density_matrix)] <- 0
for( run_this_clust in 1:max(cutree(ww_e_c127, h=cutHere)) ){
  for( q in names(cutree(ww_e_c127, h=cutHere)[cutree(ww_e_c127, h=cutHere)==run_this_clust]) ){
    if(strsplit(q, ',')[[1]][1] != "NA" && strsplit(q, ',')[[1]][2] != "NA"){
      cluster_assoc_density_matrix[ strsplit(q, ',')[[1]][1] , run_this_clust ] <- cluster_assoc_density_matrix[ strsplit(q, ',')[[1]][1] , run_this_clust ] + 1
      cluster_assoc_density_matrix[ strsplit(q, ',')[[1]][2] , run_this_clust ] <- cluster_assoc_density_matrix[ strsplit(q, ',')[[1]][2] , run_this_clust ] + 1
    }
  }
}

#save.image("20191104devel.Rdata")

sample_assoc_density_matrix <- matrix(ncol=max(cutree(ww_e_c127, h=cutHere)), 
                                      nrow=dim(E_licl)[1])

#u_CLabel_sub_.q <- unlist(u_CLabel_sub_)
#names(u_CLabel_sub_.q) <- u_ID_

colnames(sample_assoc_density_matrix) <- c(1:max(cutree(ww_e_c127, h=cutHere)))
rownames(sample_assoc_density_matrix) <- unlist(u_ID_)

sample_assoc_density_matrix[is.na(sample_assoc_density_matrix)] <- 0

for(x in rownames(sample_assoc_density_matrix)){
  for(y in colnames(sample_assoc_density_matrix)){
    sample_assoc_density_matrix[ x , y ] <- cluster_assoc_density_matrix[ unlist(u_CLabel_sub_)[ unlist(u_ID_sub_) == x ] , y ]
  }
}

#colscale = redgreen(  max(sample_assoc_density_matrix[,1][!is.na(sample_assoc_density_matrix[,1])])*2  )[38:74]
png("20191102.E_licl.tsne.COLORbyQCdensity.clust13.png", width=2000,height=500)
par(mfrow=c(1,5))
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,1]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,1]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,2]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,2]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,3]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,3]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,4]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,4]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,5]))))  )
plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,5]] )
#colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,6]))))  )
#plot( E_licl.tsne$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,6]] )
dev.off()




png("20191102.E_licl.corrected_c127.tsne.COLORbyQCdensity.clust13.png", width=2000,height=500)
par(mfrow=c(1,5))
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,1]))))  )
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,1]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,2])))))
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,2]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,3])))))
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,3]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,4])))))
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,4]] )
colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,5])))))
plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,5]] )
#colscale = redgreen(  max(as.numeric(names(table(sample_assoc_density_matrix[,6])))))
#plot( E_licl.tsne.corrected_c127$Y , pch=15, cex=0.8, col=colscale[sample_assoc_density_matrix[,6]] )
dev.off()





#############################################

primary.colors(46)[ as.numeric(unlist(u_CLabel_sub_)) ]

cold <- list()
cd <- 1
for(cq in strsplit(rownames(as.matrix(clust_par_relZ.cov.ex.dit_sam)),",")){
  if(as.numeric(cq[1]) %in% outliers || as.numeric(cq[2]) %in% outliers){
    cold[cd] <- "white"
  }else{
    cold[cd] <- "purple"
  }
  cd <- cd + 1
}
png("clust13.outlier.png", width=800,height=800)
heatmap.2( log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1), trace = "none", keysize=0.7, margins = c(10,10), col=rev(redgreen(150)),
           RowSideColors=unlist(cold)
)
dev.off()

# 2 : samples (exons)

clust_par_relZ.cov.ex.dit_cov_c127 <- vegdist(clust_par_relZ.cov.ex_c127, method = "euclidean")
png("clust_par_relZ.cov.ex.dit_cov_c127.20191102.png",width=4600,height=4600)
heatmap.2( log(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127)+1), trace = "none", keysize=0.7, margins = c(10,10), col=rev(redgreen(150))  )
dev.off()

# 1 : covariates (introns)

clust_par_relZ.cov.in.dit_cov_c127 <- vegdist(t(clust_par_relZ.cov.in_c127), method = "euclidean")
png("clust_par_relZ.cov.in.dit_cov_c127.20191102.png",width=1400,height=1400)
heatmap.2( as.matrix(clust_par_relZ.cov.in.dit_cov_c127), trace = "none", keysize=0.7, margins = c(20,20), col=rev(redgreen(150))  )
dev.off()

# 2 : samples (introns)

clust_par_relZ.cov.in.dit_sam_127 <- vegdist(clust_par_relZ.cov.in_c127, method = "euclidean")
png("clust_par_relZ.cov.in.dit_sam_127.20191102.png",width=1400,height=1400)
heatmap.2( log(as.matrix(clust_par_relZ.cov.in.dit_sam_127)+1), trace = "none", keysize=0.7, margins = c(20,20), col=rev(redgreen(150))  )
dev.off()
 

# 
# E_AB_C_absZ_.tsv: tab delimited array of cluster-pairs (rows) by covariate (columns). Each array entry indicates the absolute-value of the z-score associated with the (raw) differential-expression of that particular covariate-rank over that cluster pair.
# E_AB_K_absZ_.tsv: tab delimited array of cluster-pairs (rows) by covariate (columns). Each array entry indicates the absolute-value of the z-score associated with the (relative) differential-expression of that particular covariate-rank over that cluster pair. In this context 'relative' means corrected for the linear effect of the E gene-ranks.
# E_AB_Z_absZ_.tsv: tab delimited array of cluster-pairs (rows) by exon genes (columns). Each array entry indicates the absolute-value of the z-score associated with the (raw) differential-expression of that particular gene-rank over that cluster pair.
# E_AB_H_absZ_.tsv: tab delimited array of cluster-pairs (rows) by exon genes (columns). Each array entry indicates the absolute-value of the z-score associated with the (relative) differential-expression of that particular gene-rank over that cluster pair. In this context 'relative' means corrected for the linear effect of the covariate-rank.
# I_AB_C_absZ_.tsv: tab delimited array of cluster-pairs (rows) by covariate (columns). Each array entry indicates the absolute-value of the z-score associated with the (raw) differential-expression of that particular covariate-rank over that cluster pair.
# I_AB_K_absZ_.tsv: tab delimited array of cluster-pairs (rows) by covariate (columns). Each array entry indicates the absolute-value of the z-score associated with the (relative) differential-expression of that particular covariate-rank over that cluster pair. In this context 'relative' means corrected for the linear effect of the I gene-ranks.
# I_AB_Z_absZ_.tsv: tab delimited array of cluster-pairs (rows) by intron genes (columns). Each array entry indicates the absolute-value of the z-score associated with the (raw) differential-expression of that particular gene-rank over that cluster pair.
# I_AB_H_absZ_.tsv: tab delimited array of cluster-pairs (rows) by intron genes (columns). Each array entry indicates the absolute-value of the z-score associated with the (relative) differential-expression of that particular gene-rank over that cluster pair. In this context 'relative' means corrected for the linear effect of the covariate-rank.

pcq.cov <- prcomp(as.matrix(clust_par_relZ.cov.ex.dit_cov), scale = TRUE)
plot(pcq.cov$x, col="white", xlab=paste(apply(pcq.cov$x, 2, var)[1]/sum(apply(pcq.cov$x, 2, var)),"%"), ylab=paste(apply(pcq.cov$x, 2, var)[2]/sum(apply(pcq.cov$x, 2, var)),"%"))
text(pcq.cov$x, labels = C_VariableName_, cex=0.5)

pcq.sam <- prcomp(as.matrix(clust_par_relZ.cov.ex.dit_sam), scale = TRUE)
plot(pcq.sam$x, col="white", xlab=paste(apply(pcq.sam$x, 2, var)[1]/sum(apply(pcq.sam$x, 2, var)),"%"), ylab=paste(apply(pcq.sam$x, 2, var)[2]/sum(apply(pcq.sam$x, 2, var)),"%"))
text(pcq.sam$x, labels = E_CLabel_pair_, cex=0.5)





















##########################################################

# JUNK AND DEVEL CODE BELOW

##########################################################





load("intronANDexon.20190831.Rdata")
prefix <- "intron_and_exon"
count_cmpr <- intron_and_exon.rel

#load("intron.Rdata")
#prefix <- "intron"
#count_cmpr <- intron.counts.ctLo

#load("exon.Rdata")
#prefix <- "exon"
#count_cmpr <- exon.counts.ctLo





###########################################################################################################################################################################################

#save.image("20191102.corrected.Rdata")
#load("20191102.corrected.Rdata")

# 
# #I_corrected_.umap <- umap(I_corrected_, n_neighbors = 100, learning_rate = 0.5, init = "random")
# 
# 
# png("compare_ordination3.20191019.png", width=1400,height=1400)
# par(mfrow=c(2,2))
# plot(E_corrected_.tsne$Y, pch=16, main="TSNE: Exon", col=colorbycluster)
# plot(E_corrected_.umap$layout, pch=16, main="TSNE: Intron", col=colorbycluster)
# plot(I_corrected_.tsne$Y, pch=16, main="UMAP: Exon", col=colorbycluster)
# plot(I_corrected_.umap$layout, pch=16, main="UMAP: Intron", col=colorbycluster)
# dev.off()
#  


#write.table(I_corrected_, file = "I_corrected_20191019.tsv", sep = "\t")
#write.table(E_corrected_, file = "E_corrected_.20191019tsv", sep = "\t")

E_limp.pca <- prcomp(as.matrix(E_limp), scale = TRUE)
E_licl.pca <- prcomp(as.matrix(E_licl), scale = TRUE)
E_licr.pca <- prcomp(as.matrix(E_licr), scale = TRUE)
E_limp.corrected_c012.pca <- prcomp(as.matrix(E_limp.corrected_c012), scale = TRUE)
E_licl.corrected_c012.pca <- prcomp(as.matrix(E_licl.corrected_c012), scale = TRUE)
E_licr.corrected_c012.pca <- prcomp(as.matrix(E_licr.corrected_c012), scale = TRUE)
E_limp.corrected_c127.pca <- prcomp(as.matrix(E_limp.corrected_c127), scale = TRUE)
E_licl.corrected_c127.pca <- prcomp(as.matrix(E_licl.corrected_c127), scale = TRUE)
E_licr.corrected_c127.pca <- prcomp(as.matrix(E_licr.corrected_c127), scale = TRUE)

I_limp.pca <- prcomp(as.matrix(I_limp), scale = TRUE)
I_licl.pca <- prcomp(as.matrix(I_licl), scale = TRUE)
I_licr.pca <- prcomp(as.matrix(I_licr), scale = TRUE)
I_limp.corrected_c012.pca <- prcomp(as.matrix(I_limp.corrected_c012), scale = TRUE)
I_licl.corrected_c012.pca <- prcomp(as.matrix(I_licl.corrected_c012), scale = TRUE)
I_licr.corrected_c012.pca <- prcomp(as.matrix(I_licr.corrected_c012), scale = TRUE)
I_limp.corrected_c127.pca <- prcomp(as.matrix(I_limp.corrected_c127), scale = TRUE)
I_licl.corrected_c127.pca <- prcomp(as.matrix(I_licl.corrected_c127), scale = TRUE)
I_licr.corrected_c127.pca <- prcomp(as.matrix(I_licr.corrected_c127), scale = TRUE)

#plot(I_limp.pca$x,col=colorbycluster,pch=16)
library(cluster)


png("20191102.I_limp.pca.png", width=800,height=800)
plot(I_limp.pca$x, main="TSNE: Exon", col=colorbycluster,pch=16)
dev.off()
png("20191102.I_limp.pca.pam.png", width=800,height=800)
max <- 0
max.q <- 3
for(q in c(3:50)){
  print(q)
  I_limp.pca.pam <- pam(as.data.frame(I_limp.pca$x[,1:2]), q)
  if( I_limp.pca.pam $ silinfo $avg.width > max){
    max <- I_limp.pca.pam $ silinfo $avg.width
    max.q <- q
  }
}
I_limp.pca.pam <- pam(as.data.frame(I_limp.pca$x[,1:2]), max.q)
clusplot(I_limp.pca.pam, col.p=colorbycluster,color=T,main=max.q,metric="euclidean",shade=F,labels=0,lines=0)
dev.off()
png("20191102.I_limp.pca.pam46.png", width=800,height=800)
I_limp.pca.pam <- pam(as.data.frame(I_limp.pca$x[,1:2]), 46)
clusplot(I_limp.pca.pam, col.p=colorbycluster,color=T,metric="euclidean",shade=F,labels=0,lines=0)
dev.off()

png("20191102.E_limp.pca.png", width=800,height=800)
plot(E_limp.pca$x, col=colorbycluster,pch=16)
dev.off()
png("20191102.E_limp.pca.pam.png", width=800,height=800)
max <- 0
max.q <- 3
for(q in c(3:50)){
  print(q)
  E_limp.pca.pam <- pam(as.data.frame(I_limp.pca$x[,1:2]), q)
  if( E_limp.pca.pam $ silinfo $avg.width > max){
    max <- E_limp.pca.pam $ silinfo $avg.width
    max.q <- q
  }
}
E_limp.pca.pam <- pam(as.data.frame(E_limp.pca$x[,1:2]), max.q)
clusplot(E_limp.pca.pam, col.p=colorbycluster,color=T,main=max.q,metric="euclidean",shade=F,labels=0,lines=0)
dev.off()
png("20191102.E_limp.pca.pam46.png", width=800,height=800)
E_limp.pca.pam <- pam(as.data.frame(E_limp.pca$x[,1:2]), 46)
clusplot(E_limp.pca.pam, col.p=colorbycluster,color=T,metric="euclidean",shade=F,labels=0,lines=0)
dev.off()



png("20191102.I_licl.pca.png", width=800,height=800)
plot(I_licl.pca$x, col=colorbycluster,pch=16)
dev.off()
png("20191102.I_licl.pca.pam.png", width=800,height=800)
max <- 0
max.q <- 3
for(q in c(3:50)){
  print(q)
  I_licl.pca.pam <- pam(as.data.frame(I_licl.pca$x[,1:2]), q)
  if( I_licl.pca.pam $ silinfo $avg.width > max){
    max <- I_licl.pca.pam $ silinfo $avg.width
    max.q <- q
  }
}
I_licl.pca.pam <- pam(as.data.frame(I_licl.pca$x[,1:2]), max.q)
clusplot(I_licl.pca.pam, col.p=colorbycluster,color=T,main=max.q,metric="euclidean",shade=F,labels=0,lines=0)
dev.off()
png("20191102.I_licl.pca.pam46.png", width=800,height=800)
I_licl.pca.pam <- pam(as.data.frame(I_licl.pca$x[,1:2]), 46)
clusplot(I_licl.pca.pam, col.p=colorbycluster,color=T,metric="euclidean",shade=F,labels=0,lines=0)
dev.off()

png("20191102.E_licl.pca.png", width=800,height=800)
plot(E_licl.pca$x, col=colorbycluster,pch=16)
dev.off()
png("20191102.E_licl.pca.pam.png", width=800,height=800)
max <- 0
max.q <- 3
for(q in c(3:50)){
  print(q)
  E_licl.pca.pam <- pam(as.data.frame(E_licl.pca$x[,1:2]), q)
  if( E_licl.pca.pam $ silinfo $avg.width > max){
    max <- E_licl.pca.pam $ silinfo $avg.width
    max.q <- q
  }
}
E_licl.pca.pam <- pam(as.data.frame(E_licl.pca$x[,1:2]), max.q)
clusplot(E_licl.pca.pam, col.p=colorbycluster,color=T,main=max.q,metric="euclidean",shade=F,labels=0,lines=0)
dev.off()
png("20191102.E_licl.pca.pam46.png", width=800,height=800)
E_licl.pca.pam <- pam(as.data.frame(E_licl.pca$x[,1:2]), 46)
clusplot(E_licl.pca.pam, col.p=colorbycluster,color=T,metric="euclidean",shade=F,labels=0,lines=0)
dev.off()







p <- prcomp(as.matrix(clust_par_relZ.cov.ex.dit_cov_c127), scale=TRUE)
plot(p$x[,1:2], col="white")
text(p$x[,1:2], colnames(clust_par_relZ.cov.ex_c127),cex=0.4)
max <- 0
max.q <- 3
avgsil <- list()
for(q in c(3:50)){
  print(q)
  p.pam <- pam(as.data.frame(p$x[,1:2]), q)
  if( p.pam $ silinfo $avg.width > max){
    max <- p.pam $ silinfo $avg.width
    max.q <- q
  }
  avgsil[q-2] <- p.pam $ silinfo $avg.widt
}
names(avgsil) <- c(3:50)
colsil <- c(rep("grey",50))
colsil[names(avgsil) %in% c(42,27,13,7)] <- "red"
barplot(unlist(avgsil),horiz = TRUE, cex.names = 0.6, las=2, col=colsil)

par(mfrow=c(2,2))
for(cc in c(42,27,13,7)){
  p.pam <- pam(as.data.frame(p$x[,1:2]), cc)
  clusplot(p.pam, color=T,main=cc,metric="euclidean",shade=F,labels=0,lines=0)
}

library(ape)
library(colorRamps)

cc <- 27
par(mfrow=c(1,1))
p.pam <- pam(as.data.frame(p$x[,1:2]), cc)
clusplot(p.pam, color=T,main=cc,metric="euclidean",shade=F,labels=2,lines=0,cex=0.5)
p.pam.clustcols <- p.pam$clustering
for(r in 1:cc){
  p.pam.clustcols[p.pam.clustcols == r] <- primary.colors(cc)[r]
}

ww <-hclust(as.dist(clust_par_relZ.cov.ex.dit_cov_c127), method="complete")
hclust.clustcols <- primary.colors(max(cutree(ww, h=67)))[unname(cutree(ww, h=67))]
heatmap.2( as.matrix(clust_par_relZ.cov.ex.dit_cov_c127), trace = "none", keysize=0.7, margins = c(20,20), col=rev(redgreen(150)),
           RowSideColors = p.pam.clustcols,
           ColSideColors = hclust.clustcols
) 

#plot( as.tree(ww), col=hclust.clustcols)

