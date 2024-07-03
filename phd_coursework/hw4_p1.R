#########
S = matrix (c(50, 40, 40, 50), nrow=2, ncol=2)
eigenvalues <- eigen(S)$values
eigenvectors <- eigen(S)$vectors
PC <- as.matrix(S) %*% eigenvectors
print(round(eigenvalues/sum(eigenvalues) * 100, digits = 2))
round(cumsum(eigenvalues)/sum(eigenvalues) * 100, digits = 2)

#########

####### LOGISTIC REGRESSION

# ALL DATA
alldata <- read.csv(file="/Users/apple/Desktop/BNFO285_stat_learning/BNFO285_HW4_V3.txt",head=FALSE,sep=",")
logreg <- glm(V1 ~ ., data = alldata, family = "binomial")
results_df <-summary.glm(logreg)

# REDUCE
rownames(alldata) <- alldata[[1]]
input_set <- alldata[2:length(alldata)]
#update header
row_names <- substr(rownames(input_set),1,5)
# reduce to uniq receptor ids
uniq_set <- input_set[(row_names != 'CCR5/CXCR4'),]
row_names <- row_names[row_names != 'CCR5/CXCR4']
# and push into binary space
seq_data_bin <- as.numeric(row_names == 'CXCR4')
#regression
reduced_logreg <- glm(seq_data_bin~ ., data = uniq_set, family = "binomial")
reduced_results_df <-summary.glm(logreg)
 
####### CROSS VALIDATION
# install.packages("cvTools", dep = TRUE)
# library(cvTools)
# lrcross <- cvFit(reduced_logreg, data = input_set, y = seq_data_bin~ ., K = 10, R = 5)
folds <- cut(seq(1,nrow(uniq_set)),breaks=10,labels=FALSE)
#Perform 10 fold cross validation
for(i in 1:10){
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- uniq_set[testIndexes, ]
  trainData <- uniq_set[-testIndexes, ]
}

####### SUPPORT VECTOR MACHINE
# install.packages("e1071", dep = TRUE)
library(e1071)
reduced_svm <- svm(seq_data_bin, data = uniq_set, family = "binomial")
reduced_svm.pred <- predict(reduced_svm, uniq_set[,36])
plot(predict(reduced_svm))


