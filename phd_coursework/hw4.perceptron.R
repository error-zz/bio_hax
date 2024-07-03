
#################
# all input

setwd ("/Users/apple/Desktop/cse151")
dictionary <- read.table("hw4dictionary.txt", sep=" ", header=FALSE)
dictionary <- dictionary[,1]

training_set <- read.table("hw4train.txt", sep=" ", header=FALSE)
training_frequencies <- training_set[,1:819]
training_labels <- training_set[,820]
names(training_labels) <- rownames(training_frequencies)

test_set <- read.table("hw4test.txt", sep=" ", header=FALSE)
test_frequencies <- test_set[,1:819]
test_labels <- test_set[,820]
names(test_labels) <- rownames(test_frequencies)

# plot word appearance (y) vs. post string (x)
#plot(as.numeric(as.matrix(test_frequencies)),main="Test All Groups",col=c(rainbow(6))[test_labels])
#plot(as.numeric(as.matrix(training_frequencies)),main="Training All Groups",col=c(rainbow(6))[training_labels])

#################
# 1 and 2 

# select by content
training_set_1and2 = training_set[ which(training_set$V820==1 | training_set$V820==2), ]
#rownames(training_frequencies_1and2) <- rownames(training_set)
training_frequencies_1and2 <- training_set_1and2[,1:819]
#training_labels_1and2 <- training_set_1and2[,820]
training_labels_1and2 <- c(1,-1)[ training_set_1and2[,820] ]

# training subroutine # 
#std_train <- function(training_weights, num){
training_weights <- matrix(0,1,819)
training_weights <- unlist(as.list(training_weights))
num <- 1090
  pos <- matrix(ncol=length(training_weights), nrow=num)
  neg <- matrix(ncol=length(training_weights), nrow=num)
  pos_i <- 1
  neg_i <- 1
  for(t in 1:num){
    #yt_wtxt <- as.numeric(as.matrix(training_labels_1and2)) * ( as.numeric(as.matrix(training_weights)) %*% as.numeric(as.matrix(unlist(training_frequencies_1and2[t,]))) )
    if (pos_i <= num && (as.numeric(as.matrix(training_weights)) %*% as.numeric(as.matrix(unlist(training_frequencies_1and2[t,]))) < 0)){
      pos[pos_i,]<-as.numeric(as.matrix(unlist(training_frequencies_1and2[t,])))
      pos_i <- pos_i + 1
    }
    if (neg_i <= num && (as.numeric(as.matrix(training_weights)) %*% as.numeric(as.matrix(unlist(training_frequencies_1and2[t,]))) > 0)){
        neg[neg_i,]<-as.numeric(as.matrix(unlist(training_frequencies_1and2[t,])))
      neg_i <- neg_i + 1
    }
  }
  out <- list(pos,neg)
  names(out)<-c(-1,1)
  out
#}


# testing subroutine # 
std_test <- function(training_frequencies_1and2, trained_weights){
  pos <- matrix(ncol=length(training_weights), nrow=1)
  neg <- matrix(ncol=length(training_weights), nrow=1)
}


std_train(training_weights,1090)
