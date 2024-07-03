
setwd ("/Users/apple/Desktop/CSE280A/hw2")

plotme <- read.table("plotme.txt", sep=" ", header=FALSE)
plot(plotme, type="l", col="darkblue", lwd=10, xlab = "Generations", ylab = "Samples", main = "2b. Coalescent Rate Evaluation n=100, N=10e6, a=0.0175)")
