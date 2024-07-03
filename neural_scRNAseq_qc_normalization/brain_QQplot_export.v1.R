# install.packages("ggplot2",repos='http://cran.us.r-project.org')
library(ggplot2)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
expression <- read.table( args[1], sep="\t", header=TRUE )
prefix <- args[2]

qc <- expression[,2:30]
exp <- expression[,31:ncol(expression)]
list_of_lms <- list()

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

linear_regression  <- function(geneExpression, QC_dataframe) {
  
  count = apply( geneExpression, 2, function(x) {
    reg  <- lm( x ~ as.factor(A) + as.factor(B) + as.factor(C) + as.factor(D) + as.factor(E) + F + G + H + I + J + K + L + M + N + O + P + Q + R + S + T + U + V + W + X + Y + Z + AA + AB + AC, data=data.frame(QC_dataframe) ) 
    list_of_lms[[length(list_of_lms)+1]] <- reg
    t(coef(summary(reg)))[4,]
    
  }) 
} 

pvals = rep(NA, ncol(exp ))
names(pvals) =colnames(exp)
for(i in 1:ncol(exp))
{
  reg  <- lm( as.numeric(exp[,i]) ~ as.factor(A) + as.factor(B) + as.factor(C) + as.factor(D) + as.factor(E) + F + G + H + I + J + K + L + M + N + O + P + Q + R + S + T + U + V + W + X + Y + Z + AA + AB + AC, data=data.frame(qc) ) 
  list_of_lms[[length(list_of_lms)+1]] <- reg
  pvals[i] = t(coef(summary(reg)))[4,]
  
}

pval <- linear_regression(exp,qc)

P = pval[2,]
n = length(P)
null = sort(-log10((1:n)/(n+1)))
expect = sort(-log10(P))
d1 <- data.frame(null, expect,name =rownames(pval)[2])

P = pval[3,]
n = length(P)
null = sort(-log10((1:n)/(n+1)))
expect = sort(-log10(P))
d2 <- data.frame(null, expect,name = rownames(pval)[3])
d3 <- rbind(d1,d2)

for (name in rownames(pval)[4:nrow(pval)]) {
  P = pval[name,]
  n = length(P)
  null = sort(-log10((1:n)/(n+1)))
  expect = sort(-log10(P))
  d2 <- data.frame(null, expect,name = name)
  d3 <- rbind(d3,d2)
}

tmpTITLE <- paste(prefix, ".qq_plot.png",sep="") 
png(filename=tmpTITLE, width = 1280, height = 800)
#png(filename="qq_plot.png", width = 1280, height = 800)

tmpTITLE <- paste(prefix, " QC Metric QQ-Plot",sep="") 
plot <- ggplot(d3, aes(null, expect)) + geom_point(aes(color=name)) + ggtitle(tmpTITLE) +
  geom_abline(intercept = 0, slope = 1) + ylab("Observed -log10(P)") + xlab("Null -log10(P)")
plot

dev.off()

# save.image()
tmpTITLE <- paste(prefix, ".Rdata",sep="") 
save.image(file=tmpTITLE)

### !single sample test case!
# test = lm( exp[,1] ~ as.factor(A) + as.factor(B) + as.factor(C) + as.factor(D) + as.factor(E) + F + G + H + I + J + K + L + M + N + O + P + Q + R + S + T + U + V + W + X + Y + Z + AA + AB + AC, data=data.frame(qc) )
# test.list = list(test)
# coef(summary(test.list[[1]]))

