library(edgeR)

setwd ("/Users/apple/Desktop/BNFO285_stat_learning/hw2_20150204/")
data = read.csv('hw2-p6-data.csv',header=FALSE);

id <- NULL
A <- NULL
B <- NULL
C <- NULL
D <- NULL

for(i in 1:nrow(data)){
  id[i] = data[i,'V1']
  A[i] = data[i,'V2']
  B[i] = data[i,'V3']
  C[i] = data[i,'V4']
  D[i] = data[i,'V5']
}
all = c(A,B,C,D)
max = max(all)
median = median(all)
cutoff = 0.1*(max-median)+median

passA = 0
boolA <- NULL
passB = 0
boolB <- NULL
passC = 0
boolC <- NULL
passD = 0
boolD <- NULL
for(i in 1:nrow(data)){
  if (A[i] >= cutoff){ 
    passA = passA + 1
    boolA[i] = 1
  }else{
    boolA[i] = 0
  }
  if (B[i] >= cutoff){ 
    passB = passB + 1  
    boolB[i] = 1
  }else{
    boolB[i] = 0
  }
  if (C[i] >= cutoff){ 
    passC = passC + 1  
    boolC[i] = 1
  }else{
    boolC[i] = 0
  }
  if (D[i] >= cutoff){ 
    passD = passD + 1  
    boolD[i] = 1
  }else{
    boolD[i] = 0
  }
}
# PAIRINGS - AB AC AD BC BD CD
AB = binomTest(passA, passB, nrow(data), nrow(data), p=0.5)
AC = binomTest(passA, passC, nrow(data), nrow(data), p=0.5)
AD = binomTest(passA, passD, nrow(data), nrow(data), p=0.5)
BC = binomTest(passB, passC, nrow(data), nrow(data), p=0.5)
BD = binomTest(passB, passD, nrow(data), nrow(data), p=0.5)
CD = binomTest(passC, passD, nrow(data), nrow(data), p=0.5)

print ("AB:")
print (AB)
print ("AC:")
print (AC)
print ("AD:")
print (AD)
print ("BC:")
print (BC)
print ("BD:")
print (BD)
print ("CD:")
print (CD)



# FISHERS TEST EXCEEDS AVAILABLE MEMORY
#fisher.test(A,B, workspace = 400000, hybrid = FALSE,
#            control = list(), or = 1, alternative = "two.sided",
#            conf.int = TRUE, conf.level = 0.95,
#            simulate.p.value = FALSE, B = 2000)
#fisher.test(A,C)
#fisher.test(A,D)
#fisher.test(B,C)
#fisher.test(B,D)
#fisher.test(C,D)