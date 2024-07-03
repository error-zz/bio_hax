#################################################
#### math283 - hw7 - j.mccorrison - 11/26/14 ####
#################################################

#################
#### NOTES ######
#################

# open multiple rstudio sessions via the terminal
# open -n -a "rstudio"

# microarray data
# 9216 probes per slide, 84 slides 
# 65 samples of various types of breast cancer tumors
# 42 individuals (some before/after chemo)
# 29 samples of other issues

# our data is a SUBSET of 16 of 84 slides
# 8 patients with "liminal-like ER+ tumors"
# 4593x17 - no headers

# for this assignment, rows with any bad readings in the selected microarrays were deleted.
## Column 1 is the spot ID number on the microarray slide, between 1 and 9216.
## Columns 2 through 9 are 8 patients before treatment.
## Columns 10 through 17 are those same patients, in the same order, after treatment.
## The numbers in columns 2–17 are normalized log intensity ratios between the red and green channels.
## red = tumor, green = control

#################
#### CODE ######
#################

# Set CWD
setwd ("/Users/apple/Desktop/EXE")

# Parse input tsv
# PROD USAGE
data = read.table('hw7data.txt',header=FALSE,sep='\t');
# DEBUG : TEST CASE
# data = read.table('test250.txt',header=FALSE,sep='\t');

all <- NULL

# parse the input by row
# build list of lists (all) for parsing in next subroutine
for(r in 1:nrow(data)){
  
  # declare constants for each row
  spotID = data[r,1]
  before = as.double(data[r,2:9])
  after = as.double(data[r,10:17])
  
  # define variance
  # by row - incorrect?
  # pooled_sample_var = (((7)*var(before))+((7)*var(after)))/(14)
  beforeVAR <- apply(data[,2:9],1,var)
  afterVAR <- apply(data[,10:17],1,var)
  pooled_sample_var = (sum(beforeVAR)+sum(afterVAR))/length(all)
  
  # (1) 2-sample z test - from microarray slide 19
  # definition from slide 19 - defies instruction in hw to use single SD definition
  # using as X and Y are independent
  # sampleZ = (mean(before) - mean(after)) / (sqrt(var(before)/8 + var(after)/8))
  # OR
  # using as X and Y are dependent
  sampleZ = (mean(before) - mean(after))/(sqrt(pooled_sample_var)*sqrt(0.25))
  tmpA = pnorm(sampleZ)
  tmpB = 1 - pnorm(sampleZ)
  ztestP = 2*min(tmpA, tmpB)
  
  # (4) paired z test - from microarray slide 25
  dsum = abs(before[1]-after[1])+abs(before[2]-after[2])+abs(before[3]-after[3])+abs(before[4]-after[4])+abs(before[5]-after[5])+abs(before[6]-after[6])+abs(before[7]-after[7])+abs(before[8]-after[8])
  paired_sampleZ = (dsum/8)/((sqrt(pooled_sample_var))*sqrt(8))
  tmpA = pnorm(paired_sampleZ)
  tmpB = 1 - pnorm(paired_sampleZ)
  paired_ztestP = 2*min(tmpA, tmpB)
  #paired_ztestP = pnorm(paired_sampleZ)

  # (2) 2 sample t test
  # t-test definition from slide 20...
  # sampleT = (mean(before) - mean(after))/(pooled_sample_var*sqrt(14))
  test = t.test(before, after, var.equal=TRUE); # use this line for unpaired
  ttestP = test$p.value;
  ttestT = test$statistic;
  # paired t test
  paired_test = t.test(before, after, paired=TRUE); # use this line for paired
  paired_ttestP = paired_test$p.value;
  paired_ttestT = paired_test$statistic;
  
  # (3) mann-whitney u test
  mwtest = wilcox.test(before, after);
  mannwhitneyP = mwtest$p.value;
  mannwhitneyU = mwtest$statistic;

  # (6) w = wilcox solution
  #wtest = wilcox.test(sample, mu=sampleMEAN)
  wtest = wilcox.test(before-after,0);
  wilcoxP = wtest$p.value
  wilcoxW = wtest$statistic;
  
  # store all vars for use downstream
  # array index:
  #        1        2           3        4          5          6           7        8       9            10       11            12         13
  a <- c(spotID, sampleZ,paired_sampleZ,ttestT,paired_ttestT,mannwhitneyU,wilcoxW,ztestP,paired_ztestP,ttestP,paired_ttestP,mannwhitneyP,wilcoxP)
  all[spotID] <- list(a)

  ###############
  # !!!
  # !!! RESULTS (contents of list of lists "all") PRINTED WITH OUTPUT D (BELOW) !!!
  # !!!
  ###############  
}

##########
# SOLN B #
##########
file.remove("B_soln.tsv")
for(r in 1:length(all)){
  if(all[r][1] != "NULL"){
    tmplist = all[r][1]
    cmprID = tmplist[[1]][1]
    if ( (cmprID >= 1981) && (cmprID <= 1985) ){
    #if (cmprID >= 1981){
    #  if (cmprID <= 1985){
        tmp = toString(all[r][1])
        write(tmp, file="B_soln.tsv", append="TRUE")    
      }
    #}
  }
}

##########
# SOLN C #
##########
newall <- NULL
size_all = 0
file.remove("C_soln.tsv")
write("spotID, sampleZ,paired_sampleZ,ttestT,paired_ttestT,mannwhitneyU,wilcoxW,ztestP,paired_ztestP,ttestP,paired_ttestP,mannwhitneyP,wilcoxP", file="C_soln.tsv", append="TRUE")
for(r in 1:length(all)){
  if(all[r][1] != "NULL"){
    size_all = size_all + 1
    tmp = toString(all[r][1])
    newall[size_all] = all[r][1]
    write(tmp, file="C_soln.tsv", append="TRUE")
  }
}
# # UNIX METHODS TO PARSE OUTPUT OF ALL VALUES FOR DELIVERABLE :
# cat C_soln.tsv | tr ',' ' ' | tr 'c' ' ' | tr '(' ' ' | sort -nk 8  | awk '{print $1, $8}' | head -n 10 > C_top_ztest.txt
# cat C_soln.tsv | tr ',' ' ' | tr 'c' ' ' | tr '(' ' ' | sort -nk 9  | awk '{print $1, $9}' | head -n 10 > C_top_paired_ztestP.txt
# cat C_soln.tsv | tr ',' ' ' | tr 'c' ' ' | tr '(' ' ' | sort -nk 10  | awk '{print $1, $10}' | head -n 10 > C_top_ttestP.txt
# cat C_soln.tsv | tr ',' ' ' | tr 'c' ' ' | tr '(' ' ' | sort -nk 11  | awk '{print $1, $11}' | head -n 10 > C_top_paired_ttestP.txt
# cat C_soln.tsv | tr ',' ' ' | tr 'c' ' ' | tr '(' ' ' | sort -nk 12  | awk '{print $1, $12}' | head -n 10 > C_top_mannwhitneyP.txt
# cat C_soln.tsv | tr ',' ' ' | tr 'c' ' ' | tr '(' ' ' | tr ')' ' ' | sort -nk 13  | awk '{print $1, $13}' | head -n 10 > C_top_wilcoxP.txt


##########
# F + H  # 
##########
#         1        2           3        4          5          6           7        8       9            10       11            12         13
# a <- c(spotID, sampleZ,paired_sampleZ,ttestT,paired_ttestT,mannwhitneyU,wilcoxW,ztestP,paired_ztestP,ttestP,paired_ttestP,mannwhitneyP,wilcoxP)

################################# 8 vs. 11
Q <- NULL
Y <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][8]
  Y[size_all] = newall[[r]][11]
}
png(filename="F_a_1.png")
plot( Q, type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) 2-sample Z test [purple] vs. (Correct) Paired T-Test [black]" )
par(new=TRUE)
plot( Y, type="l", col="black", ylim=c(0, 1), )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="F_a_2.png")
plot( sort(Q), type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) SORTED 2-sample Z test [purple] vs. (Correct) Paired T-Test [black]"  )
par(new=TRUE)
plot( sort(Y), type="l", col="black", ylim=c(0, 1) )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="H_a.png")
hist( Q, main="(P-VALUES, HISTOGRAM) 2-sample Z test", breaks=50, col="purple" )
dev.off()

################################# 9 vs. 11
Q <- NULL
Y <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][9]
  Y[size_all] = newall[[r]][11]
}
png(filename="F_b_1.png")
plot( Q, type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) Paired Sample Z test [purple] vs. (Correct) Paired T-Test [black]" )
par(new=TRUE)
plot( Y, type="l", col="black", ylim=c(0, 1), )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="F_b_2.png")
plot( sort(Q), type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) SORTED Paired Sample Z test [purple] vs. (Correct) Paired T-Test [black]"  )
par(new=TRUE)
plot( sort(Y), type="l", col="black", ylim=c(0, 1) )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="H_b.png")
hist( Q, main="(P-VALUES, HISTOGRAM) Paired Sample Z test", breaks=50, col="purple" )
dev.off()
################################# 10 vs. 11'
Q <- NULL
Y <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][10]
  Y[size_all] = newall[[r]][11]
}
png(filename="F_c_1.png")
plot( Q, type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) 2-sample T test [purple] vs. (Correct) Paired T-Test [black]" )
par(new=TRUE)
plot( Y, type="l", col="black", ylim=c(0, 1), )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="F_c_2.png")
plot( sort(Q), type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) SORTED 2-sample T test [purple] vs. (Correct) Paired T-Test [black]"  )
par(new=TRUE)
plot( sort(Y), type="l", col="black", ylim=c(0, 1) )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="H_c.png")
hist( Q, main="(P-VALUES, HISTOGRAM) 2-sample T test", breaks=50, col="purple" )
dev.off()
################################# 11 vs. 11
Q <- NULL
Y <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][11]
  Y[size_all] = newall[[r]][11]
}
png(filename="F_d_1.png")
plot( Q, type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) Paired T-Test [purple] vs. (Correct) Paired T-Test [black]" )
par(new=TRUE)
plot( Y, type="l", col="black", ylim=c(0, 1), )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="F_d_2.png")
plot( sort(Q), type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) SORTED Paired T-Test [purple] vs. (Correct) Paired T-Test [black]"  )
par(new=TRUE)
plot( sort(Y), type="l", col="black", ylim=c(0, 1) )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="H_d.png")
hist( Q, main="(P-VALUES, HISTOGRAM) Paired T-Test", breaks=50, col="purple" )
dev.off()
################################# 12 vs. 11
Q <- NULL
Y <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][12]
  Y[size_all] = newall[[r]][11]
}
png(filename="F_e_1.png")
plot( Q, type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) Mann Whitney U Test [purple] vs. (Correct) Paired T-Test [black]" )
par(new=TRUE)
plot( Y, type="l", col="black", ylim=c(0, 1), )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="F_e_2.png")
plot( sort(Q), type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) SORTED Mann Whitney U Test red] vs. (Correct) Paired T-Test [black]"  )
par(new=TRUE)
plot( sort(Y), type="l", col="black", ylim=c(0, 1) )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="H_e.png")
hist( Q, main="(P-VALUES, HISTOGRAM) Mann Whitney U Test", breaks=50, col="purple" )
dev.off()
################################# 13 vs. 11
Q <- NULL
Y <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][13]
  Y[size_all] = newall[[r]][11]
}
png(filename="F_f_1.png")
plot( Q, type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) Wilcoxon test [purple] vs. (Correct) Paired T-Test [black]" )
par(new=TRUE)
plot( Y, type="l", col="black", ylim=c(0, 1), )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="F_f_2.png")
plot( sort(Q), type="l", col="purple", ylim=c(0, 1), main="(P-VALUES) SORTED Wilcoxon test [purple] vs. (Correct) Paired T-Test [black]"  )
par(new=TRUE)
plot( sort(Y), type="l", col="black", ylim=c(0, 1) )
par(new=TRUE)
abline(h=0.3, col="black")
par(new=TRUE)
abline(h=0.5, col="maroon")
dev.off()
png(filename="H_f.png")
hist( Q, main="(P-VALUES, HISTOGRAM) Wilcoxon test", breaks=50, col="purple" )
dev.off()


##########
# G  # 
##########
################################# 2
Q <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][2]
}
png(filename="G_a.png")
hist( Q, main="(TEST SCORE, HISTOGRAM) 2-sample Z score", breaks=50, col="blue" )
dev.off()
################################# 3
Q <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][3]
}
png(filename="G_b.png")
hist( Q, main="(TEST SCORE, HISTOGRAM) 2-sample Z score", breaks=50, col="blue" )
dev.off()
################################# 4
Q <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][4]
}
png(filename="G_c.png")
hist( Q, main="(TEST SCORE, HISTOGRAM) 2-sample Z score", breaks=50, col="blue" )
dev.off()
################################# 5
Q <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][5]
}
png(filename="G_d.png")
hist( Q, main="(TEST SCORE, HISTOGRAM) 2-sample Z score", breaks=50, col="blue" )
dev.off()
################################# 6
Q <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][6]
}
png(filename="G_e.png")
hist( Q, main="(TEST SCORE, HISTOGRAM) 2-sample Z score", breaks=50, col="blue" )
dev.off()
################################# 7
Q <- NULL
size_all = 0
for(r in 1:length(newall)){
  size_all = size_all + 1
  Q[size_all] = newall[[r]][7]
}
png(filename="G_f.png")
hist( Q, main="(TEST SCORE, HISTOGRAM) 2-sample Z score", breaks=50, col="blue" )
dev.off()

############
# SOLN D+A #
############
error2_z = 0
error2_paired_z = 0
error2_t = 0
error2_paired_t = 0
error2_mw = 0
error2_w = 0
passfailall <- NULL
passfailbonferroni <- NULL
passfailsidak <- NULL
# PREP FOR WRITE TO FILE
file.remove("A_cutoff_one.tsv")
file.remove("D_sidak_one.tsv")
file.remove("D_bonferroni_one.tsv")
file.remove("A_cutoff_five.tsv")
file.remove("D_sidak_five.tsv")
file.remove("D_bonferroni_five.tsv")
file.remove("A_cutoff_ten.tsv")
file.remove("D_sidak_ten.tsv")
file.remove("D_bonferroni_ten.tsv")

correct <- NULL

for(r in 1:length(newall)){
  if(all[r][1] != "NULL"){    
    #         1        2           3        4          5          6           7        8       9            10       11            12         13
    # a <- c(spotID, sampleZ,paired_sampleZ,ttestT,paired_ttestT,mannwhitneyU,wilcoxW,ztestP,paired_ztestP,ttestP,paired_ttestP,mannwhitneyP,wilcoxP)
    ztestP = all[[r]][8]
    paired_ztestP = all[[r]][9]
    ttestP = all[[r]][10]
    paired_ttestP = all[[r]][11]
    mannwhitneyP = all[[r]][12]
    wilcoxP = all[[r]][13]
    # DEBUG
    #print("ztestP", ztestP)
    #print(ztestP)
    #print("paired_ztestP", paired_ztestP)
    #print("ttestP", ttestP)
    #print("paired_ttestP", ztestP)
    #print("mannwhitneyP", ztestP)
    #print("wilcoxP", ztestP)
    
    ########
    # 0.03 #
    #correct_example
    cutoff = 0.03
    if (paired_ttestP > cutoff){ correct[r] = 0 }
    else{ correct[r] = 1 }
    
    ########
    # 0.01 #
    #boolean expression table
    onoff <- NULL
    cutoff = 0.01
    bonferroni <- cutoff/length(newall)
    onoff_bonferroni <- NULL
    sidak <- 1-(1-cutoff)^(1/length(newall))
    onoff_sidak <- NULL
    # if over cutoff, call as 
    if (ztestP > cutoff){ onoff["ztestP"] = 0 }
    else{ onoff["ztestP"] = 1 }
    if (paired_ztestP > cutoff){ onoff["paired_ztestP"] = 0 }
    else{ onoff["paired_ztestP"] = 1 }
    if (ttestP > cutoff){ onoff["ttestP"] = 0 }
    else{ onoff["ttestP"] = 1 }
    if (paired_ttestP > cutoff){ onoff["paired_ttestP"] = 0 }
    else{ onoff["paired_ttestP"] = 1 }
    if (mannwhitneyP > cutoff){ onoff["mannwhitneyP"] = 0 }
    else{ onoff["mannwhitneyP"] = 1 }
    if (wilcoxP > cutoff){ onoff["wilcoxP"] = 0 }
    else{ onoff["wilcoxP"] = 1 }
    # if over bonferroni, call as 
    if (ztestP > bonferroni){ onoff_bonferroni["ztestP"] = 0 }
    else{ onoff_bonferroni["ztestP"] = 1 }
    if (paired_ztestP > bonferroni){ onoff_bonferroni["paired_ztestP"] = 0 }
    else{ onoff_bonferroni["paired_ztestP"] = 1 }
    if (ttestP > bonferroni){ onoff_bonferroni["ttestP"] = 0 }
    else{ onoff_bonferroni["ttestP"] = 1 }
    if (paired_ttestP > bonferroni){ onoff_bonferroni["paired_ttestP"] = 0 }
    else{ onoff_bonferroni["paired_ttestP"] = 1 }
    if (mannwhitneyP > bonferroni){ onoff_bonferroni["mannwhitneyP"] = 0 }
    else{ onoff_bonferroni["mannwhitneyP"] = 1 }
    if (wilcoxP > bonferroni){ onoff_bonferroni["wilcoxP"] = 0 }
    else{ onoff_bonferroni["wilcoxP"] = 1 }
    # if over sidak, call as 
    if (ztestP > sidak){ onoff_sidak["ztestP"] = 0 }
    else{ onoff_sidak["ztestP"] = 1 }
    if (paired_ztestP > sidak){ onoff_sidak["paired_ztestP"] = 0 }
    else{ onoff_sidak["paired_ztestP"] = 1 }
    if (ttestP > sidak){ onoff_sidak["ttestP"] = 0 }
    else{ onoff_sidak["ttestP"] = 1 }
    if (paired_ttestP > sidak){ onoff_sidak["paired_ttestP"] = 0 }
    else{ onoff_sidak["paired_ttestP"] = 1 }
    if (mannwhitneyP > sidak){ onoff_sidak["mannwhitneyP"] = 0 }
    else{ onoff_sidak["mannwhitneyP"] = 1 }
    if (wilcoxP > sidak){ onoff_sidak["wilcoxP"] = 0 }
    else{ onoff_sidak["wilcoxP"] = 1 }
    passfailall["one"] <- list(onoff)
    passfailbonferroni["one"] <- list(onoff_bonferroni)
    passfailsidak["one"] <- list(onoff_sidak)
    
    #######
    # 0.10 #
    #boolean expression table
    onoff <- NULL
    cutoff = 0.1
    bonferroni <- cutoff/size_all
    onoff_bonferroni <- NULL
    sidak <- 1-(1-cutoff)^(1/size_all)
    onoff_sidak <- NULL
    # if over cutoff, call as 
    if (ztestP > cutoff){ onoff["ztestP"] = 0 }
    else{ onoff["ztestP"] = 1 }
    if (paired_ztestP > cutoff){ onoff["paired_ztestP"] = 0 }
    else{ onoff["paired_ztestP"] = 1 }
    if (ttestP > cutoff){ onoff["ttestP"] = 0 }
    else{ onoff["ttestP"] = 1 }
    if (paired_ttestP > cutoff){ onoff["paired_ttestP"] = 0 }
    else{ onoff["paired_ttestP"] = 1 }
    if (mannwhitneyP > cutoff){ onoff["mannwhitneyP"] = 0 }
    else{ onoff["mannwhitneyP"] = 1 }
    if (wilcoxP > cutoff){ onoff["wilcoxP"] = 0 }
    else{ onoff["wilcoxP"] = 1 }
    # if over bonferroni, call as 
    if (ztestP > bonferroni){ onoff_bonferroni["ztestP"] = 0 }
    else{ onoff_bonferroni["ztestP"] = 1 }
    if (paired_ztestP > bonferroni){ onoff_bonferroni["paired_ztestP"] = 0 }
    else{ onoff_bonferroni["paired_ztestP"] = 1 }
    if (ttestP > bonferroni){ onoff_bonferroni["ttestP"] = 0 }
    else{ onoff_bonferroni["ttestP"] = 1 }
    if (paired_ttestP > bonferroni){ onoff_bonferroni["paired_ttestP"] = 0 }
    else{ onoff_bonferroni["paired_ttestP"] = 1 }
    if (mannwhitneyP > bonferroni){ onoff_bonferroni["mannwhitneyP"] = 0 }
    else{ onoff_bonferroni["mannwhitneyP"] = 1 }
    if (wilcoxP > bonferroni){ onoff_bonferroni["wilcoxP"] = 0 }
    else{ onoff_bonferroni["wilcoxP"] = 1 }
    # if over sidak, call as 
    if (ztestP > sidak){ onoff_sidak["ztestP"] = 0 }
    else{ onoff_sidak["ztestP"] = 1 }
    if (paired_ztestP > sidak){ onoff_sidak["paired_ztestP"] = 0 }
    else{ onoff_sidak["paired_ztestP"] = 1 }
    if (ttestP > sidak){ onoff_sidak["ttestP"] = 0 }
    else{ onoff_sidak["ttestP"] = 1 }
    if (paired_ttestP > sidak){ onoff_sidak["paired_ttestP"] = 0 }
    else{ onoff_sidak["paired_ttestP"] = 1 }
    if (mannwhitneyP > sidak){ onoff_sidak["mannwhitneyP"] = 0 }
    else{ onoff_sidak["mannwhitneyP"] = 1 }
    if (wilcoxP > sidak){ onoff_sidak["wilcoxP"] = 0 }
    else{ onoff_sidak["wilcoxP"] = 1 }
    passfailall["ten"] <- list(onoff)
    passfailbonferroni["ten"] <- list(onoff_bonferroni)
    passfailsidak["ten"] <- list(onoff_sidak)
    
    #######
    # 0.5 #
    #boolean expression table
    onoff <- NULL
    cutoff = 0.05
    bonferroni <- cutoff/size_all
    onoff_bonferroni <- NULL
    sidak <- 1-(1-cutoff)^(1/size_all)
    onoff_sidak <- NULL
    
    # if over cutoff, call as 
    if (ztestP > cutoff){ 
      error2_z = error2_z + 1
      onoff["ztestP"] = 0
    }
    else{ onoff["ztestP"] = 1 }
    if (paired_ztestP > cutoff){ 
      error2_paired_z = error2_paired_z + 1
      onoff["paired_ztestP"] = 0
    }
    else{ onoff["paired_ztestP"] = 1 }
    if (ttestP > cutoff){ 
      error2_t = error2_t + 1
      onoff["ttestP"] = 0 
    }
    else{ onoff["ttestP"] = 1 }
    if (paired_ttestP > cutoff){ 
      error2_paired_t = error2_paired_t + 1
      onoff["paired_ttestP"] = 0 
    }
    else{ onoff["paired_ttestP"] = 1 }
    if (mannwhitneyP > cutoff){ 
      error2_mw = error2_mw + 1
      onoff["mannwhitneyP"] = 0
    }
    else{ onoff["mannwhitneyP"] = 1; }
    if (wilcoxP > cutoff){ 
      error2_w = error2_w + 1
      onoff["wilcoxP"] = 0 
    }
    else{ onoff["wilcoxP"] = 1 }   
    
    # if over bonferroni, call as 
    if (ztestP > bonferroni){ 
      error2_z = error2_z + 1
      onoff_bonferroni["ztestP"] = 0
    }
    else{ onoff_bonferroni["ztestP"] = 1 }
    if (paired_ztestP > bonferroni){ 
      error2_paired_z = error2_paired_z + 1
      onoff_bonferroni["paired_ztestP"] = 0
    }
    else{ onoff_bonferroni["paired_ztestP"] = 1 }
    if (ttestP > bonferroni){ 
      error2_t = error2_t + 1
      onoff_bonferroni["ttestP"] = 0 
    }
    else{ onoff_bonferroni["ttestP"] = 1 }
    if (paired_ttestP > bonferroni){ 
      error2_paired_t = error2_paired_t + 1
      onoff_bonferroni["paired_ttestP"] = 0 
    }
    else{ onoff_bonferroni["paired_ttestP"] = 1 }
    if (mannwhitneyP > bonferroni){ 
      error2_mw = error2_mw + 1
      onoff_bonferroni["mannwhitneyP"] = 0
    }
    else{ onoff_bonferroni["mannwhitneyP"] = 1; }
    if (wilcoxP > bonferroni){ 
      error2_w = error2_w + 1
      onoff_bonferroni["wilcoxP"] = 0 
    }
    else{ onoff_bonferroni["wilcoxP"] = 1 }
    # if over sidak, call as 
    if (ztestP > sidak){ 
      error2_z = error2_z + 1
      onoff_sidak["ztestP"] = 0
    }
    else{ onoff_sidak["ztestP"] = 1 }
    if (paired_ztestP > sidak){ 
      error2_paired_z = error2_paired_z + 1
      onoff_sidak["paired_ztestP"] = 0
    }
    else{ onoff_sidak["paired_ztestP"] = 1 }
    if (ttestP > sidak){ 
      error2_t = error2_t + 1
      onoff_sidak["ttestP"] = 0 
    }
    else{ onoff_sidak["ttestP"] = 1 }
    if (paired_ttestP > sidak){ 
      error2_paired_t = error2_paired_t + 1
      onoff_sidak["paired_ttestP"] = 0 
    }
    else{ onoff_sidak["paired_ttestP"] = 1 }
    if (mannwhitneyP > sidak){ 
      error2_mw = error2_mw + 1
      onoff_sidak["mannwhitneyP"] = 0
    }
    else{ onoff_sidak["mannwhitneyP"] = 1; }
    if (wilcoxP > sidak){ 
      error2_w = error2_w + 1
      onoff_sidak["wilcoxP"] = 0 
    }
    else{ onoff_sidak["wilcoxP"] = 1 }
    passfailall["five"] <- list(onoff)
    passfailbonferroni["five"] <- list(onoff_bonferroni)
    passfailsidak["five"] <- list(onoff_sidak)
    
    write(toString(passfailall["one"]), file="A_cutoff_one.tsv", append="TRUE")
    write(toString(passfailbonferroni["one"]), file="D_bonferroni_one.tsv", append="TRUE")
    write(toString(passfailsidak["one"]), file="D_sidak_one.tsv", append="TRUE")
    write(toString(passfailall["five"]), file="A_cutoff_five.tsv", append="TRUE")
    write(toString(passfailbonferroni["five"]), file="D_bonferroni_five.tsv", append="TRUE")
    write(toString(passfailsidak["five"]), file="D_sidak_five.tsv", append="TRUE")
    write(toString(passfailall["ten"]), file="A_cutoff_ten.tsv", append="TRUE")
    write(toString(passfailbonferroni["ten"]), file="D_bonferroni_ten.tsv", append="TRUE")
    write(toString(passfailsidak["ten"]), file="D_sidak_ten.tsv", append="TRUE")
    
    #if (correct[r] == 1){
    #  if passfailsidak["ten"])
    #}
  }
}

# # UNIX METHODS TO PARSE OUTPUT OF ALL VALUES FOR DELIVERABLE :
# cat A_cutoff_one.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $1}' | grep -c '1'
# cat A_cutoff_one.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $2}' | grep -c '1'
# cat A_cutoff_one.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $3}' | grep -c '1'
# cat A_cutoff_one.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $4}' | grep -c '1'
# cat A_cutoff_five.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $1}' | grep -c '1'
# cat A_cutoff_five.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $2}' | grep -c '1'
# cat A_cutoff_five.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $3}' | grep -c '1'
# cat A_cutoff_five.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $4}' | grep -c '1'
# cat A_cutoff_ten.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $1}' | grep -c '1'
# cat A_cutoff_ten.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $2}' | grep -c '1'
# cat A_cutoff_ten.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $3}' | grep -c '1'
# cat A_cutoff_ten.tsv | tr '(' ' ' | tr ',' ' ' | tr 'c' ' ' | tr ')' ' ' | awk '{print $4}' | grep -c '1'

##########
# SOLN E #
##########


file.remove("E_error_rates.tsv")

#################### STD CUTOFF
write(toString("STANDARD CUTOFF"), file="E_error_rates.tsv", append="TRUE")
write(toString("0.05"), file="E_error_rates.tsv", append="TRUE")
# t-test at significance = a = 0.03* is "correct"
#onoff_true3 <- NULL
#onoff_true5 <- NULL
#counter3_a = 0
#counter5_a = 0
#counter3_b = 0
#counter5_b = 0
correct_on = 0
hat <- NULL
ztest_falsepositive = 0
paired_ztest_falsepositive = 0
ttest_falsepositive = 0
paired_ttest_falsepositive = 0
mw_falsepositive = 0
wilcox_falsepositive = 0

ztest_falsenegative = 0
paired_ztest_falsenegative = 0
ttest_falsenegative = 0
paired_ttest_falsenegative = 0
mw_falsenegative = 0
wilcox_falsenegative = 0

ztest_type2 = 0
paired_ztest_type2 = 0
ttest_type2 = 0
paired_ttest_type2 = 0
mw_type2 = 0
wilcox_type2 = 0

for(r in 1:length(newall)){
  if(all[r][1] != "NULL"){   
    spotID = all[[r]][1]
    correct_bool = correct[r]
    ztestP = newall[[r]][8]
    paired_ztestP = all[[r]][9]
    ttestP = all[[r]][10]
    paired_ttestP = all[[r]][11]
    mannwhitney = all[[r]][12]
    wilcoxon = all[[r]][13]
    # a
    if (correct_bool == 0){
      correct_on = correct_on + 1
      if (ztestP <= 0.05){ ztest_falsepositive = ztest_falsepositive + 1 }
      if (paired_ztestP <= 0.05){ paired_ztest_falsepositive = paired_ztest_falsepositive + 1 }
      if (ttestP <= 0.05){ ttest_falsepositive = ttest_falsepositive + 1 }
      if (paired_ttestP <= 0.05){ paired_ttest_falsepositive = paired_ttest_falsepositive + 1 }
      if (mannwhitney <= 0.05){ mw_falsepositive = mw_falsepositive + 1 }
      if (wilcoxon <= 0.05){ wilcox_falsepositive = wilcox_falsepositive + 1 }
    }
    else{
      if (ztestP > 0.05){ ztest_falsenegative = ztest_falsenegative + 1 }
      if (paired_ztestP > 0.05){ paired_ztest_falsenegative = paired_ztest_falsenegative + 1 }
      if (ttestP > 0.05){ ttest_falsenegative = ttest_falsenegative + 1 }
      if (paired_ttestP > 0.05){ paired_ttest_falsenegative = paired_ttest_falsenegative + 1 }
      if (mannwhitney > 0.05){ mw_falsepositive = mw_falsenegative + 1 }
      if (wilcoxon > 0.05){ wilcox_falsenegative = wilcox_falsenegative + 1 }
      
      if (ztestP <= 0.05){ ztest_type2 = ztest_type2 + 1 }
      if (paired_ztestP <= 0.05){ paired_ztest_type2 = paired_ztest_type2 + 1 }
      if (ttestP <= 0.05){ ttest_falsenegative = ttest_falsenegative + 1 }
      if (paired_ttestP <= 0.05){ paired_ttest_type2 = paired_ttest_type2 + 1 }
      if (mannwhitney <= 0.05){ mw_type2 = mw_type2 + 1 }
      if (wilcoxon <= 0.05){ wilcox_type2 = wilcox_type2 + 1 }
    }
  }
}
# 
# FDR = counter5_a/size_all
# TYPE 1 error rate at 0.05 based on assumption *
# a = p(reject null given null is true)
# alpha_error1_t = abs(counter5_a-counter3_a)/size_all
# TYPE 2 error rate at 0.05 based on assumption *
# a = p(accept null given null is false)
# beta_error2_t = abs(counter5_b-counter3_b)/size_all

FDR_ztest = ztest_falsepositive/length(newall)
FDR_paired_ztest = paired_ztest_falsepositive/length(newall)
FDR_ttest = ttest_falsepositive/length(newall)
FDR_paired_ttest = paired_ttest_falsepositive/length(newall)
FDR_mw = mw_falsepositive/length(newall)
FDR_wilcox = wilcox_falsepositive/length(newall)

type1_ztest = ztest_falsenegative/length(newall)
type1_paired_ztest = paired_ztest_falsenegative/length(newall)
type1_ttest = ttest_falsenegative/length(newall)
type1_paired_ttest = paired_ttest_falsenegative/length(newall)
type1_mw = mw_falsenegative/length(newall)
type1_wilcox = wilcox_falsenegative/length(newall)

type2_ztest = ztest_type2/length(newall)
type2_paired_ztest = paired_ztest_type2/length(newall)
type2_ttest = ttest_type2/length(newall)
type2_paired_ttest = paired_ttest_type2/length(newall)
type2_mw = mw_type2/length(newall)
type2_wilcox = wilcox_type2/length(newall)
      
write(toString("FDR"), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_paired_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_paired_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_mw), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_wilcox), file="E_error_rates.tsv", append="TRUE")

write(toString("type1"), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_paired_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_paired_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_mw), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_wilcox), file="E_error_rates.tsv", append="TRUE")

write(toString("type2"), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_paired_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_paired_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_mw), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_wilcox), file="E_error_rates.tsv", append="TRUE")




#################### BONFERRONI

write(toString("BONFERRONI"), file="E_error_rates.tsv", append="TRUE")
bonferroni <- 0.05/length(newall)
write(toString(bonferroni), file="E_error_rates.tsv", append="TRUE")
# t-test at significance = a = 0.03* is "correct"
#onoff_true3 <- NULL
#onoff_true5 <- NULL
#counter3_a = 0
#counter5_a = 0
#counter3_b = 0
#counter5_b = 0
correct_on = 0
hat <- NULL
ztest_falsepositive = 0
paired_ztest_falsepositive = 0
ttest_falsepositive = 0
paired_ttest_falsepositive = 0
mw_falsepositive = 0
wilcox_falsepositive = 0

ztest_falsenegative = 0
paired_ztest_falsenegative = 0
ttest_falsenegative = 0
paired_ttest_falsenegative = 0
mw_falsenegative = 0
wilcox_falsenegative = 0

ztest_type2 = 0
paired_ztest_type2 = 0
ttest_type2 = 0
paired_ttest_type2 = 0
mw_type2 = 0
wilcox_type2 = 0

for(r in 1:length(newall)){
  if(all[r][1] != "NULL"){   
    spotID = all[[r]][1]
    correct_bool = correct[r]
    ztestP = newall[[r]][8]
    paired_ztestP = all[[r]][9]
    ttestP = all[[r]][10]
    paired_ttestP = all[[r]][11]
    mannwhitney = all[[r]][12]
    wilcoxon = all[[r]][13]
    # a
    if (correct_bool == 0){
      correct_on = correct_on + 1
      if (ztestP <= bonferroni){ ztest_falsepositive = ztest_falsepositive + 1 }
      if (paired_ztestP <= bonferroni){ paired_ztest_falsepositive = paired_ztest_falsepositive + 1 }
      if (ttestP <= bonferroni){ ttest_falsepositive = ttest_falsepositive + 1 }
      if (paired_ttestP <= bonferroni){ paired_ttest_falsepositive = paired_ttest_falsepositive + 1 }
      if (mannwhitney <= bonferroni){ mw_falsepositive = mw_falsepositive + 1 }
      if (wilcoxon <= bonferroni){ wilcox_falsepositive = wilcox_falsepositive + 1 }
    }
    else{
      if (ztestP > bonferroni){ ztest_falsenegative = ztest_falsenegative + 1 }
      if (paired_ztestP > bonferroni){ paired_ztest_falsenegative = paired_ztest_falsenegative + 1 }
      if (ttestP > bonferroni){ ttest_falsenegative = ttest_falsenegative + 1 }
      if (paired_ttestP > bonferroni){ paired_ttest_falsenegative = paired_ttest_falsenegative + 1 }
      if (mannwhitney > bonferroni){ mw_falsepositive = mw_falsenegative + 1 }
      if (wilcoxon > bonferroni){ wilcox_falsenegative = wilcox_falsenegative + 1 }
      
      if (ztestP <= bonferroni){ ztest_type2 = ztest_type2 + 1 }
      if (paired_ztestP <= bonferroni){ paired_ztest_type2 = paired_ztest_type2 + 1 }
      if (ttestP <= bonferroni){ ttest_falsenegative = ttest_falsenegative + 1 }
      if (paired_ttestP <= bonferroni){ paired_ttest_type2 = paired_ttest_type2 + 1 }
      if (mannwhitney <= bonferroni){ mw_type2 = mw_type2 + 1 }
      if (wilcoxon <= bonferroni){ wilcox_type2 = wilcox_type2 + 1 }
    }
  }
}
# 
# FDR = counter5_a/size_all
# TYPE 1 error rate at 0.05 based on assumption *
# a = p(reject null given null is true)
# alpha_error1_t = abs(counter5_a-counter3_a)/size_all
# TYPE 2 error rate at 0.05 based on assumption *
# a = p(accept null given null is false)
# beta_error2_t = abs(counter5_b-counter3_b)/size_all

FDR_ztest = ztest_falsepositive/length(newall)
FDR_paired_ztest = paired_ztest_falsepositive/length(newall)
FDR_ttest = ttest_falsepositive/length(newall)
FDR_paired_ttest = paired_ttest_falsepositive/length(newall)
FDR_mw = mw_falsepositive/length(newall)
FDR_wilcox = wilcox_falsepositive/length(newall)

type1_ztest = ztest_falsenegative/length(newall)
type1_paired_ztest = paired_ztest_falsenegative/length(newall)
type1_ttest = ttest_falsenegative/length(newall)
type1_paired_ttest = paired_ttest_falsenegative/length(newall)
type1_mw = mw_falsenegative/length(newall)
type1_wilcox = wilcox_falsenegative/length(newall)

type2_ztest = ztest_type2/length(newall)
type2_paired_ztest = paired_ztest_type2/length(newall)
type2_ttest = ttest_type2/length(newall)
type2_paired_ttest = paired_ttest_type2/length(newall)
type2_mw = mw_type2/length(newall)
type2_wilcox = wilcox_type2/length(newall)

write(toString("FDR"), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_paired_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_paired_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_mw), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_wilcox), file="E_error_rates.tsv", append="TRUE")

write(toString("type1"), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_paired_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_paired_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_mw), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_wilcox), file="E_error_rates.tsv", append="TRUE")

write(toString("type2"), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_paired_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_paired_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_mw), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_wilcox), file="E_error_rates.tsv", append="TRUE")





#################### SIDAK


write(toString("SIDAK"), file="E_error_rates.tsv", append="TRUE")
sidak <- 1-(1-(0.05))^(1/length(newall))
write(toString(sidak), file="E_error_rates.tsv", append="TRUE")
# t-test at significance = a = 0.03* is "correct"
#onoff_true3 <- NULL
#onoff_true5 <- NULL
#counter3_a = 0
#counter5_a = 0
#counter3_b = 0
#counter5_b = 0
correct_on = 0
hat <- NULL
ztest_falsepositive = 0
paired_ztest_falsepositive = 0
ttest_falsepositive = 0
paired_ttest_falsepositive = 0
mw_falsepositive = 0
wilcox_falsepositive = 0

ztest_falsenegative = 0
paired_ztest_falsenegative = 0
ttest_falsenegative = 0
paired_ttest_falsenegative = 0
mw_falsenegative = 0
wilcox_falsenegative = 0

ztest_type2 = 0
paired_ztest_type2 = 0
ttest_type2 = 0
paired_ttest_type2 = 0
mw_type2 = 0
wilcox_type2 = 0

for(r in 1:length(newall)){
  if(all[r][1] != "NULL"){   
    spotID = all[[r]][1]
    correct_bool = correct[r]
    ztestP = newall[[r]][8]
    paired_ztestP = all[[r]][9]
    ttestP = all[[r]][10]
    paired_ttestP = all[[r]][11]
    mannwhitney = all[[r]][12]
    wilcoxon = all[[r]][13]
    # a
    if (correct_bool == 0){
      correct_on = correct_on + 1
      if (ztestP <= sidak){ ztest_falsepositive = ztest_falsepositive + 1 }
      if (paired_ztestP <= sidak){ paired_ztest_falsepositive = paired_ztest_falsepositive + 1 }
      if (ttestP <= sidak){ ttest_falsepositive = ttest_falsepositive + 1 }
      if (paired_ttestP <= sidak){ paired_ttest_falsepositive = paired_ttest_falsepositive + 1 }
      if (mannwhitney <= sidak){ mw_falsepositive = mw_falsepositive + 1 }
      if (wilcoxon <= sidak){ wilcox_falsepositive = wilcox_falsepositive + 1 }
    }
    else{
      if (ztestP > sidak){ ztest_falsenegative = ztest_falsenegative + 1 }
      if (paired_ztestP > sidak){ paired_ztest_falsenegative = paired_ztest_falsenegative + 1 }
      if (ttestP > sidak){ ttest_falsenegative = ttest_falsenegative + 1 }
      if (paired_ttestP > sidak){ paired_ttest_falsenegative = paired_ttest_falsenegative + 1 }
      if (mannwhitney > sidak){ mw_falsepositive = mw_falsenegative + 1 }
      if (wilcoxon > sidak){ wilcox_falsenegative = wilcox_falsenegative + 1 }
      
      if (ztestP <= sidak){ ztest_type2 = ztest_type2 + 1 }
      if (paired_ztestP <= sidak){ paired_ztest_type2 = paired_ztest_type2 + 1 }
      if (ttestP <= sidak){ ttest_falsenegative = ttest_falsenegative + 1 }
      if (paired_ttestP <= sidak){ paired_ttest_type2 = paired_ttest_type2 + 1 }
      if (mannwhitney <= sidak){ mw_type2 = mw_type2 + 1 }
      if (wilcoxon <= sidak){ wilcox_type2 = wilcox_type2 + 1 }
    }
  }
}
# 
# FDR = counter5_a/size_all
# TYPE 1 error rate at 0.05 based on assumption *
# a = p(reject null given null is true)
# alpha_error1_t = abs(counter5_a-counter3_a)/size_all
# TYPE 2 error rate at 0.05 based on assumption *
# a = p(accept null given null is false)
# beta_error2_t = abs(counter5_b-counter3_b)/size_all

FDR_ztest = ztest_falsepositive/length(newall)
FDR_paired_ztest = paired_ztest_falsepositive/length(newall)
FDR_ttest = ttest_falsepositive/length(newall)
FDR_paired_ttest = paired_ttest_falsepositive/length(newall)
FDR_mw = mw_falsepositive/length(newall)
FDR_wilcox = wilcox_falsepositive/length(newall)

type1_ztest = ztest_falsenegative/length(newall)
type1_paired_ztest = paired_ztest_falsenegative/length(newall)
type1_ttest = ttest_falsenegative/length(newall)
type1_paired_ttest = paired_ttest_falsenegative/length(newall)
type1_mw = mw_falsenegative/length(newall)
type1_wilcox = wilcox_falsenegative/length(newall)

type2_ztest = ztest_type2/length(newall)
type2_paired_ztest = paired_ztest_type2/length(newall)
type2_ttest = ttest_type2/length(newall)
type2_paired_ttest = paired_ttest_type2/length(newall)
type2_mw = mw_type2/length(newall)
type2_wilcox = wilcox_type2/length(newall)

write(toString("FDR"), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_paired_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_paired_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_mw), file="E_error_rates.tsv", append="TRUE")
write(toString(FDR_wilcox), file="E_error_rates.tsv", append="TRUE")

write(toString("type1"), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_paired_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_paired_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_mw), file="E_error_rates.tsv", append="TRUE")
write(toString(type1_wilcox), file="E_error_rates.tsv", append="TRUE")

write(toString("type2"), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_paired_ztest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_paired_ttest), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_mw), file="E_error_rates.tsv", append="TRUE")
write(toString(type2_wilcox), file="E_error_rates.tsv", append="TRUE")



# THE END.