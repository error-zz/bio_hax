
setwd ("~/Desktop/exe.DrugLongevity")

drugList <- read.table("drugList", sep="\t", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)
geneList <- read.table("geneList", sep="\t", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)
snpList  <- read.table("snpList",  sep="\t", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)

#for(snp in snpList){
  # check for unique query indices
  # system("wget https://ldlink.nci.nih.gov/?var=rs159428&pop=CEU%2BTSI%2BFIN%2BGBR%2BIBS&r2_d=r2&tab=ldproxy")
  # https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=rs3&pop=MXL&r2_d=r2&token=faketoken123
  # "error": "Invalid or expired API token. Please register using the API Access tab: https://ldlink.nci.nih.gov/?tab=apiaccess"
  # QUERY BY HAND!
  # AND IMPORT
  # https://ldlink.nci.nih.gov/?var=XXXXXXX&pop=CEU%2BTSI%2BFIN%2BGBR%2BIBS&r2_d=r2&tab=ldproxy
  # XXXX = rsID
  # save on Download all proxy variants link
#}
snpz <- c("snp01.rs1230666.txt",
          "snp02.rs1275922.txt",
          "snp03.rs61348208.txt",
          "snp04.rs34967069.txt",
          "snp05.rs10455872.txt",
          "snp06.rs1556516.txt",
          "snp07.rs11065979.txt",
          "snp08.rs8042849.txt",
          "snp09.rs6224.txt",
          "snp10.rs12924886.txt",
          "snp11.rs142158911.txt",
          "snp12.rs429358.txt",
          "snp13.rs4970836.txt",
          "snp14.rs6744653.txt",
          "snp15.rs10211471.txt",
          "snp16.rs111333005.txt",
          "snp17.rs113160991.txt",
          "snp18.rs56179563.txt",
          "snp19.rs2519093.txt",
          "snp20.rs429358.txt",
          "snp21.rs11065979.txt",
          "snp22.rs1230666.txt",
          "snp23.rs3130507.txt",
          "snp24.rs6224.txt",
          "snp25.rs8042849.txt",
          "snp26.rs10197246.txt",
          "snp27.rs12203592.txt",
          "snp28.rs1049053.txt",
          "snp29.rs10455872.txt",
          "snp30.rs140570886.txt",
          "snp31.rs7859727.txt",
          "snp32.rs34872471.txt",
          "snp33.rs2860197.txt",
          "snp34.rs1126809.txt",
          "snp35.rs4784227.txt",
          "snp36.rs4268748.txt",
          "snp37.rs159428.txt")
snp <- as.data.frame(matrix(nrow=1, ncol = 11))
si <- 1
for( s in snpz ){
  tmp.snp <- read.table(s, sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
  if( sum(tmp.snp$R2 > 0.8)>0 ){
    nuRows <- tmp.snp[ tmp.snp$R2 > 0.8 , ]
    for(x in 1:dim(nuRows)[1]){
      snp[si,] <- c(  s , unname(unlist(nuRows[ x , ])) )
      si <- si + 1
    }
  }
}
colnames(snp) <- c( "SNPid", colnames(tmp.snp))

write.table( sort(table(snp[,1])) , "sorted.LDproxy.histogram.tsv" , sep="\t" )

write.table( snp , "snpSignifSum.tsv" , sep="\t" )
 
rsids <- list()
ri <- 1
for(s in strsplit( snp$SNPid , "\\." )){
  rsids[ri] <- s[2]
  ri <- ri + 1
}
rsids <- unlist(rsids)
snp$rsids <- rsids
 



# 

# linda <- read.table("Linda.QTL.Main_Catalog.20190109.WoNa.InternalCheckOk.Fdr0.05.eQTL_Catalog.tsv", sep="\t", stringsAsFactors = FALSE, header=FALSE, fill=TRUE)
# > sum(linda$V3 == "9:22100176")
# [1] 0
# > sum(linda$V3 == "rs1556516")
# [1] 0
# WEB QUERY GETS A RESULT:
# 9:22100176	--->	9:22084310

linda1 <- read.table("LinDA.snpsInLD.01.csv", sep=",", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
linda2 <- read.table("LinDA.snpsInLD.02.csv", sep=",", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
linda3 <- read.table("LinDA.snpsInLD.03.csv", sep=",", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
linda4 <- read.table("LinDA.snpsInLD.04.csv", sep=",", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
linda5 <- read.table("LinDA.snpsInLD.05.csv", sep=",", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
linda6 <- read.table("LinDA.snpsInLD.06.csv", sep=",", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
 







