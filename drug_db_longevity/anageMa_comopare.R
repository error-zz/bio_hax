
setwd("~/Desktop/orthoDB/")

ma.resmls      <- read.table("Ma.elife-19130-table1-data1-v2.D.resMLS.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
ma.mls         <- read.table("Ma.elife-19130-table1-data1-v2.B.MLS.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
anage          <- read.table("drugage.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
genage          <- read.table("genage.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
diopt          <- read.table("dioptscore_hm.tsv", sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)

rownames(diopt)[ rownames(diopt) %in% toupper(ma.resmls$Symbol) ]
# [1] "MMP2"
ma.resmls[ toupper(ma.resmls$Symbol) == "MMP2" , ]$p.value.all
# [1] 0.000473


rownames(diopt)[ rownames(diopt) %in% toupper(genage$symbol) ]
# [1] "APOE" "IGF2