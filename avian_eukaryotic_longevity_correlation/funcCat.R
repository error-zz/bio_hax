#funcCat.R

#############################################################################################
# EXAMPLE INPUT PREPARATION METHODS ! 

# cat ko00001.keg | grep -v '#' | grep -v '!' | tr '<' ' ' | tr '>' ' ' | grep -v '+D' | tr ']' ' ' | tr '[' '\t' > ko_sum.keg
# cat ko_sum.keg | tr "\'" "_" | tr "\"" "_" >ko_sum_clean.keg 

# grep KEGG denovo_protein_annotation.tab > KEGG.denovo_protein_annotation.tab

# qsub -P 810529 -V -e GRID_err.txt -m ea -M jmccorri@jcvi.org -N funcCat_Krona -cwd perl funcCat_Krona.pl 
# (it takes a while... output written to out.txt)
# cp out.txt per_transcript_mean_cvg.txt

#############################################################################################

setwd ("/Users/apple/Desktop/funcCat")
keg <- read.table( "ko_sum_clean.keg", sep="\t", header=FALSE, fill = TRUE, row.names = NULL )
dpa <- read.table( "KEGG.denovo_protein_annotation.tab", sep="\t", header=FALSE, fill=TRUE, quote="" )
pro <- read.table( "per_transcript_mean_cvg.txt", sep="\t", header=FALSE )
meancvg <- unlist(pro$V2)
names(meancvg) <- pro$V1
 
dpa.brief <- unique(dpa[,c(1,15)])
dpa.lol           <- matrix(nrow=length(unique(sort(dpa.brief[,1]))),ncol=999)
rownames(dpa.lol) <- unique(sort(dpa.brief[,1]))
for(k in 1:dim(dpa.brief)[1]){
  if(is.na(dpa.lol[ dpa[k,1] , 1 ])){
    strs <- strsplit(as.character(dpa.brief[k,2]),'\\|')[[1]][ grepl('KEGG', strsplit(as.character(dpa.brief[k,2]),'\\|')[[1]]) ]
    dpa.lol[ dpa[k,1] , 1:(length(strs)) ] <- strs
    print(paste('A:',dpa[k,1],':',strs))
  }else{
    strs <- strsplit(as.character(dpa.brief[k,2]),'\\|')[[1]][ grepl('KEGG', strsplit(as.character(dpa.brief[k,2]),'\\|')[[1]]) ]
    print(paste('B:',dpa[k,1],':',strs))
    if(length(strs) > 1){
      left  <- 999-sum(is.na(dpa.lol[ dpa[k,1],] ))+1
      right <- 999-sum(is.na(dpa.lol[ dpa[k,1],] ))+length(strs)
      dpa.lol[ dpa[k,1] , left:right ] <- strs
    }else{
      left  <- 999-sum(is.na(dpa.lol[ dpa[k,1],] ))+1
      dpa.lol[ dpa[k,1] , left ] <- strs
    }
    
  }
}
dpa.pop <- dpa.lol[,1:43]
  # 00531+3.2.1.50
  # ^ 
sum_matrix <- as.data.frame(matrix(ncol=6,nrow=dim(keg)[1]))
colnames(sum_matrix) <- c("Count","A","B","C","D","CD_ID")
for(k in 1:dim(keg)[1]){
#for(k in 1:4){
 s <- substring(keg[k,1], 1, 1)
 if( s == "A" ){
   cur_A <- strsplit(as.character(keg[k,1]),' ')[[1]][3]
   print(paste("A || ",cur_A))
   sum_matrix[k,] <- c(0,cur_A,"","","","")
 }else if( s == "B" ){
   if(length(strsplit(as.character(keg[k,1]),' ')[[1]]) != 1){
     cur_B <- as.character(paste( unlist(strsplit(as.character(keg[k,1]),' ')[[1]][5:(length(strsplit(as.character(keg[k,1]),' ')[[1]])-1)]), collapse=' ' ))
     print(paste("B || ",cur_A,"a->",cur_B))
     sum_matrix[k,] <- c(0,cur_A,cur_B,"","","")
   }
 }else if( s == "C" ){
   #cur_C       <-  unlist(strsplit(as.character(keg[k,1]),'    ')[[1]])[2]
   cur_C        <- substr(trimws(toupper(strsplit(as.character(keg[k,2]),':')[[1]][2])),3,nchar(trimws(toupper(strsplit(as.character(keg[k,2]),':')[[1]][2]))))
   cur_C_string <- trimws(unlist(strsplit(as.character(keg[k,1]),'    ')[[1]])[2])
   print(paste("C || ",cur_A,"a->",cur_B,"b->",cur_C,":",cur_C_string))
   sum_matrix[k,] <- c(0,cur_A,cur_B,cur_C_string,"",cur_C)
 }else if( s == "D" ){
   cur_D <- strsplit(trimws(as.character(keg[k,2])),':')[[1]][2]
   cur_D_string <- trimws(unlist(strsplit(as.character(keg[k,1]),'    ')[[1]][2]))
   print(paste("D || ",cur_A,"a->",cur_B,"b->",cur_C,"c->",cur_D,":",cur_D_string))
   sum_matrix[k,] <- c(0,cur_A,cur_B,cur_C,cur_D_string,cur_D)
 }
}
sum.matrix.pop <- sum_matrix[!is.na(sum_matrix$CD_ID),]
# 13308 of 29417

#PREPARE TO POP COUNTS BY SPLITTING THE dpa.pop (transcript ID -> KEGG ID matrix) INTO LEVEL C AND LEVEL D IDs
dpa.pop.c <- as.data.frame(matrix(nrow=dim(dpa.pop)[1],ncol=dim(dpa.pop)[2]))
dpa.pop.d <- as.data.frame(matrix(nrow=dim(dpa.pop)[1],ncol=dim(dpa.pop)[2]))
for(x in 1:dim(dpa.pop)[1]){
  for(y in 1:dim(dpa.pop)[2]){
    dpq <- dpa.pop[x,y] # "KEGG: 00970+6.1.1.16"
    if(!is.na(dpq)){
      dpq_split <- strsplit(strsplit(dpq," ")[[1]],"\\+")[[2]]
      dpa.pop.d[x,y] <- as.character(dpq_split[2])
      dpa.pop.c[x,y] <- as.character(dpq_split[1])
    }
  }
}
rownames(dpa.pop.d) <- rownames(dpa.pop)
rownames(dpa.pop.c) <- rownames(dpa.pop)
# and invert so its easier to query
uniq.dpa.pop.d.id <- unique(sort(unname(unlist(dpa.pop.d)[!is.na(unlist(dpa.pop.d))])))
uniq.dpa.pop.c.id <- unique(sort(unname(unlist(dpa.pop.c)[!is.na(unlist(dpa.pop.c))])))
#inverse.dpa.d  <- as.data.frame(matrix(nrow))
#for(d.id in uniq.dpa.pop.d.id){
#}
 

# POP COUNTS
count_matrix <- sum_matrix
#count_matrix <- count_matrix[!is.na(count_matrix$Count),]
count_matrix$Count <- as.numeric( count_matrix$Count )
unclassified_count <- 0
no_classification_match_count <- 0
count_d <- 0
count_d0 <- 0
count_c <- 0
count_n <- 0
count_u <- 0
for(g in dim(dpa.pop)[1]){
  if(sum(!is.na(dpa.pop[g,])) == 0){
    #unclassified_count += #cvg
    unclassified_count <- unclassified_count + meancvg[g]
    count_u <- count_u + 1
  }else{
    for(gs in strsplit(substr(dpa.pop[g,][!is.na(dpa.pop[g,])],7,21),"\\+")){
      #gs[1] = 5 digit kegg id
      #gs[2] = 4 entry EC id
      #if(gs[2] %in% sum_matrix$CD_ID){
      if(gs[2] %in% uniq.dpa.pop.d.id){
        # add to count for D, and also (A,B,C)
        w <- sum_matrix$CD_ID==gs[2]
        w[is.na(w)] <- FALSE
        if(sum(w)>0){
          
          # d and c level
          count_matrix[w,]$Count <- count_matrix[w,]$Count + meancvg[g]
          
          # b level
          bq <- count_matrix$B == unique(count_matrix[w,]$B)
          bq[is.na(bq)] <- FALSE
          count_matrix[rownames(count_matrix[bq,][count_matrix[bq,]$C == "",]),]$Count <- count_matrix[rownames(count_matrix[bq,][count_matrix[bq,]$C == "",]),]$Count + meancvg[g]
          
          # a level
          aq <- count_matrix$A == unique(count_matrix[w,]$A)
          aq[is.na(aq)] <- FALSE
          count_matrix[rownames(count_matrix[aq,][count_matrix[aq,]$B == "",]),]$Count <- count_matrix[rownames(count_matrix[aq,][count_matrix[aq,]$B == "",]),]$Count + meancvg[g]
          
          count_d <- count_d + 1
        }else{
          count_d0 <- count_d0 + 1
        }
      #}else if(gs[1] %in% sum_matrix$CD_ID){
      }else if(gs[1] %in% uniq.dpa.pop.c.id){
        # add to count for C, and also (A,B)
        
        print("C")
        count_c <- count_c + 1
      }else{
        # nothing
        
        no_classification_match_count <- no_classification_match_count + meancvg[g]
        count_n <- count_n + 1
        
      }
    }
  }
} 
  