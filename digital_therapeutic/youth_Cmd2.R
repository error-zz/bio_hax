load("YOUTH.incomplete_x57343.Rdata")

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
clustinput <-  #read.table( args[1], sep="\t", header=TRUE , fill=TRUE )
startx <- args[1]
endx <- args[2]
prot2clust <-    #read.table( args[2], sep="\t", header=FALSE, fill=TRUE )
#startx <- 57343
#endx <- 60000
problem_list <- list()
pi <- 1
for(z in startx:endx){
  
  if( z %% 100 == 0){ print(paste( "populating table :", z, round( (z/dim(allinput)[1])*100, digits=4 ) )) }
  
  #if(z == 57343){
  #  out_tab[z,1]  <- x_pre <- mean(pcoa.dist$vectors[ strsplit("Anxious|Concerned|Loving|Hopeful|Concerned", "\\|")[[1]] , 1 ])
  #  out_tab[z,3]  <- y_pre <- mean(pcoa.dist$vectors[ strsplit("Anxious|Concerned|Loving|Hopeful|Concerned", "\\|")[[1]] , 2 ])
  #}else{
  #}
  sss_pre <- strsplit(allinput$Pre_emotions, "\\|")[[z]]
  sss_post <- strsplit(allinput$Post_emotions, "\\|")[[z]]
  sss <- c(sss_pre, sss_post)
  
  if( sum( sss %in% unique(keep_all) ) == length(sss) ){
    out_tab[z,1]  <- x_pre <- mean(pcoa.dist$vectors[ sss_pre , 1 ])
    out_tab[z,3]  <- y_pre <- mean(pcoa.dist$vectors[ sss_pre , 2 ])
    out_tab[z,2]  <- x_post <- mean(pcoa.dist$vectors[ sss_post , 1 ])
    out_tab[z,4]  <- y_post <- mean(pcoa.dist$vectors[ sss_post , 2 ])
    
    out_tab[z,5]  <- d_effect <- sqrt( (x_post - x_pre) ^ 2 + (y_post - y_pre) ^ 2 )
    
    out_tab[z,6]  <- d_anxious_pre    <- sqrt( ( mean(pcoa.dist$vectors[anxious,1])   - x_pre) ^ 2 + (mean(pcoa.dist$vectors[anxious,2])   - y_pre) ^ 2 )
    out_tab[z,7]  <- d_depressed_pre  <- sqrt( ( mean(pcoa.dist$vectors[depressed,1]) - x_pre) ^ 2 + (mean(pcoa.dist$vectors[depressed,2]) - y_pre) ^ 2 )
    out_tab[z,8]  <- d_happy_pre      <- sqrt( ( mean(pcoa.dist$vectors[happy[happy %in% rownames(pcoa.dist$vectors)],1])     - x_pre) ^ 2 + (mean(pcoa.dist$vectors[happy[happy %in% rownames(pcoa.dist$vectors)],2])     - y_pre) ^ 2 )
    out_tab[z,9]  <- d_anxious_post   <- sqrt( ( mean(pcoa.dist$vectors[anxious,1])   - x_post) ^ 2 + (mean(pcoa.dist$vectors[anxious,2])   - y_post) ^ 2 )
    out_tab[z,10] <- d_depressed_post <- sqrt( ( mean(pcoa.dist$vectors[depressed,1]) - x_post) ^ 2 + (mean(pcoa.dist$vectors[depressed,2]) - y_post) ^ 2 )
    out_tab[z,11] <- d_happy_post     <- sqrt( ( mean(pcoa.dist$vectors[happy[happy %in% rownames(pcoa.dist$vectors)],1])     - x_post) ^ 2 + (mean(pcoa.dist$vectors[happy[happy %in% rownames(pcoa.dist$vectors)],2])     - y_post) ^ 2 )
    out_tab[z,12] <- abs(d_anxious_post) - abs(d_anxious_pre)
    out_tab[z,13] <- abs(d_depressed_post) - abs(d_depressed_pre)
    out_tab[z,14] <- abs(d_happy_post) - abs(d_happy_pre)
    
    out_tab[z,15] <- d_clinical_anxiety2_pre    <- sqrt( ( mean(pcoa.dist$vectors[clinical_anxiety2,1])   - x_pre) ^ 2 + (mean(pcoa.dist$vectors[clinical_anxiety2,2])   - y_pre) ^ 2 )
    out_tab[z,16] <- d_clinical_anxiety_pre    <- sqrt( ( mean(pcoa.dist$vectors[clinical_anxiety,1])   - x_pre) ^ 2 + (mean(pcoa.dist$vectors[clinical_anxiety,2])   - y_pre) ^ 2 )
    out_tab[z,17] <- d_clinical_depression2_pre    <- sqrt( ( mean(pcoa.dist$vectors[clinical_depression2,1])   - x_pre) ^ 2 + (mean(pcoa.dist$vectors[clinical_depression2,2])   - y_pre) ^ 2 )
    out_tab[z,18] <- d_clinical_depression_pre    <- sqrt( ( mean(pcoa.dist$vectors[clinical_depression,1])   - x_pre) ^ 2 + (mean(pcoa.dist$vectors[clinical_depression,2])   - y_pre) ^ 2 )
    out_tab[z,19] <- d_clinical_anger2_pre    <- sqrt( ( mean(pcoa.dist$vectors[clinical_anger2,1])   - x_pre) ^ 2 + (mean(pcoa.dist$vectors[clinical_anger2,2])   - y_pre) ^ 2 )
    out_tab[z,20] <- d_clinical_anger_pre    <- sqrt( ( mean(pcoa.dist$vectors[clinical_anger,1])   - x_pre) ^ 2 + (mean(pcoa.dist$vectors[clinical_anger,2])   - y_pre) ^ 2 )
    out_tab[z,21] <- d_clinical_positive2_pre    <- sqrt( ( mean(pcoa.dist$vectors[clinical_positive2,1])   - x_pre) ^ 2 + (mean(pcoa.dist$vectors[clinical_positive2,2])   - y_pre) ^ 2 )
    out_tab[z,22] <- d_clinical_positive_pre    <- sqrt( ( mean(pcoa.dist$vectors[clinical_positive,1])   - x_pre) ^ 2 + (mean(pcoa.dist$vectors[clinical_positive,2])   - y_pre) ^ 2 )
    out_tab[z,23] <- d_clinical_negative_pre    <- sqrt( ( mean(pcoa.dist$vectors[clinical_negative,1])   - x_pre) ^ 2 + (mean(pcoa.dist$vectors[clinical_negative,2])   - y_pre) ^ 2 )
    
    out_tab[z,24] <- d_clinical_anxiety2_post    <- sqrt( ( mean(pcoa.dist$vectors[clinical_anxiety2,1])   - x_post) ^ 2 + (mean(pcoa.dist$vectors[clinical_anxiety2,2])   - y_post) ^ 2 )
    out_tab[z,25] <- d_clinical_anxiety_post    <- sqrt( ( mean(pcoa.dist$vectors[clinical_anxiety,1])   - x_post) ^ 2 + (mean(pcoa.dist$vectors[clinical_anxiety,2])   - y_post) ^ 2 )
    out_tab[z,26] <- d_clinical_depression2_post    <- sqrt( ( mean(pcoa.dist$vectors[clinical_depression2,1])   - x_post) ^ 2 + (mean(pcoa.dist$vectors[clinical_depression2,2])   - y_post) ^ 2 )
    out_tab[z,27] <- d_clinical_depression_post    <- sqrt( ( mean(pcoa.dist$vectors[clinical_depression,1])   - x_post) ^ 2 + (mean(pcoa.dist$vectors[clinical_depression,2])   - y_post) ^ 2 )
    out_tab[z,28] <- d_clinical_anger2_post    <- sqrt( ( mean(pcoa.dist$vectors[clinical_anger2,1])   - x_post) ^ 2 + (mean(pcoa.dist$vectors[clinical_anger2,2])   - y_post) ^ 2 )
    out_tab[z,29] <- d_clinical_anger_post    <- sqrt( ( mean(pcoa.dist$vectors[clinical_anger,1])   - x_post) ^ 2 + (mean(pcoa.dist$vectors[clinical_anger,2])   - y_post) ^ 2 )
    out_tab[z,30] <- d_clinical_positive2_post    <- sqrt( ( mean(pcoa.dist$vectors[clinical_positive2,1])   - x_post) ^ 2 + (mean(pcoa.dist$vectors[clinical_positive2,2])   - y_post) ^ 2 )
    out_tab[z,31] <- d_clinical_positive_post    <- sqrt( ( mean(pcoa.dist$vectors[clinical_positive,1])   - x_post) ^ 2 + (mean(pcoa.dist$vectors[clinical_positive,2])   - y_post) ^ 2 )
    out_tab[z,32] <- d_clinical_negative_post    <- sqrt( ( mean(pcoa.dist$vectors[clinical_negative,1])   - x_post) ^ 2 + (mean(pcoa.dist$vectors[clinical_negative,2])   - y_post) ^ 2 )
    
    out_tab[z,33] <- abs(d_clinical_anxiety2_post) - abs(d_clinical_anxiety2_pre)
    out_tab[z,34] <- abs(d_clinical_anxiety_post) - abs(d_clinical_anxiety_pre)
    out_tab[z,35] <- abs(d_clinical_depression2_post) - abs(d_clinical_depression2_pre)
    out_tab[z,36] <- abs(d_clinical_depression_post) - abs(d_clinical_depression_pre)
    out_tab[z,37] <- abs(d_clinical_anger2_post) - abs(d_clinical_anger2_pre)
    out_tab[z,38] <- abs(d_clinical_anger_post) - abs(d_clinical_anger_pre)
    out_tab[z,39] <- abs(d_clinical_positive2_post) - abs(d_clinical_positive2_pre)
    out_tab[z,40] <- abs(d_clinical_positive_post) - abs(d_clinical_positive_pre)
    out_tab[z,41] <- abs(d_clinical_negative_post) - abs(d_clinical_negative_pre)
    
    
    #### DEVEL ####
    # closest cluster pre/post
    best_pre <- 0
    max_pre <- 99
    best_post <- 0
    max_post <- 99
    for(r in 1:k.best){
      
      w <- sqrt( (x_pre - -pk$medoids[r,1]) ^ 2 + (y_pre - -pk$medoids[r,2]) ^ 2 )
      if(w < max_pre){
        best_pre <- r
        max_pre <- w
      }
      # print(paste("pre",r,w,best_pre))
      
      w <- sqrt( (x_post - -pk$medoids[r,1]) ^ 2 + (y_post - -pk$medoids[r,2]) ^ 2 )
      if(w < max_post){
        best_post <- r
        max_post <- w
      }
      # print(paste("post",r,w,best_post))
      
    }
    # plot(pcoa.dist$vectors[,1:2],
    #      xlab=paste("PCoA.1(",round(100*pcoa.dist$values[1,3], digits=2),"%)"), 
    #      ylab=paste("PCoA.2(",round(100*pcoa.dist$values[2,3], digits=2),"%)"),
    #      main=paste0("Emotion Scores from ", dim(allinput)[1]," Samplings\nBray-Curtis Distance - PCoA\n"),
    #      cex=2,
    #      pch=19,
    #      col=unlist(grouped_cols)
    # ) 
    # points(x_pre,y_pre,col="darkorange4")
    # points(x_post,y_post,pch=8,col="darkorange1")
    # points(-pk$medoids[best_pre,1],-pk$medoids[best_pre,2],col="darkorange4")
    # points(-pk$medoids[best_post,1],-pk$medoids[best_post,2],pch=8,col="darkorange1")
    
    
    out_tab[z,42] <- best_pre
    out_tab[z,43] <- best_post
    out_tab[z,44] <- max_pre
    out_tab[z,45] <- max_post
    
    
  }else{
    problem_list[pi] <- z
    print(paste("PROBLEM: ",pi,z))
    pi <- pi + 1
  }
  
  #points(-pk$medoids[3,1],-pk$medoids[3,2],pch=4)
  #points(-pk$medoids[3,1],-pk$medoids[3,2],pch=4)
  
}
write.table(unlist(problem_list), file = paste0("problem_list.",startx,"to",endx,".tsv"), sep="\t", row.names = FALSE, col.names = FALSE)
write.table(out_tab[startx:endx,], file = paste0("out_tab.",startx,"to",endx,".tsv"), sep="\t")
save.image(file=paste0("YOUTH.incomplete_x",endx,".Rdata"))
