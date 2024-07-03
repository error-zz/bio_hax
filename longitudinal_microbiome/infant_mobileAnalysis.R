library(vegan)
library(ggplot2)
library(directlabels)
library(lme4)
library(Hmisc)
library(gplots)

######
# x  #
######

setwd ("/Users/jmccorri/Desktop/infant")
# load("Genus_eq_Weight_DaysAfter_RunNose_SoreThroat_AnyAntibiotChronPast_given_Weight_RunNose_SoreThroat_AnyAntibiotChronPast_groupBy_InfantID.Rdata")

##################################


load(file="megatable_20181007.Rdata")

#############
# FIX ALL_* #
#############


megatable.sortByInfant[ rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RITM054" , ])[ grepl( "NP1", rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RITM054" , ]) ) ] , ]$compare_id <- "RITM054_NP1"
megatable.sortByInfant[ rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RITM054" , ])[ grepl( "NP2", rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RITM054" , ]) ) ] , ]$compare_id <- "RITM054_NP2"
megatable.sortByInfant[ rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RITM054" , ])[ grepl( "NP3", rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RITM054" , ]) ) ] , ]$compare_id <- "RITM054_NP3"
megatable.sortByInfant[ rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RITM054" , ])[ grepl( "NP4", rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RITM054" , ]) ) ] , ]$compare_id <- "RITM054_NP4"

megatable.sortByInfant[ rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RMPRU063" , ])[ grepl( "NP1", rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RMPRU063" , ]) ) ] , ]$compare_id <- "RMPRU063_NP1"
megatable.sortByInfant[ rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RMPRU063" , ])[ grepl( "NP2", rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RMPRU063" , ]) ) ] , ]$compare_id <- "RMPRU063_NP2"
megatable.sortByInfant[ rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RMPRU063" , ])[ grepl( "NP3", rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RMPRU063" , ]) ) ] , ]$compare_id <- "RMPRU063_NP3"
megatable.sortByInfant[ rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RMPRU063" , ])[ grepl( "NP4", rownames(megatable.sortByInfant[ megatable.sortByInfant$compare_id == "RMPRU063" , ]) ) ] , ]$compare_id <- "RMPRU063_NP4"


# SWITCH TO MEGTABLE.Q!

megtable.q <- megatable.sortByInfant[ order(megatable.sortByInfant$compare_id) , ]

#########
# fever #
#########
megtable.q$fever_1_today    <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$fever_2_past     <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$fever_3_cumPast  <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$fever_4_daysPast <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$fever_5_chronPast <- rep(NA,dim(megatable.sortByInfant)[1])
# today
megtable.q[ !megatable.sortByInfant$all_fever %in% c(" 1", "1"," 2", "2"),  ]$fever_1_today <- NA
megtable.q[ megatable.sortByInfant$all_fever %in% c(" 1", "1"),  ]$fever_1_today <- 1
megtable.q[ megatable.sortByInfant$all_fever %in% c(" 2", "2"),  ]$fever_1_today <- 0
# past
cur_id <- "qqqqqq"
#for(x in 1:dim(megtable.q)[1]){
happened_days_after <- 0
for(x in 1:4654){
  if(!is.na(megtable.q$compare_id[x])){  
    
    if(megtable.q$compare_id[x] != cur_id){
      
      cur_id <- megtable.q$compare_id[x]
      print(cur_id)
      prev <- 0
      
    }
    
    if(prev == 0){
      
      if(is.na( megtable.q[x,"fever_1_today"] )){
        
        megtable.q[x,"fever_2_past"]      <- 0
        megtable.q[x,"fever_3_cumPast"]   <- 0
        megtable.q[x,"fever_4_daysPast"]  <- 0
        megtable.q[x,"fever_5_chronPast"] <- 0
        prev <- 0
        
      }else if(megtable.q[x,"fever_1_today"] == 0){
        
        megtable.q[x,"fever_2_past"]      <- 0
        megtable.q[x,"fever_3_cumPast"]   <- 0
        megtable.q[x,"fever_4_daysPast"]  <- 0
        megtable.q[x,"fever_5_chronPast"] <- 0
        prev <- 0
        
      }else{
        
        megtable.q[x,"fever_2_past"] <- 1
        megtable.q[x,"fever_3_cumPast"] <- megtable.q[x-1,"fever_3_cumPast"] + 1
        megtable.q[x,"fever_4_daysPast"] <- 1
        megtable.q[x,"fever_5_chronPast"] <- 1
        happened_days_after <- megtable.q[x,"all_days_after"]
        prev <- 1
        
      }
      
    }else{
      
      if(is.na( megtable.q[x,"fever_1_today"] )){
        
        megtable.q[x,"fever_2_past"]     <- 1
        megtable.q[x,"fever_3_cumPast"]  <- megtable.q[x-1,"fever_3_cumPast"]
        megtable.q[x,"fever_4_daysPast"] <- megtable.q[x-1,"fever_4_daysPast"]+1
        megtable.q[x,"fever_5_chronPast"] <- megtable.q[x-1,"fever_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- megtable.q[x-1,"fever_2_past"]
        
      }
      else if(megtable.q[x,"fever_1_today"] == 1 ){
        
        megtable.q[x,"fever_2_past"] <- 1
        megtable.q[x,"fever_3_cumPast"] <- megtable.q[x-1,"fever_3_cumPast"]
        megtable.q[x,"fever_4_daysPast"] <- megtable.q[x-1,"fever_4_daysPast"]+1
        megtable.q[x,"fever_5_chronPast"] <- megtable.q[x-1,"fever_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- 1
        
      }else{
        
        megtable.q[x,"fever_2_past"]       <- 0
        megtable.q[x,"fever_3_cumPast"]    <- megtable.q[x-1,"fever_3_cumPast"]
        megtable.q[x,"fever_4_daysPast"]   <- 0 
        # megtable.q[x-1,"fever_4_daysPast"] 
        megtable.q[x,"fever_5_chronPast"]  <- 0 
        prev <- 0
        
      }
    }
    
  }
}

png(filename="fever.observation_debug.png", width=500,height=750)
plot.new()
par( mfrow=c(5,1) )
plot(megtable.q$fever_1_today,     cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="fever : observed")
plot(megtable.q$fever_2_past,      cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="fever : observed in past")
#col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)])
plot(megtable.q$fever_3_cumPast,   cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="fever : cumulative past observation unique events")
plot(megtable.q$fever_4_daysPast,  cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="fever : observations since start of unique event")
plot(megtable.q$fever_5_chronPast, cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="fever : days since start of unique event")
dev.off()







#########
# cough #
#########
megtable.q$cough_1_today    <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$cough_2_past     <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$cough_3_cumPast  <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$cough_4_daysPast <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$cough_5_chronPast <- rep(NA,dim(megatable.sortByInfant)[1])
# today
megtable.q[ !megatable.sortByInfant$all_cough %in% c(" 1", "1"," 2", "2"),  ]$cough_1_today <- NA
megtable.q[ megatable.sortByInfant$all_cough %in% c(" 1", "1"),  ]$cough_1_today <- 1
megtable.q[ megatable.sortByInfant$all_cough %in% c(" 2", "2"),  ]$cough_1_today <- 0
# past
#megtable.q <- megatable.sortByInfant[ order(megatable.sortByInfant$compare_id) , ]
cur_id <- "qqqqqq"
#for(x in 1:dim(megtable.q)[1]){
happened_days_after <- 0
for(x in 1:4654){
  if(!is.na(megtable.q$compare_id[x])){  
    
    if(megtable.q$compare_id[x] != cur_id){
      
      cur_id <- megtable.q$compare_id[x]
      print(cur_id)
      prev <- 0
      
    }
    
    if(prev == 0){
      
      if(is.na( megtable.q[x,"cough_1_today"] )){
        
        megtable.q[x,"cough_2_past"]      <- 0
        megtable.q[x,"cough_3_cumPast"]   <- 0
        megtable.q[x,"cough_4_daysPast"]  <- 0
        megtable.q[x,"cough_5_chronPast"] <- 0
        prev <- 0
        
      }else if(megtable.q[x,"cough_1_today"] == 0){
        
        megtable.q[x,"cough_2_past"]      <- 0
        megtable.q[x,"cough_3_cumPast"]   <- 0
        megtable.q[x,"cough_4_daysPast"]  <- 0
        megtable.q[x,"cough_5_chronPast"] <- 0
        prev <- 0
        
      }else{
        
        megtable.q[x,"cough_2_past"] <- 1
        megtable.q[x,"cough_3_cumPast"] <- megtable.q[x-1,"cough_3_cumPast"] + 1
        megtable.q[x,"cough_4_daysPast"] <- 1
        megtable.q[x,"cough_5_chronPast"] <- 1
        happened_days_after <- megtable.q[x,"all_days_after"]
        prev <- 1
        
      }
      
    }else{
      
      if(is.na( megtable.q[x,"cough_1_today"] )){
        
        megtable.q[x,"cough_2_past"]     <- 1
        megtable.q[x,"cough_3_cumPast"]  <- megtable.q[x-1,"cough_3_cumPast"]
        megtable.q[x,"cough_4_daysPast"] <- megtable.q[x-1,"cough_4_daysPast"]+1
        megtable.q[x,"cough_5_chronPast"] <- megtable.q[x-1,"cough_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- megtable.q[x-1,"cough_2_past"]
        
      }
      else if(megtable.q[x,"cough_1_today"] == 1 ){
        
        megtable.q[x,"cough_2_past"] <- 1
        megtable.q[x,"cough_3_cumPast"] <- megtable.q[x-1,"cough_3_cumPast"]
        megtable.q[x,"cough_4_daysPast"] <- megtable.q[x-1,"cough_4_daysPast"]+1
        megtable.q[x,"cough_5_chronPast"] <- megtable.q[x-1,"cough_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- 1
        
      }else{
        
        megtable.q[x,"cough_2_past"]       <- 0
        megtable.q[x,"cough_3_cumPast"]    <- megtable.q[x-1,"cough_3_cumPast"]
        megtable.q[x,"cough_4_daysPast"]   <- 0 
        # megtable.q[x-1,"cough_4_daysPast"] 
        megtable.q[x,"cough_5_chronPast"]  <- 0 
        prev <- 0
        
      }
    }
    
  }
}

png(filename="cough.observation_debug.png", width=500,height=750)
plot.new()
par( mfrow=c(5,1) )
plot(megtable.q$cough_1_today,     cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="cough : observed")
plot(megtable.q$cough_2_past,      cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="cough : observed in past")
#col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)])
plot(megtable.q$cough_3_cumPast,   cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="cough : cumulative past observation unique events")
plot(megtable.q$cough_4_daysPast,  cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="cough : observations since start of unique event")
plot(megtable.q$cough_5_chronPast, cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="cough : days since start of unique event")
dev.off()



#########
# runnose #
#########
megtable.q$runnose_1_today    <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$runnose_2_past     <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$runnose_3_cumPast  <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$runnose_4_daysPast <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$runnose_5_chronPast <- rep(NA,dim(megatable.sortByInfant)[1])
# today
megtable.q[ !megatable.sortByInfant$all_runnose %in% c(" 1", "1"," 2", "2"),  ]$runnose_1_today <- NA
megtable.q[ megatable.sortByInfant$all_runnose %in% c(" 1", "1"),  ]$runnose_1_today <- 1
megtable.q[ megatable.sortByInfant$all_runnose %in% c(" 2", "2"),  ]$runnose_1_today <- 0
# past
#megtable.q <- megatable.sortByInfant[ order(megatable.sortByInfant$compare_id) , ]
cur_id <- "qqqqqq"
#for(x in 1:dim(megtable.q)[1]){
happened_days_after <- 0
for(x in 1:4654){
  if(!is.na(megtable.q$compare_id[x])){  
    
    if(megtable.q$compare_id[x] != cur_id){
      
      cur_id <- megtable.q$compare_id[x]
      print(cur_id)
      prev <- 0
      
    }
    
    if(prev == 0){
      
      if(is.na( megtable.q[x,"runnose_1_today"] )){
        
        megtable.q[x,"runnose_2_past"]      <- 0
        megtable.q[x,"runnose_3_cumPast"]   <- 0
        megtable.q[x,"runnose_4_daysPast"]  <- 0
        megtable.q[x,"runnose_5_chronPast"] <- 0
        prev <- 0
        
      }else if(megtable.q[x,"runnose_1_today"] == 0){
        
        megtable.q[x,"runnose_2_past"]      <- 0
        megtable.q[x,"runnose_3_cumPast"]   <- 0
        megtable.q[x,"runnose_4_daysPast"]  <- 0
        megtable.q[x,"runnose_5_chronPast"] <- 0
        prev <- 0
        
      }else{
        
        megtable.q[x,"runnose_2_past"] <- 1
        megtable.q[x,"runnose_3_cumPast"] <- megtable.q[x-1,"runnose_3_cumPast"] + 1
        megtable.q[x,"runnose_4_daysPast"] <- 1
        megtable.q[x,"runnose_5_chronPast"] <- 1
        happened_days_after <- megtable.q[x,"all_days_after"]
        prev <- 1
        
      }
      
    }else{
      
      if(is.na( megtable.q[x,"runnose_1_today"] )){
        
        megtable.q[x,"runnose_2_past"]     <- 1
        megtable.q[x,"runnose_3_cumPast"]  <- megtable.q[x-1,"runnose_3_cumPast"]
        megtable.q[x,"runnose_4_daysPast"] <- megtable.q[x-1,"runnose_4_daysPast"]+1
        megtable.q[x,"runnose_5_chronPast"] <- megtable.q[x-1,"runnose_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- megtable.q[x-1,"runnose_2_past"]
        
      }
      else if(megtable.q[x,"runnose_1_today"] == 1 ){
        
        megtable.q[x,"runnose_2_past"] <- 1
        megtable.q[x,"runnose_3_cumPast"] <- megtable.q[x-1,"runnose_3_cumPast"]
        megtable.q[x,"runnose_4_daysPast"] <- megtable.q[x-1,"runnose_4_daysPast"]+1
        megtable.q[x,"runnose_5_chronPast"] <- megtable.q[x-1,"runnose_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- 1
        
      }else{
        
        megtable.q[x,"runnose_2_past"]       <- 0
        megtable.q[x,"runnose_3_cumPast"]    <- megtable.q[x-1,"runnose_3_cumPast"]
        megtable.q[x,"runnose_4_daysPast"]   <- 0 
        # megtable.q[x-1,"runnose_4_daysPast"] 
        megtable.q[x,"runnose_5_chronPast"]  <- 0 
        prev <- 0
        
      }
    }
    
  }
}

png(filename="runnose.observation_debug.png", width=500,height=750)
plot.new()
par( mfrow=c(5,1) )
plot(megtable.q$runnose_1_today,     cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="runnose : observed")
plot(megtable.q$runnose_2_past,      cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="runnose : observed in past")
#col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)])
plot(megtable.q$runnose_3_cumPast,   cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="runnose : cumulative past observation unique events")
plot(megtable.q$runnose_4_daysPast,  cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="runnose : observations since start of unique event")
plot(megtable.q$runnose_5_chronPast, cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="runnose : days since start of unique event")
dev.off()

# TEMP, SORETHROAT, IMPRESSION #

megtable.q[,"all_temp"]         <- rep(NA,dim(megtable.q)[1])
megtable.q[,"all_sorethroat"]   <- rep(NA,dim(megtable.q)[1])
megtable.q[,"all_impression"]   <- rep(NA,dim(megtable.q)[1])
tmplist <- list()
for(g in 1:dim(megtable.q)[1]){
  
  if(!is.na( megtable.q[g,"Temp"]        )){   megtable.q[g,"all_temp"] <- megtable.q[g,"Temp"]   }
  else if(!is.na( megtable.q[g,"ACTEMP"] )){   megtable.q[g,"all_temp"] <- megtable.q[g,"ACTEMP"]
  }else{                                                   megtable.q[g,"all_temp"] <- NA                                   }
  
  if(!is.na( megtable.q[g,"Exam_pharyn"]   )){  megtable.q[g,"all_sorethroat"] <- megtable.q[g,"Exam_pharyn"]   }
  else if(!is.na( megtable.q[g,"SORETHRT"] )){  megtable.q[g,"all_sorethroat"] <- megtable.q[g,"SORETHRT"]
  }else{                                                    megtable.q[g,"all_sorethroat"] <- NA                                   }
  
  if(!is.na( megtable.q[g,"Impress_diag"]    )){  megtable.q[g,"all_impression"] <- megtable.q[g,"Impress_diag"]   }
  else if(!is.na( megtable.q[g,"Impression"] )){  megtable.q[g,"all_impression"] <- megtable.q[g,"Impression"] }
  else if(!is.na( megtable.q[g,"DIAG2SE1"] )){    megtable.q[g,"all_impression"] <- megtable.q[g,"DIAG2SE1"]
  }else{                                                      megtable.q[g,"all_impression"] <- NA                                   }
  
}

# temp
qr            <- as.numeric(megtable.q$all_temp) > 50
qr[is.na(qr)] <- FALSE
megtable.q[qr,"all_temp"] <- NA
plot(megtable.q$all_temp, col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)], main="temperature")

# sorethroat

megtable.q[ ! megtable.q$all_sorethroat %in% c(" 1", "1"," 2", "2"),  ]$all_sorethroat <- NA
megtable.q[ megtable.q$all_sorethroat %in% c(" 1", "1"),  ]$all_sorethroat <- 1
megtable.q[ megtable.q$all_sorethroat %in% c(" 2", "2"),  ]$all_sorethroat <- 0
plot(megtable.q$all_sorethroat, col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)], main="sorethroat",ylim=c(-0.1,1.1))

#########
# sorethroat #
#########
megtable.q$sorethroat_1_today    <- megtable.q$all_sorethroat #rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$sorethroat_2_past     <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$sorethroat_3_cumPast  <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$sorethroat_4_daysPast <- rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$sorethroat_5_chronPast <- rep(NA,dim(megatable.sortByInfant)[1])
# today
megtable.q[ !megtable.q$all_sorethroat %in% c(" 1", "1"," 2", "2"),  ]$sorethroat_1_today <- NA
megtable.q[ megtable.q$all_sorethroat %in% c(" 1", "1"),  ]$sorethroat_1_today <- 1
megtable.q[ megtable.q$all_sorethroat %in% c(" 2", "2"),  ]$sorethroat_1_today <- 0
# past
#megtable.q <- megatable.sortByInfant[ order(megatable.sortByInfant$compare_id) , ]
cur_id <- "qqqqqq"
#for(x in 1:dim(megtable.q)[1]){
happened_days_after <- 0
for(x in 1:4654){
  if(!is.na(megtable.q$compare_id[x])){  
    
    if(megtable.q$compare_id[x] != cur_id){
      
      cur_id <- megtable.q$compare_id[x]
      print(cur_id)
      prev <- 0
      
    }
    
    if(prev == 0){
      
      if(is.na( megtable.q[x,"sorethroat_1_today"] )){
        
        megtable.q[x,"sorethroat_2_past"]      <- 0
        megtable.q[x,"sorethroat_3_cumPast"]   <- 0
        megtable.q[x,"sorethroat_4_daysPast"]  <- 0
        megtable.q[x,"sorethroat_5_chronPast"] <- 0
        prev <- 0
        
      }else if(megtable.q[x,"sorethroat_1_today"] == 0){
        
        megtable.q[x,"sorethroat_2_past"]      <- 0
        megtable.q[x,"sorethroat_3_cumPast"]   <- 0
        megtable.q[x,"sorethroat_4_daysPast"]  <- 0
        megtable.q[x,"sorethroat_5_chronPast"] <- 0
        prev <- 0
        
      }else{
        
        megtable.q[x,"sorethroat_2_past"] <- 1
        megtable.q[x,"sorethroat_3_cumPast"] <- megtable.q[x-1,"sorethroat_3_cumPast"] + 1
        megtable.q[x,"sorethroat_4_daysPast"] <- 1
        megtable.q[x,"sorethroat_5_chronPast"] <- 1
        happened_days_after <- megtable.q[x,"all_days_after"]
        prev <- 1
        
      }
      
    }else{
      
      if(is.na( megtable.q[x,"sorethroat_1_today"] )){
        
        megtable.q[x,"sorethroat_2_past"]     <- 1
        megtable.q[x,"sorethroat_3_cumPast"]  <- megtable.q[x-1,"sorethroat_3_cumPast"]
        megtable.q[x,"sorethroat_4_daysPast"] <- megtable.q[x-1,"sorethroat_4_daysPast"]+1
        megtable.q[x,"sorethroat_5_chronPast"] <- megtable.q[x-1,"sorethroat_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- megtable.q[x-1,"sorethroat_2_past"]
        
      }
      else if(megtable.q[x,"sorethroat_1_today"] == 1 ){
        
        megtable.q[x,"sorethroat_2_past"] <- 1
        megtable.q[x,"sorethroat_3_cumPast"] <- megtable.q[x-1,"sorethroat_3_cumPast"]
        megtable.q[x,"sorethroat_4_daysPast"] <- megtable.q[x-1,"sorethroat_4_daysPast"]+1
        megtable.q[x,"sorethroat_5_chronPast"] <- megtable.q[x-1,"sorethroat_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- 1
        
      }else{
        
        megtable.q[x,"sorethroat_2_past"]       <- 0
        megtable.q[x,"sorethroat_3_cumPast"]    <- megtable.q[x-1,"sorethroat_3_cumPast"]
        megtable.q[x,"sorethroat_4_daysPast"]   <- 0 
        # megtable.q[x-1,"sorethroat_4_daysPast"] 
        megtable.q[x,"sorethroat_5_chronPast"]  <- 0 
        prev <- 0
        
      }
    }
    
  }
}


png(filename="sorethroat.observation_debug.png", width=500,height=750)
plot.new()
par( mfrow=c(5,1) )
plot(megtable.q$sorethroat_1_today,     cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="sorethroat : observed")
plot(megtable.q$sorethroat_2_past,      cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="sorethroat : observed in past")
#col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)])
plot(megtable.q$sorethroat_3_cumPast,   cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="sorethroat : cumulative past observation unique events")
plot(megtable.q$sorethroat_4_daysPast,  cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="sorethroat : observations since start of unique event")
plot(megtable.q$sorethroat_5_chronPast, cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="sorethroat : days since start of unique event")
dev.off()




#######
# PCV # 
#######

plot.new()
par(mfrow=c(3,1))
plot(megtable.q$PCV1DATE_days_from)
plot(megtable.q$PCV2DATE_days_from)
plot(megtable.q$PCV3DATE_days_from)

megtable.q$pcv_1_today     <- rep(NA,dim(megtable.q)[1]) #rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$pcv_2_past      <- rep(NA,dim(megtable.q)[1])
megtable.q$pcv_3_cumPast   <- rep(NA,dim(megtable.q)[1])
megtable.q$pcv_4_daysPast  <- rep(NA,dim(megtable.q)[1])
megtable.q$pcv_5_chronPast <- rep(NA,dim(megtable.q)[1])

plot.new()
par(mfrow=c(5,1))

megtable.q[ ! megtable.q$PCV %in% c(" 1", "1"," 2", "2"),  ]$pcv_1_today <- NA
megtable.q[ megtable.q$PCV %in% c(" 1", "1"),  ]$pcv_1_today <- 1
megtable.q[ megtable.q$PCV %in% c(" 2", "2"),  ]$pcv_1_today <- 0
plot(megtable.q$pcv_1_today, col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)], main="PCV",ylim=c(-0.1,1.1))

#megtable.q[ ! megtable.q$PCVboost %in% c(" 1", "1"," 2", "2"),  ]$pcv_1_today <- NA
megtable.q[ megtable.q$PCVboost %in% c(" 1", "1"),  ]$pcv_1_today <- 1
megtable.q[ megtable.q$PCVboost %in% c(" 2", "2"),  ]$pcv_1_today <- 0
plot(megtable.q$pcv_1_today, col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)], main="PCV+BOOST",ylim=c(-0.1,1.1))

for(q in 1:dim(megtable.q)[1]){
  if( (!is.na(megtable.q$PCV1DATE[q]) && megtable.q$PCV1DATE[q] != megtable.q$PCV1DATE[q-1]) || (!is.na(megtable.q$PCV1DATE[q]) && megtable.q$compare_id[q] != megtable.q$compare_id[q-1]) ){
    megtable.q[q,"pcv_1_today"] <- 1
  }
  else if( (!is.na(megtable.q$PCV2DATE[q]) && megtable.q$PCV2DATE[q] != megtable.q$PCV2DATE[q-1]) || (!is.na(megtable.q$PCV2DATE[q]) && megtable.q$compare_id[q] != megtable.q$compare_id[q-1]) ){
    megtable.q[q,"pcv_1_today"] <- 1
  }
  else if( (!is.na(megtable.q$PCV3DATE[q]) && megtable.q$PCV3DATE[q] != megtable.q$PCV3DATE[q-1]) || (!is.na(megtable.q$PCV1DATE[q]) && megtable.q$compare_id[q] != megtable.q$compare_id[q-1]) ){
    megtable.q[q,"pcv_1_today"] <- 1
  }else if( grepl("RMPRU", megtable.q$compare_id[q]) ){
    megtable.q[q,"pcv_1_today"] <- 0
  }
}
#megtable.q[ grepl("RMPRU", megtable.q$compare_id) , ]$pcv_1_today
plot(megtable.q$pcv_1_today, col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)], main="PCV+BOOST+RMPRU",ylim=c(-0.1,1.1))


cur_id <- "qqqqqq"
#for(x in 1:dim(megtable.q)[1]){
happened_days_after <- 0
for(x in 1:4654){
  if(!is.na(megtable.q$compare_id[x])){  
    
    if(megtable.q$compare_id[x] != cur_id){
      
      cur_id <- megtable.q$compare_id[x]
      print(cur_id)
      prev <- 0
      
    }
    
    if(prev == 0){
      
      if(is.na( megtable.q[x,"pcv_1_today"] )){
        
        megtable.q[x,"pcv_2_past"]      <- 0
        megtable.q[x,"pcv_3_cumPast"]   <- 0
        megtable.q[x,"pcv_4_daysPast"]  <- 0
        megtable.q[x,"pcv_5_chronPast"] <- 0
        prev <- 0
        
      }else if(megtable.q[x,"pcv_1_today"] == 0){
        
        megtable.q[x,"pcv_2_past"]      <- 0
        megtable.q[x,"pcv_3_cumPast"]   <- 0
        megtable.q[x,"pcv_4_daysPast"]  <- 0
        megtable.q[x,"pcv_5_chronPast"] <- 0
        prev <- 0
        
      }else{
        
        megtable.q[x,"pcv_2_past"] <- 1
        megtable.q[x,"pcv_3_cumPast"] <- megtable.q[x-1,"pcv_3_cumPast"] + 1
        megtable.q[x,"pcv_4_daysPast"] <- 1
        megtable.q[x,"pcv_5_chronPast"] <- 1
        happened_days_after <- megtable.q[x,"all_days_after"]
        prev <- 1
        
      }
      
    }else{
      
      if(is.na( megtable.q[x,"pcv_1_today"] )){
        
        megtable.q[x,"pcv_2_past"]     <- 1
        megtable.q[x,"pcv_3_cumPast"]  <- megtable.q[x-1,"pcv_3_cumPast"]
        megtable.q[x,"pcv_4_daysPast"] <- megtable.q[x-1,"pcv_4_daysPast"]+1
        megtable.q[x,"pcv_5_chronPast"] <- megtable.q[x-1,"pcv_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- megtable.q[x-1,"pcv_2_past"]
        
      }
      else if(megtable.q[x,"pcv_1_today"] == 1 ){
        
        megtable.q[x,"pcv_2_past"] <- 1
        megtable.q[x,"pcv_3_cumPast"] <- megtable.q[x-1,"pcv_3_cumPast"]
        megtable.q[x,"pcv_4_daysPast"] <- megtable.q[x-1,"pcv_4_daysPast"]+1
        megtable.q[x,"pcv_5_chronPast"] <- megtable.q[x-1,"pcv_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- 1
        
      }else{
        
        megtable.q[x,"pcv_2_past"]       <- 0
        megtable.q[x,"pcv_3_cumPast"]    <- megtable.q[x-1,"pcv_3_cumPast"]
        #megtable.q[x,"pcv_4_daysPast"]   <- 0 
        megtable.q[x,"pcv_4_daysPast"] <- megtable.q[x-1,"pcv_4_daysPast"]+1
        megtable.q[x,"pcv_5_chronPast"] <- megtable.q[x-1,"pcv_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        #megtable.q[x,"pcv_5_chronPast"]  <- 0 
        prev <- 0
        
      }
    }
    
  }
}


png(filename="pcv.observation_debug.png", width=500,height=750)
plot.new()
par( mfrow=c(5,1) )
plot(megtable.q$pcv_1_today,     cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="pcv : vaccination administered")
plot(megtable.q$pcv_2_past,      cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="pcv : vaccination administered in past")
#col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)])
plot(megtable.q$pcv_3_cumPast,   cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="pcv : cumulative past vaccination unique events")
plot(megtable.q$pcv_4_daysPast,  cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="pcv : vaccinations since start of unique vaccination event")
plot(megtable.q$pcv_5_chronPast, cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="pcv : days since start of unique vaccination event")
dev.off()




##############
# ANTIBIOTIC #
##############

megtable.q$antibioticAll_1_today     <- rep(NA,dim(megtable.q)[1]) #rep(NA,dim(megatable.sortByInfant)[1])
megtable.q$antibioticAll_2_past      <- rep(NA,dim(megtable.q)[1])
megtable.q$antibioticAll_3_cumPast   <- rep(NA,dim(megtable.q)[1])
megtable.q$antibioticAll_4_daysPast  <- rep(NA,dim(megtable.q)[1])
megtable.q$antibioticAll_5_chronPast <- rep(NA,dim(megtable.q)[1])

plot.new()
par(mfrow=c(5,1))

megtable.q[ ! megtable.q$AntibiotTaken %in% c(" 1", "1"," 2", "2"),  ]$antibioticAll_1_today <- NA
megtable.q[ megtable.q$AntibiotTaken %in% c(" 1", "1"),  ]$antibioticAll_1_today <- 1
megtable.q[ megtable.q$AntibiotTaken %in% c(" 2", "2"),  ]$antibioticAll_1_today <- 0

#megtable.q[ ! megtable.q$ANTIBIOT %in% c(" 1", "1"," 2", "2"),  ]$antibioticAll_1_today <- NA
megtable.q[ megtable.q$ANTIBIOT %in% c(" 1", "1"),  ]$antibioticAll_1_today <- 1
megtable.q[ megtable.q$ANTIBIOT %in% c(" 2", "2"),  ]$antibioticAll_1_today <- 0

plot(megtable.q$antibioticAll_1_today, col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)], main="ANTIBIOTIC",ylim=c(-0.1,1.1))


cur_id <- "qqqqqq"
#for(x in 1:dim(megtable.q)[1]){
happened_days_after <- 0
for(x in 1:4654){
  if(!is.na(megtable.q$compare_id[x])){  
    
    if(megtable.q$compare_id[x] != cur_id){
      
      cur_id <- megtable.q$compare_id[x]
      print(cur_id)
      prev <- 0
      
    }
    
    if(prev == 0){
      
      if(is.na( megtable.q[x,"antibioticAll_1_today"] )){
        
        megtable.q[x,"antibioticAll_2_past"]      <- 0
        megtable.q[x,"antibioticAll_3_cumPast"]   <- 0
        megtable.q[x,"antibioticAll_4_daysPast"]  <- 0
        megtable.q[x,"antibioticAll_5_chronPast"] <- 0
        prev <- 0
        
      }else if(megtable.q[x,"antibioticAll_1_today"] == 0){
        
        megtable.q[x,"antibioticAll_2_past"]      <- 0
        megtable.q[x,"antibioticAll_3_cumPast"]   <- 0
        megtable.q[x,"antibioticAll_4_daysPast"]  <- 0
        megtable.q[x,"antibioticAll_5_chronPast"] <- 0
        prev <- 0
        
      }else{
        
        megtable.q[x,"antibioticAll_2_past"] <- 1
        megtable.q[x,"antibioticAll_3_cumPast"] <- megtable.q[x-1,"antibioticAll_3_cumPast"] + 1
        megtable.q[x,"antibioticAll_4_daysPast"] <- 1
        megtable.q[x,"antibioticAll_5_chronPast"] <- 1
        happened_days_after <- megtable.q[x,"all_days_after"]
        prev <- 1
        
      }
      
    }else{
      
      if(is.na( megtable.q[x,"antibioticAll_1_today"] )){
        
        megtable.q[x,"antibioticAll_2_past"]     <- 1
        megtable.q[x,"antibioticAll_3_cumPast"]  <- megtable.q[x-1,"antibioticAll_3_cumPast"]
        megtable.q[x,"antibioticAll_4_daysPast"] <- megtable.q[x-1,"antibioticAll_4_daysPast"]+1
        megtable.q[x,"antibioticAll_5_chronPast"] <- megtable.q[x-1,"antibioticAll_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- megtable.q[x-1,"antibioticAll_2_past"]
        
      }
      else if(megtable.q[x,"antibioticAll_1_today"] == 1 ){
        
        megtable.q[x,"antibioticAll_2_past"] <- 1
        megtable.q[x,"antibioticAll_3_cumPast"] <- megtable.q[x-1,"antibioticAll_3_cumPast"]
        megtable.q[x,"antibioticAll_4_daysPast"] <- megtable.q[x-1,"antibioticAll_4_daysPast"]+1
        megtable.q[x,"antibioticAll_5_chronPast"] <- megtable.q[x-1,"antibioticAll_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        prev <- 1
        
      }else{
        
        megtable.q[x,"antibioticAll_2_past"]       <- 0
        megtable.q[x,"antibioticAll_3_cumPast"]    <- megtable.q[x-1,"antibioticAll_3_cumPast"]
        #megtable.q[x,"antibioticAll_4_daysPast"]   <- 0 
        megtable.q[x,"antibioticAll_4_daysPast"] <- megtable.q[x-1,"antibioticAll_4_daysPast"]+1
        megtable.q[x,"antibioticAll_5_chronPast"] <- megtable.q[x-1,"antibioticAll_5_chronPast"] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
        #megtable.q[x,"antibioticAll_5_chronPast"]  <- 0 
        prev <- 0
        
      }
    }
    
  }
}

png(filename="antibioticAll.observation_debug.png", width=500,height=750)
plot.new()
par( mfrow=c(5,1) )
plot(megtable.q$antibioticAll_1_today,     cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="antibioticAll : vaccination administered")
plot(megtable.q$antibioticAll_2_past,      cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="antibioticAll : vaccination administered in past")
#col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)])
plot(megtable.q$antibioticAll_3_cumPast,   cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="antibioticAll : cumulative past vaccination unique events")
plot(megtable.q$antibioticAll_4_daysPast,  cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="antibioticAll : vaccinations since start of unique vaccination event")
plot(megtable.q$antibioticAll_5_chronPast, cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main="antibioticAll : days since start of unique vaccination event")
dev.off()


# Update unique labels
uniqAntibiot       <- as.data.frame(table( c( as.character(megatable.sortByInfant$ANTINAME1), as.character( megatable.sortByInfant$Antibiotname1 ), as.character( megatable.sortByInfant$Antibiotname2 ), as.character( megatable.sortByInfant$Antibiotname3 )) ) )
uniqAntibiot$label <- uniqAntibiot$Var1
uniqAntibiot$label <- as.character(uniqAntibiot$label)
uniqAntibiot$label[as.character(uniqAntibiot$label) %in% c("ANTIBIOTIC", "MEDS UNKNOWN","NAME UNKNOWN","null","UNKNOWN","Unrecalled","Unrecalled (Cloxacillin?)","Unrecalled antibiotics")] <- "Unknown"
uniqAntibiot$label[uniqAntibiot$label %in% c("Amoxicillin","AMOXICILLIN","AMOXIL","AMOXILLIN","AMOXYCILLIN")] <- "Amoxicillin"
uniqAntibiot$label[uniqAntibiot$label %in% c("Ampicillin IV")] <- "Ampicillin"
uniqAntibiot$label[uniqAntibiot$label %in% c("Cefalexin","Cefalexin 125mg","Cefalexin 250mg")] <- "Cefalexin"
uniqAntibiot$label[uniqAntibiot$label %in% c("Co amoxyclav","Co-amoxiclav","Co-Amoxiclav","CO-AMOXICLAV","Coamoxiclav")] <- "Co-amoxiclav"
uniqAntibiot$label[uniqAntibiot$label %in% c("Gentamicin IV")] <- "Gentamicin"
uniqAntibiot$label[uniqAntibiot$label %in% c("Ilosone (Erythromycin)")] <- "Erythromycin"
uniqAntibiot$label[uniqAntibiot$label %in% c("IV","IV antibiotic","IV Antibiotics (Unrecalled)","IV antibiotics","IV Antibiotic","IV Antibiotics","Unknown intravenous antibiotics")] <- "IV"
uniqAntibiot$label[uniqAntibiot$label %in% c("Permethrin lotion")] <- "Permethrin"
uniqAntibiot$label[uniqAntibiot$label == ""] <- NA
uniqAntibiot.noNA <- uniqAntibiot[-1,]
uniqAntibiot.noNA[ uniqAntibiot.noNA$label == "Cefalexin" , ]$label = "Cephalexin"
uniq_counts        <- rep(0, 31)
names(uniq_counts) <- unique(uniqAntibiot.noNA$label)[ !is.na(unique(uniqAntibiot.noNA$label)) ]
for(u in names(uniq_counts)){
  uniq_counts[u] <- sum(uniqAntibiot.noNA[uniqAntibiot.noNA$label == u,"Freq"])
}
par(mar=c(4,12,1,1))
plot.new()
par(mfrow=c(1,1))
barplot(log(sort(uniq_counts)),horiz=TRUE,las=2,main="Antibiotic") 

megtable.q[,"uniqAntiobiotLabel"] <- rep(NA,dim(megtable.q)[1])
for(x in 1:dim(megtable.q)[1]){
  if(!is.na(megtable.q[x,"ANTINAME1"])){
    #print(paste( x, megtable.q[x,"ANTINAME1"] , uniqAntibiot.noNA[ as.character(uniqAntibiot.noNA$Var1)==as.character(megtable.q[x,"ANTINAME1"]), ]$label ))
    megtable.q[x,"uniqAntiobiotLabel"] <- uniqAntibiot.noNA[ as.character(uniqAntibiot.noNA$Var1)==as.character(megtable.q[x,"ANTINAME1"]), ]$label
  }
  if(!is.na(megtable.q[x,"Antibiotname1"]) && megtable.q[x,"Antibiotname1"] != ""){
    megtable.q[x,"uniqAntiobiotLabel"] <- uniqAntibiot.noNA[ as.character(uniqAntibiot.noNA$Var1)==as.character(megtable.q[x,"Antibiotname1"]), ]$label
  }
} 

allr <- list()
rr <- 1
for(u in unique(megtable.q[,"uniqAntiobiotLabel"])){
  lll <- megtable.q$uniqAntiobiotLabel == u
  lll[is.na(lll)] <- FALSE
  mmm <- megtable.q[ lll , ] 
  print(u)
  #table(mmm$SiteID)
  print( table(megtable.q[lll,]$compare_id) )
  for( r in  names(table(megtable.q[lll,]$compare_id)) ){
    allr[rr] <- r
    rr <- rr + 1
  }
}
allr.uniq <- sort(unlist(unique(allr)))
    
allr.abund <- matrix( nrow=length(allr.uniq), ncol=length(unique(megtable.q[,"uniqAntiobiotLabel"])) )
rownames(allr.abund) <- allr.uniq
colnames(allr.abund) <- unique(megtable.q[,"uniqAntiobiotLabel"])
allr.abund[is.na(allr.abund)] <- 0
for(u in unique(megtable.q[,"uniqAntiobiotLabel"])){
  lll <- megtable.q$uniqAntiobiotLabel == u
  lll[is.na(lll)] <- FALSE
  mmm <- megtable.q[ lll , ] 
  print(u)
  for( r in  names(table(megtable.q[lll,]$compare_id)) ){
    allr.abund[r,u] <- table(megtable.q[lll,]$compare_id)[r]
  }
}
heatmap.2(allr.abund,trace='none',col=c("white","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue",
                                        "purple","purple","purple","purple","purple","purple","purple","purple","purple","purple"),
          RowSideColors = c("red","green")[ as.factor(grepl("RMPRU",rownames(allr.abund))) ]
        )
 
allr.counts <- megtable.q[ megtable.q$compare_id %in% rownames(allr.abund) , colnames(genus.i)  ]

# ITERATE THROUGH ALL ANTIBIOTICS AND PARSE TO BASIC RESPONSE CURVES

for(colm in unique( uniqAntibiot.noNA$label )){
  
  print(colm)
  
  colm1 <- paste0("antibiotic.",colm,".1_today")
  colm2 <- paste0("antibiotic.",colm,".2_past")
  colm3 <- paste0("antibiotic.",colm,".3_cumPast")
  colm4 <- paste0("antibiotic.",colm,".4_daysPast")
  colm5 <- paste0("antibiotic.",colm,".5_chronPast")
  megtable.q[,colm1]     <- rep(NA,dim(megtable.q)[1]) #rep(NA,dim(megatable.sortByInfant)[1])
  megtable.q[,colm2]      <- rep(NA,dim(megtable.q)[1])
  megtable.q[,colm3]      <- rep(NA,dim(megtable.q)[1])
  megtable.q[,colm4]      <- rep(NA,dim(megtable.q)[1])
  megtable.q[,colm5]      <- rep(NA,dim(megtable.q)[1])
  
  #megtable.q[ !megtable.q$AntibiotTaken %in% c(" 1", "1"," 2", "2"), colm1 ] <- NA
  #megtable.q[ megtable.q$AntibiotTaken %in% c(" 1", "1"), colm1 ] <- 1
  #megtable.q[ megtable.q$AntibiotTaken %in% c(" 2", "2"), colm1 ] <- 0
  megtable.q[ ! megtable.q$Antibiotname1 %in% as.character( uniqAntibiot.noNA[ uniqAntibiot.noNA$label == colm , 1 ] ) , colm1 ] <- 0
  megtable.q[ ! megtable.q$ANTINAME1 %in% as.character( uniqAntibiot.noNA[ uniqAntibiot.noNA$label == colm , 1 ] ) , colm1 ] <- 0
  megtable.q[ megtable.q$Antibiotname1 %in% as.character( uniqAntibiot.noNA[ uniqAntibiot.noNA$label == colm , 1 ] ) , colm1 ] <- 1
  megtable.q[ megtable.q$ANTINAME1 %in% as.character( uniqAntibiot.noNA[ uniqAntibiot.noNA$label == colm , 1 ] ) , colm1 ] <- 1
  
  #megtable.q[ megtable.q$ANTIBIOT %in% c(" 1", "1"), colm1 ] <- 1
  #megtable.q[ megtable.q$ANTIBIOT %in% c(" 2", "2"), colm1 ] <- 0
  
  #plot(megtable.q$colm1, col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)], main=colm ,xlim=c(1,dim(megtable.q)[1]),ylim=c(-0.1,1.1))
    
  cur_id <- "qqqqqq"
  #for(x in 1:dim(megtable.q)[1]){
  happened_days_after <- 0
  for(x in 1:4654){
    if(!is.na(megtable.q$compare_id[x])){  
      
      if(megtable.q$compare_id[x] != cur_id){
        
        cur_id <- megtable.q$compare_id[x]
        print(cur_id)
        prev <- 0
        
      }
      
      if(prev == 0){
        
        if(is.na( megtable.q[x,colm1] )){
          
          megtable.q[x,colm2] <- 0
          megtable.q[x,colm3] <- 0
          megtable.q[x,colm4] <- 0
          megtable.q[x,colm5] <- 0
          prev <- 0
          
        }else if(megtable.q[x,colm1] == 0){
          
          megtable.q[x,colm2] <- 0
          megtable.q[x,colm3] <- 0
          megtable.q[x,colm4] <- 0
          megtable.q[x,colm5] <- 0
          prev <- 0
          
        }else{
          
          megtable.q[x,colm2] <- 1
          if( x > 2 ){
            megtable.q[x,colm3] <- megtable.q[x-1,colm3] + 1
          }else{
            megtable.q[x,colm3] <- 1
          }
          megtable.q[x,colm4] <- 1
          megtable.q[x,colm5] <- 1
          happened_days_after <- megtable.q[x,"all_days_after"]
          prev <- 1
          
        }
        
      }else{
        
        if(is.na( megtable.q[x,colm1] )){
          
          megtable.q[x,colm2]     <- 1
          megtable.q[x,colm3]  <- megtable.q[x-1,colm3]
          megtable.q[x,colm4] <- megtable.q[x-1,colm4]+1
          megtable.q[x,colm5] <- megtable.q[x-1,colm5] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
          prev <- megtable.q[x-1,colm2]
          
        }
        else if(megtable.q[x,colm1] == 1 ){
          
          megtable.q[x,colm2] <- 1
          megtable.q[x,colm3] <- megtable.q[x-1,colm3]
          megtable.q[x,colm4] <- megtable.q[x-1,colm4]+1
          megtable.q[x,colm5] <- megtable.q[x-1,colm5] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
          prev <- 1
          
        }else{
          
          megtable.q[x,colm2]       <- 0
          megtable.q[x,colm3]    <- megtable.q[x-1,colm3]
          #megtable.q[x,"antibioticAll_4_daysPast"]   <- 0 
          megtable.q[x,colm4] <- megtable.q[x-1,colm4]+1
          megtable.q[x,colm5] <- megtable.q[x-1,colm5] + (as.numeric(megtable.q[x,]$all_days_after) - as.numeric(happened_days_after))
          #megtable.q[x,"antibioticAll_5_chronPast"]  <- 0 
          prev <- 0
          
        }
      }
      
    }
  }
  
  plot.new()
  par(mfrow=c(5,1))
  png(filename=paste0("antibiotic_indicator_variables.",colm,".observation_debug.png"), width=500,height=750)
  plot.new()
  par( mfrow=c(5,1) )
  plot(megtable.q[ , colm1 ], cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main=paste0(colm," : antibiotic administered"))
  plot(megtable.q[ , colm2  ],      cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main=paste0(colm," : antibiotic administered in past"))
  #col=rainbow(length(unique(megtable.q$compare_id)))[as.factor(megtable.q$compare_id)])
  plot(megtable.q[ , colm3 ],   cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main=paste0(colm," : cumulative past antibiotic unique events"))
  plot(megtable.q[ , colm4 ],  cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main=paste0(colm," : antibiotics since start of unique antibiotic event"))
  plot(megtable.q[ , colm5 ], cex=0.3,col=rep(c("red","darkblue"), 263)[as.factor(megtable.q$compare_id)], main=paste0(colm," : days since start of unique antibiotic event"))
  dev.off()
  
  
  sum_matrix <- matrix(nrow=length(colnames(genus.i)),ncol=3)
  rownames(sum_matrix) <- colnames(genus.i)
  comparable_set <- megtable.q[ !is.na(megtable.q$all_days_after), ]
  comparable_set <- comparable_set[ !is.na(comparable_set$all_weight), ]
  comparable_set <- comparable_set[ !is.na(comparable_set$all_height), ]
  comparable_set <- comparable_set[ !is.na(comparable_set$all_runnose), ]
  comparable_set <- comparable_set[ !is.na(comparable_set$all_sorethroat), ]
  
  colm_pval_matrix     <- matrix(nrow=dim(genus.i)[2],ncol=8)
  colm_pval_matrix.fdr <- matrix(nrow=dim(genus.i)[2],ncol=8)
  rownames(colm_pval_matrix) <- rownames(colm_pval_matrix.fdr) <- colnames( genus.i )
  colnames(colm_pval_matrix) <- colnames(colm_pval_matrix.fdr) <- c("intercept","all_days_after","all_weight","runnose_5_chronPast","sorethroat_5_chronPast","fever_5_chronPast","cough_5_chronPast","colm")
  gg <- 1
  for(g in colnames(genus.i)){
    
    if( gg %% 100 == 0){
      print( paste(g, (gg / length(colnames(genus.i))) * 100, "%" ) )
    }
    
    greg <- lmer(   as.numeric(comparable_set[,g]) ~ as.numeric(comparable_set$all_days_after) + (as.numeric(comparable_set$all_weight) * as.numeric(comparable_set$all_weight)) + as.numeric(comparable_set$runnose_5_chronPast) + as.numeric(comparable_set$sorethroat_5_chronPast) + as.numeric(comparable_set$fever_5_chronPast) + as.numeric(comparable_set$cough_5_chronPast) 
                    + as.numeric(comparable_set[,colm5])   
                    + (  + as.numeric(comparable_set$runnose_5_chronPast) + as.numeric(comparable_set$sorethroat_5_chronPast) + as.numeric(comparable_set$fever_5_chronPast) + as.numeric(comparable_set$cough_5_chronPast) 
                         + as.numeric(comparable_set[,colm5]) 
                         + ( as.numeric(comparable_set$all_weight) * as.numeric(comparable_set$all_weight)) | comparable_set$compare_id )   )
  
    coefs <- data.frame(coef(summary(greg)))
    coefs.p <- 2 * (1 - pnorm(abs(coefs$t.value)))
    nupnames <- list()
    np <- 1
    pnames <- strsplit( unlist( strsplit( unlist(rownames(as.matrix(coefs))) , ')' ) ) , '$' , fixed=TRUE)
    for(n in pnames){
      nupnames[np] <- n[2]
      np <- np + 1
    }
    pnames.short <- list()
    pi <- 1
    for(p in pnames){
      pnames.short[pi] <- p[2]
      pi <- pi + 1
    }
    pnames.short <- unlist(pnames.short)
    pnames.short[1] <- "intercept"
    pnames.short[8] <- "colm"
    names(coefs.p) <- pnames.short
    coefs.p.fdr <- p.adjust(coefs.p, method="fdr")
    #ppp[qi,names(coefs.p)] <- coefs.p
    #ppp.fdr[qi,names(coefs.p.fdr)] <- coefs.p.fdr
                                                                                                                                            
    gg <- gg + 1
    
  }
  
}
 







