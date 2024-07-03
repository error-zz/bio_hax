######################################################
#### statlearning - hw2-4 - j.mccorrison - 2/2/15 ####
######################################################

setwd ("/Users/apple/Desktop/BNFO285_stat_learning/hw2_20150204/")
data = read.csv('hw2-p4-data.csv',header=FALSE);

euro_x <- NULL
euro_y <- NULL
afri_x <- NULL
afri_y <- NULL
query_x <- NULL
query_y <- NULL
euro_c = 0
afri_c = 0
query_c = 0

euro_centroid_x = 15.168
euro_centroid_y = 49.869
afri_centroid_x = 20.844
afri_centroid_y = 7.648
all_euro_dist <- NULL
e_x = 0
all_afri_dist <- NULL
a_x = 0
for(i in 1:nrow(data)){
  if (data[i,'V3'] == "European"){
    euro_x[euro_c]=data[i,'V1']
    euro_y[euro_c]=data[i,'V2']
    all_euro_dist[e_x] = sqrt( (abs(euro_x[euro_c]-euro_centroid_x))^2 + (abs(euro_y[euro_c]-euro_centroid_y))^2 )
    e_x = e_x + 1
    euro_c = euro_c + 1
  }
  else if (data[i,'V3'] == "African"){
    afri_x[afri_c]=data[i,'V1']
    afri_y[afri_c]=data[i,'V2']
    all_afri_dist[e_x] = sqrt( (abs(afri_x[afri_c]-afri_centroid_x))^2 + (abs(afri_y[euro_c]-afri_centroid_y))^2 )
    a_x = a_x + 1
    afri_c = afri_c + 1
  }
  else{
    query_x[query_c]=data[i,'V1']
    query_y[query_c]=data[i,'V2']
    query_c = query_c + 1
  }
}

setwd ("/Users/apple/Desktop/BNFO285_stat_learning/hw2_20150204/")
data = read.table('perl_output_sorted.txt',header=FALSE);

# probability test
query_id <- NULL
euro_dist <- NULL
afri_dist <- NULL
for(i in 1:nrow(data)){
  query_id[i] = data[i,'V1']
  euro_dist[i] = data[i,'V2']
  afri_dist[i] = data[i,'V3'] 
}

all_dist = c(euro_dist, afri_dist)
mean_dist = mean(all_dist)
sd_dist = sd(all_dist)

# EURO AFRI
for(i in 1:nrow(data)){
  t = (abs(euro_dist[i] - mean_dist))/(sd_dist/sqrt(nrow(data)))
  prob = 2*pt(-abs(t),df=nrow(data)-1)
  print(query_id[i])
  print(prob)
}

for(i in 1:nrow(data)){
  t = (abs(afri_dist[i] - mean_dist))/(sd_dist/sqrt(nrow(data)))
  prob = 2*pt(-abs(t),df=nrow(data)-1)
  print(query_id[i])
  print(prob)
}

