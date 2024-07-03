#######
# 1.1 #
#######

plot.new()

plot(0,0,xlim=c(-1,5),ylim=c(-1,3),bg="red",pch=22)
par(new=TRUE)
plot(2,2,xlim=c(-1,5),ylim=c(-1,3),bg="green",pch=22)
par(new=TRUE)
plot(4,0,xlim=c(-1,5),ylim=c(-1,3),bg="blue",pch=22)
x <- seq(-1,3,by=0.1)
par(new=TRUE)
abline(v=2)
par(new=TRUE)
curve(x-2,xlim=c(-1,5),ylim=c(-1,3))
par(new=TRUE)
curve(2-x,xlim=c(-1,5),ylim=c(-1,3))


# 1.2 # 

plot.new()

plot(0,0,xlim=c(-1,5),ylim=c(-1,3),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(2,2,xlim=c(-1,5),ylim=c(-1,3),bg="green",pch=22,cex=3.0)
par(new=TRUE)
plot(4,0,xlim=c(-1,5),ylim=c(-1,3),bg="blue",pch=22,cex=3.0)
x <- seq(-1,3,by=0.1)
par(new=TRUE)
curve((x/2)-0.5,xlim=c(-1,5),ylim=c(-1,3))
par(new=TRUE)
curve(2.5-(x/2),xlim=c(-1,5),ylim=c(-1,3),col="purple")
par(new=TRUE)
abline(v=2)



# 1.3.1 #


plot.new()

plot(0,0,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(1,1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(-1,1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="blue",pch=22,cex=3.0)
x <- seq(-1,3,by=0.1)
par(new=TRUE)
curve(x+1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5))
par(new=TRUE)
curve(-x+1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5))
par(new=TRUE)
abline(v=0)




# 1.3.2 #
plot.new()
plot(0,0,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(1,1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(-1,1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="blue",pch=22,cex=3.0)
x <- seq(-1,3,by=0.1)
par(new=TRUE)
curve((x/2)+0.25,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5))
par(new=TRUE)
curve(-(x/2)+0.25,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5))
par(new=TRUE)
abline(v=0)






# 1.3.2 B #
plot.new()
plot(0,0,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(1,1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(-1,1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="blue",pch=22,cex=3.0)
x <- seq(-1,3,by=0.1)
par(new=TRUE)
curve((x/2)-1/2,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5))
par(new=TRUE)
curve(-(x/2)-1/25,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5))
par(new=TRUE)
abline(v=0)









5/2-x/2

plot.new()

plot(0,0,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(1,1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(-1,1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5),bg="blue",pch=22,cex=3.0)
x <- seq(-1,3,by=0.1)
par(new=TRUE)
curve(x+1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5))
par(new=TRUE)
curve(-x+1,xlim=c(-1.5,1.5),ylim=c(-0.5,1.5))
par(new=TRUE)
abline(v=0)


plot.new()

plot(0,0,xlim=c(-1.5,6),ylim=c(-0.5,2),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(1,1,xlim=c(-1.5,6),ylim=c(-0.5,2),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(-1,1,xlim=c(-1.5,6),ylim=c(-0.5,2),bg="blue",pch=22,cex=3.0)
x <- seq(-2,7,by=0.1)







par(new=TRUE)
curve(x-0.25*x^2+1,xlim=c(-1.5,6),ylim=c(-0.5,2))
par(new=TRUE)
curve(-0.5*x^2-x+1,xlim=c(-1.5,6),ylim=c(-0.5,2))
par(new=TRUE)
curve(x+1,xlim=c(-1.5,6),ylim=c(-0.5,2))
par(new=TRUE)
curve(-x+1,xlim=c(-1.5,6),ylim=c(-0.5,2))
par(new=TRUE)
abline(v=3+2*sqrt(2))





plot.new()

plot(0,0,xlim=c(-1.5,6),ylim=c(-0.5,2),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(1,1,xlim=c(-1.5,6),ylim=c(-0.5,2),bg="red",pch=22,cex=3.0)
par(new=TRUE)
plot(-1,1,xlim=c(-1.5,6),ylim=c(-0.5,2),bg="blue",pch=22,cex=3.0)
x <- seq(-2,7,by=0.1)
par(new=TRUE)
curve(x-0.25*x^2+1,xlim=c(-1.5,6),ylim=c(-0.5,2))
par(new=TRUE)
curve(-0.5*x^2-x+1,xlim=c(-1.5,6),ylim=c(-0.5,2))
par(new=TRUE)
curve(x+1,xlim=c(-1.5,6),ylim=c(-0.5,2))
par(new=TRUE)
curve(-x+1,xlim=c(-1.5,6),ylim=c(-0.5,2))
par(new=TRUE)
abline(v=3+2*sqrt(2))









#par(new=TRUE)
#curve((1/8)*(8*x-x^2),xlim=c(-1,5),ylim=c(-1,3))

#par(new=TRUE)
#curve(((1/8)*(16-x^2)),xlim=c(-1,5),ylim=c(-1,3))

par(new=TRUE)
abline(v=2)
par(new=TRUE)
curve(x-2,xlim=c(-1,5),ylim=c(-1,3))
par(new=TRUE)
curve(2-x,xlim=c(-1,5),ylim=c(-1,3))


par(new=TRUE)











#par(new=TRUE)
#curve(x,xlim=c(-1,5),ylim=c(-1,3))

x <- seq(-1,3,by=0.1)
y <- seq(-1,5,by=0.1)
color_assign <- matrix(ncol=length(y),nrow = length(x))
colnames(color_assign) <- y
rownames(color_assign) <- x
for(xx in x){
  for (yy in y){
    color_assign[xx][yy] <- "gray"
  }
}
for(xx in x){
  for (yy in y){
    xx <- as.numeric(xx)
    if (xx < 2){
      if(yy<=xx){
        color_assign[xx,yy] = "red"
        print(paste("red : ",xx, yy))
        par(new=TRUE)
        plot(xx,yy,xlim=c(-1,5),ylim=c(-1,3),col="red",pch=22,cex=0.7)
      }
    }
    if (xx > 2){
      if(yy >= 4-xx){
        color_assign[xx,yy] = "blue"
        print(paste("blue : ",xx, yy))
        par(new=TRUE)
        plot(xx,yy,xlim=c(-1,5),ylim=c(-1,3),col="blue",pch=22,cex=0.7)
      }else{
        
      }
    }
    if(xx < yy){
      if(yy < 4-xx){
        color_assign[xx,yy] = "green"
        print(paste("green : ",xx, yy))
        par(new=TRUE)
        plot(xx,yy,xlim=c(-1,5),ylim=c(-1,3),col="green",pch=22,cex=0.7)
      }
    }
  }
}

plot(color_assign)

