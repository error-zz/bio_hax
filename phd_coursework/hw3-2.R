
p<-function(r){
  c( (1-r)^2, (2*r)*(1-r), r^2 )
}

r.est<-function(N,eps=10^(-10)){
  n<-sum(N)
  r<-.1
  dif=1
  while(dif>eps){
    r.old<-r
    phi<-(r^2)/((1-r)^2+r^2)
    r<-1/(2*n)*(2*(N[3,1]+N[1,3])+
                  N[3,2]+N[2,3]+N[2,1]+N[1,2]+
                  2*phi*N[2,2])
    dif<-abs(r-r.old)
  }
  r
}

n= 1029
r= 0.4247

N<-rmultinom(1,n,p(r))
r<-r.est(matrix(N,3,3))
LR<-2*sum(N*log(p(r))-N*log(p(1/2))) 
X2<- pchisq(LR,1,lower=FALSE)



