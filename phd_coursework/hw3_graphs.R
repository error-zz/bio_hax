# problem 3.3.a

f_n <- function(i){
  eval.parent(substitute(i <- i + ( (2^(i-1))*i )))
} 

f_b <- function (n) {
  sum <- 0
  sum <- sum + f_n(n)
  return(sum)
}

f_z <- function (z,a,b) ((1-a)/(1+a))*(a^abs(z-f_b(b)))
  
f_z_individual <- function(z) f_z(z,a,b)

f_z_sum <- function(z){
  abs( f_z(z,a,b) ) + abs( f_z_sum(z-1) )
}

for (b in seq(1,6,by=1)){
  plot(f_b, 0, 10, main="f(B) = sum[i:n] ( (2^(i-1))*i )", lwd=3, col="purple")
  for (a in seq(0.05,0.95,by=0.3)){
    plot(f_z_individual, 1, 400, main=paste("((1-a)/(1+a))*(a^abs(z-f_b(b))) a=", a," b=", b), lwd=3, col="red")
    nums <- seq(1,400,by=1)#c(1:400)
    #cur_sum <- array(0,dim=c(1,400))
    #cur_sum[0] <- 0
    for (q in 2:length(nums)+1){
      q <- as.numeric(q)
      cur_sum[q] <- cur_sum[q-1] + f_z_individual(q)
    }
    plot(cur_sum, main=paste("SUM : ((1-a)/(1+a))*(a^abs(z-f_b(b))) a=", a," b=", b), lwd=3, col="blue")
  }
}   




