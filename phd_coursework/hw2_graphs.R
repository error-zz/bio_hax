
r22b <- function (x) 1+(-(0.5)^x)
r22b2 <- function (x) 1+(-(0.5)^x) - 1+(-(0.5)^x)-1
plot(r22b, 0,10,main="1+(-(0.5)^x)", lwd=3, col="purple")



r1 <- function(x) 1/(1+2.71828182846^(-x))
r2 <- function(x) ( 1/(1+2.71828182846^(-x)) ) * ( 1/(1+2.71828182846^(x)) )
r3 <- function(x) ( 1/(1+2.71828182846^(-x)) ) + ( 1/(1+2.71828182846^(x)) )
r4 <- function(x) log ( ( 1/(1+2.71828182846^(-x)))/( 1-( 1/(1+2.71828182846^(-x)))  ))
r4validate <- function (x) x
plot (r1, -10, 10, col="purple", lwd=3)
plot (r2, -10, 10, col="red", lwd=3, add=TRUE)
plot (r3, -10, 10, col="blue", lwd=3, add=TRUE)
plot (r4, -10, 10, col="darkorange", lwd=3, add=TRUE)
plot (r4validate, -10, 10, col="yellow", pch=22, lty=2, lwd=3, add=TRUE)
abline(0,0)
abline(,,0,0)

