


setwd ("/Users/apple/Desktop/BNFO285_stat_learning/hw3/")
x = read.csv('BNFO285_HW3_P3.csv',header=TRUE);
x <- sort(x[,2])
plot(data, pch=20, col="coral3")
# AIC = −2 log L( O|y ) + 2k
# AIC = bias + variance
binom <- fitdistr(x, "logistic")
binom_aic = AIC(binom)
poiss <- fitdistr(x, "Poisson")
poiss_aic = AIC(poiss)
geome <- fitdistr(x, "geometric")
geome_aic = AIC(geome)
negbi <- fitdistr(x, "negative binomial")
negbi_aic = AIC(negbi)
unifo <- fitdistr(x, "normal")
unifo_aic = AIC(unifo)

print(binom_aic)
print(poiss_aic)
print(geome_aic)
print(negbi_aic)
print(unifo_aic)