library(NMF)
library(car)


load("BRCA21.RData")
#Simulate poisson and NB data from the BRCA21 dataset
data <- t(V)
NMF = nmf(data, rank = 6, nrun = 10, method = "KL")
W = basis(NMF)
H = coef(NMF)

V_pois_obs <- matrix(rpois(21*96, lambda = W%*%H), nrow = 21)
V_NB_obs <- matrix(rnbinom(21*96, size = 50, prob = W%*%H/(50+W%*%H)), nrow = 21)

#Estimate W and H in both bases
NMF_pois <- nmf(V_pois_obs, rank = 6, nrun = 10, method = "KL")
W_pois = basis(NMF_pois)
H_pois = coef(NMF_pois)
NMF_NB <- nmf(V_NB_obs, rank = 6, nrun = 10, method = "KL")
W_NB = basis(NMF_NB)
H_NB = coef(NMF_NB)

#Summary function sum of squares deviation for V/WH
A_pois_obs <- V_pois_obs/(W_pois%*%H_pois)
SSA_pois_obs <- sum((A_pois_obs-1)^2)/length(A_pois_obs)
A_NB_obs <- V_NB_obs/(W_NB%*%H_NB)
SSA_NB_obs <- sum((A_NB_obs-1)^2)/length(A_NB_obs)


id = var(as.vector(A_pois_obs))/mean(as.vector(A_pois_obs))
pchisq(id, df = length(as.vector(A_pois_obs))-1)/(length(as.vector(A_pois_obs))-1)

library(car)
anscombeTrans = function(x) 2*sqrt(x + 3/8)
plot(anscombeTrans(as.vector(data/(W%*%H))))
plot(anscombeTrans(as.vector(data)), anscombeTrans(as.vector(W%*%H)) - anscombeTrans(as.vector(data)), xlim = c(0,50))
shapiro.test(anscombeTrans(as.vector(data/(W%*%H))))
shapiro.test(anscombeTrans(as.vector(W%*%H)) - anscombeTrans(as.vector(data)))

hist(anscombeTrans(as.vector(W%*%H)) - anscombeTrans(as.vector(data)), prob = T, breaks = 50)
curve(dnorm(x, mean = 0, sd = sqrt(1/4)), add = T)
curve(dnorm(x, mean = 0, sd = sd(anscombeTrans(as.vector(W%*%H)) - anscombeTrans(as.vector(data)))), add = T)

qqPlot(anscombeTrans(as.vector(W%*%H)) - anscombeTrans(as.vector(data)))

qqPlot(anscombeTrans(as.vector(data/(W%*%H))))
ihs <- function(x) {
  y <- log(x + sqrt(x^2+1))
  return(y)
}

inverse = function (f, lower = -100, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}
square_inverse = inverse(function (x) x^2, 0.1, 100)
sinh_inverse = inverse(function(x) sinh(x), 0.1, 100)

anscombeTrans_NB = function(x, alpha) sqrt(x-1/2)*ihs(sqrt((alpha + 3/8)/(x-3/4)))
anscombeTrans_NB1 = function(x, alpha) sqrt(x-1/2)*sinh_inverse(sqrt((alpha + 3/8)/(x-3/4)))$root
plot(anscombeTrans_NB(as.vector(data/(W%*%H)), 10))
plot(anscombeTrans_NB(as.vector(data), 10), anscombeTrans_NB(W%*%H, 10) - anscombeTrans_NB(as.vector(data),10))
plot(anscombeTrans_NB1(as.vector(data), 10), anscombeTrans_NB(W%*%H, 10) - anscombeTrans_NB(as.vector(data),10))

NBdata = rnbinom(100, size = 10, prob = runif(100))

m = lm(NBdata ~ 1)
plot(NBdata, residuals(m))
plot(anscombeTrans_NB(NBdata, 10))


################################################################################
#
#Anscombe Poisson transformation
#
################################################################################

#
#INITIAL TESTING USING GENERATED VALUES
#

anscombeTrans = function(x) 2*sqrt(x + 3/8)

m = 20
Pdata = rpois(5000, m)
plot(anscombeTrans(Pdata))

#sd theoretically = 1 + O(1/mean^2)
sd(anscombeTrans(Pdata))
qqPlot(anscombeTrans(Pdata))
shapiro.test(anscombeTrans(Pdata))

hist(anscombeTrans(Pdata), prob = T)
curve(dnorm(x,mean = 2*sqrt(m + 3/8)-1/(4*sqrt(m)), sd = 1), add = T, col = "red")

resi = (anscombeTrans(Pdata) - 2*sqrt(m + 3/8)-1/(4*sqrt(m)))
hist(resi)

#
#TESTING IN NMF SETTING
#

V_pois_obs <- matrix(rpois(21*96, lambda = W%*%H), nrow = 21)
NMF_pois <- nmf(V_pois_obs, rank = 6, nrun = 10, method = "KL")
W_pois = basis(NMF_pois)
H_pois = coef(NMF_pois)
################################################################################
#
#Anscombe Negative Binomial transformation
#
################################################################################

source("NMFNBMMsquarem.R")
inverse = function (f, lower = -100, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}
sinh_inverse = inverse(function(x) sinh(x), 0.1, 100)
ihs <- function(x) {
  y <- log(x + sqrt(x^2+1))
  return(y)
}

anscombeTrans_NB = function(x, alpha, c = 3/8) sqrt(x-1/2)*ihs(sqrt((alpha + c)/(x-c)))
anscombeTrans_NB_alt = function(x, alpha) log(alpha + x/2)

k = 10
NBdata = rnbinom(5000, size = k, 0.3)
hist(NBdata)
hist(anscombeTrans_NB(NBdata, k), prob = T, breaks = 50)
curve(dnorm(x, mean = mean(anscombeTrans_NB(NBdata, k)), sd = trigamma(k)), add = T, col = "red")


hist(anscombeTrans_NB_alt(NBdata, k), prob = T, breaks = 50)
curve(dnorm(x, mean = mean(anscombeTrans_NB_alt(NBdata, k)), sd = trigamma(k)), add = T, col = "red")


#
#
#

k=10
V_NB_obs <- matrix(rnbinom(21*96, size = k, prob = W%*%H/(k+W%*%H)), nrow = 21)
NMF_NB <- NMFNBMMsquarem(V_NB_obs, 6, k)
W_NB = NMF_NB$E
H_NB = NMF_NB$P

res = (anscombeTrans_NB_alt(as.vector(V_NB_obs), k) - anscombeTrans_NB_alt(as.vector(H_NB%*%W_NB), k))/trigamma(k)
hist(res, breaks = 50, prob = T)
curve(dnorm(x), add = T, col = "red")

res = (anscombeTrans_NB_alt(as.vector(V_NB_obs), k) - anscombeTrans_NB_alt(as.vector(H_NB%*%W_NB), k))/(1/(k-1/2))
hist(res, breaks = 50, prob = T)
curve(dnorm(x), add = T, col = "red")


res = (anscombeTrans_NB(as.vector(V_NB_obs), k) - anscombeTrans_NB(as.vector(H_NB%*%W_NB), k))/trigamma(k)
hist(res, breaks = 50, prob = T)
curve(dnorm(x), add = T, col = "red")


res = (anscombeTrans_NB(as.vector(V_NB_obs), k, c = 0.3) - anscombeTrans_NB(as.vector(H_NB%*%W_NB), k))/(1/(k-1/2))
hist(res, breaks = 50, prob = T)
curve(dnorm(x), add = T, col = "red")

