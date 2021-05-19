library(NMF)
library(car)
source("NMFNBMMsquarem.R")


load("BRCA21.RData")
load("data/patients.rda")
#Simulate poisson and NB data from the BRCA21 dataset
#data <- t(V)
data = patients
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

m = lm(anscombeTrans_NB(NBdata, 10) ~ 1)
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
#plot(anscombeTrans(Pdata))

#sd theoretically = 1 + O(1/mean^2)
sd(anscombeTrans(Pdata), na.rm = T)
#qqPlot(anscombeTrans(Pdata))
#shapiro.test(anscombeTrans(Pdata))

hist(anscombeTrans(Pdata), prob = T, breaks = 25)
curve(dnorm(x,mean = 2*sqrt(m + 3/8)-1/(4*sqrt(m)), sd = 1), add = T, col = "red")

#
#MIXTURE OF LAMBDAS
#

Pdata = c(rpois(1000, 10),
          rpois(1000, 20),
          rpois(1000, 50),
          rpois(1000, 100),
          rpois(1000, 200)
          )
hist(anscombeTrans(Pdata), prob = T, breaks = 25)

#
#TESTING IN NMF SETTING
#

V_pois_obs <- matrix(rpois(21*96, lambda = W%*%H), nrow = 21)
NMF_pois <- nmf(V_pois_obs, rank = 6, nrun = 10, method = "KL")
W_pois = basis(NMF_pois)
H_pois = coef(NMF_pois)


resi = anscombeTrans(V_pois_obs) - anscombeTrans(W_pois%*%H_pois)
hist(resi, prob = T, breaks = 30)
curve(dnorm(x,mean = 0, sd = 1), add = T, col = "red")


#
#TESTING ON DATA
#

NMF_pois <- nmf(data, rank = 6, nrun = 10, method = "KL")
W_pois = basis(NMF_pois)
H_pois = coef(NMF_pois)


resi = anscombeTrans(data) - anscombeTrans(W_pois%*%H_pois)
hist(resi, prob = T, breaks = 30)
curve(dnorm(x,mean = 0, sd = 1), add = T, col = "red")
curve(dnorm(x, sd = sd(resi, na.rm = T)), add = T, col = "blue")
shapiro.test(resi)


################################################################################
#
#Anscombe Negative Binomial transformation
#
################################################################################

#
#INITIAL TESTING
#

anscombeTrans_NB = function(x, alpha, c = 3/8) sqrt(alpha-1/2)*asinh(sqrt((x + c)/(alpha-2*c)))
anscombeTrans_NB_alt = function(x, alpha) log(x + alpha/2)
anscombeTrans_NB2 = function(x, alpha, c = 3/8) 2*asinh(sqrt((x + c)/(alpha-2*c)))

k = 40
NBdata = rnbinom(5000, size = k, 0.3)
#hist(NBdata)
hist(anscombeTrans_NB(NBdata, k), prob = T, breaks = 50)
curve(dnorm(x, mean = mean(anscombeTrans_NB(NBdata, k)), sd = trigamma(k)), add = T, col = "red")
curve(dnorm(x, mean = mean(anscombeTrans_NB(NBdata, k)), sd = sd(anscombeTrans_NB(NBdata, k))), add = T, col = "blue")


hist(anscombeTrans_NB2(NBdata, k), prob = T, breaks = 50)
curve(dnorm(x, mean = mean(anscombeTrans_NB2(NBdata, k)), sd = trigamma(k)), add = T, col = "red")
curve(dnorm(x, mean = mean(anscombeTrans_NB2(NBdata, k)), sd = sd(anscombeTrans_NB2(NBdata, k))), add = T, col = "blue")


hist(anscombeTrans_NB_alt(NBdata, k), prob = T, breaks = 50)
curve(dnorm(x, mean = mean(anscombeTrans_NB_alt(NBdata, k)), sd = trigamma(k)), add = T, col = "red")
curve(dnorm(x, mean = mean(anscombeTrans_NB_alt(NBdata, k)), sd = sd(anscombeTrans_NB_alt(NBdata, k))), add = T, col = "blue")

#
#MIXTURE OF ALPHAS
#

NBdata = c(rnbinom(1000, size = 10, 0.3), 
           rnbinom(1000, size = 20, 0.3),
           rnbinom(1000, size = 50, 0.3),
           rnbinom(1000, size = 100, 0.3),
           rnbinom(1000, size = 200, 0.3)
           )

hist(NBdata)
par(mfrow=c(3,1))
hist(anscombeTrans_NB(NBdata, (10 + 20 + 50 + 100 + 200)/5), prob = T, breaks = 50)
hist(anscombeTrans_NB(NBdata, 10), prob = T, breaks = 50)
hist(anscombeTrans_NB(NBdata, 500), prob = T, breaks = 50)

hist(anscombeTrans_NB2(NBdata, (10 + 20 + 50 + 100 + 200)/5), prob = T, breaks = 50)
hist(anscombeTrans_NB2(NBdata, 10), prob = T, breaks = 50)
hist(anscombeTrans_NB2(NBdata, 500), prob = T, breaks = 50)

hist(anscombeTrans_NB_alt(NBdata, (10 + 20 + 50 + 100 + 200)/5), prob = T, breaks = 50)
hist(anscombeTrans_NB_alt(NBdata, 10), prob = T, breaks = 50)
hist(anscombeTrans_NB_alt(NBdata, 500), prob = T, breaks = 50)

U = runif(1000)
NBdata = c(rnbinom(1000, size = 10, U), 
           rnbinom(1000, size = 20, U),
           rnbinom(1000, size = 50, U),
           rnbinom(1000, size = 100, U),
           rnbinom(1000, size = 200, U)
)
hist(anscombeTrans_NB(NBdata, (10 + 20 + 50 + 100 + 200)/5), prob = T, breaks = 50)
hist(anscombeTrans_NB(NBdata, 10), prob = T, breaks = 50)
hist(anscombeTrans_NB(NBdata, 500), prob = T, breaks = 50)

hist(anscombeTrans_NB2(NBdata, (10 + 20 + 50 + 100 + 200)/5), prob = T, breaks = 50)
hist(anscombeTrans_NB2(NBdata, 10), prob = T, breaks = 50)
hist(anscombeTrans_NB2(NBdata, 500), prob = T, breaks = 50)

hist(anscombeTrans_NB_alt(NBdata, (10 + 20 + 50 + 100 + 200)/5), prob = T, breaks = 50)
hist(anscombeTrans_NB_alt(NBdata, 10), prob = T, breaks = 50)
hist(anscombeTrans_NB_alt(NBdata, 500), prob = T, breaks = 50)


#
#TESTING IN NMF SETTING
#


k=50
V_NB_obs <- matrix(rnbinom(21*96, size = k, prob = W%*%H/(k+W%*%H)), nrow = 21)
NMF_NB <- NMFNBMMsquarem(V_NB_obs, 6, k)
W_NB = NMF_NB$E
H_NB = NMF_NB$P


res = (anscombeTrans_NB_alt(as.vector(V_NB_obs), k) - anscombeTrans_NB_alt(as.vector(H_NB%*%W_NB), k))
hist(res, breaks = 50, prob = T)
curve(dnorm(x), add = T, col = "red")
curve(dnorm(x, sd = sd(res, na.rm = T)), add = T, col = "blue")

#res = (anscombeTrans_NB_alt(as.vector(V_NB_obs), k) - anscombeTrans_NB_alt(as.vector(H_NB%*%W_NB), k))/(1/(k-1/2))
#hist(res, breaks = 50, prob = T)
#curve(dnorm(x), add = T, col = "red")
#curve(dnorm(x, sd = sd(res, na.rm = T)), add = T, col = "blue")


res = (anscombeTrans_NB(as.vector(V_NB_obs), k) - anscombeTrans_NB(as.vector(H_NB%*%W_NB), k))
hist(res, breaks = 50, prob = T)
curve(dnorm(x), add = T, col = "red")
curve(dnorm(x, sd = sd(res, na.rm = T)), add = T, col = "blue")


#res = (anscombeTrans_NB(as.vector(V_NB_obs), k, c = 0.3) - anscombeTrans_NB(as.vector(H_NB%*%W_NB), k))/(1/(k-1/2))
#hist(res, breaks = 50, prob = T)
#curve(dnorm(x), add = T, col = "red")
#curve(dnorm(x, sd = sd(res, na.rm = T)), add = T, col = "blue")


res = (anscombeTrans_NB2(as.vector(V_NB_obs), k) - anscombeTrans_NB2(as.vector(H_NB%*%W_NB), k))
hist(res, breaks = 50, prob = T)
curve(dnorm(x), add = T, col = "red")
curve(dnorm(x, sd = sd(res, na.rm = T)), add = T, col = "blue")


#
#TESTING ON DATA
#


k=500
NMF_NB <- NMFNBMMsquarem(data, 6, k)
W_NB = NMF_NB$E
H_NB = NMF_NB$P

WH = H_NB%*%W_NB
alpha = k
k <- optimize(function(alpha) -sum(dnbinom(x = data, size = alpha, prob = WH/(alpha + WH), log = T)), interval = c(0,100))$minimum

NMF_NB <- NMFNBMMsquarem(data, 6, k)
W_NB = NMF_NB$E
H_NB = NMF_NB$P

res = (anscombeTrans_NB_alt(as.vector(data), k) - anscombeTrans_NB_alt(as.vector(H_NB%*%W_NB), k))
hist(res, breaks = 50, prob = T)
curve(dnorm(x), add = T, col = "red")
curve(dnorm(x, sd = sd(res, na.rm = T)), add = T, col = "blue")
shapiro.test(res)

res = (anscombeTrans_NB(as.vector(data), k) - anscombeTrans_NB(as.vector(H_NB%*%W_NB), k))
hist(res, breaks = 50, prob = T)
curve(dnorm(x), add = T, col = "red")
curve(dnorm(x, sd = sd(res, na.rm = T)), add = T, col = "blue")

res = (anscombeTrans_NB2(as.vector(data), k) - anscombeTrans_NB2(as.vector(H_NB%*%W_NB), k))
hist(res, breaks = 50, prob = T)
curve(dnorm(x), add = T, col = "red")
curve(dnorm(x, sd = sd(res, na.rm = T)), add = T, col = "blue")

mean(data/(H_NB%*%W_NB))
sd(data/(H_NB%*%W_NB))




################################################################################
#
#Profile likelihood
#
################################################################################


alpha = c(1,25,50,75,100,125,150,175,200,250,300,500)
profileDF = NA
for (k in 2:10) {
  for (i in 1:length(alpha)) {
    model = NMFNBMMsquarem(V, k, alpha = alpha[i])
    W_NB = model$E
    H_NB = model$P
    WH = H_NB%*%W_NB
    loglik = sum(dnbinom(x = V, size = alpha[i], prob = WH/(alpha[i] + WH), log = T))
    if (is.na(profileDF)){
      profileDF <- c(k, alpha[i], loglik)
    }
    
    else{
      profileDF <- rbind(profileDF, c(k, alpha[i], loglik))
    }
  }
}
colnames(profileDF) <- c("K", "alpha", "loglik")

profileDF = profileDF %>% data.frame()
p = ggplot(profileDF) + geom_line(aes(x = alpha, y = loglik, colour = factor(K)))
#abline(v = alpha[which.max(loglik)])


q = profileDF %>% group_by(K) %>% summarize(value = max(loglik))
profileDF %>% group_by(alpha) %>% summarize(value = max(loglik)) %>% arrange(desc(value))
p + geom_hline(yintercept = q$value)


profileDFfinal = NA
for (i in 2:10) {
  index = which.max(profileDF[profileDF$K == i,"loglik"])
  bestAlpha = profileDF[profileDF$K == i,"alpha"][index]
  value = profileDF[profileDF$K == i,"loglik"][index]
  if (is.na(profileDF)){
    profileDFfinal <- c(i, bestAlpha, value)
  }
  
  else{
    profileDFfinal <- rbind(profileDFfinal, c(i, bestAlpha, value))
  }
}
colnames(profileDFfinal) = c("K", "Alpha", "loglik")
profileDFfinal %>% data.frame() %>% arrange(desc(loglik))



max(profileDF[profileDF$K == 7, "loglik"])
