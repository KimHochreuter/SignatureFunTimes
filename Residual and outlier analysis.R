library(NMF)

load("DATA/patients.rda")
M <- patients
#M = t(V)
M <- M[rowSums(M)>1000,]

######################################################
##
##RESIDUAL PLOT, POISSON
##
######################################################

NMF_final = nmf(M, rank = 6, nrun = 10, method = "KL")

H_po = basis(NMF_final)
W_po = coef(NMF_final)

diff = M - H_po%*%W_po
plot(M, diff, xlim = c(0,1000), ylim=c(-500, 500))
plot(M, diff, xlim = c(0,10000), ylim=c(-1000, 1000))
plot(M, diff)

plot(M, diff, xlim = c(0,1000), ylim=c(-500, 500))
plot(M, diff)
curve(5*sqrt(x), add = T, col = "red")
curve(-5*sqrt(x), add = T, col = "red")
curve(2*sqrt(x + (x^2)/10), add = T, col = "green")
curve(-2*(sqrt(x + (x^2)/10)), add = T, col = "green")

A_po = M/(H_po%*%W_po) #Color patients, mutation type

hist(A_po, breaks = 50, probability = T)
curve(dgamma(x,shape=100,rate = 100), add = T)
mean(A_po)
var(as.numeric(A_po))


######################################################
##
##RESIDUAL PLOT, POISSON
##
######################################################
source("NMFNBMMsquarem.R")

NBNMF = NMFNBMMsquarem(M, 6, 10)
H = NBNMF$P
W = NBNMF$E

diff_nb = M - H%*%W
plot(M,diff_nb)
curve(5*sqrt(x), add = T, col = "red")
curve(-5*sqrt(x), add = T, col = "red")
curve(2*sqrt(x + (x^2)/10), add = T, col = "green")
curve(-2*(sqrt(x + (x^2)/10)), add = T, col = "green")

par(mfrow=c(1,2))
plot(0)
