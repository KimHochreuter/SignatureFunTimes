load("DATA/patients.rda")
M <- patients
#M = t(V)
M <- M[rowSums(M)>1000,]

######################################################
#RESIDUAL PLOT
######################################################

NMF_final = nmf(M, rank = 6, nrun = 10, method = "KL")

H = basis(NMF_final)
W = coef(NMF_final)

diff = M - H%*%W
plot(M, diff, xlim = c(0,1000), ylim=c(-500, 500))
plot(M, diff, xlim = c(0,10000), ylim=c(-1000, 1000))
plot(M, diff)

plot(M, diff, xlim = c(0,1000), ylim=c(-500, 500))
plot(M, diff)
curve(5*sqrt(x), add = T, col = "red")
curve(-5*sqrt(x), add = T, col = "red")
curve(2*sqrt(x + (x^2)/10), add = T, col = "green")
curve(-2*(sqrt(x + (x^2)/10)), add = T, col = "green")

A = M/(H%*%W) #Color patients, mutation type

hist(A, breaks = 50, probability = T)
curve(dgamma(x,shape=100,rate = 100), add = T)
mean(A)
var(as.numeric(A))




