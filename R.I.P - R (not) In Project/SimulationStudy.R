library(NMF)
load("BRCA21.RData")

###############################################################################
#
#    POISSON
#
###############################################################################


D = rpois(dim(V)[1]*dim(V)[2], lambda = V)

pois_generated_data = matrix(D, ncol = 21, nrow = 96)
colnames(pois_generated_data) = colnames(V)
rownames(pois_generated_data) = rownames(V)

NMF_ = nmf(pois_generated_data, rank = 6, method = "KL", nrun = 10)

W = t(coef(NMF_))
H = t(basis(NMF_))

diff = V - t(H)%*%t(W)
plot(V, diff, xlim = c(0,1000), ylim=c(-500, 500))
plot(V, diff, xlim = c(0,10000), ylim=c(-1000, 1000))
plot(V, diff)


A = V/(t(H)%*%t(W))

hist(A, breaks = 50, probability = T)
curve(dgamma(x, shape=71,rate = 71), add = T)
mean(A)
var(as.numeric(A))


Nboot = 1000
res = numeric(Nboot)
res2 = numeric(Nboot)

for (i in 1:Nboot) {
  D = rpois(dim(V)[1]*dim(V)[2], lambda = as.vector(V))
  pois_generated_data = matrix(D, ncol = 21, nrow = 96)
  colnames(pois_generated_data) = colnames(V)
  rownames(pois_generated_data) = rownames(V)
  
  NMF_ = nmf(pois_generated_data, rank = 6, method = "KL")
  
  W = t(coef(NMF_))
  H = t(basis(NMF_))
  
  A = V/(t(H)%*%t(W))
  res2[i] = sum(A-1)^2
  res[i] = mean(A)
  if (i %% 100 == 0) {
    message(i/Nboot)
  }
}
hist(res[res<=2], breaks = 50, xlim = c(0,2))
res

###############################################################################
#
#    NEGATIVE BINOMIAL
#
###############################################################################

P = matrix(data = rep(0,21*96),ncol = 21, nrow = 96)

for (i in 1:dim(V)[1]) {
  for (j in 1:dim(V)[2]) {
    P[i,j] = rnbinom(n = 1, size = 50, prob = V[i,j]/(50 + V[i,j]))
  }
}

P

NB_generated_data = P
colnames(NB_generated_data) = colnames(V)
rownames(NB_generated_data) = rownames(V)

NMF_ = nmf(NB_generated_data, rank = 6, method = "KL")

W = t(coef(NMF_))
H = t(basis(NMF_))

A = V/(t(H)%*%t(W))

hist(A, breaks = 50, probability = T)
curve(dgamma(x, shape=100,rate = 100), add = T)
mean(A)