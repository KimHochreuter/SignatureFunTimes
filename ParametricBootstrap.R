library(NMF)
#Bootstrap analysis of distribution

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
SSA_pois_obs <- sum((A_pois_obs-1)^2)

A_NB_obs <- V_NB_obs/(W_NB%*%H_NB)
SSA_NB_obs <- sum((A_NB_obs-1)^2)

N_boot = 10000
SSA_pois <- rep(NULL, N_boot)
SSA_NB <- rep(NULL, N_boot)
for(n in 1:N_boot){
  V_pois <- rpois(21*96, lambda = as.vector(W_pois%*%H_pois))
  V_NB <- rnbinom(21*96, size = 50, prob = W_NB%*%H_NB/(50+W_NB%*%H_NB))
  A_pois <- V_pois/(W_pois%*%H_pois)
  A_NB <- V_NB/(W_NB%*%H_NB)
  SSA_pois[n] <- sum((A_pois-1)^2)
  SSA_NB[n] <- sum((A_NB-1)^2)
}

#Tester om hhv pois og nb ligner fordelingen af de simulerede poisson
hist(SSA_pois)
abline(v=SSA_pois_obs, col="red")
mean(SSA_pois>SSA_pois_obs)


hist(SSA_pois)
abline(v=SSA_NB_obs, col="red")
mean(SSA_pois>SSA_NB_obs)
