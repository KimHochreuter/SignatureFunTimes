D = rnbinom(48000, 100, 0.9)
D = matrix(D, ncol = 96)

NMF_NB = NMFNBMMsquarem(D, 4, alpha = 500)
alpha_NB = t(NMF_NB$E)
beta_NB = t(NMF_NB$P)
log_lik_NB <- function(alpha){
  WH <- t(alpha_NB%*%beta_NB)
  #M <- round(M_CV_NB)
  return(-sum(dnbinom(x = D, size = alpha, prob = WH/(alpha + WH), log = T)))
}
alpha <- optimize(log_lik_NB, interval = c(0,100))$minimum


load("BRCA21.RData")
V = t(V)
NMF_NB = NMFNBMMsquarem(V, 14, alpha = 100)
alpha_NB = t(NMF_NB$E)
beta_NB = t(NMF_NB$P)

D = rnbinom(2016, size = 100, prob = alpha_NB%*%beta_NB/(alpha_NB%*%beta_NB + 100))
D = matrix(D, ncol = 96)


log_lik_NB <- function(alpha){
  WH <- t(alpha_NB%*%beta_NB)
  #M <- round(M_CV_NB)
  return(-sum(dnbinom(x = D, size = alpha, prob = WH/(alpha + WH), log = T)))
}
alpha <- optimize(log_lik_NB, interval = c(0,5000))$minimum


D_CV <- CVNB(D, K=15)

D_CV <- as.data.frame(D_CV)
D_CV_ = D_CV[D_CV$n_update == 5,]
g = D_CV_ %>%
  group_by(K) %>% {.}
#summarise(medMSE = median(MSE), medBIC = median(BIC))
ggplot(g) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = MSE)) + ylim(c(0,1000000))
ggplot(g) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = alpha))
