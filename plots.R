source("CV_function_version.R")

load("BRCA21.RData")
load("data/patients.rda")

porun = CVPO(patients, K = 20)
nbrun = CVNB(patients, K = 20)
porun1 = CVPO(V, K = 20)
nbrun1 = CVNB(V, K = 20)


################################################################################
##
##  POISSON MSE/BIC PLOT
##
################################################################################
poruna = data.frame(porun)
#poruna = data.frame(porun1)
po_plot_df = poruna[poruna$n_update == 5,]

d = po_plot_df %>%
  group_by(K) %>% {.}
  #summarise(medMSE = median(MSE), medBIC = median(BIC))
ggplot(d) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = MSE)) + ylim(c(0,50000))
ggplot(d) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = BIC))



################################################################################
##
##  NEGATIVE BINOMIAL MSE/BIC PLOT
##
################################################################################
nbruna = data.frame(nbrun)
#nbruna = data.frame(nbruna1)
nb_plot_df = nbruna[nbruna$n_update == 5,]
g = nb_plot_df %>%
  group_by(K) %>% {.}
  #summarise(medMSE = median(MSE), medBIC = median(BIC))
ggplot(g) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = MSE)) + ylim(c(0,10000))
ggplot(g) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = BIC))



################################################################################
##
##  RESIDUAL PLOTS
##
################################################################################
data = t(V)
poNMF = nmf(data, rank = 6, nrun = 10, method = "KL")
H_po = basis(poNMF)
W_po = coef(poNMF)
alpha = 56
nbNMF = NMFNBMMsquarem(data, 6, alpha)
H_nb = nbNMF$P
W_nb = nbNMF$E

diff_po = data - H_po%*%W_po
diff_nb = data - H_nb%*%W_nb
z = data.frame(data = as.vector(data), diff_po = as.vector(diff_po), as.vector(diff_nb))
a = (ggplot(z) + geom_point(aes(data, diff_po)) 
          + xlim(c(0,500))
          + ylim(c(-200,200))
          + geom_function(fun = function(x) 2*sqrt(x), aes(colour = "Poisson"))
          + geom_function(fun = function(x) -2*sqrt(x), aes(colour = "Poisson"))
          + geom_function(fun = function(x) 2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial"))
          + geom_function(fun = function(x) -2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial")))
a
b = (ggplot(z) + geom_point(aes(data, diff_nb)) 
  + xlim(c(0,500)) 
  + ylim(c(-200,200))
  + geom_function(fun = function(x) 2*sqrt(x), aes(colour = "Poisson"))
  + geom_function(fun = function(x) -2*sqrt(x), aes(colour = "Poisson"))
  + geom_function(fun = function(x) 2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial"))
  + geom_function(fun = function(x) -2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial")))
b



################################################################################
##
##  WHEN IS ALPHA = INFTY?
##
################################################################################
nbr = rnbinom(1000, 100, 0.9)
hist(nbr, prob = T)
points(dpois(0:1200, lambda = mean(nbr)), pch = 16)

por = rpois(1000, lambda = mean(nbr))
hist(por, prob = T)
points(dpois(0:300, lambda = mean(nbr)), pch = 16)
