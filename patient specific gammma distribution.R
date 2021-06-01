load("~/Uni/Speciale/SignatureFunTimes/Data/BRCA21.RData")
source("NMFNBMMsquarem.R")
library(ggplot2)
library(ggpubr)
library(ggforce)
library(lsa)

NMF_pois <- nmf(V, rank = 4, nrun = 10, method = "KL")
W_pois = basis(NMF_pois)
H_pois = coef(NMF_pois)

A_pois <- V/W_pois%*%H_pois
hist(A_pois[,1], breaks = 15)
plot(density(A_pois[,2]))
for (i in 2:ncol(A_pois)){

}

par(mfrow = c(1,5))

hist(A_pois[,1], breaks = 17, probability = T)
param <- mean(A_pois[,1])/var(A_pois[,1])
curve(dgamma(x,shape=param,rate = param), add = T)

hist(A_pois[,5], breaks = 17, probability = T)
param <- mean(A_pois[,5])/var(A_pois[,5])
curve(dgamma(x,shape=param,rate = param), add = T)

hist(A_pois[,12], breaks = 17, probability = T)
param <- mean(A_pois[,12])/var(A_pois[,12])
curve(dgamma(x,shape=param,rate = param), add = T)

hist(A_pois[,17], breaks = 17, probability = T)
param <- mean(A_pois[,17])/var(A_pois[,17])
curve(dgamma(x,shape=param,rate = param), add = T)


hist(A_pois[,21], breaks = 17, probability = T)
param <- mean(A_pois[,21])/var(A_pois[,21])
curve(dgamma(x,shape=param,rate = param), add = T)

hist1 <- ggplot(as.data.frame(A_pois[,21]), aes(x = A_pois[, 21])) +
  geom_histogram(aes(y = ..density..), bins = 17, fill = "grey", col = "black") + 
  theme_bw() + 
  xlab("") +
  stat_function(fun=dgamma, args=list(shape=params[5,2], rate=params[5,2]), col = "red")

hist2 <- ggplot(as.data.frame(A_pois[,1]), aes(x = A_pois[, 1])) +
  geom_histogram(aes(y = ..density..),bins = 17, fill = "grey", col = "black") + 
  theme_bw() + 
  xlab("") +
  stat_function(fun=dgamma, args=list(shape=params[1,2], rate=params[1,2]), col = "red")


param <- mean(A_pois[,21])/var(A_pois[,21])
curve(dgamma(x,shape=param,rate = param), add = T)

gammaplot_BRCA21 <- data.frame(A_pois[,c(1,5,12,17,21)])

param <- function(x){return(mean(x)/var(x))}
params <- data.frame(cbind(colnames(gammaplot_BRCA21), 
                           apply(gammaplot_BRCA21, MARGIN = 2,param), colSums(V[,c(1,5,12,17,21)])))

colnames(params) <- c("Patient", "param", "Count")
params[,2] <- as.numeric(params[,2])

gammaplot <- ggplot(data.frame(x = c(0 , 2.5)), aes(x = x)) + 
  stat_function(fun=dgamma, args=list(shape=params[1,2], rate=params[1,2]), aes(colour = paste(params[1,3], " (", round(params[1,2],2), ")")))+
  stat_function(fun=dgamma, args=list(shape=params[2,2], rate=params[2,2]), aes(colour = paste(params[2,3], " (", round(params[2,2],2), ")")))+
  stat_function(fun=dgamma, args=list(shape=params[3,2], rate=params[3,2]), aes(colour = paste(params[3,3], " (", round(params[3,2],2), ")")))+
  stat_function(fun=dgamma, args=list(shape=params[4,2], rate=params[4,2]), aes(colour = paste(params[4,3], " (", round(params[4,2],2), ")")))+
  stat_function(fun=dgamma, args=list(shape=params[5,2], rate=params[5,2]), aes(colour = paste(params[5,3], " (", round(params[5,2],2), ")")))+
  ylab("Gamma density") + xlab(" ") +
  labs(color=expression(paste("Total counts (", alpha, ")"))) +
  theme_bw()



hist_gamma <- (ggdraw() + draw_plot(hist1, x = 0, y = 0, width = 0.3, height = 0.5) 
          + draw_plot(hist2, x = 0, y = 0.5, width = 0.3, height = 0.45)
          + draw_plot(gammaplot, x = 0.3, y = 0, width = 0.7, height = 0.95)
          + draw_label("Fit of Gamma variable in Poisson-Gamma mixture ", x = 0.05, hjust = 0, y = 0.97, size = 13))

ggsave(plot = hist_gamma, file = "pictures/gammaplot.png", width = 200, height = 105.83332, units = "mm")

