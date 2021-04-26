library(NMF)
library(tidyverse)
library(stringr)

start_time <- Sys.time()

load("DATA/patients.rda")
M <- patients
#M = t(V)
M <- M[rowSums(M)>1000,]#Removing patients with less than 1000 mutations
#M = head(M, n = 50)
G <- nrow(M) #Number of patients
m <- ncol(M) #Number of mutation types

K_range <- 2:10 #Number of signatures to investigate

MSE <- list(0) #find better solution
MSE_CV <- matrix(0, ncol = 10, nrow = 5)

AIC_collect = list(0)
BIC_collect = list(0)
KL_collect = list(0)

AIC_pois = matrix(0, ncol = 10, nrow = 5)
BIC_pois = matrix(0, ncol = 10, nrow = 5)
KL_pois = matrix(0, ncol = 10, nrow = 5)

CV_idx <- replicate(10, list(as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*0.01),replace=FALSE),])))

for (k in 2:max(K_range)){
    #MSE_CV <- 0
    for (i in 1:10){
      M_CV <- M
      M_CV[CV_idx[[i]]] <- 0 #Replace approximately 1% of cells in matrix by 0.
      for (s in 1:5){
        NMF <- nmf(M_CV, rank = k, nrun = 10, method = "KL")
        alpha = t(coef(NMF))
        beta = t(basis(NMF))
        #we replace the zero cells with the estimate from the CV. This is run 5 times, with updated M
        M_CV[CV_idx[[i]]] <- (t(alpha%*%beta))[CV_idx[[i]]]
        MSE_CV[s,i] <- mean((as.vector(M[CV_idx[[i]]])-(t(alpha%*%beta))[CV_idx[[i]]])^2)
        AIC_pois[s,i] = -2*sum(as.vector(M)*log(as.vector(alpha%*%beta))-as.vector(alpha%*%beta)) + 2*(k + ncol(M) + nrow(M)) 
        BIC_pois[s,i] = -2*sum(as.vector(M)*log(as.vector(alpha%*%beta))-as.vector(alpha%*%beta)) + 2*(k + ncol(M) + nrow(M))*log(ncol(M) * nrow(M)) 
        KL_pois[s,i] = sum(as.vector(M)*log(as.vector(as.vector(M)/as.vector(alpha%*%beta))) - as.vector(M) + as.vector(alpha%*%beta))
        #print(MSE_CV)
        
      }
    }
  AIC_collect[[k]] = AIC_pois
  BIC_collect[[k]] = BIC_pois
  KL_collect[[k]] = KL_pois
  MSE[[k]] = MSE_CV
  #print(MSE)
  message(100*k/(max(K_range)), "%")
}
plot(0,0,xlim = c(min(K_range), max(K_range)), ylim = c(0, 10000))
for(k in K_range){
  points(rep(k, 10), MSE[[k]][5,], xlim = c(min(K_range), max(K_range)))
  points(k, median(MSE[[k]][5,]), col = 2, pch = 16)
}

end_time <- Sys.time()
end_time - start_time


##############################################
#SIGNATURE PLOTS
##############################################

NMF_final = nmf(M, rank = 6, nrun = 10)

W = coef(NMF_final)
W_df = data.frame(t(W))
colnames(W_df) = c("s1", "s2", "s3", "s4", "s5", "s6")
W_df$MutationType = rownames(W_df)
W_df$muta2 = str_sub(W_df$MutationType, 3, -3)


( ggplot(W_df) 
    + geom_col(aes(x = MutationType, y = s1, fill = muta2))
    #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
    + theme_bw() 
    + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
    + facet_grid(~muta2,scales="free_x") )

( ggplot(W_df) 
  + geom_col(aes(x = MutationType, y = s2, fill = muta2))
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
  + facet_grid(~muta2,scales="free_x") )



################################################################################
##
##  CV FOR POISSON
##
################################################################################
source("NMFNBMMsquarem.R")

start_time <- Sys.time()

load("DATA/patients.rda")
M <- patients
#M = t(V)
M <- M[rowSums(M)>1000,]#Removing patients with less than 1000 mutations
#M = head(M, n = 50)
G <- nrow(M) #Number of patients
m <- ncol(M) #Number of mutation types

K_range <- 2:10 #Number of signatures to investigate

MSE_NB <- list(0) #find better solution
MSE_NB2 <- list(0)
MSE_CV_NB <- matrix(0, ncol = 10, nrow = 5)

MSE_po <- list(0) #find better solution
MSE_CV_po <- matrix(0, ncol = 10, nrow = 5)

CV_idx <- replicate(10, list(as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*0.01),replace=FALSE),])))

for (k in 2:max(K_range)){
  #MSE_CV <- 0
  for (i in 1:10){
    M_CV_po <- M
    M_CV_po[CV_idx[[i]]] <- 0 #Replace approximately 1% of cells in matrix by 0.
    for (s in 1:5){
      NMF_po <- nmf(M_CV_po, rank = k, method = "KL")
      alpha_po = t(coef(NMF_po))
      beta_po = t(basis(NMF_po))
      
      #we replace the zero cells with the estimate from the CV. This is run 5 times, with updated M
      M_CV_po[CV_idx[[i]]] <- (t(alpha_po%*%beta_po))[CV_idx[[i]]]
      
      MSE_CV_po[s,i] <- mean((as.vector(M[CV_idx[[i]]])-(t(alpha_po%*%beta_po))[CV_idx[[i]]])^2)
    }
  }
  MSE_po[[k]] = MSE_CV_po
  message(round(100*(k-1)/(max(K_range)-1), digits = 2), "%")
}

end_time <- Sys.time()
end_time - start_time


################################################################################
##
##  CV FOR POISSON PT. II (WITHOUT LISTS)
##
################################################################################
source("NMFNBMMsquarem.R")

start_time <- Sys.time()

load("DATA/patients.rda")
M <- patients
#M = t(V)
M <- M[rowSums(M)>1000,]#Removing patients with less than 1000 mutations
#M = head(M, n = 50)
G <- nrow(M) #Number of patients
m <- ncol(M) #Number of mutation types

K_range <- 2:15 #Number of signatures to investigate

MSE_CV_po <- NA

CV_idx <- replicate(10, list(as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*0.01),replace=FALSE),])))

for (k in 2:max(K_range)){
  #MSE_CV <- 0
  for (i in 1:10){
    M_CV_po <- M
    M_CV_po[CV_idx[[i]]] <- 0 #Replace approximately 1% of cells in matrix by 0.
    for (s in 1:5){
      NMF_po <- nmf(M_CV_po, rank = k, method = "KL")
      alpha_po = t(coef(NMF_po))
      beta_po = t(basis(NMF_po))
      
      #we replace the zero cells with the estimate from the CV. This is run 5 times, with updated M
      M_CV_po[CV_idx[[i]]] <- (t(alpha_po%*%beta_po))[CV_idx[[i]]]
      
      if (is.na(MSE_CV_po)){
        MSE_CV_po <- c(k, i, s, mean((as.vector(M[CV_idx[[i]]])-(t(alpha_po%*%beta_po))[CV_idx[[i]]])^2))
        
      }
      else{
        MSE_CV_po <- rbind(MSE_CV_po,c(k, i, s, mean((as.vector(M[CV_idx[[i]]])-(t(alpha_po%*%beta_po))[CV_idx[[i]]])^2)))
      }
    }
  }
  message(round(100*(k-1)/(max(K_range)-1), digits = 2), "%")
}

end_time <- Sys.time()
end_time - start_time

colnames(MSE_CV_po) <- c("K", "CV_set", "n_update", "MSE")

################################################################################
##
##  CV FOR NB (alpha as CV param)
##
################################################################################

start_time <- Sys.time()

K_range <- 2:10 #Number of signatures to investigate
alpha =  c(1,10,20,50,100)#seq(1,100, length.out = 5)


#MSE_NB <- list(0) #find better solution
#MSE_NB2 <- list(0)
#MSE_CV_NB <- matrix(0, ncol = 10, nrow = 5) # the five runs for the 10 CV sets
MSE_CV_NB <- NA

CV_idx <- replicate(10, list(as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*0.01),replace=FALSE),])))

for (k in 2:max(K_range)){
  for (z in 1:length(alpha)){
    for (i in 1:10){
    M_CV_NB <- M
    M_CV_NB[CV_idx[[i]]] <- 0 #Replace approximately 1% of cells in matrix by 0.
      for (s in 1:5) {
        NMF_NB = NMFNBMMsquarem(M_CV_NB, k, alpha = alpha[z])
        alpha_NB = t(NMF_NB$E)
        beta_NB = t(NMF_NB$P)
        #MSE_CV_NB[s,z] = mean((as.vector(M[CV_idx[[i]]])-(t(alpha_NB%*%beta_NB))[CV_idx[[i]]])^2)
        MSE <-  mean((as.vector(M[CV_idx[[i]]])-(t(alpha_NB%*%beta_NB))[CV_idx[[i]]])^2)
        if (is.na(MSE_CV_NB)){
          MSE_CV_NB <- c(k, alpha[z], i, s,MSE)
        }
        else{
          MSE_CV_NB <- rbind(MSE_CV_NB, c(k, alpha[z], i, s,MSE))
        }
      }
      
      #we replace the zero cells with the estimate from the CV. This is run 5 times, with updated M
      M_CV_NB[CV_idx[[i]]] <- (t(alpha_NB%*%beta_NB))[CV_idx[[i]]]
      
      
    }
    #MSE_NB2[[i]] = MSE_CV_NB
  }
  #MSE_NB[[k]] = MSE_NB2
  message(round(100*(k-1)/(max(K_range)-1), digits = 2), "%")
}

end_time <- Sys.time()
end_time - start_time
colnames(MSE_CV_NB) <- c("K", "alpha", "CV_set", "n_update", "MSE")

################################################################################
##
##  CV FOR NB (alpha updated in run)
##
################################################################################

start_time <- Sys.time()

K_range <- 2:10 #Number of signatures to investigate
start_alpha <- 10

log_lik_NB <- function(alp, x, wh){
  x <- as.vector(x)
  wh <- as.vector(wh)
  print(c(gamma(alp+x[91]),gamma(alp+x[93]), gamma(alp+x[111]), gamma(alp+x[121])))
  return(log(gamma(x+alp)) - log(gamma(alp))) + x*log(wh/(alp+wh)) + alp*log(1-wh/(wh+alp))
}

MSE_CV_NB <- NA
CV_idx <- replicate(10, list(as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*0.01),replace=FALSE),])))

for (k in 2:max(K_range)){
  for (i in 1:10){
    alpha <- start_alpha
    M_CV_NB <- M
    M_CV_NB[CV_idx[[i]]] <- 0 #Replace approximately 1% of cells in matrix by 0.
    for (s in 1:5) {
      NMF_NB = NMFNBMMsquarem(M_CV_NB, k, alpha = alpha)
      alpha_NB = t(NMF_NB$E)
      beta_NB = t(NMF_NB$P)
      alpha <- optim(par = alpha, log_lik_NB(alp = 10, x = M_CV_NB, wh = alpha_NB%*%beta_NB))
      #MSE_CV_NB[s,z] = mean((as.vector(M[CV_idx[[i]]])-(t(alpha_NB%*%beta_NB))[CV_idx[[i]]])^2)
      MSE <-  mean((as.vector(M[CV_idx[[i]]])-(t(alpha_NB%*%beta_NB))[CV_idx[[i]]])^2)
      if (is.na(MSE_CV_NB)){
        MSE_CV_NB <- c(k, alpha[z], i, s,MSE)
      }
      else{
        MSE_CV_NB <- rbind(MSE_CV_NB, c(k, alpha[z], i, s,MSE))
      }
    }
    
    #we replace the zero cells with the estimate from the CV. This is run 5 times, with updated M
    M_CV_NB[CV_idx[[i]]] <- (t(alpha_NB%*%beta_NB))[CV_idx[[i]]]
    
  }
  message(round(100*(k-1)/(max(K_range)-1), digits = 2), "%")
}

end_time <- Sys.time()
end_time - start_time
colnames(MSE_CV_NB) <- c("K", "alpha", "CV_set", "n_update", "MSE")

################################################################################
##
##  Plots
##
################################################################################

plot(0,0,xlim = c(min(K_range), max(K_range)), ylim = c(0, 10000))
for(k in K_range){
  points(rep(k, 10), MSE[[k]][1,], xlim = c(min(K_range), max(K_range)))
  points(k, median(MSE[[k]][1,]), col = 2, pch = 16)
}
abline(h = 2000)



plot(0,0,xlim = c(min(K_range), max(K_range)), ylim = c(0, 10000))
for(k in K_range){
  points(rep(k, 10), MSE_po[[k]][5,], xlim = c(min(K_range), max(K_range)))
  points(k, median(MSE_po[[k]][5,]), col = 2, pch = 16)
}
abline(h = 2000)
abline(h = 1000)

plot(0,0,xlim = c(min(K_range), max(K_range)), ylim = c(0, 10000))
for (j in 1:length(alpha)) {
  for(k in K_range){
    #points(rep(k, 10), MSE_NB[[k]][[j]][5,], xlim = c(min(K_range), max(K_range)))
    points(k, median(MSE_NB[[k]][[j]][5,]), col = j, pch = 16, )
  }
}
abline(h = 2000)
abline(h = 1000)

plot(0,0,xlim = c(1, 100), ylim = c(0, 10000))
for (j in seq(1,100, length.out = 10)) {
  for(k in K_range){
    #points(rep(k, 10), MSE_NB[[k]][[j]][5,], xlim = c(min(K_range), max(K_range)))
    points(j, median(MSE_NB[[k]][[j]][5,]), col = j, pch = 16, )
  }
}
abline(h = 1000)

plot(0,0,xlim = c(1, 100), ylim = c(0, 20))
for (j in 1:length(alpha)) {
  for(k in K_range){
    #points(rep(k, 10), MSE_NB[[k]][[j]][5,], xlim = c(min(K_range), max(K_range)))
    points(alpha[j], k, cex = MSE_NB[[k]][[j]][5,]/10000,  col = j, pch = 16)
  }
}
abline(h = 1000)

library(ggplot2)
library(tidyverse)


#plot cv for nb

par(mfrow=c(2,1))

plot_CV_NB <- data.frame(subset(MSE_CV_NB, MSE_CV_NB[,4] == 5)) %>%
  mutate(K = K-1)%>%
  mutate(K = replace_na(K,10))%>%
  group_by(K, alpha) %>%
  summarise(medMSE = median(MSE))
ggplot(data = plot_CV_NB, aes(x = K, y = medMSE, colour = factor(alpha))) + 
  geom_point()+geom_line(aes(group = alpha))
ggplot(data = plot_CV_NB, aes(x = factor(alpha), y = medMSE, colour = factor(K))) + 
  geom_point()+geom_line(aes(group = K))

#plot cv for pois
plot_CV_po <- data.frame(subset(MSE_CV_po, MSE_CV_po[,3] == 5)) %>%
  group_by(K) %>%
  summarise(medMSE = median(MSE))
ggplot(data = plot_CV_po, aes(x = K, y = medMSE)) + geom_point()+geom_line()


fill <- "skyblue2"
line <- "steelblue"
plot_CV_po_BOX <- data.frame(subset(MSE_CV_po, MSE_CV_po[,3] == 5)) 
ggplot(data = plot_CV_po_BOX, aes(x = factor(K), y = MSE, group = K)) + 
  geom_boxplot(fill = fill, alpha = 0.7, outlier.shape = 20) + 
  scale_x_discrete(name = "Number of mutational signatures") +
  scale_y_continuous(name = "Mean squared error", limits=c(0, 0.75*10^5))+ theme_bw()

       