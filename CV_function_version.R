source("NMFNBMMsquarem.R")
library(NMF)
library(tidyverse)
library(stringr)
CVNB = function(M, K = 10, start_alpha = 10, n_mutatypes = 96){
  start_time <- Sys.time()
  if (dim(M)[2] != n_mutatypes) {
    M = t(M)
  }
  G <- nrow(M) #Number of patients
  m <- ncol(M) #Number of mutation types
  K_range <- 2:K #Number of signatures to investigate
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
        #log_lik_NB <- function(alpha){
        #  WH <- t(alpha_NB%*%beta_NB)
        #  M <- round(M_CV_NB)
        #  return(-sum(dnbinom(x = M, size = alpha, prob = WH/(alpha + WH), log = T)))
        #}
        log_lik_NB <- function(alpha){
          WH <- as.vector(t(alpha_NB%*%beta_NB))
          M <- as.vector(M_CV_NB)
          return(-sum(lgamma(M + alpha) - lgamma(alpha) 
                      + M*(log(WH) - log(WH + alpha)) + alpha*log(1-WH/(WH+alpha))))
        }
        alpha <- optimize(log_lik_NB, interval = c(0,5000))$minimum
        MSE <-  mean((as.vector(M[CV_idx[[i]]])-(t(alpha_NB%*%beta_NB))[CV_idx[[i]]])^2)
        params <- (dim(alpha_NB)[1]*dim(alpha_NB)[2] + dim(beta_NB)[1]*dim(beta_NB)[2] + 1)
        nobs <- ncol(M)*nrow(M)
        BIC <- -2*log_lik_NB(alpha) + params*log(nobs)
        if (is.na(MSE_CV_NB)){
          MSE_CV_NB <- c(k, alpha, i, s, MSE, BIC)
        }
        else{
          MSE_CV_NB <- rbind(MSE_CV_NB, c(k, alpha, i, s, MSE, BIC))
        }
        M_CV_NB[CV_idx[[i]]] <- (t(alpha_NB%*%beta_NB))[CV_idx[[i]]]
      }
    }
    message(round(100*(k-1)/(max(K_range)-1), digits = 2), "%")
  }
  end_time <- Sys.time()
  colnames(MSE_CV_NB) <- c("K", "alpha", "CV_set", "n_update", "MSE", "BIC")
  print(end_time - start_time)
  return(MSE_CV_NB)
}
CVPO = function(M, K = 10, n_mutatypes = 96){
  start_time <- Sys.time()
  if (dim(M)[2] != n_mutatypes) {
    M = t(M)
  }
  G <- nrow(M) #Number of patients
  m <- ncol(M) #Number of mutation types
  K_range <- 2:K #Number of signatures to investigate
  MSE_CV_po <- NA
  CV_idx <- replicate(10, list(as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*0.01),replace=FALSE),])))
  for (k in 2:max(K_range)){
    for (i in 1:10){
      M_CV_po <- M
      M_CV_po[CV_idx[[i]]] <- 0 #Replace approximately 1% of cells in matrix by 0.
      for (s in 1:5){
        NMF_po <- nmf(M_CV_po, rank = k, method = "KL")
        alpha_po = t(coef(NMF_po))
        beta_po = t(basis(NMF_po))
        #we replace the zero cells with the estimate from the CV. This is run 5 times, with updated M
        M_CV_po[CV_idx[[i]]] <- (t(alpha_po%*%beta_po))[CV_idx[[i]]]
        params <- (dim(alpha_po)[1]*dim(alpha_po)[2] + dim(beta_po)[1]*dim(beta_po)[2])
        nobs <- ncol(M)*nrow(M)
        log_lik_M <- sum(as.vector(M_CV_po)*log(as.vector(t(alpha_po%*%beta_po)))-as.vector(t(alpha_po%*%beta_po)))
        BIC <- -2*log_lik_M + params*log(nobs) 
        if (is.na(MSE_CV_po)){
          MSE_CV_po <- c(k, i, s, mean((as.vector(M[CV_idx[[i]]])-(t(alpha_po%*%beta_po))[CV_idx[[i]]])^2), BIC)
        }
        else{
          MSE_CV_po <- rbind(MSE_CV_po,c(k, i, s, mean((as.vector(M[CV_idx[[i]]])-(t(alpha_po%*%beta_po))[CV_idx[[i]]])^2), BIC))
        }
      }
    }
    message(round(100*(k-1)/(max(K_range)-1), digits = 2), "%")
  }
  end_time <- Sys.time()
  colnames(MSE_CV_po) <- c("K", "CV_set", "n_update", "MSE", "BIC")
  print(end_time - start_time)
  return(MSE_CV_po)
}