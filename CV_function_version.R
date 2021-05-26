source("NMFNBMMsquarem.R")
library(NMF)
library(tidyverse)
library(stringr)
library(lsa)
CVNB = function(M, K = 10, start_alpha = 10, n_mutatypes = 96, n_cv_sets = 10, n_updates = 5, CVentries_percent = 0.01){
  start_time <- Sys.time()
  if (dim(M)[2] != n_mutatypes) {
    M = t(M)
  }
  G <- nrow(M) #Number of patients
  m <- ncol(M) #Number of mutation types
  K_range <- 2:K #Number of signatures to investigate
  MSE_CV_NB <- NA
  CV_idx <- replicate(n_cv_sets, list(as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*CVentries_percent),replace=FALSE),])))
  for (k in 2:max(K_range)){
    for (i in 1:n_cv_sets){
      alpha <- start_alpha
      M_CV_NB <- M
      M_CV_NB[CV_idx[[i]]] <- 0 #Replace approximately 1% of cells in matrix by 0.
      for (s in 1:n_updates) {
        NMF_NB = NMFNBMMsquarem(M_CV_NB, k, alpha = alpha, arrange = FALSE)
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
        alpha <- optimize(log_lik_NB, interval = c(0,100000))$minimum
        MSE <-  mean((as.vector(M[CV_idx[[i]]])-(t(alpha_NB%*%beta_NB))[CV_idx[[i]]])^2)
        params <- (dim(alpha_NB)[1]*dim(alpha_NB)[2] + dim(beta_NB)[1]*dim(beta_NB)[2] + 1)
        nobs <- ncol(M)*nrow(M)
        BIC <- 2*log_lik_NB(alpha) + params*log(nobs)
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
CVPO = function(M, K = 10, n_mutatypes = 96, n_cv_sets = 10, n_updates = 5, CVentries_percent = 0.01){
  start_time <- Sys.time()
  if (dim(M)[2] != n_mutatypes) {
    M = t(M)
  }
  G <- nrow(M) #Number of patients
  m <- ncol(M) #Number of mutation types
  K_range <- 2:K #Number of signatures to investigate
  MSE_CV_po <- NA
  CV_idx <- replicate(n_cv_sets, list(as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*CVentries_percent),replace=FALSE),])))
  for (k in 2:max(K_range)){
    for (i in 1:n_cv_sets){
      M_CV_po <- M
      M_CV_po[CV_idx[[i]]] <- 0 #Replace approximately 1% of cells in matrix by 0.
      for (s in 1:n_updates){
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

SignaturePairing = function(Nsig, H_NB, H_PO){
  NamingOrder = c()
  for (i in 1:Nsig) {
    CosineSim = cosine(cbind(H_PO,H_NB[,i]))
    CosineSim[CosineSim == 1] = 0
    NamingOrder[i] = names(which.max(CosineSim[,length(CosineSim[1,])]))
    #Fixes NB signature order and finds naming order for PO signatures
  }
  return(NamingOrder)
}

NB_PO_sig_comparison = function(Nsig, H_NB, H_PO){
  H_PO = H_PO[, order(colnames(H_PO))]
  H_NB = H_NB[, order(colnames(H_NB))]
  Signature <- colnames(H_NB)
  CosineSim = c()
  for (i in 1:Nsig) {
    CosineSim[i] = cosine(H_NB[,i], H_PO[,i])  
  }
  CosineSim = as.numeric(round(CosineSim, digits = 4))
  res = data.frame(cbind(Signature, CosineSim))
  res[,2] = as.numeric(res[,2])
  return(res)
}

Cosmic_comparison = function(Nsig, H_NB, H_PO, COSMIC){
  COSMIC = COSMIC[match(rownames(H), COSMIC$Type),]
  H_PO = H_PO[, order(colnames(H_PO))]
  H_NB = H_NB[, order(colnames(H_NB))]
  Signature <- colnames(H_NB)
  CosineSim = c(0,0,0,0)
  for (i in 1:Nsig){
    cosine_po = (cosine(as.matrix(cbind(H_PO[,i], COSMIC[4:75])))[,1])[-1]
    cosine_nb = (cosine(as.matrix(cbind(H_NB[,i], COSMIC[4:75])))[,1])[-1]
    match_po = names(which.max(cosine_po))
    match_nb = names(which.max(cosine_nb))
    CosineSim = rbind(CosineSim ,c("PO", paste("s",i,sep=""), max(cosine_po), match_po))
    CosineSim = rbind(CosineSim, c("NB", paste("s",i,sep=""), max(cosine_nb), match_nb))
  }
  CosineSim = as.data.frame(CosineSim[-1,])
  colnames(CosineSim) = c("Distribution", "Signature", "Similarity", "Cosmic Signature")
  CosineSim[,3] = as.numeric(CosineSim[,3])
  return(CosineSim)
}
