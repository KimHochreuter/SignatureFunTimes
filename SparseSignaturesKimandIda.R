library(NMF)
library(nnls)
library(Biobase)
library(nnlasso)
#THIS SCRIPT ASSUMES THAT K'S LOWER VALUE ALWAYS IS 2!!!!

###########################################################
# Step 1
###########################################################

load("patients.rda")
load("background.rda")
M <- patients

###########################################################
# Step 2
###########################################################

M <- M[rowSums(M)>1000,] #Removing patients with less than 1000 mutations
G <- nrow(M) #Number of patients
m <- ncol(M) #Number of mutation types

###########################################################
# Step 3
###########################################################

K_range <- 2:10 #Number of signatures to investigate
lambda_range <- seq(0,1,length.out = 5) #Level of sparsity to investigate

###########################################################
# Step 4
###########################################################

#Estimation of alpha_0, i.e. background exposure.
#Beta_0 is the background signature, which is given beforehand ind the file "background".

#background exposure
alpha_0 <- matrix(0,nrow=G,ncol=1)
background = as.matrix(background)
for(i in 1:G) {
  alpha_0[i,] <- nnls(background,as.matrix(M[i,]))$x
}

#Initialize beta using repeated NMF
starting_beta <- array(list(),c(length(K_range),1))
curr_M <- M - (alpha_0 %*% t(background)) #alt det der mangler at blive forklaret af af de ægte signaturer
curr_M[curr_M<0] <- 0 # negative entries are set to zero - problematic?

#For each value of K we perform repeated NMF to gain an initial value of beta, which is the estimated signatures.
k_idx <- 0
for(k in K_range) {
  k_idx <- k_idx + 1
  NMF <- nmf(t(curr_M),rank=(k-1),nrun=10)
  beta <- basis(NMF)
  beta <- t(beta)
  beta <- rbind(t(background),beta)
  
  # normalize and save beta for the current value of K
  rownames(beta) <- c("Background",paste0("S",1:(nrow(beta)-1)))
  colnames(beta) <- colnames(M)
  beta <- beta / rowSums(beta)
  starting_beta[[k_idx,1]] <- beta
  
  alpha <- t(coef(NMF))
  rownames(alpha) <- 1:nrow(alpha)
  colnames(alpha) <- rownames(beta[[k]])[-1]
  starting_alpha[[k_idx]] <- cbind(alpha_0,alpha)
}
beta <- starting_beta
alpha <- starting_alpha

#MSE <- matrix(rep(0, length(K_range)*length(lambda_range)*3), ncol = 3)
MSE <- c(0, 0, 10^7) #find better solution
for (k in 1:length(K_range)){
  for (lambda in lambda_range) {
    MSE_CV <- 0
    for (i in 1:10){
      M_CV <- M
      CV_idx <- as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*0.1),replace=FALSE),])
      M_CV[CV_idx] <- 0
      MSE_CV <- 0
      for (s in 1:5){
        for (l in 1:20){
          for (j in 1:m){ #beta can only be updated col-wise, due to the function nnlasso
            target <- M_CV[,j] - alpha[[k]][,1]*background[j] # vi får negative target values?
            #target[target<0] <- 0
            beta[[k]][2:(k+1),j] <- as.vector(nnlasso(x = alpha[[k]][,2:(k+1), drop=FALSE], 
                                                      y = target, 
                                                      family = "normal", 
                                                      lambda = lambda, 
                                                      intercept = FALSE, 
                                                      normalize = FALSE, 
                                                      tol = 1e-05, 
                                                      maxiter = 10000, 
                                                      path = FALSE)$coef[2,])
          }
          for(i in 1:G){ #update alphas
            alpha[[k]][i,] <- nnls(t(as.matrix(beta[[k]])),as.matrix(M_CV[i,]))$x
          }
        }
        #we replace the zero cells with the estimate from the CV. This is run 5 times, with updated M
        M_CV[CV_idx] <- (alpha[[k]]%*%beta[[k]])[CV_idx]
      }
      MSE_CV <- MSE_CV +  mean((as.vector(M[CV_idx,])-(alpha[[k]]%*%beta[[k]])[CV_idx,])^2)
    }
    MSE <- rbind(MSE, c(k, lambda, MSE_CV))
  }
}
colnames(MSE) <- c("K", "lambda", "MSE")
MSE[which.min(MSE[,3]),]
K_best <- MSE[which.min(MSE[,3]),1]
lambda_best <- MSE[which.min(MSE[,3]),2]

alpha <- starting_alpha[[K_best]]
beta <- starting_beta[[K_best]]
for (l in 1:20){
  for (j in 1:m){ #beta can only be updated col-wise, due to the function nnlasso
    target <- M[,j] - alpha[,1]*background[j]
    beta[2:(K_best+1),j] <- as.vector(nnlasso(x = alpha[, 2:(K_best+1),drop=FALSE], 
                                              y = target, 
                                              family = "normal", 
                                              lambda = lambda_best, 
                                              intercept = FALSE, 
                                              normalize = FALSE, 
                                              tol = 1e-05, 
                                              maxiter = 10000, 
                                              path = FALSE)$coef[2,])
  }
  for(i in 1:G){ #update alphas
    alpha[i,] <- nnls(t(as.matrix(beta)),as.vector(M[i,]))$x
  }
}
par(mfrow=c(2,2))
for(i in 1:(K_best+1)){
  barplot(beta[i,])
}


mean((cbind(as.vector(M) - as.vector(as.matrix(alpha)%*%as.matrix(beta[[K_best]]))))^2)
as.vector(as.matrix(alpha)%*%as.matrix(beta[[K_best]]))
dim(alpha)
dim(beta)
