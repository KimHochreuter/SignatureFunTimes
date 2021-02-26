MSE <- c(0, 0, 10^7) #find better solution
for (k in 1:length(K_range)){
  for (lambda in lambda_range) {
    MSE_CV <- 0
    for (i in 1:10){
      M_CV <- M
      CV_idx <- as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*0.1),replace=FALSE),])
      M_CV[CV_idx] <- 0 #Replace approximately 1% of cells in matrix by 0.
      MSE_CV <- 0
      for (s in 1:5){ #As suggested by paper, step 5c is repeated 5 times.
        for (l in 1:20){ #Repeated estimation to ensure convergence.
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