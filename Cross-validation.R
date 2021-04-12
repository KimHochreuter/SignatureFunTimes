library(NMF)

load("patients.rda")
M <- patients
M <- M[rowSums(M)>1000,] #Removing patients with less than 1000 mutations
G <- nrow(M) #Number of patients
m <- ncol(M) #Number of mutation types

K_range <- 2:10 #Number of signatures to investigate

MSE <- list(0) #find better solution
MSE_CV <- matrix(0, ncol = 10, nrow = 5)
for (k in 2:max(K_range)){
    #MSE_CV <- 0
    for (i in 1:10){
      M_CV <- M
      CV_idx <- as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*0.01),replace=FALSE),])
      M_CV[CV_idx] <- 0 #Replace approximately 1% of cells in matrix by 0.
      for (s in 1:5){
        NMF <- nmf(M_CV, rank = k, nrun = 10)
        alpha = t(coef(NMF))
        #we replace the zero cells with the estimate from the CV. This is run 5 times, with updated M
        M_CV[CV_idx] <- (t(alpha%*%beta))[CV_idx]
        MSE_CV[s,i] <- mean((as.vector(M[CV_idx])-(t(alpha%*%beta))[CV_idx])^2)
        #print(MSE_CV)
      }
    }
  MSE[[k]] = MSE_CV
  print(MSE)
  message(100*k/(max(K_range)), "%")
}
plot(0,0,xlim = c(min(K_range), max(K_range)), ylim = c(0, 60000))
for(k in K_range){
  points(rep(k, 10), MSE[[k]][5,], xlim = c(min(K_range), max(K_range)))
  points(k, median(MSE[[k]][5,]), col = 2, pch = 16)
}
