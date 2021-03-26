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

K_range <- 2:15 #Number of signatures to investigate

MSE <- list(0) #find better solution
MSE_CV <- matrix(0, ncol = 10, nrow = 5)

AIC_collect = list(0)
BIC_collect = list(0)

AIC_normal = matrix(0, ncol = 10, nrow = 5)
BIC_normal = matrix(0, ncol = 10, nrow = 5)

CV_idx <- replicate(10, list(as.matrix(expand.grid(1:G, 1:m)[sample(1:(G*m), ceiling(G*m*0.01),replace=FALSE),])))

for (k in 2:max(K_range)){
    #MSE_CV <- 0
    for (i in 1:10){
      M_CV <- M
      M_CV[CV_idx[[i]]] <- 0 #Replace approximately 1% of cells in matrix by 0.
      for (s in 1:5){
        NMF <- nmf(M_CV, rank = k, nrun = 10)
        alpha = t(coef(NMF))
        beta = t(basis(NMF))
        #we replace the zero cells with the estimate from the CV. This is run 5 times, with updated M
        M_CV[CV_idx[[i]]] <- (t(alpha%*%beta))[CV_idx[[i]]]
        MSE_CV[s,i] <- mean((as.vector(M[CV_idx[[i]]])-(t(alpha%*%beta))[CV_idx[[i]]])^2)
        AIC_normal[s,i] = 96*k + MSE_CV[s,i]
        BIC_normal[s,i] = k*log(m*G) + MSE_CV[s,i]
        #print(MSE_CV)
        
      }
    }
  AIC_collect[[k]] = AIC_normal
  BIC_collect[[k]] = BIC_normal
  MSE[[k]] = MSE_CV
  #print(MSE)
  message(100*k/(max(K_range)), "%")
}
plot(0,0,xlim = c(min(K_range), max(K_range)), ylim = c(0, 50000))
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
