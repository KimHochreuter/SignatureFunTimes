# SPARSE EXPOSURES??
library(nnls)

Sigs <- as.matrix(COSMIC_Mutational_Signatures_v3_1[,4:75])
rownames(Sigs) <- c(COSMIC_Mutational_Signatures_v3_1[,3]$TYPE)

#initialize W
K <- ncol(Sigs)
V <- patients
M <- nrow(V)
N <- ncol(V)
H <- t(Sigs)

W <- matrix(0,nrow=M,ncol=K)
rownames(W) <- 1:nrow(W)
colnames(W) <- rownames(H)
for(i in 1:M) {
  W[i,] <- nnls(t(H),as.vector(V[i,]))$x
}
table(rowSums(W) == 0)
table(colSums(W) == 0)


#Pick out the k largest sigs. k determined by CV. Regularise the remaining signatures.
klargestsigs <- function(vec, k = 10){
  index <- which(vec >= sort(vec, decreasing=T)[k], arr.ind=TRUE)
  #return(colnames(W)[index])
  return(index)
}
largestsigs <- apply(W, MARGIN = 1, FUN = klargestsigs)

#likelihood
logL <- function(w,lam,H,largestsigs, v){
  regterm <- 0
  regks <- setdiff(c(1:K)[w > 0], largestsigs)
  for(k in 1:K){
    if(w[k]<0) return(10^6)
  }
  for(k in regks){
    ## Prior (regularization) for second weight Gamma distribution 
    regterm <- regterm -lam*w[k]
  }
  ## Likelihood: Poisson distribution
  term3 <- sum( -w%*%H+v*log(w%*%H) )
  ## Return negative regularized log-likelihood
  res <- -regterm-term3
  print(res)
  return(res)
}

#test for 4'th genome
logL(as.numeric(W[4,]), lam = 0.1, H, largestsigs[,4], V[4,])

fit1 <- nlm(logL,as.numeric(W[4,]),lam = 0.1, H = as.matrix(H), largestsigs = largestsigs[,4], v = V[4,])
fit1$estimate
GKL1 <- fit1$minimum
barplot(V[4,])
points(1:96, W[4,]%*%H, col = "red", pch = 16)
