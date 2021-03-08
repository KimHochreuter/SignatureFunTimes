# SPARSE EXPOSURES??
library(nnls)

Sigs <- as.matrix(COSMIC_Mutational_Signatures_v3_1[,4:75])
rownames(Sigs) <- c(COSMIC_Mutational_Signatures_v3_1[,3]$TYPE)

#initialize W and H - can also just use NMF? Maybe ignore
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

klargestsigs <- function(vec, k = 10){
  index <- which(vec >= sort(vec, decreasing=T)[k], arr.ind=TRUE)
  #return(colnames(W)[index])
  return(index)
}
largestsigs <- apply(W, MARGIN = 1, FUN = klargestsigs)
