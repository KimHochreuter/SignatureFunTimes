# SPARSE EXPOSURES??
'''This function/script takes a predetermined K as input. 
Combine with the Cross-validation. Maybe???
In this we use the assumption V = W%*%t(H)'''

#initialize W and H - can also just use NMF? Maybe ignore
K <- 5
V <- patients
M <- nrow(V)
N <- ncol(V)
W <- matrix(runif(M*K, 0, 100), ncol = K, nrow = M)
H <- matrix(runif(K*N, 0, 100), ncol = N, nrow = K)

#initialize W and H
NMF <- nmf(t(V),rank=(K),nrun=10)
H <- basis(NMF)

colnames(H) <- c(paste0("S",1:(ncol(H))))
rownames(H) <- colnames(V)

W <- t(coef(NMF))
rownames(W) <- 1:nrow(W)
colnames(W) <- colnames(H)

#update cols in H 
updateColInH <- function(hi, Ri, wi, delta, i){
  v <- (t(Ri)%*%wi + delta*hi)/(norm(wi, type = "2") + delta)
  v <- sort(v, decreasing = T)
  alpha <- 1
  for (j in 1:n){
    alpha <- alpha - v[i]
    C = alpha/j
    if (j == n || C<= -v[(j+1)])
      break
  }
  h[(j+1):n] <- 0
  h(j) <- h
  return(h)
}

delta <- 10^(-4)
R <- V-W%*%t(H)
