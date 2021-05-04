library("SQUAREM")
library("Rcpp")
NMFNBMMsquarem = function(M, N, alpha, seed = sample(1:100,1), arrange = TRUE, tol = 1e-5){
  K <- dim(M)[1]  # mutations
  G <- dim(M)[2]  # patients
  
  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of P matrices
  Elist <- list()            # list of E matrices
  reslist <- list()
  
  NB.em = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
    
    P <- P * ( ( ( M / (P%*%E) ) %*% t(E) ) / ( ( (alpha + M) / (alpha + P%*%E) ) %*% t(E) ) )     # update of signatures
    
    E <- E * ( (t(P) %*% ( M / (P%*%E) ) ) / (t(P) %*% ( (alpha + M) / (alpha + P%*%E) ) ) )     # update of exposures
    
    par = c(as.vector(P),as.vector(E))
    par[par <= 0] = 1e-10
    return(log(par))
  }
  
  NBobj = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
    
    y = as.vector(M)
    mu = as.vector(P%*%E)
    r = mu 
    p = which(y > 0)
    r[p] = (y * (log(y)- log(mu)) - (alpha + y) * log((alpha + y)/(alpha + mu)))[p]
    
    obj = sum(r) # euclidean distance
    
    return(obj)
  }
  
  for(i in 1:length(seed)){ 
    set.seed(seed[i])
    
    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize E
    
    init = log(c(as.vector(P),as.vector(E)))
    sres = squarem(init, fixptfn = NB.em, objfn = NBobj, control = list(tol = tol))
    
    P = matrix(exp(sres$par[1:(K*N)]), nrow = K, ncol = N)
    E = matrix(exp(sres$par[-c(1:(K*N))]), nrow = N, ncol = G)
    E = diag(colSums(P)) %*% E # normalizing 
    P = P %*% diag(1/colSums(P))
    
    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div[i] <- NBobj(sres$par) # final generalized Kullback-Leibler divergence
    reslist[[i]] = sres
  }
  
  best <- which.min(div) # Smallest GKLD value
  P = Plist[[best]]
  E = Elist[[best]]
  
  if(arrange == TRUE){
    idx = order(rowSums(E),decreasing = TRUE)
    P = P[,idx]
    E = E[idx,]
  }
  
  Output <- list()
  Output$P <-  P
  Output$E <-  E
  Output$obj <- div[best]
  Output$results <- reslist
  
  return(Output)
}