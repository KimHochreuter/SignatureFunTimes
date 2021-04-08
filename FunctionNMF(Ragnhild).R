
library("SQUAREM")
library("Rcpp")

#############################################
##  CONE PROJECTION FOR THE NLS PROBLEM    ##
#############################################
##
## Date: June 2018
## Authors: Hobolth, Guo, Kousholt, Jensen
##
## The algorithm finds the non-negative vector hhat
## that minimizes the residual sum of squares
## RSS(h)=(v-W h)'(v-W h)
##
## The solution to the problem is based on a modification
## of the cone projection algorithm in Meyer (2013).
##
## Name: ConeProjNonNegRegr
## 
## Input:
## v: vector of length n with positive entries
## W: n times k matrix with positive entries
## 
## Output:
## x: minimizer of the RSS 
##    (entries constrained to be non-negative)
## objfct: final value of the RSS
## time: time for the projection
##
##############################################################
ConeProjNonNegRegr <- function(v,W){
  tol=1e-8 ; n=length(v) ; k=length(W)/n
  start=Sys.time()
  ##-------------------
  ## Initial step
  ##-------------------
  if(det(t(W)%*%W) == 0){print(W)}
  a=solve(t(W)%*%W)%*%t(W)%*%v
  if(min(a)>=0){
    h = a ; chck=0         ## Finish and return h=a
  } 
  if(max(a)<tol){
    h=rep(0,k) ; chck=0    ## Finish and return h=0
  }else{
    J0=which.max(a) ; chck=1
  }
  while(chck==1){
    ##------------------
    ## Step 1
    ##------------------
    W0 = W[,J0,drop=FALSE]
    if(det(t(W0)%*%W0) == 0){print(W0)}
    a = solve( t(W0)%*%W0 )%*%t(W0)%*%v
    cc = t(W)%*%(v-W0%*%a)
    if(max(cc)<tol){
      h=rep(0,k) ; h[J0]=a ; chck=0  ## Finish and return h=a
    }else{
      j1=which.max(cc) ; J1=sort(c(J0,j1))
      hw=rep(0,k) ; hw[J0]=a ; h1=hw[J1]
      chck=2
      while(chck==2){
        ##------------------
        ## Step 2
        ##------------------
        W1=W[,J1,drop=FALSE]
        if(det(t(W1)%*%W1) == 0){print(W1)}
        d1=solve(t(W1)%*%W1)%*%t(W1)%*%v
        if(min(d1)<(-tol)){
          xhat=min( h1[d1<0]/(h1[d1<0]-d1[d1<0]) )
          h1hat=xhat*d1+(1-xhat)*h1
          J1=J1[h1hat>tol]
          h1=h1hat[h1hat>tol]
        }else{
          J0=J1[d1>tol] ; chck=1
        }
      }
    }
  }
  ##---------------------------------
  time=Sys.time()-start
  RSS=sum((v-W%*%h)^2)
  return(list(x=h,objfct=RSS,time=time))
}

###############################################################
## The function runs NMF using the Cone projection algorithm, 
## where it is assumed that the data is normally distributed. 
##
## Input: V    - A non-negative matrix to factorize
##        K    - A small dimension for the new factorized matrices
##        seed - A vector of different seeds to initialize the algorithm
##        tol  - The change in the frobenius norm for when to stop
##        
## Output: W   - A non-negative matrix of dimension row(V) x K
##         H   - A non-negative matrix of dimension K x col(V), where rows sum to 1
##         frob- The value of the frobrinius norm
###############################################################
NMFNormCone <- function(V,K,tol,seed, initialH = NULL){
  N <- dim(V)[1]  # patients
  M <- dim(V)[2]  # mutations
  
  frob <- rep(0,length(seed))# vector of different frobenius values
  Wlist <- list()            # list of W matrices
  Hlist <- list()            # list of H matrices
  
  for(k in 1:length(seed)){
    set.seed(seed[k])
    if(is.null(initialH)){
    repeat{
      W <- matrix(runif(K*N)+1, nrow = N, ncol = K) # Initialize W
      H <- matrix(runif(K*M)+1, nrow = K, ncol = M) # Initialize H
    if(det(H %*% t(H)) != 0 && det(t(W) %*% W) != 0 ){break}
    }
    }else{
      W <- matrix(0, nrow = N, ncol = K)
      H <- initialH
    }
    
    RSSprev <- 0 
    
    repeat{
      for(i in 1:N){
        ConeW <-  ConeProjNonNegRegr(V[i,],t(H)) # Updating rows in W
        W[i,] <- ConeW$x
      }
      for(j in 1:M){
        ConeH <- ConeProjNonNegRegr(V[,j],W)     # Updating columns in H
        H[,j] <- ConeH$x
      }
      RSS <- sqrt(sum((V-W%*%H)^2)) # Frobenius Norm
      print(RSS)
      if(abs(RSSprev-RSS) < tol){break}
      RSSprev <- RSS
    }
    
    Wlist[[k]] <- W # final value of W
    Hlist[[k]] <- H # final value of H
    frob[k] <- RSS  # final Frobenius Norm
  }
  best <- which.min(frob)
  W <- Wlist[[best]]
  H <- Hlist[[best]]
  
  finalW <- W%*%diag(rowSums(H))   # Making sure the rows of H
  finalH <- diag(1/rowSums(H))%*%H # all sum to one.
  
  Output <- list()
  Output$W <- finalW
  Output$H <- finalH
  Output$frob <- frob[best]
  return(Output)
}

#########################################################
## Function for finding the boundaries for the matrix A,
## when there are only 2 signatures changing
## 
## Input:  W - matrix with two columns
##         H - matrix with two rows
##
## Output: Boundaries for a and b in A. 
#########################################################
Bound2Sig <- function(W,H){
  h1 <- H[1,]
  h2 <- H[2,]
  hdiff <- h1-h2
  hpos <- (hdiff > 0)
  hneg <- (hdiff < 0)
  
  rmin <- max(-h2[hpos]/hdiff[hpos])
  rmax <- min(-h2[hneg]/hdiff[hneg])
  
  w1 <- W[,1]
  w2 <- W[,2]
  zero = ((w1 == 0)*(w2 == 0))
  wfrac <- w1/(w1+w2) 
  qmin <- min(wfrac[!zero])
  qmax <- max(wfrac[!zero])
  
  Output <- c(rmax,qmax,qmin,rmin)
  return(Output)
}

#######################################################################
## Function for calculating the generalised Kullback Leibler divergence
#######################################################################
gkl.dev <- function(y, mu){
  r <- mu
  p <- which(y > 0)
  r[p] <- (y * (log(y)- log(mu)) - y + mu)[p]
  return(sum(r))
}

## C++ implementation
cppFunction('double gkldev(NumericVector y, NumericVector mu){
  double sum = 0;
  int ny = y.size();
  for(int i = 0; i < ny; i++){
    if (y[i] > 0){
      sum +=  y[i] * log(y[i]/mu[i]) - y[i] + mu[i];
    }else{
      sum += mu[i];
    }
  }
  return sum;
}')
#############################################################
## Function factorizing M into two matrices P and E of
## dimension ncol(M) x N and N x nrow(M). 
## The objective function is the generalized Kullbach-Leibler divergence(GKLD).
##
## Input: M     - non-negative data matrix of size
##        N     - Small dimension of the two new matrices
##        tol   - Change of GKLD when algorithm is stopped 
##        seed  - Vector of random seeds to initialize the matrices
##
## Output: P   - Non-negative matrix of dimension ncol(V) x K, with columns summing to one
##         E   - Non-negative matrix of dimension K x nrow(V), where rows sum to one  
##         gkl - Smallest Value of the generalized kullbach leibler
###########################################################
NMFPoisEM <- function(M,N,tol,seed, arrange = TRUE){
  K <- dim(M)[1]  # patients
  G <- dim(M)[2]  # mutations
  
  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of P matrices
  Elist <- list()            # list of E matrices
  
  
  for(i in 1:length(seed)){ 
    set.seed(seed[i])
    
    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize E
    
    GKLold <- 0
    
    repeat{
      PE <- P%*%E
      P <- P * ((M/PE) %*% t(E))      # update of signatures
      P <- P %*% diag(1/colSums(P))   # make sure the columns sum to one
      
      PE <- P%*%E
      E <- E * (t(P) %*% (M/PE))      # update of exposures
      
      
      GKL <- gkl.dev(as.vector(M),as.vector(P%*%E)) # GKLD value
      
      print(GKL)                                    # print GKLD value
      
      if(abs(GKLold - GKL) < tol){break}            # stop iterating if GKLD change less than tol
      GKLold <- GKL
    }
    
    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div[i] <- GKL   # final generalized Kullback-Leibler divergence
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
  Output$gkl <- div[best]
  
  return(Output)
}

#############################################################
## Function factorizing M into two matrices P and E of
## dimension ncol(M) x N and N x nrow(M). 
## The objective function is the generalized Kullbach-Leibler divergence(GKLD).
##
## Input: M     - non-negative data matrix of size
##        N     - Small dimension of the two new matrices
##        tol   - Maximum change of P and E when stopping 
##        seed  - Vector of random seeds to initialize the matrices
##
## Output: P   - Non-negative matrix of dimension ncol(V) x K, with columns summing to one
##         E   - Non-negative matrix of dimension K x nrow(V), where rows sum to one  
##         gkl - Smallest Value of the generalized kullbach leibler
###########################################################
NMFPoisEMsquarem = function(M,N,seed, arrange = TRUE, tol = 1e-5){
  K <- dim(M)[1]  # mutations
  G <- dim(M)[2]  # patients
  
  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of P matrices
  Elist <- list()            # list of E matrices
  reslist <- list()
  
  poisson.em = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
    
    PE <- P%*%E
    P <- P * ((M/PE) %*% t(E))      # update of signatures
    P <- P %*% diag(1/colSums(P))   # make sure the columns sum to one
    
    PE <- P%*%E
    E <- E * (t(P) %*% (M/PE))      # update of exposures
    
    par = c(as.vector(P),as.vector(E))
    par[par <= 0] = 1e-10
    return(log(par))
  }
  
  gklobj = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
    
    GKL <- gkl.dev(as.vector(M),as.vector(P%*%E)) # GKLD value
    
    return(GKL)
  }
  
  for(i in 1:length(seed)){ 
    set.seed(seed[i])
    
    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize E
    
    init = log(c(as.vector(P),as.vector(E)))
    sres = squarem(init, fixptfn = poisson.em, objfn = gklobj, control = list(tol = tol))
    
    P = matrix(exp(sres$par[1:(K*N)]), nrow = K, ncol = N)
    E = matrix(exp(sres$par[-c(1:(K*N))]), nrow = N, ncol = G)
    
    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div[i] <- gklobj(sres$par)   # final generalized Kullback-Leibler divergence
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
  Output$gkl <- div[best]
  Output$results <- reslist
  
  return(Output)
}

#############################################################
## Function factorizing M into two matrices P and E of
## dimension ncol(M) x N and N x nrow(M). 
## The objective function is the generalized Kullbach-Leibler divergence(GKLD).
## updates come from lee and seung.
##
## Input: M     - non-negative data matrix of size
##        N     - Small dimension of the two new matrices
##        tol   - Change of GKLD when algorithm is stopped 
##        seed  - Vector of random seeds to initialize the matrices
##
## Output: P   - Non-negative matrix of dimension ncol(V) x K, with columns summing to one
##         E   - Non-negative matrix of dimension K x nrow(V), where rows sum to one  
##         gkl - Smallest Value of the generalized kullbach leibler
###########################################################
NMFPoisMMsquarem = function(M,N,seed, arrange = TRUE, tol = 1e-5){
  K <- dim(M)[1]  # mutations
  G <- dim(M)[2]  # patients
  
  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of P matrices
  Elist <- list()            # list of E matrices
  reslist <- list()
  
  poisson.em = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
    
    PE <- P%*%E
    P <- P * ((M/PE) %*% t(E))      # update of signatures
    P <- P %*% diag(1/rowSums(E))   
    
    PE <- P%*%E
    E <- E * (t(P) %*% (M/PE))      # update of exposures
    E <- diag(1/colSums(P)) %*% E
    
    par = c(as.vector(P),as.vector(E))
    par[par <= 0] = 1e-10
    return(log(par))
  }
  
  gklobj = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
    
    GKL <- gkl.dev(as.vector(M),as.vector(P%*%E)) # GKLD value
    
    return(GKL)
  }
  
  for(i in 1:length(seed)){ 
    set.seed(seed[i])
    
    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize E
    
    init = log(c(as.vector(P),as.vector(E)))
    sres = squarem(init, fixptfn = poisson.em, objfn = gklobj, control = list(tol = tol))
    
    P = matrix(exp(sres$par[1:(K*N)]), nrow = K, ncol = N)
    E = matrix(exp(sres$par[-c(1:(K*N))]), nrow = N, ncol = G)
    E = diag(colSums(P)) %*% E # normalizing 
    P = P %*% diag(1/colSums(P))
    
    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div[i] <- gklobj(sres$par)   # final generalized Kullback-Leibler divergence
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
  Output$gkl <- div[best]
  Output$results <- reslist
  
  return(Output)
}

#############################################################
## Function factorizing M into two matrices P and E of
## dimension ncol(M) x N and N x nrow(M) with the use of opportunities.
## The objective function is the generalized Kullbach-Leibler divergence(GKLD).
##
## Input: M     - non-negative data matrix of size
##        N     - Small dimension of the two new matrices
##        tol   - Change of GKLD when algorithm is stopped 
##        seed  - Vector of random seeds to initialize the matrices
##        opp   - The opportunities for the different mutation types
##
## Output: P   - Non-negative matrix of dimension ncol(V) x K, with columns summing to one
##         E   - Non-negative matrix of dimension K x nrow(V), where rows sum to one  
##         gkl - Smallest Value of the generalized kullbach leibler
###########################################################
NMFPoisEMopp <- function(M,N,tol,seed, arrange = TRUE,opp){
  K <- dim(M)[1]  # patients
  G <- dim(M)[2]  # mutations
  ODfrac = diag(1/opp)
  OD = diag(opp)
  
  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of W matrices
  Elist <- list()            # list of H matrices
  
  
  for(i in 1:length(seed)){ 
    set.seed(seed[i])
    
    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize W
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize H
    
    GKLold <- 0
    
    repeat{
      PE <- P%*%E
      P <-  P * ((M/PE) %*% t(E))                      # update of signatures
      P <-  (P %*% diag(1/colSums(P)))                 # make sure the columns sum to one
      P <- ODfrac %*% P                                # add effects of opportunities
      
      PE <- P%*%E
      E <- E * (t(P) %*% (M/PE))                       # update of exposures
      
      
      GKL <- gkl.dev(as.vector(M),as.vector(OD%*%P%*%E)) # GKLD value
      
      print(GKL)                                    # print GKLD value
      
      if(abs(GKLold - GKL) < tol){break}            # stop iterating if GKLD change less than tol
      GKLold <- GKL
    }
    
    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div <- GKL      # final generalized Kullback-Leibler divergence
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
  Output$gkl <- div[best]
  
  return(Output)
}


#############################################################
## A function for factorizing V into two matrices W and H of
## dimension col(V) x K and K x row(V), where we assume 
## the data matrix V is normal distributed. 
## The objective function used is the Frobrinius norm 
######################################################
NMFNormEM <- function(V,K,tol,seed){
  N <- dim(V)[1]  # patients
  M <- dim(V)[2]  # mutations
  
  div <- rep(0,length(seed)) # vector of different GKL values
  Wlist <- list()            # list of W matrices
  Hlist <- list()            # list of H matrices

  for(i in 1:length(seed)){
    set.seed(seed[i])
    
    W <- matrix(runif(K*N), nrow = N, ncol = K)  # Initialize W
    H <- matrix(runif(K*M), nrow = K, ncol = M)  # Initialize H
    
    Frobold <- 0
    
    repeat{
      top <- t(W)%*%V
      bottom <- t(W)%*%W%*%H
      H <- H * top/bottom # update of H
      
      top <- V%*%t(H)
      bottom <- W%*%H%*%t(H)
      W <- W * top/bottom # update of W
      
      Frob <- sum((V-(W%*%H))^2)
      print(Frob)
      
      if(abs(Frobold - Frob) < tol){break} 
      Frobold <- Frob
    }
    
    Wlist[[i]] <- W # final W
    Hlist[[i]] <- H # final H
    div <- Frob     # final generalised Kullback-Leibler divergence
  }
  
  best <- which.min(div) # Smallest generalised Kullback-Leibler divergence
  
  Output <- list()
  Output$W <- Wlist[[best]]
  Output$H <- Hlist[[best]]
  Output$gkl <- div[best]
  
  return(Output)
}

#############################################################
## Function factorizing M into two matrices P and E of
## dimension ncol(M) x N and N x nrow(M).
## The objective function is the euclidian distance also called Frobrinious norm.
## updates are from lee and seung.
##
## Input: M     - non-negative data matrix of size
##        N     - Small dimension of the two new matrices
##        tol   - Change of GKLD when algorithm is stopped 
##        seed  - Vector of random seeds to initialize the matrices
##
## Output: P   - Non-negative matrix of dimension ncol(V) x K, with columns summing to one
##         E   - Non-negative matrix of dimension K x nrow(V), where rows sum to one  
##         gkl - Smallest Value of the generalized kullbach leibler
###########################################################
NMFNormMMsquarem = function(M,N,seed, arrange = TRUE, tol = 1e-5){
  K <- dim(M)[1]  # mutations
  G <- dim(M)[2]  # patients
  
  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of P matrices
  Elist <- list()            # list of E matrices
  reslist <- list()
  
  norm.em = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
    
    P <- P * ( (M %*% t(E)) / (P%*%E%*%t(E)) )     # update of signatures
    
    E <- E * ( (t(P) %*% M) / (t(P)%*%P%*%E) )    # update of exposures
    
    par = c(as.vector(P),as.vector(E))
    par[par <= 0] = 1e-10
    return(log(par))
  }
  
  frobobj = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
    
    frob <- sum((as.vector(M) - as.vector(P%*%E))^2) # euclidian distance
    
    return(frob)
  }
  
  for(i in 1:length(seed)){ 
    set.seed(seed[i])
    
    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize E
    
    init = log(c(as.vector(P),as.vector(E)))
    sres = squarem(init, fixptfn = norm.em, objfn = frobobj, control = list(tol = tol))
    
    P = matrix(exp(sres$par[1:(K*N)]), nrow = K, ncol = N)
    E = matrix(exp(sres$par[-c(1:(K*N))]), nrow = N, ncol = G)
    E = diag(colSums(P)) %*% E # normalizing 
    P = P %*% diag(1/colSums(P))
    
    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div[i] <- frobobj(sres$par)   # final generalized Kullback-Leibler divergence
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
  Output$gkl <- div[best]
  Output$results <- reslist
  
  return(Output)
}

#############################################################
## Function factorizing M into two matrices P and E of
## dimension ncol(M) x N and N x nrow(M).
## The objective function is likelihood of the negative binomial distribution
##
## Input: M     - non-negative data matrix of size
##        N     - Small dimension of the two new matrices
##        tol   - Change of GKLD when algorithm is stopped 
##        seed  - Vector of random seeds to initialize the matrices
##
## Output: P   - Non-negative matrix of dimension ncol(V) x K, with columns summing to one
##         E   - Non-negative matrix of dimension K x nrow(V), where rows sum to one  
##         obj - Smallest object value
###########################################################
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

#####################################################
## Multiplicative model updates to factorize a non-negative V(NxM) into 
## two non-negative matrices W(NxK) and H(KxM).
## 
## Input : V    - Data matrix of dimension NxM
##         K    - Number of signatures
##         seed - Seed for initialization of W and H
##         tol  - Stopping tolerance of the change of GKL
##         X    - Design matrix used to update H.
##
## Output: W    - Non-negative matrix of dimension NxK
##         H    - Non-negative matrix of dimension KxM
##         gkl  - Generalized kullback leibler value for final data
###########################################################################
NMFglm.Rglm <- function(V,K,seed,tol,X){
  N <- dim(V)[1]  # patients
  M <- dim(V)[2]  # mutations
  
  set.seed(seed)
  
  W <- matrix(runif(K*N)+1, nrow = N, ncol = K) # Initialize W
  W <- W %*% diag(1/colSums(W))                 # Make sure the column sum to one
  
  H <- matrix(runif(K*M)+1, nrow = K, ncol = M) # Initialize H
  
  GKLold <- 0
  
  repeat{
    WH <- W%*%H
    y <- as.vector(t(H * (t(W) %*% (V/WH))))                  # output
    fit <- glm(y ~ 0+X, family = quasipoisson(link = "log"))  # fit                                   
    H <- matrix(exp(predict(fit)), nrow = K, byrow = T)       # update of H
    
    WH <- W%*%H
    W <- W * ((V/WH) %*% t(H))      # update of W
    W <- W %*% diag(1/colSums(W))   # make sure the column sum to one
    
    GKL <- gkl.dev(as.vector(V),as.vector(W%*%H))
    print(GKL)
    
    if(abs(GKLold - GKL) < tol){break} 
    GKLold <- GKL
  }
  
  Output <- list()
  Output$W <- W
  Output$H <- H
  Output$gkl <- GKL
  
  return(Output)
}

###########################################################
## Function updates H via GLM instead of poisson update
##
## Input : y - Original update via poisson
## 
## Output: eta - New update with glm dependency
###########################################################
glm.update <- function(y,X, maxiter = 25, epsilon = 1e-8){
  mu <- y + 0.1            # starting value
  dev.old <- gkl.dev(y,mu) # start deviance
  for(d in 1:maxiter){
    z <- log(mu) + (y-mu)/mu
    w <- sqrt(mu)
    fit <- .Call(stats:::C_Cdqrls, X*w, z*w, epsilon, FALSE) # solving least squares
    coef <- fit$coefficients                                 # coefficients
    eta <- as.vector(X%*%coef)
    mu <- exp(eta)
    dev.new <- gkl.dev(y,mu)
    if(2*abs(dev.old - dev.new)/(0.1 + abs(2*dev.new)) < epsilon){break}
    dev.old <- dev.new
  } 
  return(mu)
}

#####################################################
## Multiplicative model updates to factorize a non-negative V(NxM) into 
## two non-negative matrices W(NxK) and H(KxM).
## 
## Input : V    - Data matrix of dimension NxM
##         K    - Number of signatures
##         seed - Seed for initialization of W and H
##         tol  - Stopping tolerance of the change of GKL
##         X    - Design matrix used to update H.
##
## Output: W    - Non-negative matrix of dimension NxK
##         H    - Non-negative matrix of dimension KxM
##         gkl  - Generalized kullback leibler value for final data
###########################################################################
NMFglm <- function(V,K,seed,tol,X, W.init = NULL, H.init = NULL){
  N <- dim(V)[1]  # patients
  M <- dim(V)[2]  # mutations
  
  set.seed(seed)
  
  if(is.null(W.init)){
    W <- matrix(runif(K*N)+1, nrow = N, ncol = K) # Initialize W
  }else{
    W <- W.init
  }
  W <- W %*% diag(1/colSums(W))                   # make sure the column sum to one
  
  
  if(is.null(H.init)){
    H <- matrix(runif(K*M)+1, nrow = K, ncol = M)  # Initialize H
  }else{
    H <- H.init
  }
  
  
  GKLold <- 0
  
  repeat{
    WH <- W%*%H
    y <- as.vector(t(H * (t(W) %*% (V/WH))))                # output
    H <- matrix(glm.update(y,X), nrow = K, byrow = T)       # update of H
    
    
    WH <- W%*%H
    W <- W * ((V/WH) %*% t(H))      # update of W
    W <- W %*% diag(1/colSums(W))   # make sure the column sum to one
    
    GKL <- gkl.dev(as.vector(V),as.vector(W%*%H))
    print(GKL)
    
    if(abs(GKLold - GKL) < tol){break} 
    GKLold <- GKL
  }
  
  Output <- list()
  Output$W <- W
  Output$H <- H
  Output$gkl <- GKL
  
  return(Output)
}

#########################################################
## Function for calculating the angle between two vectors 
## Input: x,y - two vectors of same dimensions 
## 
## Output: The angle in radians
#########################################################
angle <- function(x,y){
  dot.prod <- sum(x*y) 
  norm.x <- sqrt(sum(x^2))
  norm.y <- sqrt(sum(y^2))
  frac <- dot.prod / (norm.x * norm.y)
  theta <- ifelse(frac > 1, 0, acos(frac))
  return(as.numeric(theta))
}

##################################################
## Cosine similarity between a vector x and y
###############################################
similarity <- function(x,y){
  dot.prod <- sum(x*y) 
  norm.x <- sqrt(sum(x^2))
  norm.y <- sqrt(sum(y^2))
  frac <- dot.prod / (norm.x * norm.y)
  return(as.numeric(frac))
}


#############################################################
## Distance between the two matrices H1 and H2, by summing the distance 
## between signatures with the smallest distance. 
frob_dist <- function(H1,H2){
  K <- nrow(H1)
  d <- numeric(K)
  m <- numeric(K)
  dist <- sapply(1:K, function(y) sapply(1:K,function(x) sum((H1[x,] - H2[y,])^2)))
  residual <- 0
  for(s in 1:K){
    min.dist <- min(dist)
    remove = which(dist == min.dist, arr.ind = TRUE)
    dist[remove[1,1],] <- Inf
    dist[,remove[1,2]] <- Inf
    d[remove[1,1]] <- min.dist
    m[remove[1,1]] <- remove[1,2]
    residual = residual + min.dist
  }
  
  Output <- list()
  Output$total <- residual # The total distance
  Output$sig   <- d        # Individual distances numbered after signatures in H1
  Output$match <- m
  return(Output)
}

#############################################################
## Distance between the two matrices H1 and H2, by taking the 
## average cosine similarity between matrices
cos.sim <- function(H1,H2){
  K <- nrow(H1)
  d <- numeric(K)
  m <- numeric(K)
  dist <- sapply(1:K, function(y) sapply(1:K,function(x) similarity(H1[x,],H2[y,])))
  dist <- as.matrix(dist)
  residual <- 0
  for(s in 1:K){
    max.dist <- max(dist)
    remove = which(dist == max.dist, arr.ind = TRUE)
    dist[remove[1,1],] <- 0
    dist[,remove[1,2]] <- 0
    d[remove[1,1]] <- max.dist
    m[remove[1,1]] <- remove[1,2]
    residual = residual + max.dist
  }
  
  Output <- list()
  Output$total <- residual/K # The average cosine similarity
  Output$sig   <- d          # Individual cosine similarity numbered after signatures in H1
  Output$match <- m
  return(Output)
}

##################################
## Function for calculating the Kullback Leibler divergence
##################################
kl.dev <- function(y, mu){
  logy = logmu = rep(-100,length(y))
  py <- which(y > 0)
  pmu <- which(mu > 0)
  
  logy[py] = log(y)[py]
  logmu[pmu] = log(mu)[pmu]
  
  r <- (y * (logy - logmu))
  return(sum(r))
}

#############################################
## Hellinger distance between two vectors
#############################################
hellinger = function(y,mu){
  1/sqrt(2)*sqrt(sum((sqrt(y)-sqrt(mu))^2))
}

#############################################################
## Kullback-Leibler dist between the two matrices P1 and P2, by summing the distance 
## between signatures with the smallest distance. 
kl_dist <- function(P1,P2){
  N <- ncol(P1)
  d <- numeric(N)
  m <- numeric(N)
  dist <- sapply(1:N, function(y) sapply(1:N,function(x) kl.dev(P1[,x],P2[,y])))
  residual <- 0
  for(s in 1:N){
    min.dist <- min(dist)
    remove = which(dist == min.dist, arr.ind = TRUE)
    dist[remove[1,1],] <- Inf
    dist[,remove[1,2]] <- Inf
    d[remove[1,1]] <- min.dist
    m[remove[1,1]] <- remove[1,2]
    residual = residual + min.dist
  }
  
  Output <- list()
  Output$total <- residual # The total distance
  Output$sig   <- d        # Individual distances numbered after signatures in P1
  Output$match <- m
  return(Output)
}
#############################################################
## Hellinger dist between the two matrices P1 and P2, by summing the distance 
## between signatures with the smallest distance. 
h_dist = function(P1,P2){
  N <- ncol(P1)
  d <- numeric(N)
  m <- numeric(N)
  dist <- sapply(1:N, function(y) sapply(1:N,function(x) hellinger(P1[,x],P2[,y])))
  residual <- 0
  for(s in 1:N){
    min.dist <- min(dist)
    remove = which(dist == min.dist, arr.ind = TRUE)
    dist[remove[1,1],] <- Inf
    dist[,remove[1,2]] <- Inf
    d[remove[1,1]] <- min.dist
    m[remove[1,1]] <- remove[1,2]
    residual = residual + min.dist
  }
  
  Output <- list()
  Output$total <- residual # The total distance
  Output$sig   <- d        # Individual distances numbered after signatures in P1
  Output$match <- m
  return(Output)
}

##########################################################################
## Tranformation matrix of one signature with another
#########################################################################
Tmat = function(lambda,s,s.mix){
  A = diag(N)
  A[c(s,s.mix),s] = c(lambda,1-lambda)
  return(A)
}

#########################################################
## Bound of one signature w.r.t. another 
#########################################################
lambda.range = function(P,E){
  if(ncol(P) > 2 | nrow(E) > 2){
    print("Wrong dimension of input")
  }
  p1 = P[,1]
  p2 = P[,2]
  pdiff = p1 - p2
  pos = (pdiff < 0)
  
  lmin = max(p1[pos]/pdiff[pos])
  
  e1 = E[1,]
  e2 = E[2,]
  zero = ((e1 == 0)*(e2 == 0))
  efrac <- e2/(e1+e2)
  lmax = min(efrac[!zero])
  res = c(lmin,lmax)
  
  if(lmin > lmax){cat("Lambda range weird:",res,"\n")}
  
  return(res)
}


## C++ implementation of SFS function
sourceCpp("C:/Users/au543194/Desktop/Kodning/data/SFS_c.cpp")
#########################################################
## Find range of feasible solutions by simulation
## 
## Input :  P     - the resulting signatures
##          E     - the resulting exposures
##          iter  - number of iterations in simulation
##
## Output: 
#########################################################
sfs = function(P,E, max.iter = 10^5, check = 1000, beta.par = 1, save = TRUE, eps = 1e-8){
  start_time = Sys.time()
  N = ncol(P)
  K = nrow(P)
  G = ncol(E)
  
  sig = rep(c(1:N), check + 2)
  avgDev = numeric(max.iter/check)
  
  if(save){
  expos = matrix(0, nrow = (max.iter+1)*N, ncol = G)
  probs = matrix(0, nrow = (max.iter+1)*N, ncol = K)
  
  probs[1:N,] = t(P) # Inserting the estimated solution as a feasible solution
  expos[1:N,] = E
  }else{
    probs = matrix(P, nrow = (check+2)*N, ncol = K, byrow = T)
  }
  
  diff.res = 0;
  
  ## Transformation matrix of one signature with another
  Tmat = function(lambda,s,s.mix){
    A = diag(N)
    A[c(s,s.mix),s] = c(1-lambda,lambda)
    return(A)
  }
  
  Pmat = P
  Emat = E
  
  ## Setting very small entries to zero
  zeroP = which(abs(Pmat) < 1e-10, arr.ind = TRUE)                      
  zeroE = which(abs(Emat) < 1e-10, arr.ind = TRUE)
  np = nrow(zeroP)
  ne = nrow(zeroE)
  if(np>0){
    for(i in 1:np){
      Pmat[zeroP[i,1],zeroP[i,2]] = 1e-10
    }
  }
  if(ne>0){
    for(i in 1:ne){
      Emat[zeroE[i,1],zeroE[i,2]] = 1e-10
    }
  }
  
  ## Running the algorithm
  for(i in 1:max.iter){
    for(s in 1:N){
        
        # Change signature s, with another signature 
        s.mix = ifelse(N > 2, sample(x = c(1:N)[-s], size = 1),   # sample another signature to mix with
                                     c(1:N)[-s])
        bound = lambda.range(Pmat[,c(s,s.mix)],Emat[c(s,s.mix),]) # the range for that signature
        lr = bound[2]-bound[1]
        if(bound[2] > 0.5) cat("Bound is large", bound, "for signature", s, "mixed with",s.mix, "\n")
        if(abs(lr) > 1e-5){
          #if(lr > 0.5) cat("Range is", bound, "for signature", s, "mixed with",s.mix, "\n")
        b = rbeta(1,shape1 = beta.par, shape2 = beta.par)
        lambda = b*bound[1]+ (1-b)*bound[2]                       # choose random number in range from beta dist
        
        mat = Tmat(lambda,s,s.mix)                                # transformation matrix
        Pmat = Pmat%*%mat                                         # new signatures
        Emat = solve(mat)%*%Emat                                  # new exposures
        
        # set very small entrances to zero
        zeroP = which(abs(Pmat[,s]) < 1e-10)                      
        zeroE = which(abs(Emat[s.mix,]) < 1e-10)
        np = length(zeroP)
        ne = length(zeroE)
        if(np>0){
          Pmat[zeroP,s] = 1e-10
        }
        if(ne>0){
          Emat[s.mix,zeroE] = 1e-10
        }
      }
    }
    
    if((i %% (max.iter/10)) == 0){print(i)} # print iteration number sometimes
    
    if(sum(Pmat < 0) > 0){print(Pmat);print(i);print(s);break}      # stop if any entry is negative
    if(sum(Emat < 0) > 0){print(Emat);print(i);print(s.mix);break}  # stop if any entry is negative
    
    # save the results
    if(save){
      probs[(1+N*i):(N+N*i),] = t(Pmat)
      expos[(1+N*i):(N+N*i),] = Emat 
    }else{
      iter = ((i-1) %% check)
      probs[(1+N*(iter+2)):(N+N*(iter+2)),] = t(Pmat)
    }
    
    # check whether to stop the algorithm
    if((i %% check) == 0){
      if(save){
        sig = rep(c(1:N), i+1)
        prob.diff = sapply(1:K, function(y) tapply(probs[1:(N+N*i),y], sig, function(r) diff(range(r))))
      }else{
        prob.diff           = sapply(1:K, function(y) tapply(probs[,y], sig, function(r) diff(range(r))))
        probs[1:N,]         = sapply(1:K, function(y) tapply(probs[,y], sig, min))
        probs[(1+N):(2*N),] = sapply(1:K, function(y) tapply(probs[,y], sig, max)) 
      }
      
      diff.ny = mean(prob.diff)
      avgDev[i/check] = diff.ny
      if((diff.ny - diff.res) < eps){
        print(i)
        break
      }else{
        diff.res = diff.ny
      }
    }
  }
  
  time = Sys.time() - start_time
  Output = list()
  if(save){
    Output$probs = probs[1:(i*N),] # include all solutions
    Output$expos = expos[1:(i*N),] # include all solutions
  }
  Output$avgDev = avgDev
  Output$range = probs[1:(2*N),] # only include minimum and maximum probabilities 
  Output$diff = diff.res
  Output$n.iter = i
  Output$time = time
  return(Output)
}

###############################################################################
## Change because of variance. Illustrated by poisson parametric bootstrapping
##
## Input :  P     - the resulting signatures
##          E     - the resulting exposures
##          iter  - number of iterations in bootstrap
##
## Output: 

boot_pois = function(P,E, iter, same.init = FALSE){
  N = ncol(P)
  K = nrow(P)
  G = ncol(E)
  
  P.probs = matrix(0, nrow = iter*N, ncol = K)
  E.expos = matrix(0, nrow = iter*N, ncol = G)
  
  PE = P%*%E
  no.seed = c(1,200,400,600,800,1000,1200,1400,1600,1800)

  for(i in 1:iter){
    set.seed(i)
    M.sim = matrix(rpois(K*G, lambda = PE), nrow = K) # compute random sample
    res = NMFPoisEMsquarem(M.sim, N, seed = ifelse(same.init,no.seed,i*c(1,200,400,600,800,1000,1200,1400,1600,1800)), tol = 1e-3)
    dist = cos.sim(t(P),t(res$P))
    P.probs[(1+N*(i-1)):(N+N*(i-1)),] = t(res$P[,dist$match])
    E.expos[(1+N*(i-1)):(N+N*(i-1)),] = res$E[dist$match,]
  }
  
  Output = list()
  
  Output$P = P.probs
  Output$E = E.expos
  return(Output)
}


###############################################################################
## Multinomial bootstraps of original genome matrix
###############################################################################
boot_multi = function(M, N, iter, same.init = FALSE){
  K = nrow(M)
  G = ncol(M)
  
  P.probs = matrix(0, nrow = iter*N, ncol = K)
  E.expos = matrix(0, nrow = iter*N, ncol = G)

  no.seed = c(1,200,400,600,800,1000,1200,1400,1600,1800)
  
  for(i in 1:iter){
    set.seed(i)
    M.boot = apply(M, 2, function(x) rmultinom(1,size = sum(x), prob = x))
    res = NMFPoisEMsquarem(M.boot, N, seed = ifelse(same.init,no.seed,i*c(1,200,400,600,800,1000,1200,1400,1600,1800)), tol = 1e-3)
    P.probs[(1+N*(i-1)):(N+N*(i-1)),] = t(res$P)
    E.expos[(1+N*(i-1)):(N+N*(i-1)),] = res$E
  }
  
  Output = list()
  
  Output$P = P.probs
  Output$E = E.expos
  return(Output)
}

