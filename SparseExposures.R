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
    if(w[k]<0) return(10^8)
  }
  #for(k in regks){
    ## Prior (regularization) for second weight Gamma distribution 
  #  regterm <- regterm + lam*w[k]
  #}
  ## Likelihood: Poisson distribution
  #term3 <- sum(v*log(w%*%H)-w%*%H)
  ## Return negative regularized log-likelihood
  res <- -sum(v*log(w%*%H)-w%*%H)+lam*sum(w[regks])
  return(res)
}

#test for 4'th genome
logL(as.numeric(W[4,]), lam = 0.1, H, largestsigs[,4], V[4,])

fit1 <- nlm(logL,as.numeric(W[4,]),lam = 0.8, H = as.matrix(H), largestsigs = largestsigs[,4], v = V[4,])
fit1$estimate
GKL1 <- fit1$minimum
barplot(V[4,])
points(1:96, W[4,]%*%H, col = "red", pch = 16)



set.seed(4)
n <- 100
p <- 10
x <- matrix(rnorm(n*p), n, p)
X <- model.matrix(~x)
y <- rnorm(n)


## i: index of which parameters aren't tuned
## l: tuning parameter
LASSO.OLS.LOVECHILD <- function(i, l) {
  penLogLik <- function(b) {
    -sum(dnorm(y - b%*%t(X), log=TRUE)) + l*sum(abs(b[!i]))
  }
  b <- vector('numeric', p+1)
  val <- nlm(penLogLik, b)
  val$estimate
}

## ordinary least squares 1 way: all unpenalized betas
LASSO.OLS.LOVECHILD(i=rep(TRUE, p), l=20)
lm(y ~ x)

## ordinary least squares another way (no tuning)
LASSO.OLS.LOVECHILD(i=rep(FALSE, p), l=0)
lm(y ~ x)

## unconstrained first half, heavy tuning forces last 5 to 0
## (not too heavy or else algorithm doesn't converge to anything...)
round(LASSO.OLS.LOVECHILD(i=1:10 %in% 1:5, l=20), 2)
round(lm(y ~ x[, 1:5])$coef, 2)
