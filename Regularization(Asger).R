##------------------------------------------------------------
## Regularization experiment 
##------------------------------------------------------------
## 20 mutation types and 2 signatures 
H1 <- c(rep(4,10),rep(1,10))/50
H2 <- c(rep(1,10),rep(4,10))/50
plot(1:20,H1,pch=19,col="black",ylim=c(0,0.5),
     xlab="Mutation type",ylab="Probability")
points(1:20,H2,pch=19,col="red")
legend("topright",c("H1","H2"),col=c("black","red"),
       lty=1,lwd=2,bty="n",cex=2)
##-----------------------------------------------------------
## True weight of signatures and simulation from Poisson
mn.v <- 10*H1+100*H2
plot(1:20,mn.v,xlab="Mutation type",ylab="Number of mutations",
     ylim=c(0,12),col="black",cex=1.8,pch=4)
set.seed(5)
v <- rpois(n=20,mn.v)
cat("Total number of mutations:",sum(v),"\n")
points(1:20,v,col="purple",pch=3,cex=1.5)
##------------------------------------------------------------
## Negative regularized log-likelihood (for minimization)
##------------------------------------------------------------
lnL <- function(w,lam,a,b,H1,H2,v){
  ## Weights must be positive
  if(w[1]<0.001) return(10^5)
  if(w[2]<0.001) return(10^5)
  ## Prior (shrinkage) for first weight: Exponential distribution
  term1 <- -lam*w[1]
  ## Prior (regularization) for second weight Gamma distribution 
  term2 <- -b*w[2]+(a-1)*log(w[2])
  ## Likelihood: Poisson distribution
  term3 <- sum( -w[1]*H1-w[2]*H2+v*log(w[1]*H1+w[2]*H2) )
  ## Return negative regularized log-likelihood
  return(-term1-term2-term3)
}
##-------------------------------------------------------------------------
## Test the function
lnL(c(10,100),lam=0,a=1,b=0,H1=H1,H2=H2,v=v)
## No regularization 
fit1 <- nlm(lnL,c(10,100),lam=0,a=1,b=0,H1=H1,H2=H2,v=v)
fit1$estimate
GKL1 <- fit1$minimum
points(1:20,fit1$estimate[1]*H1+fit1$estimate[2]*H2,pch=19,col="red")
## Regularize first weight
fit2 <- nlm(lnL,c(10,100),lam=0.1,a=1,b=0,H1=H1,H2=H2,v=v)
fit2$estimate
GKL2 <- lnL(fit2$estimate,lam=0,a=1,b=0,H1=H1,H2=H2,v=v)
points(1:20,fit2$estimate[1]*H1+fit2$estimate[2]*H2,pch=19,col="blue")
## Regularize both weights
fit3 <- nlm(lnL,c(10,100),lam=0.1,a=1,b=0.1,H1=H1,H2=H2,v=v)
fit3$estimate
GKL3 <- lnL(fit3$estimate,lam=0,a=1,b=0,H1=H1,H2=H2,v=v)
points(1:20,fit3$estimate[1]*H1+fit3$estimate[2]*H2,pch=19,col="green")
## Severe regularization of first weight
fit4 <- nlm(lnL,c(10,100),lam=0.5,a=1,b=0,H1=H1,H2=H2,v=v)
fit4$estimate
GKL4 <- lnL(fit4$estimate,lam=0,a=1,b=0,H1=H1,H2=H2,v=v)
points(1:20,fit4$estimate[1]*H1+fit4$estimate[2]*H2,pch=19,col="cyan")
## Add legend
legend("topleft",c("True Poisson rates","Observed mutation counts",
                   "No regularization",
                   "Regularization on first weight",
                   "Regularization on both weights",
                   "Severe regularization on first weight"),
       col=c("black","purple","red","blue","green","cyan"),
       lty=1,lwd=3,bty="n",cex=1.1)
##-----------------------------------------------------------------------
## Summary
c(GKL1,GKL2,GKL3,GKL4)
c(fit1$estimate[1],fit2$estimate[1],fit3$estimate[1],fit4$estimate[1])
c(fit1$estimate[2],fit2$estimate[2],fit3$estimate[2],fit4$estimate[2])
cat("No regularization. Expected number of mutations:",
    sum(fit1$estimate[1]*H1+fit1$estimate[2]*H2),"\n")
cat("First weight reg. Expected number of mutations:",
    sum(fit2$estimate[1]*H1+fit2$estimate[2]*H2),"\n")
cat("Both weights reg. Expected number of mutations:",
    sum(fit3$estimate[1]*H1+fit3$estimate[2]*H2),"\n")
cat("First weight hard reg. Expected number of mutations:",
    sum(fit4$estimate[1]*H1+fit4$estimate[2]*H2),"\n")
##-----------------------------------------------------------------------

