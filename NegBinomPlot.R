par(mfrow = c(1,2)) 

alph <- c(0.1,1,10,100)
#x <- seq(0,10, length.out = 1000)
df1 <- cbind(1:10, rep(0, 10),dpois(1:10, lambda = 5))
for(i in 1:4){
  df1 <- rbind(df1, cbind(1:10, rep(alph[i],10),dnbinom(1:10, size = alph[i], prob = (1-5/(5+alph[i])))))
}
colnames(df1) <- c("x", "alpha", "val")
df1 <- as.data.frame(df1)%>%
  mutate(alpha = as.character(alpha))
greeks <- list(bquote(alpha %->% infinity) ,bquote(alpha==.(0.1)), 
               bquote(alpha==.(1)), bquote(alpha==.(10)), bquote(alpha==.(100)))
NBp <- ggplot(data = df1, aes(x=x, y=val, color = alpha))+geom_line() +geom_point() +
  scale_colour_manual(values = 1:5, labels=greeks) +
  theme(legend.title=element_blank(), legend.position="top", legend.text = element_text(size = 14),
        text = element_text(size=14)) + 
  labs(x = "y",y = expression(paste("NB( y | ", alpha, ")")))


D_alpha <- function(alp, y = 3, l = 5){
  return(y*log(y/l) - (alp + y)*log((alp+y)/(alp+l)))
}
D_KL <- function(y = 3, l= 5){
  return(y*log(y/l)- y + l)
}

seq <- seq(0,4, length.out = 10000)
df <- cbind(seq, rep(0, 10000),D_KL(l = seq))
for(i in 1:4){
  df <- rbind(df, cbind(seq, rep(alph[i],10000),D_alpha(alp = alph[i], l = seq)))
}
colnames(df) <- c("x", "alpha", "val")
df <- as.data.frame(df)%>%
  mutate(alpha = as.character(alpha))
#greeks <- list(bquote(alpha %->% infinity) ,bquote(alpha==.(0.1)), bquote(alpha==.(1)), bquote(alpha==.(10)), bquote(alpha==.(100)))
DP <- ggplot(data = df, aes(x=x, y=val, color = alpha)) + ylim(0,2.5) + geom_line() +
  scale_colour_manual(values = 1:5, labels=greeks) +
  theme(legend.title=element_blank(), legend.position="top", text = element_text(size=14)) + 
  labs(x = expression(lambda), y = expression(paste(D[alpha], "(y = 3 | ", lambda, ")"))) 


#install.packages("ggpubr")
library(ggpubr)
ggarrange(NBp, DP + theme(legend.position = "none"), 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

