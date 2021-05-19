source("CV_function_version.R")
load("BRCA21.RData")
load("data/patients.rda")
################################################################################
##
##  Initialize dataset
##
################################################################################

dataset = V #Pick dataset
if(dim(dataset)[1] != 96){
  dataset = t(dataset)
}

################################################################################
##
##  Cross - Validation
##
################################################################################
#porun = CVPO(patients, K = 20)
#nbrun = CVNB(patients, K = 20)
porun1 = CVPO(dataset, K = 20)
nbrun1 = CVNB(dataset, K = 20)
porun2 = CVPO(dataset, K = 20, n_cv_sets = 20, n_updates = 10)
nbrun2 = CVNB(dataset, K = 20, n_cv_sets = 20, n_updates = 10)


##------------------------------------------------------------------------------
## POISSON MSE/BIC PLOT
#poruna = data.frame(porun)
poruna = data.frame(porun2)
po_plot_df = poruna[poruna$n_update == 10,]

d = po_plot_df %>%
  group_by(K) %>% 
  #summarise(medMSE = median(MSE), medBIC = median(BIC))
  {.}
ggplot(d) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = MSE))
ggplot(d) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = BIC))

##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL MSE/BIC PLOT
#nbruna = data.frame(nbrun)
nbruna = data.frame(nbrun2)
nb_plot_df = nbruna[nbruna$n_update == 10,]
g = nb_plot_df %>%
  group_by(K) %>% {.}
#summarise(medMSE = median(MSE), medBIC = median(BIC))
ggplot(g) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = MSE)) + ylim(0,10000)
ggplot(g) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = BIC))
ggplot(g) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = alpha))

##------------------------------------------------------------------------------


################################################################################
##
##  Residuals
##
################################################################################

data = t(dataset)
poNMF = nmf(data, rank = 4, nrun = 10, method = "KL")
H_po = coef(poNMF) # coef = H, in X = WH
W_po = basis(poNMF) # basis = W
alpha = 56
nbNMF = NMFNBMMsquarem(data, 6, alpha)
H_nb = nbNMF$P
W_nb = nbNMF$E

diff_po = data - H_po%*%W_po
diff_nb = data - H_nb%*%W_nb
z = data.frame(data = as.vector(data), diff_po = as.vector(diff_po), as.vector(diff_nb))


##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL MSE/BIC PLOT
p1 = (ggplot(z) + geom_point(aes(data, diff_po)) 
     + xlim(c(0,500))
     + ylim(c(-200,200))
     + geom_function(fun = function(x) 2*sqrt(x), aes(colour = "Poisson"))
     + geom_function(fun = function(x) -2*sqrt(x), aes(colour = "Poisson"))
     + geom_function(fun = function(x) 2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial"))
     + geom_function(fun = function(x) -2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial"))
     + ggtitle("Poisson residuals"))
p1

##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL RESIDUAL PLOT
p2 = (ggplot(z) + geom_point(aes(data, diff_nb)) 
     + xlim(c(0,500)) 
     + ylim(c(-200,200))
     + geom_function(fun = function(x) 2*sqrt(x), aes(colour = "Poisson"))
     + geom_function(fun = function(x) -2*sqrt(x), aes(colour = "Poisson"))
     + geom_function(fun = function(x) 2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial"))
     + geom_function(fun = function(x) -2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial"))
     + ggtitle("Negative binomial residuals"))
p2
##------------------------------------------------------------------------------




################################################################################
##
##  Signature/Exposure plots
##
################################################################################

best_k = 4

NMF_final = nmf(dataset, rank = best_k, nrun = 10)
NMF_final_scaled = scale(NMF_final)


##------------------------------------------------------------------------------
## POISSON SIGNATURES

W = basis(NMF_final_scaled)
W_df = data.frame(W)
colnames(W_df) = paste("s",1:best_k,sep="")
W_df$MutationType = rownames(W_df)
W_df$muta2 = str_sub(W_df$MutationType, 3, -3)

W_df = pivot_longer(W_df, s1:s4)
colnames(W_df)[c(3,4)] = c("Signature", "Intensity")
( ggplot(W_df) 
  + geom_col(aes(x = MutationType, y = Intensity, fill = muta2))
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
  + facet_grid(vars(Signature), vars(muta2), scales="free_x") )
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
( ggplot(W_df) 
  + geom_col(aes(x = MutationType, y = s3, fill = muta2))
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
  + facet_grid(~muta2,scales="free_x") )

( ggplot(W_df) 
  + geom_col(aes(x = MutationType, y = s4, fill = muta2))
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
  + facet_grid(~muta2,scales="free_x") )


##------------------------------------------------------------------------------
## POISSON EXPOSURES
H = coef(NMF_final_scaled)
H_df = data.frame(H)
#colnames(H_df) = paste("p",1:dim(H)[2],sep="")
H_df$Signature = paste("s",rownames(H_df),sep="")
H_df = pivot_longer(H_df, colnames(H_df)[1]:colnames(H_df)[dim(H)[2]])
colnames(H_df)[c(2,3)] = c("Patient", "Exposure")



( ggplot(H_df) 
  + geom_bar(aes(x = Patient, y = Exposure, fill = Signature), position = "fill", stat = "identity")
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90), legend.position = "none"))


##------------------------------------------------------------------------------
## POISSON MUTATIONAL PROFILE'

muta_profile = colSums(dataset)
muta_profile = data.frame(muta_profile)
muta_profile$Patient = rownames(muta_profile)
colnames(muta_profile)[1] = "real_count"

H_df_total_count = H_df %>% group_by(Patient) %>% summarize(count = sum(Exposure))
H_df_total_count = merge(H_df_total_count, muta_profile, by = "Patient")


( ggplot(H_df_total_count) 
  + geom_bar(aes(x = Patient, y = count), stat = "identity")
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
  + geom_point(aes(x = Patient, y = real_count, color = "red"), shape = 4, size = 5, stroke = 2) )



##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL SIGNATURES



##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL EXPOSURES



##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL MUTATIONAL PROFILE




################################################################################
##
##  Rodebunke
##
################################################################################


##------------------------------------------------------------------------------
## 
plot(apply(dataset, 2, mean), apply(dataset, 2, sd), xlim = c(0,150), ylim = c(0,200))
abline(a=0,b=1)

plot(apply(patients, 2, mean), apply(patients, 2, sd), xlim = c(0,250), ylim = c(0,500))
abline(a=0,b=1)

##------------------------------------------------------------------------------


