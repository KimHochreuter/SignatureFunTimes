source("CV_function_version.R")
load("data/BRCA21.RData")
load("data/patients.rda")
load("data/Liver326.RData")
################################################################################
##
##  Initialize dataset
##
################################################################################

dataset = Liver #Pick dataset
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
ggplot(g) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = -BIC))
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

diff_po = data - W_po%*%H_po
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

NMF_NB_final = NMFNBMMsquarem(dataset, 4, median(g[g$K == 4,]$alpha), arrange = F)


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
#H_df$Signature = OMEGANAME


H_df = pivot_longer(H_df, colnames(H)[1]:colnames(H)[dim(H)[2]])
colnames(H_df)[c(2,3)] = c("Patient", "Exposure")



( ggplot(H_df) 
  + geom_bar(aes(x = Patient, y = Exposure, fill = Signature), position = "fill", stat = "identity")
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90)) )


##------------------------------------------------------------------------------
## POISSON MUTATIONAL PROFILE

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

H_NB = NMF_NB_final$P
rownames(H_NB) = rownames(dataset)
colnames(H_NB) = paste("s",1:best_k,sep="")


H_NB_df = data.frame(H_NB)
H_NB_df$MutationType = rownames(H_NB_df)
H_NB_df$muta2 = str_sub(H_NB_df$MutationType, 3, -3)

H_NB_df = pivot_longer(H_NB_df, colnames(H_NB)[1]:colnames(H_NB)[length(colnames(H_NB))])
colnames(H_NB_df)[c(3,4)] = c("Signature", "Intensity")
( ggplot(H_NB_df) 
  + geom_col(aes(x = MutationType, y = Intensity, fill = muta2))
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
  + facet_grid(vars(Signature), vars(muta2), scales="free_x") )


##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL EXPOSURES
W_NB = NMF_NB_final$E
W_NB_df = data.frame(W_NB)
#colnames(H_df) = paste("p",1:dim(H)[2],sep="")
colnames(W_NB) = colnames(dataset)
colnames(W_NB_df) = colnames(dataset)
W_NB_df$Signature = paste("s",rownames(W_NB_df),sep="")
W_NB_df = pivot_longer(W_NB_df, colnames(W_NB_df)[1]:colnames(W_NB_df)[dim(W_NB)[2]])
colnames(W_NB_df)[c(2,3)] = c("Patient", "Exposure")



( ggplot(W_NB_df) 
  + geom_bar(aes(x = Patient, y = Exposure, fill = Signature), position = "fill", stat = "identity")
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90)) )


##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL MUTATIONAL PROFILE
muta_profile = colSums(dataset)
muta_profile = data.frame(muta_profile)
muta_profile$Patient = rownames(muta_profile)
colnames(muta_profile)[1] = "real_count"

W_NB_df_total_count$Patient = colnames(dataset)
W_NB_df_total_count = W_NB_df %>% group_by(Patient) %>% summarize(count = sum(Exposure))
W_NB_df_total_count = merge(W_NB_df_total_count, muta_profile, by = "Patient")


( ggplot(W_NB_df_total_count) 
  + geom_bar(aes(x = Patient, y = count), stat = "identity")
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
  + geom_point(aes(x = Patient, y = real_count, color = "red"), shape = 4, size = 5, stroke = 2) )



##------------------------------------------------------------------------------



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



library(lsa)
asd = cosine(cbind(H_NB,W[,2]))
asd[asd == 1] = 0
names(which.max(asd[,length(asd[1,])]))

OMEGANAME = c()
for (i in 1:4) {
  asd = cosine(cbind(H_NB,W[,i]))
  asd[asd == 1] = 0
  OMEGANAME[i] = names(which.max(asd[,length(asd[1,])]))
}
OMEGANAME
colnames(W) = OMEGANAME
W
rownames(H) = OMEGANAME


library(readxl)
COSMIC_Mutational_Signatures_v3_1 <- read_excel("DATA/COSMIC_Mutational_Signatures_v3.1.xlsx")
View(COSMIC_Mutational_Signatures_v3_1)
cosmic = COSMIC_Mutational_Signatures_v3_1

library(readr)
COSMIC_v3_2_SBS_GRCh38 <- read_table2("DATA/COSMIC_v3.2_SBS_GRCh38.txt")
View(COSMIC_v3_2_SBS_GRCh38)
cosmic = COSMIC_v3_2_SBS_GRCh38

dumoggrim = as.matrix(cbind(W[,4], cosmic[4:75]))
dummy = cosine(dumoggrim)
barplot(dummy[1,])
abline(a=0.8,b=0)



sigmixture = (cosmic$SBS2 + cosmic$SBS13)/sum(cosmic$SBS2 + cosmic$SBS13)
cosine(W[,3], sigmixture)

par(mfrow = c(1,1))
barplot(cosmic$SBS2)
barplot(cosmic$SBS13)
barplot(W[,3])

dummy1 = data.frame(W)
dummy1$type = rownames(W)
dummy1 = dummy1 %>% arrange(type)

dummy2 = cosmic %>% arrange(Type)

dummy1$type == dummy2$Type


sigmixture = (dummy2$SBS2 + dummy2$SBS13)/sum(dummy2$SBS2 + dummy2$SBS13)
cosine(dummy1[,3], sigmixture)

par(mfrow = c(2,2))
barplot(dummy2$SBS2)
barplot(dummy2$SBS13)
barplot(dummy1[,3])
barplot(sigmixture)
