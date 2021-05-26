library(ggpubr)
library(ggforce)
library(lsa)
source("CV_function_version.R")
load("data/BRCA21.RData")
load("data/patients.rda")
load("data/Liver326.RData")
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
porun = CVPO(V, K = 30)
nbrun = CVNB(V, K = 30)
porun1 = CVPO(Liver, K = 30)
nbrun1 = CVNB(Liver, K = 30)
#porun2 = CVPO(dataset, K = 20, n_cv_sets = 20, n_updates = 10)
#nbrun2 = CVNB(dataset, K = 20, n_cv_sets = 20, n_updates = 10)

#dataset = V
#porun3 = CVPO(dataset, K = 20)
#nbrun3 = CVNB(dataset, K = 20)

##------------------------------------------------------------------------------
## POISSON MSE/BIC PLOT
#poruna = data.frame(porun)
poruna = data.frame(porund1_ida)
po_plot_df = poruna[poruna$n_update == 5,]

d = po_plot_df %>%
  group_by(K) %>% 
  #summarise(medMSE = median(MSE), medBIC = median(BIC))
  {.}
p_MSE_PO = ( ggplot(d) 
       + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = MSE)) 
       + ggtitle("BRCA21 Poisson MSE")
       + xlab("Number of mutational signatures")
       + theme_bw() )
p_DKL_PO = ( ggplot(d) 
       + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = DKL)) 
       + ggtitle("BRCA21 Poisson MSE")
       + xlab("Number of mutational signatures")
       + theme_bw() )
p_BIC_PO = (ggplot(d) 
      + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = BIC)) 
      + ggtitle("BRCA21 Poisson BIC")
      + xlab("Number of mutational signatures")
      + theme_bw() )

##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL MSE/BIC PLOT
#nbruna = data.frame(nbrun)
nbruna = data.frame(nbrund1_ida)
nb_plot_df = nbruna[nbruna$n_update == 5,]
g = nb_plot_df %>%
  group_by(K) %>% {.}
#summarise(medMSE = median(MSE), medBIC = median(BIC))
p_MSE_NB = ( ggplot(g) 
             + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = MSE))
             + xlab("Number of mutational signatures")
             + ggtitle("BRCA21 Negative Binomial MSE")
             + theme_bw()
             #+ ylim(0,6000)
              )
p_Da_NB = ( ggplot(g) 
             + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = D_alpha))
             + xlab("Number of mutational signatures")
             + ggtitle("BRCA21 Negative Binomial MSE")
             + theme_bw()
             #+ ylim(0,6000)
)
p_BIC_NB = (ggplot(g) + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = BIC)) 
            + xlab("Number of mutational signatures")
            + ggtitle("BRCA21 Negative Binomial BIC")
            + theme_bw()
            #+ ylim(-2050013, 0)
            )
p_ALPHA_NB = (ggplot(g) 
              + geom_boxplot(fill = "skyblue2", aes(x = factor(K), y = alpha))
              +theme_bw()
              + xlab("Number of mutational signatures")
              + ylab(expression(alpha))
              + ggtitle("BRCA21 Negative Binomial Dispersion Parameter"))

##------------------------------------------------------------------------------

mse = ggarrange(p_MSE_PO, p_MSE_NB)
dkl = ggarrange(p_DKL_PO, p_Da_NB)
bic = ggarrange(p_BIC_PO, p_BIC_NB)
alpha = p_ALPHA_NB# + ylim(0,5000)
ggsave(plot = mse,file = "pictures/BRCA21mse.png", width = 200, height = 105.83332, units = "mm")
ggsave(plot = bic,file = "pictures/BRCA21bic.png", width = 200, height = 105.83332, units = "mm")
ggsave(plot = alpha,file = "pictures/BRCA21alpha.png", width = 132.29165, height = 105.83332, units = "mm")
#ggsave(plot = mse,file = "pictures/PCAWGmse.png", width = 200, height = 105.83332, units = "mm")
#ggsave(plot = bic,file = "pictures/PCAWGbic.png", width = 200, height = 105.83332, units = "mm")
#ggsave(plot = alpha,file = "pictures/PCAWGalpha.png", width = 132.29165, height = 105.83332, units = "mm")

################################################################################
##
##  Residuals
##
################################################################################

#data = t(dataset)
Nsig = 4
poNMF = nmf(t(V), rank = Nsig, nrun = 10, method = "KL")
H_po = coef(poNMF) # coef = H, in X = WH
W_po = basis(poNMF) # basis = W

alpha = median(g[g$K == Nsig,]$alpha)
nbNMF = NMFNBMMsquarem(V, Nsig, alpha)
H_nb = nbNMF$P
W_nb = nbNMF$E

diff_po = t(V) - H_po%*%W_po
diff_nb = V - H_nb%*%W_nb
z = data.frame(data = as.vector(V), diff_po = as.vector(diff_po), as.vector(diff_nb))


##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL MSE/BIC PLOT
p1 = (ggplot(z) + geom_point(aes(data, diff_po)) 
     + xlim(c(0,500))
     + ylim(c(-200,200))
     + geom_function(fun = function(x) 2*sqrt(x), aes(colour = "Poisson"))
     + geom_function(fun = function(x) -2*sqrt(x), aes(colour = "Poisson"))
     + geom_function(fun = function(x) 2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial"))
     + geom_function(fun = function(x) -2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial"))
     + ggtitle("Poisson residuals")
     + ylab("Residuals")
     + xlab("Observed values")
     + theme_bw()
     + theme(legend.title = element_blank()) 
     )

##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL RESIDUAL PLOT
p2 = (ggplot(z) + geom_point(aes(data, diff_nb)) 
     + xlim(c(0,500)) 
     + ylim(c(-200,200))
     + geom_function(fun = function(x) 2*sqrt(x), aes(colour = "Poisson"))
     + geom_function(fun = function(x) -2*sqrt(x), aes(colour = "Poisson"))
     + geom_function(fun = function(x) 2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial"))
     + geom_function(fun = function(x) -2*sqrt(x + (x^2)/alpha), aes(colour = "Negative Binomial"))
     + ggtitle("Negative binomial residuals")
     + ylab("Residuals")
     + xlab("Observed values")
     + theme_bw()
     + theme(legend.title = element_blank()) 
     )

##------------------------------------------------------------------------------


( resi = ggarrange(p1, p2, common.legend = T) )
ggsave(plot = resi,file = "pictures/BRCA21residuals.png", width = 200, height = 105.83332, units = "mm")
#ggsave(plot = resi,file = "pictures/PCAWGresiduals.png", width = 200, height = 105.83332, units = "mm")

################################################################################
##
##  Signature/Exposure plots
##
################################################################################

Nsig = 4

NMF_final = nmf(V, rank = Nsig, nrun = 10)
NMF_final_scaled = scale(NMF_final)

NMF_NB_final = NMFNBMMsquarem(V, Nsig, median(g[g$K == Nsig,]$alpha), arrange = F)


##------------------------------------------------------------------------------
## POISSON SIGNATURES

H = basis(NMF_final_scaled)
colnames(H) = paste("s",1:Nsig,sep="")
H_df = data.frame(H)
colnames(H_df) = paste("s",1:Nsig,sep="")
H_df$MutationType = rownames(H_df)
H_df$muta2 = str_sub(H_df$MutationType, 3, -3)

H_df = pivot_longer(H_df, colnames(H_df)[1]:colnames(H_df)[Nsig])
colnames(H_df)[c(3,4)] = c("Signature", "Intensity")
(POsig = ( ggplot(H_df) 
  + geom_col(aes(x = MutationType, y = Intensity, fill = muta2))
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90, size = 5), legend.position = "none")
  + facet_grid(vars(Signature), vars(muta2), scales="free_x") ))

ggsave(plot = POsig,file = "pictures/BRCA21poSIG.png", width = 200, height = 105.83332, units = "mm")
#ggsave(plot = POsig,file = "pictures/PCAWGpoSIG.png", width = 200, height = 105.83332, units = "mm")

##------------------------------------------------------------------------------
## POISSON EXPOSURES
W = coef(NMF_final_scaled)
W_df = data.frame(W)
#colnames(H_df) = paste("p",1:dim(H)[2],sep="")
W_df$Signature = paste("s",rownames(W_df),sep="")
#H_df$Signature = OMEGANAME
W_df = pivot_longer(W_df, colnames(W)[1]:colnames(W)[dim(W)[2]])
colnames(W_df)[c(2,3)] = c("Patient", "Exposure")



(POexpo = ( ggplot(W_df) 
  + geom_bar(aes(x = Patient, y = Exposure, fill = Signature), position = "fill", stat = "identity")
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90)) ))

##------------------------------------------------------------------------------
## POISSON MUTATIONAL PROFILE

muta_profile = colSums(V)
muta_profile = data.frame(muta_profile)
muta_profile$Patient = rownames(muta_profile)
colnames(muta_profile)[1] = "real_count"

H_df_total_count = W_df %>% group_by(Patient) %>% summarize(count = sum(Exposure))
H_df_total_count = merge(H_df_total_count, muta_profile, by = "Patient")


(POcount = ( ggplot(H_df_total_count) 
  + geom_bar(aes(x = Patient, y = count), stat = "identity")
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90))
  + geom_point(aes(x = Patient, y = real_count, color = "Estimated count"), shape = 4, size = 2, stroke = 2)
  + theme(legend.title = element_blank()) ))



##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL SIGNATURES

H_NB = NMF_NB_final$P
rownames(H_NB) = rownames(V)
colnames(H_NB) = paste("s",1:Nsig,sep="")
colnames(H_NB) = SignaturePairing(Nsig, H_NB, H)

H_NB_df = data.frame(H_NB)
H_NB_df$MutationType = rownames(H_NB_df)
H_NB_df$muta2 = str_sub(H_NB_df$MutationType, 3, -3)

H_NB_df = pivot_longer(H_NB_df, colnames(H_NB)[1]:colnames(H_NB)[Nsig])
colnames(H_NB_df)[c(3,4)] = c("Signature", "Intensity")



(NBsig = ( ggplot(H_NB_df) 
  + geom_col(aes(x = MutationType, y = Intensity, fill = muta2))
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
  + facet_grid(vars(Signature), vars(muta2), scales="free_x") ))
ggsave(plot = NBsig,file = "pictures/BRCA21nbSIG.png", width = 200, height = 105.83332, units = "mm")

##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL EXPOSURES
W_NB = NMF_NB_final$E
W_NB_df = data.frame(W_NB)
#colnames(H_df) = paste("p",1:dim(H)[2],sep="")
colnames(W_NB) = colnames(V)
colnames(W_NB_df) = colnames(V)
W_NB_df$Signature = paste("s",rownames(W_NB_df),sep="")
W_NB_df = pivot_longer(W_NB_df, colnames(W_NB_df)[1]:colnames(W_NB_df)[dim(W_NB)[2]])
colnames(W_NB_df)[c(2,3)] = c("Patient", "Exposure")



(NBexpo = ( ggplot(W_NB_df) 
  + geom_bar(aes(x = Patient, y = Exposure, fill = Signature), position = "fill", stat = "identity")
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90)) )
  + theme(legend.title = element_blank()))


##------------------------------------------------------------------------------
## NEGATIVE BINOMIAL MUTATIONAL PROFILE
muta_profile = colSums(V)
muta_profile = data.frame(muta_profile)
muta_profile$Patient = rownames(muta_profile)
colnames(muta_profile)[1] = "real_count"


W_NB_df_total_count = W_NB_df %>% group_by(Patient) %>% summarize(count = sum(Exposure))
#W_NB_df_total_count$Patient = colnames(V)
W_NB_df_total_count = merge(W_NB_df_total_count, muta_profile, by = "Patient")


(NBcount = ( ggplot(W_NB_df_total_count) 
  + geom_bar(aes(x = Patient, y = count), stat = "identity")
  #+ geom_col(aes(x = MutationType, y = s2, fill = muta2))
  + theme_bw() 
  + theme(axis.text.x = element_text(angle = 90))
  + geom_point(aes(x = Patient, y = real_count, color = "Estimated count"), shape = 4, size = 2, stroke = 2) ))



##------------------------------------------------------------------------------

( BRCA21expo = ggarrange(POexpo, NBexpo, common.legend = TRUE) )
ggsave(plot = BRCA21expo, file = "pictures/BRCA21expo.png", width = 200, height = 105.83332, units = "mm")
#ggsave(plot = BRCA21expo, file = "pictures/PCAWGpoSIG.png", width = 200, height = 105.83332, units = "mm")

( BRCA21count = ggarrange(POcount, NBcount, common.legend = T) )
ggsave(plot = BRCA21count, file = "pictures/BRCA21count.png", width = 200, height = 105.83332, units = "mm")
#ggsave(plot = BRCA21expo, file = "pictures/PCAWGpoSIG.png", width = 200, height = 105.83332, units = "mm")

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



colnames(H) = SignaturePairing(4,H_NB,H)


NB_PO_sig_comparison(Nsig, H_NB, H)
ggplot(NB_PO_sig_comparison(Nsig, H_NB, H)) + geom_bar(aes(x = Signature, y = CosineSim), stat = "identity") + ylim(c(0,1))

Cosmic_comparison(Nsig, H_NB, H, cosmic)
(ggplot(Cosmic_comparison(Nsig, H_NB, H, cosmic)) 
  + geom_bar(aes(x =`Cosmic Signature`, y = Similarity, fill = Distribution), stat = "identity", position = position_dodge()) 
  + facet_grid(cols = vars(Signature),  scales="free_x") + ylim(c(0,1))
  #+ geom_text(aes(label=Similarity))
  )
  







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

COSMIC$SBS2 == cosmic$SBS2
par(mfrow = c(1,2))
barplot(cosmic$SBS2)
barplot(COSMIC$SBS2)
barplot(H[,2])

sigmixture = (COSMIC$SBS2 + COSMIC$SBS13)/sum(COSMIC$SBS2 + COSMIC$SBS13)
cosine(H[,2], sigmixture)
barplot(sigmixture)
par(mfrow = c(1,1))
barplot(cosmic$SBS2)
barplot(cosmic$SBS13)
barplot(H[,2])

dummy1 = data.frame(H)
dummy1$type = rownames(H)
dummy1 = dummy1 %>% arrange(type)

dummy2 = cosmic %>% arrange(Type)

dummy1$type == dummy2$Type


sigmixture = (dummy2$SBS2 + dummy2$SBS13)/sum(dummy2$SBS2 + dummy2$SBS13)
cosine(dummy1[,2], sigmixture)

par(mfrow = c(2,2))
barplot(dummy2$SBS2)
barplot(dummy2$SBS13)
barplot(dummy1[,3])
barplot(sigmixture)
