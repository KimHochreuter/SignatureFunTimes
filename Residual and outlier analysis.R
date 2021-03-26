######################################################
#RESIDUAL PLOT
######################################################
NMF_final = nmf(M, rank = 6, nrun = 10)

H = basis(NMF_final)
W = coef(NMF_final)

diff = M - H%*%W
plot(M, diff, xlim = c(0,1000), ylim=c(-500, 500))