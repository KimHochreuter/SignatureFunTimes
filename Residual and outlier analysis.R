######################################################
#RESIDUAL PLOT
######################################################
NMF_final = nmf(M, rank = 6, nrun = 10)

H = basis(NMF_final)

W = coef(NMF_final)

diff = M - H%*%W
diff
plot(M[1,],diff[1,])N
plot(M[3,], diff[3,])
plot(H%*%W,diff, xlim = c(0,1000),  ylim=c(-500, 500))