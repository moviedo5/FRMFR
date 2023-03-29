# ver en  wavelets pkg si sale la misma indep de la numero de nivel!!!
# existen diferencias entre n y 2^J
# el n.levels no va bien, hacer el idwt? o coger el mra y modifcar para acepte matrices
# waveslim:::dwt.C==wavelets:::dwt.forward

library(fda.usc)
#X1 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,.7,.9,0,.3)
#X2 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,-.7,.9,0,.3)
# newX <- fdata(rbind(X1,X2))
mra.fdata <- fda.usc.devel:::mra.fdata 
#source("mra.fdata.R")
# combine them and compute MRA
par(mfrow=c(1,1))

tt<- seq(0,1,len=64)
newX <- rproc2fdata(100,t=tt,mu=sin(2*pi*tt))
nlev=5
mra.out <- mra.fdata(newX, n.levels=nlev)
class(mra.out)
names(mra.out)
plot(newX,col=1)
lines(mra.out$fdata.est,col=nlev+1)

#####################
data(aemet)
x <- aemet$temp
res2 <- mra.fdata(x, n.levels=4)

res3 <- mra.fdata(x, n.levels=3)

plot(x[1:5],col=1,lwd=2)
lines(res2$fdata.est[1:5],col=2,lwd=2)
lines(res3$fdata.est[1:5],col=3,lwd=2)

# #####################
# X1 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,.7,.9)
# X2 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,-.7,.9)
# newX <- fdata(rbind(X1,X2))
# #mra.out <- mra.fdata(newX, n.levels=4)
# x<-1:length(X1)/length(X1)
# y<-X1
# x0<-approx(x,y,n=16)
# fx<-fdata(x0$y,x0$x)
# fx0<-fdata(y,x)
# plot(fx)
# lines(fx0,col=2)
# 
# plot(newX)
# #serie completa (da igual el n.level)
# #lines(mra.out$S$S1+mra.out$D$D1,col=4:5)
# lines(mra.out$S$S2+mra.out$D$D2,col=4:5)
# 
# #library(fda.usc)
# xx<-seq(0,1,len=64)
# x1<-sin(2*pi*xx)+rnorm(xx,sd=.2)
# x2<-cos(2*pi*xx)+rnorm(xx,sd=.2)
# newX<-fdata(rbind(x1,x2))
# mra.out <- mra.fdata(newX, n.levels=3)
# plot(newX)
# lines(mra.out$S$S2+mra.out$D$D2,col=4:5)
# 
# plot(newX-mra.out$S$S1-mra.out$D$D1)
# plot(newX-mra.out$S$S2-mra.out$D$D2-mra.out$D$D3)
# 
# 
# #fit<-y$S$S3+y$S$S2+y$S$S1
# y<- mra(x12,"d4", 3)
# fit<-+y$D$D1+y$D$D2+y$D$D3
# #fit<-y$S$S3+y$S$S2+y$S$S1
# matplot(x12,type="l")
# matplot(fit,type="l",add=TRUE,col=3:4,lwd=2)
# 
# 
# y<- mra(matrix(x,1), "haar", 3)
# y<- mra(matrix(x,ncol=1), "haar", 3)
# 
# # obtain the two series listed in Percival and Walden (2000), page 42
# X1 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,.7,.9,0,.3)
# X2 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,-.7,.9,0,.3)
# 
# # combine them and compute MRA
# newX <- cbind(X1,X2)
# mra.out <- mra(newX, n.levels=3)
# matplot(newX,type="l")
# #aa<-matrix(mra.out@S$S1,ncol=2,nrow=16)
# matplot(mra.out$S$S1,type="l",add=TRUE,col=5:6,lwd=4)
# 
# 
# 
# # obtain the two series listed in Percival and Walden (2000), page 42
# X1 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,.7,.9,0,.3)
# X2 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,-.7,.9,0,.3)
# # combine them and compute MRA
# newX <- cbind(X1,X2)
# mra.out <- mra(newX, n.levels=3,fast=FALSE)
# matplot(newX,type="l")
# aa<-matrix(mra.out@S$S1,ncol=2,nrow=16)
# matplot(aa,type="l",add=TRUE,col=5:6,lwd=4)

