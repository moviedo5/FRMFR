
#library(goffda)
data(ontario)
# data("aemet_temp")
# aemet_temp$df


library(fda.usc.devel)
library(refund)
library(FRegSigCom)
df=ontario$df
ldat <- ldata(df,y=ontario$elec,x=fdata(ontario$temp$data[,25:49],0:24))
ldat$x1 <- fdata.deriv(ldat$x)
tt <- ldat$x$argvals
tty <- ldat$y$argvals
nbasis.x=11; nbasis.y=11
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis.y=create.bspline.basis(rangeval=range(tty),nbasis=nbasis.y)
basis.x=list("x"=basis1,"x1"=basis1)
npcx <- 5;npcy <- 2

s.n.basis <- nbasis.x;t.n.basis <- nbasis.y;x.n.basis <- nbasis.x

nk <- nbasis.y; nkk <- c(nbasis.x,nbasis.y)

nrep  <- 25
rr <- matrix(NA,nrep,10)
i <- 1
for (i in 1:nrep){
 print(i)
 set.seed(i*42586)
 ii <-sample(nrow(df),floor(nrow(df)*.5))
 ldat2<-ldat[ii,row=T]
 ldat3<-ldat[-ii,row=T]
 rest <- sum(norm.fdata(ldat3$y-func.mean(ldat2$y))^2)

# res0=fregre.lm.fr(y~x,ldat2
#                   basis.x=basis.x,  #ojo problema con basis.x spline
#                   basis.y=basis.y)
# basis.x=list("x"=basis1)
#basis.y=create.bspline.basis(rangeval=range(tty),nbasis=nbasis.y)
 basis.x <- list("x"=create.pc.basis(ldat2$x,1:4))
 basis.y <- create.pc.basis(ldat2$y,1:4)
 #~summary(basis.x$x); summary(basis.y)
 
 # res0=fregre.lm.fr(y~x,ldat2,  basis.x=basis.x )#ok
 # pred0<-predict.fregre.lm.fr(res0,ldat3)
 # sum(norm.fdata(ldat3$y-pred0)^2)/rest
 # 
 # bsp=create.bspline.basis(rangeval=range(tty),nbasis=nbasis.y)
 # res0=fregre.lm.fr(y~x,ldat2,  basis.x=basis.x, basis.y=bsp) #OK
 # 
 # pred0<-predict.fregre.lm.fr(res0,ldat3)
 # sum(norm.fdata(ldat3$y-pred0)^2)/rest
 # plot(pred0,col=1); lines(pred1,col=2)
 # basis.y=create.pc.basis(ldat2$x,1:10)
 # res0=fregre.lm.fr(y~x,ldat2,  basis.x=basis.x, basis.y=basis.y)
 # pred0<-predict(res0,ldat3)
 # sum(norm.fdata(ldat3$y-pred0)^2)/rest
 # plot(pred0,col=1); lines(pred1,col=2)
 # res0=fregre.lm.fr(y~x,ldat2,  basis.x=basis.x
 #                   #, basis.y=basis.y
 #                  )
# pred0<-predict(res0,ldat3)
# sum(norm.fdata(ldat3$y-pred0)^2)/rest;sum(norm.fdata(ldat3$y-pred1)^2)/rest
# plot(pred0,col=1);lines(pred1,col=2)

res1=fregre.mlm.fr(y~x,ldat2,  basis.x=basis.x,
                   basis.y=basis.y)
pred1<-predict(res1,ldat3)
# sum(norm.fdata(ldat3$y-pred1)^2)/rest

res2=fregre.sam.fr(y~s(x),ldat2,basis.x=basis.x,
                   basis.y=basis.y)
pred2<-predict(res2,ldat3)
rr[i,1] <- sum(norm.fdata(ldat3$y-pred1)^2)/rest
rr[i,2] <- sum(norm.fdata(ldat3$y-pred2)^2)/rest

basis.x<-list("x"=create.pc.basis(ldat2$x,1:npcx))
basis.y<-create.pc.basis(ldat2$y,1:npcy)

res1=fregre.mlm.fr(y~x,ldat2,  basis.x=basis.x,basis.y=basis.y)
pred1<-predict(res1,ldat3)
res2=fregre.sam.fr(y~s(x),ldat2,basis.x=basis.x,basis.y=basis.y)
pred2<-predict(res2,ldat3)
rr[i,3] <- sum(norm.fdata(ldat3$y-pred1)^2)/rest
rr[i,4] <- sum(norm.fdata(ldat3$y-pred2)^2)/rest

res1=fregre.kam.fr(y~x,ldat2)
pred1<-predict(res1,ldat3)
rr[i,5] <- sum(norm.fdata(ldat3$y-pred1)^2)/rest

#par.np=list("x"=list("type.S"="S.KNN","Ker"=AKer.unif,"h"=7))
par.np<-list()
par.metric=list(x=list(metric=metric.lp,lp=1))
res1=fregre.kam.fr(y~x,ldat2,par.np=par.np,par.metric=par.metric)
pred1<-predict(res1,ldat3)
rr[i,6] <- sum(norm.fdata(ldat3$y-pred1)^2)/rest

ydat<-ldat2$y$data
xdat <-ldat2$x$data
T<-ncol(ydat)
S<-ncol(xdat)
tj<-seq(1,T,len=T)
si<-seq(1,S,len=S)

respff1=pffr(ydat~ sff(xdat,xind=si,
      splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk)),
             yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
             #bs.int=list(bs="ps",k=21,m=c(2,1))
             )
prff1=fdata(predict(respff1,list(xdat=ldat3$x$data)),argvals=tty)
rr[i,7] <- sum(norm.fdata(ldat3$y-prff1)^2)/rest

respff1pc=pffr(ydat~ffpc(xdat,xind=si,npc.max=npcx),yind=tj)
prff1.pc=fdata(predict(respff1pc,list(xdat=ldat3$x$data)),argvals=tty)
rr[i,8] <- sum(norm.fdata(ldat3$y-prff1.pc)^2)/rest

#install.packages("D:/Users/moviedo/OneDrive - Universidade da Coruña/fda.usc/Morteza/FRegSigCom_0.3.0.tar.gz", repos = NULL, type = "source")
#library(FRegSigCom)
XX=list(X1=xdat)
XXnew=list(X1=ldat3$x$data)
t.x.list=list(X1=si)

#----- ---- ---- ----     cv.nonlinear
resSC.NL1=cv.nonlinear(XX, ydat, t.x.list, 
                       tj,s.n.basis=s.n.basis,
                       x.n.basis=x.n.basis,
                       t.n.basis=t.n.basis)
prSCNL1=fdata(pred.nonlinear(resSC.NL1,XXnew),tty)
rr[i,9] <- sum(norm.fdata(ldat3$y-prSCNL1)^2)/rest
#----- ---- ---- ----     cv.sigcom
resSC.Lin1=cv.sigcom(XX, ydat, t.x.list,
                     tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
prSCLin1=fdata(pred.sigcom(resSC.Lin1,XXnew),tty)
rr[i,10] <- sum(norm.fdata(ldat3$y-prSCLin1)^2)/rest
#--
print(rr[i,])
print(which.min(rr[i,]))
}


par(mfrow=c(1,1))
#colnames(rr)<-c("LM_bspx21bspy11","SAM_bspx21bspy11",
colnames(rr)<-c("LM_PCx10PCy4","SAM_PCx10PCy4",
                "LM_PCx5PCy2","SAM_PCx5PCy2",
                "KAM_L2","KAM_L1",
                "pffr31_11","pffr_PC5",
                "cv.nonlinear11_11","cv.sigcom11_11"
                )
boxplot(rr)
range(ldat2$df$date)
100-round(colMeans(rr)*100,1)
100-round(apply(rr,2,median,na.rm=T)*100,1)
round(apply(rr,2,sd),4)

write.table(rr,"out.txt")
[1] 0.1971852 0.1545276 0.2129331 0.1495876 0.1589951 0.1570534
[7] 0.1646810 0.1955769 0.1630982 0.1916253
