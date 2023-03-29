library(fda.usc.devel)
library(refund)
library(FRegSigCom)
names(omel)
df=omel$df
nn<-nrow(df)
i1 <-which(omel$df$fechas=="2011-01-01")
i2 <-which(omel$df$fechas=="2012-12-31")
i1 <-which(omel$df$fechas==omel$df$fechas[1])
i1 <-which(omel$df$fechas=="2010-01-01")
i1 <-which(omel$df$fechas=="2013-01-01")
i2 <-which(omel$df$fechas==max(omel$df$fechas))-1
i1<-1
i2 <-which(omel$df$fechas=="2010-01-01")-1
#i2 <-500
length(i1:i2)
ldat <- ldata(df[i1:i2,],y=omel$precio[i1:i2],x=omel$energia[i1:i2])
par(mfrow=c(1,2))
ts.plot(ldat$x$data[,1])
ts.plot(ldat$y$data[,1])
#ldat$x1 <- fdata.deriv(ldat$x)
tt <- ldat$x$argvals
tty <- ldat$y$argvals
nbasis.x=11; nbasis.y=11
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis.y=create.bspline.basis(rangeval=range(tty),nbasis=nbasis.y)
basis.x=list("x"=basis1,"x1"=basis1)
npcx <- 5;npcy <- 3

s.n.basis <- nbasis.x;t.n.basis <- nbasis.y;x.n.basis <- nbasis.x

nk <- nbasis.y; nkk <- c(nbasis.x,nbasis.y)

nrep  <- 100
rr <- matrix(NA,nrep,10)
i <- 1
for (i in 1:nrep){
  print(i)
  set.seed(i*42586)
  ii <-sample(nrow(ldat$df),floor(nrow(ldat$df)*.75))
  ldat2<-ldat[ii,row=T]
  ldat3<-ldat[-ii,row=T]
  rest <- sum(norm.fdata(ldat3$y-func.mean(ldat2$y))^2)
  
  # res0=fregre.lm.fr(y~x,ldat2
  #                   basis.x=basis.x,  #ojo problema con basis.x spline
  #                   basis.y=basis.y)
  # basis.x=list("x"=basis1)
  #basis.y=create.bspline.basis(rangeval=range(tty),nbasis=nbasis.y)
  basis.x<-list("x"=create.pc.basis(ldat2$x,1:7))
  basis.y<-create.pc.basis(ldat2$y,1:7)
  
  res1=fregre.mlm.fr(y~x,ldat2,  basis.x=basis.x,
                     basis.y=basis.y)
  pred1<-predict(res1,ldat3)
  res2=fregre.sam.fr(y~s(x),ldat2,basis.x=basis.x,
                     basis.y=basis.y)
  pred2<-predict(res2,ldat3)
  rr[i,1] <- sum(norm.fdata(ldat3$y-pred1)^2)/rest
  rr[i,2] <- sum(norm.fdata(ldat3$y-pred2)^2)/rest
  
  # basis.x<-list("x"=create.pc.basis(ldat2$x,1:npcx))
  # basis.y<-create.pc.basis(ldat2$y,1:npcy)
  #basis.x=list("x"=basis1)
  # basis.y=create.bspline.basis(rangeval=range(tty),nbasis=nbasis.y)
  basis.x<-list("x"=create.pc.basis(ldat2$x,1:3))
  basis.y<-create.pc.basis(ldat2$y,1:3)
  #summary(basis.x$x,biplot=F);  summary(basis.y,biplot=F)
  
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
colnames(rr)<-c("LM_PCx7PCy7","SAM_PCx7PCy7",
                "LM_PCx3PCy3","SAM_PCx3PCy3",
                "KAM_L2","KAM_L1",
                "pffr31_11","pffr_PC3",
                "cv.nonlinear11_11","cv.sigcom11_11"
)
boxplot(rr)

nrow(ldat$df);range(ldat$df$fechas);nrep
100-round(colMeans(rr)*100,1)

100-round(apply(rr,2,median,na.rm=T)*100,1)
round(apply(rr,2,sd),4)

write.table(rr,"out.txt")
[1] 0.1971852 0.1545276 0.2129331 0.1495876 0.1589951 0.1570534
[7] 0.1646810 0.1955769 0.1630982 0.1916253

> nrow(ldat$df);range(ldat$df$fechas);nrep
[1] 500
[1] "2008-04-01" "2009-08-13"
[1] 100
> 100-round(colMeans(rr)*100,1)
LM_PCx7PCy7      SAM_PCx7PCy7       LM_PCx3PCy3 
48.2              63.2              19.5 
SAM_PCx3PCy3            KAM_L2            KAM_L1 
41.0              70.7              68.1 
pffr31_11          pffr_PC3 cv.nonlinear11_11 
43.7              22.2              56.8 
cv.sigcom11_11 
47.2
