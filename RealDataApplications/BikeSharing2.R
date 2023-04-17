library(fda.usc.devel)
library(refund)
library(FRegSigCom)

data("BikeSharing")
names(BikeSharing)

ddays=unique(BikeSharing$df$date)
lesp=which(ddays %in% as.Date(c("2011-09-17","2012-01-21","2012-06-09")))

# pdf(file="bikesharing.pdf",width=10.67,height=6)
m <- cbind(c(1,3,5),c(2,4,5))
layout(m)

x1 <- BikeSharing$temp
x2 <- BikeSharing$humidity
x3 <- BikeSharing$windspeed
x4 <- BikeSharing$feeltemp
y <- BikeSharing$logNBCR


plot(x1,col=gray(.5))
lines(x1[lesp],lwd=2)
plot(x2,col=gray(.5))
lines(x2[lesp],lwd=2)
plot(x3,col=gray(.5))
lines(x3[lesp],lwd=2)
plot(x4,col=gray(.5))
lines(x4[lesp],lwd=2)
plot(y,col=gray(.5))
lines(y[lesp],lwd=2)
#dev.off()

# The plot of Bike-sharing variables over 112 Saturdays from January 1st, 2011 
# to December 1st, 2012. The black, red, and green lines correspond to dates 
# September 17th, 2011, January 21st, 2012, and January 9th, 2012, respectively.


##### Refund parameters 
kk=c(9,9,9) # Model has more coefficients than data
kk=c(7,9,9) # Model has more coefficients than data
mm <- c(2,2,2)
#kk <- 7
k0 <- 7
m0 <- c(2,1)
kpc <-  11

#SigCom parameters 
snb <- 15
tnb <- 15

# fda.usc parameters
npcy <- 6
npc1 <- 2
npc2 <- 4
npc3 <- 6
npc4 <- 2
tj <- seq(0.5,23.5)

nrep=20
# rr=matrix(NA,nrep,7)
# colnames(rr)=c("FLMFR","FSAMFR","FKAMFR","PFR","FAMM","LSC","DISC")


rr1 <- matrix(NA,ncol=7,nrow=nrep)
colnames(rr1) <- c("FLMFR","FSAMFR","FKAMFR","PFFR","FAMM","LSIG","DISC")
rr3 <- rr2 <- rr1

compModel <- rep(TRUE,len=7)
names(compModel) <- colnames(mtrain1)
compModel["FAMM"] <- FALSE # Time consuming
compModel <- !compModel  # Only FAMM

compModel <- !compModel  # Only FAMM

compModel[1:7]<-F
compModel["FAMM"]<-T
kk=c(7,9,9) # Model has more coefficients than data

set.seed(20030101)
for (i in 1:nrep){ # i<-1
  index=sample(1:102,82)
  
  yt=y[index]
  X1=x1[index]
  X2=x2[index]
  X3=x3[index]
  X4=x4[index]
  
  ytnew=y[-index]
  X1new=x1[-index]
  X2new=x2[-index]
  X3new=x3[-index]
  X4new=x4[-index]
  
  ldatos=ldata(X1=X1,X2=X2,X3=X3,X4=X4,yt=yt)
  ldatnew=ldata(X1=X1new,X2=X2new,X3=X3new,X4=X4new,yt=ytnew)
  
  rest=sum(norm.fdata(ldatnew$yt-func.mean(ldatos$yt))^2)
  
  ####### fda.usc.devel
  if (any(compModel[c("FLMFR","FSAMFR")])){
    b.x=list(X1=create.pc.basis(X1,1:npc1),X2=create.pc.basis(X2,1:npc2),
             X3=create.pc.basis(X3,1:npc3),X4=create.pc.basis(X4,1:npc4))
    b.y=create.pc.basis(yt,1:npcy)
    }
  
  
  if (compModel["FLMFR"]){
    form <- yt~X4
    reslin=fregre.mlm.fr(form,data=ldatos,basis.y=b.y,basis.x=b.x)
    prlin=predict(reslin,ldatnew)
    rr1[i,1]=1-sum(norm.fdata(ldatnew$yt-prlin)^2)/rest
    
    form <- yt~X4+X2
    reslin=fregre.mlm.fr(form,data=ldatos,basis.y=b.y,basis.x=b.x)
    prlin=predict(reslin,ldatnew)
    rr2[i,1]=1-sum(norm.fdata(ldatnew$yt-prlin)^2)/rest
    
    form <- yt~X4+X2+X3
    reslin=fregre.mlm.fr(form,data=ldatos,basis.y=b.y,basis.x=b.x)
    prlin=predict(reslin,ldatnew)
    rr3[i,1]=1-sum(norm.fdata(ldatnew$yt-prlin)^2)/rest
    }
  if (compModel["FSAMFR"]){
    forms <- yt~s(X4)
    ressam=fregre.sam.fr(forms,data=ldatos,basis.y=b.y,basis.x=b.x)
    prsam=predict(ressam,ldatnew)
    rr1[i,2]=1-sum(norm.fdata(ldatnew$yt-prsam)^2)/rest
    
    forms <- yt~s(X4)+s(X2,k=5)
    ressam=fregre.sam.fr(forms,data=ldatos,basis.y=b.y,basis.x=b.x)
    prsam=predict(ressam,ldatnew)
    rr2[i,2]=1-sum(norm.fdata(ldatnew$yt-prsam)^2)/rest
    
    forms <- yt~s(X4,k=5)+s(X2,k=5)+s(X3,k=5)
    ressam=fregre.sam.fr(forms,data=ldatos,basis.y=b.y,basis.x=b.x)
    prsam=predict(ressam,ldatnew)
    rr3[i,2]=1-sum(norm.fdata(ldatnew$yt-prsam)^2)/rest
    }
  if (compModel["FKAMFR"]){
    pmetric=list(df=data.frame(idx=1:nrow(yt)),yt=list(metric=metric.lp,lp=2), 
                 X1=list(metric=metric.lp,lp=2),X2=list(metric=metric.lp,lp=2),
                 X3=list(metric=metric.lp,lp=2),X4=list(metric=metric.lp,lp=2))
    p.np=list(X1=list(Ker=AKer.norm),X2=list(Ker=AKer.norm),
              X3=list(Ker=AKer.norm),X4=list(Ker=AKer.norm))
    
    form <- yt~X4
    reskam=fregre.kam.fr(form,data=ldatos,par.metric=pmetric)
    prkam=predict(reskam,ldatnew)
    rr1[i,3]=1-sum(norm.fdata(ldatnew$yt-prkam)^2)/rest
    
    form <- yt~X4+X2
    reskam=fregre.kam.fr(form,data=ldatos,par.metric=pmetric)
    prkam=predict(reskam,ldatnew)
    rr1[i,3]=1-sum(norm.fdata(ldatnew$yt-prkam)^2)/rest
    
    form <- yt~X4+X2+X3
    reskam=fregre.kam.fr(form,data=ldatos,par.metric=pmetric)
    prkam=predict(reskam,ldatnew)
    rr1[i,3]=1-sum(norm.fdata(ldatnew$yt-prkam)^2)/rest
  }
  
  Y1D=yt$data;X1D=X1$data;X2D=X2$data;X3D=X3$data;X4D=X4$data
  if (compModel["PFFR"]){
  respff1pc=pffr(Y1D~ffpc(X4D,xind=tj,npc.max=npc4,splinepars=list(k=kpc))
                         ,yind=tj)
  respff1pc.fit=fdata(matrix(respff1pc$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)
  prff1.pc=fdata(predict(respff1pc,list(X4D=X4new$data)),argvals=tj)
  rr1[i,4]=1-sum(norm.fdata(ldatnew$yt-prff1.pc)^2)/rest
  
  respff1pc=pffr(Y1D~ffpc(X4D,xind=tj,npc.max=npc4,splinepars=list(k=kpc))
                                   +ffpc(X2D,xind=tj,npc.max=npc2,splinepars=list(k=kpc))
                 ,yind=tj)
  respff1pc.fit=fdata(matrix(respff1pc$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)
  prff1.pc=fdata(predict(respff1pc,list(X4D=X4new$data,X2D=X2new$data)),argvals=tj)
  rr2[i,4]=1-sum(norm.fdata(ldatnew$yt-prff1.pc)^2)/rest
  
  respff1pc=pffr(Y1D~ffpc(X4D,xind=tj,npc.max=npc4,splinepars=list(k=kpc))
                                   +ffpc(X2D,xind=tj,npc.max=npc2,splinepars=list(k=kpc))
                                   +ffpc(X3D,xind=tj,npc.max=npc3,splinepars=list(k=kpc))
                 ,yind=tj)
  respff1pc.fit=fdata(matrix(respff1pc$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)
  prff1.pc=fdata(predict(respff1pc,list(X4D=X4new$data,X2D=X2new$data,X3D=X3new$data)),argvals=tj)
  rr1[i,4]=1-sum(norm.fdata(ldatnew$yt-prff1.pc)^2)/rest
  }
  
  if (compModel["FAMM"]){
    respff1=pffr(Y1D~sff(X4D,xind=tj,splinepars=list(bs="ps",m=mm,k=kk))
                   ,yind=tj,bs.yindex=list(bs="ps",m=m0,k=k0))
    
    respff1.fit=fdata(matrix(respff1$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)
    prff1=fdata(predict(respff1,list(X4D=X4new$data)),argvals=tj)
    rr1[i,5]=1-sum(norm.fdata(ldatnew$yt-prff1)^2)/rest
    
    respff1=pffr(Y1D~sff(X4D,xind=tj,splinepars=list(bs="ps",m=mm,k=kk))
                 #  ojo estaba así       +sff(X3D,xind=tj,splinepars=list(bs="ps",m=mm,k=kk))
                 +sff(X2D,xind=tj,splinepars=list(bs="ps",m=mm,k=kk))
                 ,yind=tj,bs.yindex=list(bs="ps",m=m0,k=k0))
    
    respff1.fit=fdata(matrix(respff1$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)
    prff1=fdata(predict(respff1,list(X4D=X4new$data,X2D=X2new$data)),argvals=tj)
    rr2[i,5]=1-sum(norm.fdata(ldatnew$yt-prff1)^2)/rest
    
    respff1=pffr(Y1D~sff(X4D,xind=tj,splinepars=list(bs="ps",m=mm,k=kk))
                   +sff(X2D,xind=tj,splinepars=list(bs="ps",m=mm,k=kk))
                          +sff(X3D,xind=tj,splinepars=list(bs="ps",m=mm,k=kk))
                 ,yind=tj,bs.yindex=list(bs="ps",m=m0,k=k0))
    
    respff1.fit=fdata(matrix(respff1$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)
    prff1=fdata(predict(respff1,list(X4D=X4new$data,X2D=X2new$data,X3D=X3new$data)),argvals=tj)
    rr3[i,5]=1-sum(norm.fdata(ldatnew$yt-prff1)^2)/rest
  }
  ########## FRegSigCom 
  if (compModel["LSIG"]){
    #----- ---- ---- ----     cv.sigcom
    XX=list(X4=X4$data)
    XXnew=list(X4=X4new$data)
    t.x.list=list(X4=tj)
    resSC.Lin1=cv.sigcom(XX, Y1D, t.x.list, tj,s.n.basis=snb,t.n.basis=tnb,
                         upper.comp=max(npc1,npc2,npc3,npc4))
    resSC.Lin1.fit=fdata(pred.sigcom(resSC.Lin1, XX),tj)
    prSCLin1=fdata(pred.sigcom(resSC.Lin1,XXnew),tj)
    rr1[i,6]=1-sum(norm.fdata(ldatnew$yt-prSCLin1)^2)/rest
    
    XX=list(X4=X4$data,X2=X2$data)
    XXnew=list(X4=X4new$data,X2=X2new$data)
    t.x.list=list(X4=tj,X2=tj)
    resSC.Lin1=cv.sigcom(XX, Y1D, t.x.list, tj,
                         s.n.basis=snb,t.n.basis=tnb,upper.comp=max(npc1,npc2,npc3,npc4))
    resSC.Lin1.fit=fdata(pred.sigcom(resSC.Lin1, XX),tj)
    prSCLin1=fdata(pred.sigcom(resSC.Lin1,XXnew),tj)
    rr2[i,6]=1-sum(norm.fdata(ldatnew$yt-prSCLin1)^2)/rest
    
    XX=list(X4=X4$data,X2=X2$data,X3=X3$data)
    XXnew=list(X4=X4new$data,X2=X2new$data,X3=X3new$data)
    t.x.list=list(X4=tj,X2=tj,X3=tj)
    resSC.Lin1=cv.sigcom(XX, Y1D, t.x.list, tj,s.n.basis=15,
                         t.n.basis=15,upper.comp=max(npc1,npc2,npc3,npc4))
    resSC.Lin1.fit=fdata(pred.sigcom(resSC.Lin1, XX),tj)
    prSCLin1=fdata(pred.sigcom(resSC.Lin1,XXnew),tj)
    rr3[i,6]=1-sum(norm.fdata(ldatnew$yt-prSCLin1)^2)/rest
  }
  if (compModel["DISC"]){
    #----- ---- ---- ----     cv.nonlinear
    XX=list(X4=X4$data)
    XXnew=list(X4=X4new$data)
    t.x.list=list(X4=tj)
    resSC.NL1=cv.nonlinear(XX, Y1D, t.x.list, tj,
                           s.n.basis=15,x.n.basis=11,t.n.basis=15,
                           upper.comp=max(npc1,npc2,npc3,npc4))
    resSC.NL1.fit=fdata(pred.nonlinear(resSC.NL1, XX),tj)
    prSCNL1=fdata(pred.nonlinear(resSC.NL1,XXnew),tj)
    rr1[i,7]=1-sum(norm.fdata(ldatnew$yt-prSCNL1)^2)/rest
    
    XX=list(X4=X4$data,X2=X2$data)
    XXnew=list(X4=X4new$data,X2=X2new$data)
    t.x.list=list(X4=tj,X2=tj)
    resSC.NL1=cv.nonlinear(XX, Y1D, t.x.list, tj,
                           s.n.basis=15,x.n.basis=11,t.n.basis=15,
                           upper.comp=max(npc1,npc2,npc3,npc4))
    resSC.NL1.fit=fdata(pred.nonlinear(resSC.NL1, XX),tj)
    prSCNL1=fdata(pred.nonlinear(resSC.NL1,XXnew),tj)
    rr2[i,7]=1-sum(norm.fdata(ldatnew$yt-prSCNL1)^2)/rest
    
    XX=list(X4=X4$data,X2=X2$data,X3=X3$data)
    XXnew=list(X4=X4new$data,X2=X2new$data,X3=X3new$data)
    t.x.list=list(X4=tj,X2=tj,X3=tj)
    resSC.NL1=cv.nonlinear(XX, Y1D, t.x.list, tj,
                           s.n.basis=15,x.n.basis=11,t.n.basis=15,
                           upper.comp=max(npc1,npc2,npc3,npc4))
    resSC.NL1.fit=fdata(pred.nonlinear(resSC.NL1, XX),tj)
    prSCNL1=fdata(pred.nonlinear(resSC.NL1,XXnew),tj)
    rr3[i,7]=1-sum(norm.fdata(ldatnew$yt-prSCNL1)^2)/rest
  }
  print(i)
  print(round(rr1[i,],2))
  print(round(rr2[i,],2))
  print(round(rr3[i,],2))
  
}
par(mfrow=c(1,1))
boxplot(rr1)
round(apply(rr1,2,median,na.rm=TRUE),3)
round(apply(rr2,2,median,na.rm=TRUE),3)
round(apply(rr3,2,median,na.rm=TRUE),3)

# # PC-MOF
# FLMFR FSAMFR FKAMFR    PFR   FAMM    LSC   DISC 
# 0.551  0.622  0.607  0.544  0.581  0.544  0.592
                                        #   OJO el FAMM da 0.581!!!
# #PAPPER
# covariates FLMFR FSAMFR FKAMFR PFR FAMM LSC DISC
# FT 0.551 0.622 0.607 0.544 0.373 0.544 0.592
# FT, H 0.502 0.608 0.636 0.512 0.373 0.486 0.636
# FT, H, WS 0.522 0.560 0.644 0.507 0.529 0.500 0.652
# Table 3: Median of R^2_p for the methods (including the covariates in order of
# importance) for predicting log(NCBR+1). The highest values per row are
# printed in bold.
