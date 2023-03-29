
library(fda.usc.devel)
library(FRegSigCom)
library(refund)

data("omel2018")
npx=4
npy=4
ls()
ldat <- omel2018

ttx=ldat$energia$argvals
tty=ldat$precio$argvals
nbasis.x=11;nbasis.y=11
s.n.basis=11;t.n.basis=11;x.n.basis=11
nk=nbasis.y;nkk=c(nbasis.x,nbasis.y)

set.seed(20030101)
nrep=10
rr=matrix(NA,nrep,7)
colnames(rr)=c("FLMFR","FSAMFR","FKAMFR","PFR","FAMM","LSC","DISC")
for (i in 1:nrep){
  ii=sample(nrow(ldat$df),floor(nrow(ldat$df)*.75))
  ldatest=ldat[ii,row=T]
  ldatpred=ldat[-ii,row=T]
  b.x=list("energia"=create.pc.basis(ldatest$energia,1:npx))
  b.y=create.pc.basis(ldatest$precio,1:npy)
  rest=sum(norm.fdata(ldatpred$precio-func.mean(ldatest$precio))^2)
  reslm=fregre.mlm.fr(precio~energia,ldatest,basis.x=b.x,basis.y=b.y)
  predlm=predict(reslm,ldatpred)
  rr[i,1]=1-sum(norm.fdata(ldatpred$precio-predlm)^2)/rest
  
  ressam=fregre.sam.fr(precio~s(energia),ldatest,basis.x=b.x,basis.y=b.y)
  predsam=predict(ressam,ldatpred)
  rr[i,2]=1-sum(norm.fdata(ldatpred$precio-predsam)^2)/rest
  
  par.np=list()
  par.metric=list(energia=list(metric=metric.lp,lp=2))
  reskam=fregre.kam.fr(precio~energia,ldatest,par.np=par.np,par.metric=par.metric)
  predkam=predict(reskam,ldatpred)
  rr[i,3]=1-sum(norm.fdata(ldatpred$precio-predkam)^2)/rest
  
  ydat<-ldatest$precio$data
  xdat <-ldatest$energia$data
  tj<-ldatest$precio$argvals
  si<-ldatest$energia$argvals
  
  respffpc=pffr(ydat~ffpc(xdat,xind=si,npc.max=npx),yind=tj)
  predffpc=fdata(predict(respffpc,list(xdat=ldatpred$energia$data)),argvals=tty)
  rr[i,4]=1-sum(norm.fdata(ldatpred$precio-predffpc)^2)/rest
  
  respffnl=pffr(ydat~ sff(xdat,xind=si,
                          splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk)),
                yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
                #bs.int=list(bs="ps",k=21,m=c(2,1))
  )
  predffnl=fdata(predict(respffnl,list(xdat=ldatpred$energia$data)),argvals=tty)
  rr[i,5]=1-sum(norm.fdata(ldatpred$precio-predffnl)^2)/rest
  XX=list(Xen=ldatest$energia$data)
  XXnew=list(Xen=ldatpred$energia$data)
  t.x.list=list(Xen=si)
  resSCLin=cv.sigcom(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
  predSCLin=fdata(pred.sigcom(resSCLin,XXnew),tty)
  rr[i,6]=1-sum(norm.fdata(ldatpred$precio-predSCLin)^2)/rest
  resSCNL=cv.nonlinear(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,
                       x.n.basis=x.n.basis,t.n.basis=t.n.basis)
  predSCNL=fdata(pred.nonlinear(resSCNL,XXnew),tty)  
  rr[i,7]=1-sum(norm.fdata(ldatpred$precio-predSCNL)^2)/rest
}

# Simulación varias variables
load("omel2008-09.rda")
npx=4
npy=4

nn=nrow(ldat$energia)
nlag=7
nl=(nlag+1):nn

ldatm=ldata(df=ldat$df[nl,],precio=ldat$precio[nl],ener=ldat$energia[nl],ener1=ldat$energia[nl-1],
            ener7=ldat$energia[nl-nlag],prec1=ldat$precio[nl-1],prec7=ldat$precio[nl-nlag])


ttx=ldatm$ener$argvals
tty=ldatm$precio$argvals
nbasis.x=11;nbasis.y=11
s.n.basis=11;t.n.basis=11;x.n.basis=11
nk=nbasis.y;nkk=c(nbasis.x,nbasis.y)

set.seed(20030101)
nrep=10
rr=matrix(NA,nrep,7)
rtim=matrix(NA,nrep,7)
colnames(rr)=c("FLMFR","FSAMFR","FKAMFR","PFR","FAMM","LSC","DISC")
for (i in 1:nrep){
  print(paste("Repetition:",i))
  ii=sample(nrow(ldatm$df),floor(nrow(ldatm$df)*.75))
  ldatest=ldatm[ii,row=T]
  ldatpred=ldatm[-ii,row=T]
  itime=Sys.time()
  b.x=list("ener"=create.pc.basis(ldatest$ener,1:npx),"ener1"=create.pc.basis(ldatest$ener1,1:npx),
           "ener7"=create.pc.basis(ldatest$ener7,1:npx),"prec1"=create.pc.basis(ldatest$prec1,1:npy),
           "prec7"=create.pc.basis(ldatest$prec7,1:npy))
  b.y=create.pc.basis(ldatest$precio,1:npy)
  rest=sum(norm.fdata(ldatpred$precio-func.mean(ldatest$precio))^2)
  #reslm=fregre.mlm.fr(precio~ener1+ener7,ldatest,basis.x=b.x,basis.y=b.y)
  reslm=fregre.mlm.fr(precio~prec1+prec7,ldatest,basis.x=b.x,basis.y=b.y)
  predlm=predict(reslm,ldatpred)
  rr[i,1]=1-sum(norm.fdata(ldatpred$precio-predlm)^2)/rest
  itime2=Sys.time()
  print(paste("FLMFR:",round(difftime(itime2,itime,units="mins"),2)))
  rtim[i,1]=difftime(itime2,itime,units="mins")
  #ressam=fregre.sam.fr(precio~s(ener1)+s(ener7),ldatest,basis.x=b.x,basis.y=b.y)
  ressam=fregre.sam.fr(precio~s(prec1)+s(prec7),ldatest,basis.x=b.x,basis.y=b.y)
  predsam=predict(ressam,ldatpred)
  rr[i,2]=1-sum(norm.fdata(ldatpred$precio-predsam)^2)/rest
  itime3=Sys.time()
  print(paste("FSAMFR:",round(difftime(itime3,itime2,units="mins"),2)))
  rtim[i,2]=difftime(itime3,itime2,units="mins")
  #par.np=list(ener=list(Ker=AKer.norm,type.S="S.NW"),ener1=list(Ker=AKer.norm,type.S="S.NW"),ener7=list(Ker=AKer.norm,type.S="S.NW"))
  #par.metric=list(ener=list(metric=metric.lp,lp=2),ener1=list(metric=metric.lp,lp=2),ener7=list(metric=metric.lp,lp=2))
  #reskam=fregre.kam.fr(precio~ener1+ener7,ldatest,par.np=par.np,par.metric=par.metric)
  par.np=list(precio=list(Ker=AKer.norm,type.S="S.NW"),prec1=list(Ker=AKer.norm,type.S="S.NW"),prec7=list(Ker=AKer.norm,type.S="S.NW"))
  par.metric=list(precio=list(metric=metric.lp,lp=2),prec1=list(metric=metric.lp,lp=2),prec7=list(metric=metric.lp,lp=2))
  reskam=fregre.kam.fr(precio~prec1+prec7,ldatest,par.np=par.np,par.metric=par.metric)
  predkam=predict(reskam,ldatpred)
  rr[i,3]=1-sum(norm.fdata(ldatpred$precio-predkam)^2)/rest
  itime4=Sys.time()
  print(paste("FKAMFR:",round(difftime(itime4,itime3,units="mins"),2)))
  rtim[i,3]=difftime(itime4,itime3,units="mins")
  ydat<-ldatest$precio$data
  xdat <-ldatest$ener$data
  xdat1<-ldatest$ener1$data
  xdat7<-ldatest$ener7$data
  ydat1<-ldatest$prec1$data
  ydat7<-ldatest$prec7$data
  tj<-ldatest$precio$argvals
  si<-ldatest$ener$argvals
  
  #respffpc=pffr(ydat~ffpc(xdat1,xind=si,npc.max=npx)+ffpc(xdat7,xind=si,npc.max=npx),yind=tj)
  #predffpc=fdata(predict(respffpc,list(xdat1=ldatpred$ener1$data,xdat7=ldatpred$ener7$data)),argvals=tty)
  respffpc=pffr(ydat~ffpc(ydat1,xind=tj,npc.max=npy)+ffpc(ydat7,xind=tj,npc.max=npy),yind=tj)
  predffpc=fdata(predict(respffpc,list(ydat1=ldatpred$prec1$data,ydat7=ldatpred$prec7$data)),argvals=tty)
  rr[i,4]=1-sum(norm.fdata(ldatpred$precio-predffpc)^2)/rest
  itime5=Sys.time()
  print(paste("PFR:",round(difftime(itime5,itime4,units="mins"),2)))
  rtim[i,4]=difftime(itime5,itime4,units="mins")
  
  #respffnl=pffr(ydat~ sff(xdat1,xind=si,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))+
  #					sff(xdat7,xind=si,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk)),
  #               yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
  #bs.int=list(bs="ps",k=21,m=c(2,1))
  #			  )
  #predffnl=fdata(predict(respffnl,list(xdat1=ldatpred$ener1$data,xdat7=ldatpred$ener7$data)),argvals=tty)
  respffnl=pffr(ydat~ sff(ydat1,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))+
                  sff(ydat7,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk)),
                yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
                #bs.int=list(bs="ps",k=21,m=c(2,1))
  )
  predffnl=fdata(predict(respffnl,list(ydat1=ldatpred$prec1$data,ydat7=ldatpred$prec7$data)),argvals=tty)
  rr[i,5]=1-sum(norm.fdata(ldatpred$precio-predffnl)^2)/rest
  itime6=Sys.time()
  print(paste("FAMM:",round(difftime(itime6,itime5,units="mins"),2)))
  rtim[i,5]=difftime(itime6,itime5,units="mins")
  #XX=list(Xen1=ldatest$ener1$data,Xen7=ldatest$ener7$data)
  #XXnew=list(Xen1=ldatpred$ener1$data,Xen7=ldatpred$ener7$data)
  #t.x.list=list(Xen1=si,Xen7=si)
  XX=list(Yen1=ldatest$prec1$data,Yen7=ldatest$prec7$data)
  XXnew=list(Yen1=ldatpred$prec1$data,Yen7=ldatpred$prec7$data)
  t.x.list=list(Yen1=tj,Yen7=tj)
  resSCLin=cv.sigcom(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
  predSCLin=fdata(pred.sigcom(resSCLin,XXnew),tty)
  rr[i,6]=1-sum(norm.fdata(ldatpred$precio-predSCLin)^2)/rest
  itime7=Sys.time()
  print(paste("LSC:",round(difftime(itime7,itime6,units="mins"),2)))
  rtim[i,6]=difftime(itime7,itime6,units="mins")
  
  #resSCNL=cv.nonlinear(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,
  #				 x.n.basis=x.n.basis,t.n.basis=t.n.basis)
  #predSCNL=fdata(pred.nonlinear(resSCNL,XXnew),tty)  
  resSCNL=cv.nonlinear(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,
                       x.n.basis=x.n.basis,t.n.basis=t.n.basis)
  predSCNL=fdata(pred.nonlinear(resSCNL,XXnew),tty)  
  rr[i,7]=1-sum(norm.fdata(ldatpred$precio-predSCNL)^2)/rest
  itime8=Sys.time()
  print(paste("DISC:",round(difftime(itime8,itime7,units="mins"),2)))
  rtim[i,7]=difftime(itime8,itime7,units="mins")
  print(apply(rr,2,mean,na.rm=TRUE))
  print(apply(rtim,2,mean,na.rm=TRUE))
}




