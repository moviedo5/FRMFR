
library(fda.usc.devel)
library(FRegSigCom)
library(refund)
# Pr y En--->fda.usc.devel no poner añudas
# @noRd o 
#' @keywords internal ponerlo en las funciones auxiliares
#' Price por Pr
#' Energia por En
#' logNBCR --> log(NBCR+1) # incluir 1 
#' Exportar Minverse-->hacer ayuda
#' fregre.gsam.cl/con, classif.bootstrap
#' 
#data("omel2018")
#ldat <- omel2018

load("omel2018-19.rda")
npx=4
npy=4
ttx=ldat$energia$argvals
tty=ldat$precio$argvals
nbasis.x=11;nbasis.y=11
s.n.basis=11;t.n.basis=11;x.n.basis=11
nk=nbasis.y;nkk=c(nbasis.y,nbasis.x,nbasis.y)
spars=list(bs = "ps",m = c(2, 2, 2), k=nkk)

nrep=20
rr <- matrix(NA,ncol=7,nrow=nrep)
colnames(rr) <- c("FLMFR","FSAMFR","FKAMFR","PFFR","FAMM","LSIG","DISC")

compModel <- rep(TRUE,len=7)
names(compModel) <- colnames(rr)
compModel["FAMM"] <- FALSE # All except FAMM Time consuming
#compModel <- !compModel  # Only FAMM
#compModel[1:7]<-T
#compModel[c(1,2,4,6,7)] <- TRUE
compModel

n <- nrow(ldat$df)

set.seed(20030101)

for (i in 1:nrep){
  ii=sample(n,floor(n*.75))
  ldatest=ldat[ii,row=T]
  ldatpred=ldat[-ii,row=T]
  rest=sum(norm.fdata(ldatpred$precio-func.mean(ldatest$precio))^2)
  
  if (any(compModel[c("FLMFR","FSAMFR")])){
    b.x=list("energia"=create.pc.basis(ldatest$energia,1:npx))
    b.y=create.pc.basis(ldatest$precio,1:npy)
  }  
  
  if (compModel["FLMFR"]){
    reslm <- fregre.mlm.fr(precio~energia,ldatest,basis.x=b.x,basis.y=b.y)
    predlm <- predict(reslm,ldatpred)
    rr[i,1] <- 1-sum(norm.fdata(ldatpred$precio-predlm)^2)/rest
  }
  
  if (compModel["FSAMFR"]){
    ressam <- fregre.sam.fr(precio~s(energia),ldatest,basis.x=b.x,basis.y=b.y)
    predsam <- predict(ressam,ldatpred)
    rr[i,2] <- 1-sum(norm.fdata(ldatpred$precio-predsam)^2)/rest
  }
  
  if (compModel["FKAMFR"]){
    par.np=list()
    par.metric=list(energia=list(metric=metric.lp,lp=2))
    reskam <- fregre.kam.fr(precio~energia,ldatest,par.np=par.np,par.metric=par.metric)
    predkam <- predict(reskam,ldatpred)
    rr[i,3] <- 1-sum(norm.fdata(ldatpred$precio-predkam)^2)/rest
  }
  
  
  ydat<-ldatest$precio$data
  xdat <-ldatest$energia$data
  tj<-ldatest$precio$argvals
  si<-ldatest$energia$argvals
  if (compModel["PFFR"]){
    respffpc=pffr(ydat~ffpc(xdat,xind=si,npc.max=npx),yind=tj)
    predffpc=fdata(predict(respffpc,list(xdat=ldatpred$energia$data)),argvals=tty)
    rr[i,4] <- 1-sum(norm.fdata(ldatpred$precio-predffpc)^2)/rest
  }
  
  if (compModel["FAMM"]){
    respffnl=pffr(ydat~ sff(xdat,xind=si, splinepars = spars),
                  yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
                  #bs.int=list(bs="ps",k=21,m=c(2,1))
    )
    predffnl=fdata(predict(respffnl,list(xdat=ldatpred$energia$data)),argvals=tty)
    rr[i,5] <- 1-sum(norm.fdata(ldatpred$precio-predffnl)^2)/rest
  }
  if (compModel["LSIG"]){
    XX=list(Xen=ldatest$energia$data)
    XXnew=list(Xen=ldatpred$energia$data)
    t.x.list=list(Xen=si)
    resSCLin=cv.sigcom(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
    predSCLin=fdata(pred.sigcom(resSCLin,XXnew),tty)
    rr[i,6] <- 1-sum(norm.fdata(ldatpred$precio-predSCLin)^2)/rest
  }
  if (compModel["DISC"]){
    resSCNL=cv.nonlinear(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,
                         x.n.basis=x.n.basis,t.n.basis=t.n.basis)
    predSCNL=fdata(pred.nonlinear(resSCNL,XXnew),tty)  
    rr[i,7] <- 1-sum(norm.fdata(ldatpred$precio-predSCNL)^2)/rest
  }
  print(i)
  print(round(rr[i,],2))
}
rr2018 <- rr
colMeans(rr2018)
#####################################


#data("omel2008")
#ldat <- omel2008
load("omel2008-09.rda")
npx=4
npy=4
ttx=ldat$energia$argvals
tty=ldat$precio$argvals
nbasis.x=11;nbasis.y=11
s.n.basis=11;t.n.basis=11;x.n.basis=11
nk=nbasis.y;nkk=c(nbasis.y,nbasis.x,nbasis.y)
spars=list(bs = "ps",m = c(2, 2, 2), k=nkk)

nrep=20
rr <- matrix(NA,ncol=7,nrow=nrep)
colnames(rr) <- c("FLMFR","FSAMFR","FKAMFR","PFFR","FAMM","LSIG","DISC")

compModel <- rep(TRUE,len=7)
names(compModel) <- colnames(rr)
compModel["FAMM"] <- FALSE # All except FAMM Time consuming
#compModel <- !compModel  # Only FAMM
#compModel[1:7]<-T
#compModel[c(1,2,4,6,7)] <- TRUE
compModel

n <- nrow(ldat$df)

set.seed(20030101)

for (i in 1:nrep){
  ii=sample(n,floor(n*.75))
  ldatest=ldat[ii,row=T]
  ldatpred=ldat[-ii,row=T]
  rest=sum(norm.fdata(ldatpred$precio-func.mean(ldatest$precio))^2)
  
  if (any(compModel[c("FLMFR","FSAMFR")])){
    b.x=list("energia"=create.pc.basis(ldatest$energia,1:npx))
    b.y=create.pc.basis(ldatest$precio,1:npy)
  }  
  
  if (compModel["FLMFR"]){
    reslm <- fregre.mlm.fr(precio~energia,ldatest,basis.x=b.x,basis.y=b.y)
    predlm <- predict(reslm,ldatpred)
    rr[i,1] <- 1-sum(norm.fdata(ldatpred$precio-predlm)^2)/rest
  }
  
  if (compModel["FSAMFR"]){
    ressam <- fregre.sam.fr(precio~s(energia),ldatest,basis.x=b.x,basis.y=b.y)
    predsam <- predict(ressam,ldatpred)
    rr[i,2] <- 1-sum(norm.fdata(ldatpred$precio-predsam)^2)/rest
  }
  
  if (compModel["FKAMFR"]){
    par.np=list()
    par.metric=list(energia=list(metric=metric.lp,lp=2))
    reskam <- fregre.kam.fr(precio~energia,ldatest,par.np=par.np,par.metric=par.metric)
    predkam <- predict(reskam,ldatpred)
    rr[i,3] <- 1-sum(norm.fdata(ldatpred$precio-predkam)^2)/rest
  }
  
  
  ydat<-ldatest$precio$data
  xdat <-ldatest$energia$data
  tj<-ldatest$precio$argvals
  si<-ldatest$energia$argvals
  if (compModel["PFFR"]){
    respffpc=pffr(ydat~ffpc(xdat,xind=si,npc.max=npx),yind=tj)
    predffpc=fdata(predict(respffpc,list(xdat=ldatpred$energia$data)),argvals=tty)
    rr[i,4] <- 1-sum(norm.fdata(ldatpred$precio-predffpc)^2)/rest
  }
  
  if (compModel["FAMM"]){
    respffnl=pffr(ydat~ sff(xdat,xind=si, splinepars = spars),
                  yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
                  #bs.int=list(bs="ps",k=21,m=c(2,1))
    )
    predffnl=fdata(predict(respffnl,list(xdat=ldatpred$energia$data)),argvals=tty)
    rr[i,5] <- 1-sum(norm.fdata(ldatpred$precio-predffnl)^2)/rest
  }
  if (compModel["LSIG"]){
    XX=list(Xen=ldatest$energia$data)
    XXnew=list(Xen=ldatpred$energia$data)
    t.x.list=list(Xen=si)
    resSCLin=cv.sigcom(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
    predSCLin=fdata(pred.sigcom(resSCLin,XXnew),tty)
    rr[i,6] <- 1-sum(norm.fdata(ldatpred$precio-predSCLin)^2)/rest
  }
  if (compModel["DISC"]){
    resSCNL=cv.nonlinear(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,
                         x.n.basis=x.n.basis,t.n.basis=t.n.basis)
    predSCNL=fdata(pred.nonlinear(resSCNL,XXnew),tty)  
    rr[i,7] <- 1-sum(norm.fdata(ldatpred$precio-predSCNL)^2)/rest
  }
  print(i)
  print(round(rr[i,],2))
}
rr2008 <- rr
colMeans(rr2008)
#####################################
# Plots included in the paper
#####################################

pdf(file="omel.pdf",width=10.67,height=6)
load("omel2018-19.rda")
ldat18=ldat
load("omel2008-09.rda")
ldat$energia[[4]]$main="Electricity Demand 2008-09"
ldat18$energia[[4]]$main="Electricity Demand 2018-19"

ldat$precio[[4]]$main="Electricity Price 2008-09"
ldat18$precio[[4]]$main="Electricity Price 2018-19"
par(mfrow=c(2,2))

#m=rbind(c(1,2,5),c(3,4,5))
#layout(m)
rre=range(ldat$energia$data,ldat18$energia$data)
rrp=range(ldat$precio$data,ldat18$precio$data)
plot(ldat$energia,col="gray50",ylim=rre)
lines(func.mean(ldat$energia),lwd=2)
plot(ldat18$energia,col="gray50",ylim=rre)
lines(func.mean(ldat18$energia),lwd=2)

plot(ldat$precio,col="gray50",ylim=rrp)
lines(func.mean(ldat$precio),lwd=2)
plot(ldat18$precio,col="gray50",ylim=rrp)
lines(func.mean(ldat18$precio),lwd=2)

dev.off()


pdf(file="omel2.pdf",width=10.67,height=6)

par(mfrow=c(1,2))

#m=rbind(c(1,2,5),c(3,4,5))
#layout(m)
ran=range(rr2008,rr2018)
boxplot(rr2008,cex.axis=0.75,ylim=ran,main=expression( R[p]^2: textstyle(2008-2009.)))
boxplot(rr2018,cex.axis=0.75,ylim=ran,main=expression( R[p]^2: textstyle(2018-2019.)))

dev.off()

#####################################
# Results for table 5 using the Electricity Demand of the previous day and the previous week as covariates
#load("omel2008-09.rda")
library(microbenchmark)

npx <- 4
npy <- 4
nn <- nrow(ldat$energia)
nlag <- 7
nl <- (nlag+1):nn
names(ldat)

ldatm <- ldata(df=ldat$df[nl,],precio=ldat$precio[nl],ener=ldat$energia[nl],ener1=ldat$energia[nl-1],
               ener7=ldat$energia[nl-nlag],prec1=ldat$precio[nl-1],prec7=ldat$precio[nl-nlag])

ttx <- ldatm$ener$argvals
tty <- ldatm$precio$argvals
nbasis.x <- 11;nbasis.y <- 11
s.n.basis <- 11;t.n.basis <- 11;x.n.basis <- 11
nk=nbasis.y;nkk=c(nbasis.y,nbasis.x,nbasis.y)
spars=list(bs = "ps",m = c(2, 2, 2), k=nkk)


nrep <- 20
rr <- matrix(NA,nrep,7)
rtim <- matrix(NA,nrep,7)
colnames(rr) <- c("FLMFR","FSAMFR","FKAMFR","PFFR","FAMM","LSIG","DISC")
colnames(rtim) <- paste0("t",colnames(rr))

compModel <- rep(TRUE,len=7)
names(compModel) <- colnames(rr)
compModel["FAMM"] <- FALSE # All except FAMM Time consuming
#compModel <- !compModel  # Only FAMM
#compModel[1:7]<-T
#compModel[c(1,2,4,6,7)] <- TRUE
compModel

set.seed(20030101)
for (i in 1:nrep){
  print(paste("Repetition:",i))
  ii <- sample(nrow(ldatm$df),floor(nrow(ldatm$df)*.75))
  ldatest <- ldatm[ii,row=T]
  ldatpred <- ldatm[-ii,row=T]
  rest <- sum(norm.fdata(ldatpred$precio-func.mean(ldatest$precio))^2)
  
  if (any(compModel[c("FLMFR","FSAMFR")])){
    b.x <- list("ener"=create.pc.basis(ldatest$ener,1:npx),"ener1"=create.pc.basis(ldatest$ener1,1:npx),
                "ener7"=create.pc.basis(ldatest$ener7,1:npx),"prec1"=create.pc.basis(ldatest$prec1,1:npy),
                "prec7"=create.pc.basis(ldatest$prec7,1:npy))
    b.y <- create.pc.basis(ldatest$precio,1:npy)
  }
  if (compModel["FLMFR"]){
    aux=microbenchmark({
    reslm <- fregre.mlm.fr(precio~ener1+ener7,ldatest,basis.x=b.x,basis.y=b.y)
    predlm <- predict(reslm,ldatpred)
    rr[i,1] <- 1-sum(norm.fdata(ldatpred$precio-predlm)^2)/rest
    },times=1)
#    print(paste("FLMFR:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,1] <- aux$time/10^9
  }
  if (compModel["FSAMFR"]){
    aux=microbenchmark({
    ressam <- fregre.sam.fr(precio~s(ener1)+s(ener7),ldatest,basis.x=b.x,basis.y=b.y)
    predsam <- predict(ressam,ldatpred)
    rr[i,2] <- 1-sum(norm.fdata(ldatpred$precio-predsam)^2)/rest
    },times=1)
#    print(paste("FSAMFR:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,2] <- aux$time/10^9
  }
  if (compModel["FKAMFR"]){
    aux=microbenchmark({
    par.np=list(ener=list(Ker=AKer.norm,type.S="S.NW"),ener1=list(Ker=AKer.norm,type.S="S.NW"),ener7=list(Ker=AKer.norm,type.S="S.NW"))
    par.metric=list(ener=list(metric=metric.lp,lp=2),ener1=list(metric=metric.lp,lp=2),ener7=list(metric=metric.lp,lp=2))
    reskam <- fregre.kam.fr(precio~ener1+ener7,ldatest,par.np=par.np,par.metric=par.metric)
    predkam <- predict(reskam,ldatpred)
    rr[i,3] <- 1-sum(norm.fdata(ldatpred$precio-predkam)^2)/rest
    },times=1)
#   print(paste("FKAMFR:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,3] <- aux$time/10^9
  }
  ydat<-ldatest$precio$data
  xdat <-ldatest$ener$data
  xdat1<-ldatest$ener1$data
  xdat7<-ldatest$ener7$data
  ydat1<-ldatest$prec1$data
  ydat7<-ldatest$prec7$data
  tj<-ldatest$precio$argvals
  si<-ldatest$ener$argvals
  
  if (compModel["PFFR"]){
    aux=microbenchmark({
    respffpc=pffr(ydat~ffpc(xdat1,xind=si,npc.max=npx)+ffpc(xdat7,xind=si,npc.max=npx),yind=tj)
    predffpc=fdata(predict(respffpc,list(xdat1=ldatpred$ener1$data,xdat7=ldatpred$ener7$data)),argvals=tty)
    rr[i,4] <-1-sum(norm.fdata(ldatpred$precio-predffpc)^2)/rest
    },times=1)
#    print(paste("PFR:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,4] <- aux$time/10^9
  }
  
  if (compModel["FAMM"]){
    aux=microbenchmark({
    respffnl=pffr(ydat~ sff(xdat1,xind=si,splinepars = spars)+
    					sff(xdat7,xind=si,splinepars = spars),
                   yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
    			  )
    predffnl=fdata(predict(respffnl,list(xdat1=ldatpred$ener1$data,xdat7=ldatpred$ener7$data)),argvals=tty)
    rr[i,5] <- 1-sum(norm.fdata(ldatpred$precio-predffnl)^2)/rest
    },times=1)
#    print(paste("FAMM:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,5] <- aux$time/10^9
  }
  
  if (compModel["LSIG"]){
    XX=list(Xen1=ldatest$ener1$data,Xen7=ldatest$ener7$data)
    XXnew=list(Xen1=ldatpred$ener1$data,Xen7=ldatpred$ener7$data)
    t.x.list=list(Xen1=si,Xen7=si)
    aux=microbenchmark({
    resSCLin <- cv.sigcom(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
    predSCLin <- fdata(pred.sigcom(resSCLin,XXnew),tty)
    rr[i,6] <- 1-sum(norm.fdata(ldatpred$precio-predSCLin)^2)/rest
    },times=1)
#    print(paste("LSC:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,6] <- aux$time/10^9
  }
  if (compModel["DISC"]){
    aux=microbenchmark({
    resSCNL=cv.nonlinear(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,
    				 x.n.basis=x.n.basis,t.n.basis=t.n.basis)
    predSCNL=fdata(pred.nonlinear(resSCNL,XXnew),tty)  
    rr[i,7] <- 1-sum(norm.fdata(ldatpred$precio-predSCNL)^2)/rest
    },times=1)
#    print(paste("DISC:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,7] <- aux$time/10^9
  }
  print(round(apply(rr,2,mean,na.rm=TRUE),3))
  print(round(apply(rtim,2,mean,na.rm=TRUE),3))
}
rr_En <- rr
rtim_En <- rtim
round(colMeans(rr_En),3)
round(colMeans(rtim_En),3)


#Table 5: Average of R-square for prediction models with response Prd. Period: 2008-09.
# Rsquare  FLMFR FSAMFR FKAMFR   PFR  FAMM   LSC  DISC
#  Tabla 5 0.336  0.456  0.631 0.324 0.564  0.600  0.655
# MOF PC   0.336  0.456  0.631  0.324  0.564  0.601  0.655 
#####################################
# Results for table 5 using the price of the previous day and the previous week as covariates
#load("omel2008-09.rda")
data("omel2008")
ldat <- omel2008

npx <- 4
npy <- 4
nn <- nrow(ldat$energia)
nlag <- 7
nl <- (nlag+1):nn
names(ldat)

ldatm <- ldata(df=ldat$df[nl,],precio=ldat$precio[nl],ener=ldat$energia[nl],ener1=ldat$energia[nl-1],
            ener7=ldat$energia[nl-nlag],prec1=ldat$precio[nl-1],prec7=ldat$precio[nl-nlag])


ttx <- ldatm$ener$argvals
tty <- ldatm$precio$argvals
nbasis.x <- 11;nbasis.y <- 11
s.n.basis <- 11;t.n.basis <- 11;x.n.basis <- 11
nk=nbasis.y;nkk=c(nbasis.y,nbasis.x,nbasis.y)
spars=list(bs = "ps",m = c(2, 2, 2), k=nkk)

nrep <- 20
rr <- matrix(NA,nrep,7)
rtim <- matrix(NA,nrep,7)
colnames(rr) <- c("FLMFR","FSAMFR","FKAMFR","PFFR","FAMM","LSIG","DISC")
colnames(rtim) <- paste0("t",colnames(rr))

compModel <- rep(TRUE,len=7)
names(compModel) <- colnames(rr)
compModel["FAMM"] <- FALSE # All except FAMM Time consuming
#compModel <- !compModel  # Only FAMM
#compModel[1:7]<-T
#compModel[c(1,2,4,6,7)] <- TRUE
compModel


set.seed(20030101)
for (i in 1:nrep){
  print(paste("Repetition:",i))
  ii <- sample(nrow(ldatm$df),floor(nrow(ldatm$df)*.75))
  ldatest <- ldatm[ii,row=T]
  ldatpred <- ldatm[-ii,row=T]
  rest <- sum(norm.fdata(ldatpred$precio-func.mean(ldatest$precio))^2)
  
  if (any(compModel[c("FLMFR","FSAMFR")])){
    b.x <- list("ener"=create.pc.basis(ldatest$ener,1:npx),"ener1"=create.pc.basis(ldatest$ener1,1:npx),
           "ener7"=create.pc.basis(ldatest$ener7,1:npx),"prec1"=create.pc.basis(ldatest$prec1,1:npy),
           "prec7"=create.pc.basis(ldatest$prec7,1:npy))
    b.y <- create.pc.basis(ldatest$precio,1:npy)
  }
  if (compModel["FLMFR"]){
    aux=microbenchmark({
    reslm <- fregre.mlm.fr(precio~prec1+prec7,ldatest,basis.x=b.x,basis.y=b.y)
    predlm <- predict(reslm,ldatpred)
    rr[i,1] <- 1-sum(norm.fdata(ldatpred$precio-predlm)^2)/rest
    },times=1)
#    print(paste("FLMFR:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,1] <- aux$time/10^9 
  }
  if (compModel["FSAMFR"]){
    aux=microbenchmark({
    ressam <- fregre.sam.fr(precio~s(prec1)+s(prec7),ldatest,basis.x=b.x,basis.y=b.y)
    predsam <- predict(ressam,ldatpred)
    rr[i,2] <- 1-sum(norm.fdata(ldatpred$precio-predsam)^2)/rest
    },times=1)
#    print(paste("FSAMFR:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,2] <- aux$time/10^9 
  }
  if (compModel["FKAMFR"]){
    aux=microbenchmark({
    par.np <- list(precio=list(Ker=AKer.norm,type.S="S.NW"),prec1=list(Ker=AKer.norm,type.S="S.NW"),prec7=list(Ker=AKer.norm,type.S="S.NW"))
    par.metric <- list(precio=list(metric=metric.lp,lp=2),prec1=list(metric=metric.lp,lp=2),prec7=list(metric=metric.lp,lp=2))
    reskam <- fregre.kam.fr(precio~prec1+prec7,ldatest,par.np=par.np,par.metric=par.metric)
    predkam <- predict(reskam,ldatpred)
    rr[i,3] <- 1-sum(norm.fdata(ldatpred$precio-predkam)^2)/rest
    },times=1)
#    print(paste("FKAMFR:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,3] <-aux$time/10^9 
}
  ydat<-ldatest$precio$data
  xdat <-ldatest$ener$data
  xdat1<-ldatest$ener1$data
  xdat7<-ldatest$ener7$data
  ydat1<-ldatest$prec1$data
  ydat7<-ldatest$prec7$data
  tj<-ldatest$precio$argvals
  si<-ldatest$ener$argvals
  
  if (compModel["PFFR"]){
    aux=microbenchmark({
    respffpc <- pffr(ydat~ffpc(ydat1,xind=tj,npc.max=npy)+ffpc(ydat7,xind=tj,npc.max=npy),yind=tj)
    predffpc <- fdata(predict(respffpc,list(ydat1=ldatpred$prec1$data,ydat7=ldatpred$prec7$data)),argvals=tty)
    rr[i,4] <-1-sum(norm.fdata(ldatpred$precio-predffpc)^2)/rest
    },times=1)
#    print(paste("PFR:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,4] <- aux$time/10^9 
  }
  
  if (compModel["FAMM"]){
    aux=microbenchmark({
    respffnl <- pffr(ydat~ sff(ydat1,xind=tj,splinepars = spars)+
                    sff(ydat7,xind=tj,splinepars = spars),
                  yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
    )
    predffnl <- fdata(predict(respffnl,list(ydat1=ldatpred$prec1$data,ydat7=ldatpred$prec7$data)),argvals=tty)
    rr[i,5] <- 1-sum(norm.fdata(ldatpred$precio-predffnl)^2)/rest
    },times=1)
#    print(paste("FAMM:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,5] <- aux$time/10^9 
  }
 
  if (compModel["LSIG"]){

    XX <- list(Yen1=ldatest$prec1$data,Yen7=ldatest$prec7$data)
    XXnew <- list(Yen1=ldatpred$prec1$data,Yen7=ldatpred$prec7$data)
    t.x.list <- list(Yen1=tj,Yen7=tj)
    aux=microbenchmark({
    resSCLin <- cv.sigcom(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
    predSCLin <- fdata(pred.sigcom(resSCLin,XXnew),tty)
    rr[i,6] <- 1-sum(norm.fdata(ldatpred$precio-predSCLin)^2)/rest
    },times=1)
#    print(paste("LSC:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,6] <- aux$time/10^9 
    }
  if (compModel["DISC"]){
    aux=microbenchmark({
    resSCNL <- cv.nonlinear(XX,ydat,t.x.list,tj,s.n.basis=s.n.basis,
                         x.n.basis=x.n.basis,t.n.basis=t.n.basis)
    predSCNL <- fdata(pred.nonlinear(resSCNL,XXnew),tty)  
    rr[i,7] <- 1-sum(norm.fdata(ldatpred$precio-predSCNL)^2)/rest
    },times=1)

#    print(paste("DISC:",round(difftime(itime2,itime,units="mins"),2)))
    rtim[i,7] <- aux$time/10^9 
    }
  print(round(apply(rr,2,mean,na.rm=TRUE),3))
  print(round(apply(rtim,2,mean,na.rm=TRUE),3))
}

rr_Pr <- rr
rtim_Pr <- rtim
round(colMeans(rr_Pr),3)
round(colMeans(rtim_Pr),3)

#Table 5: Average of R-square for prediction models with response Prd. Period: 2008-09.
# Rsquare          FLMFR FSAMFR FKAMFR   PFR  FAMM   LSC  DISC
# Tabla 5 results  0.900  0.881  0.886 0.886 0.859 0.898 0.881
# MOF PC results   0.890  0.881  0.885 0.886 0.859 0.898 0.881       

               