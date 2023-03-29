#' Nonlinear Estimation of Additive Functional
#' Regression Models with Functional Response
#' Authors: Morteza Amini, Mohammad Darbalaei,
#' Manuel Febrero–Bande and Manuel Oviedo–de la Fuente
#' 
#' 3.2.1 Bike-sharing data
#' Source code
rm(list=ls())
#'


library(fda.usc.devel)
library(refund)
library(FRegSigCom)

data("BikeSharing")
names(BikeSharing)

draw <- FALSE
x1 <- BikeSharing$temp
x2 <- BikeSharing$humidity
x3 <- BikeSharing$windspeed
x4 <- BikeSharing$feeltemp
y <- BikeSharing$logNBCR
if (draw){
pdf(file="bikesharing.pdf",width=10.67,height=6)
  m=cbind(c(1,3,5),c(2,4,5))
  layout(m)
  lesp <- c(35, 53, 73)
  #par(mfrow=c(2,2))
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
  dev.off()
}
#' Figure 5: The Plot of Bike-sharing Variables over 112 Saturdays from Jan
#' 1st, 2011 to Dec 1st, 2012 . The Black, Red, and Green Lines Correspond
#' to Dates Sep 17th, 2011, Jan 21th, 2012, and Jun 9th, 2012

#' Models in fda.usc package: 
#' fregre.lm.fr: Functional Linear Model with Functional Response (FLMFR)
#' fregre.sam.fr: Functional Spectral Additive Model with Functional Response (FSAMFR)
#' fregre.kam.fr: Functional Kernel Additive Model with Functional Response (FKAMFR)
#' 
#' Models in refund package: 
#' pffr with ffpc: Penalized flexible functional regression with linear effect
#' (PFR) with PC-based function-on-function regression term.
#' pffr with sff : Penalized flexible functional regression with nonlinear
#' effect (FAMM) with smooth function-on-function regression term.
#' 
#' Models in FRegSigCom (obsolete) package: 
#' cv.sigcom: Cross-validation for linear function-on-function regression (LSC)
#' cv.nonlinear: Cross-validation for nonlinear function-on-function regression (DISC)

#' Repeat procedure for the following covariates (including the covariates
#' in order of importance):
#' i) feeltemp (FT),
#' ii) feeltemp (FT), humidity (H) 
#' iii) feeltemp (FT), humidity (H), windspeed (WS)
#' 
tj <- seq(0.5,23.5)

nrep=20
rr=matrix(NA,nrep,7)
colnames(rr)=c("FLMFR","FSAMFR","FKAMFR","PFR","FAMM","LSC","DISC")
set.seed(20030101)
for (i in 1:nrep){
  print(i)
index=sample(1:102,82)

npcy=6;npc1=2;npc2=4;npc3=6;npc4=2
yt=y[index];X1=x1[index]; X2=x2[index]; X3=x3[index]; X4=x4[index]
 
ytnew=y[-index]; X1new=x1[-index]; X2new=x2[-index];
X3new=x3[-index]; X4new=x4[-index]
ldatos=ldata(X1=X1,X2=X2,X3=X3,X4=X4,yt=yt)
ldatnew=ldata(X1=X1new,X2=X2new,X3=X3new,X4=X4new,yt=ytnew)

# ldatos=BikeSharing[index,row=T]
# ldatnew=BikeSharing[-index,row=T]
# names(ldatnew)[-1]<-names(ldatos)[-1]<-c("yt","X1","X2","X3","X4")
# X1 <- ldatos$X1;X2 <- ldatos$X2;X3 <- ldatos$X3;X4 <- ldatos$X4
# yt <- ldatos$yt

rest=sum(norm.fdata(ldatnew$yt-func.mean(ldatos$yt))^2)

##### fda.usc.devel
b.x=list(X1=create.pc.basis(X1,1:npc1),X2=create.pc.basis(X2,1:npc2),
X3=create.pc.basis(X3,1:npc3),X4=create.pc.basis(X4,1:npc4))
b.y=create.pc.basis(yt,1:npcy)
#
pmetric=list(df=data.frame(idx=1:nrow(yt)),yt=list(metric=metric.lp,lp=2), 
X1=list(metric=metric.lp,lp=2),X2=list(metric=metric.lp,lp=2),
X3=list(metric=metric.lp,lp=2),X4=list(metric=metric.lp,lp=2))
p.np=list(X1=list(Ker=AKer.norm),X2=list(Ker=AKer.norm),
		  X3=list(Ker=AKer.norm),X4=list(Ker=AKer.norm))

#reslin=fregre.mlm.fr(yt~X4+X2+X3,data=ldatos,basis.y=b.y,basis.x=b.x)
#ressam=fregre.sam.fr(yt~s(X4,k=5)+s(X2,k=5)+s(X3,k=5),data=ldatos,basis.y=b.y,basis.x=b.x)
#reskam=fregre.kam.fr(yt~X4+X2+X3,data=ldatos,par.metric=pmetric)
reslin=fregre.mlm.fr(yt~X4+X2,data=ldatos,basis.y=b.y,basis.x=b.x)
ressam=fregre.sam.fr(yt~s(X4)+s(X2,k=5),data=ldatos,basis.y=b.y,basis.x=b.x)
ressam2=fregre.sam.fr(yt~s(X4,k=3)+s(X2,k=3),data=ldatos,basis.x=b.x)

reskam=fregre.kam.fr(yt~X4+X2,data=ldatos,par.metric=pmetric)
# reslin=fregre.mlm.fr(yt~X4,data=ldatos,basis.y=b.y,basis.x=b.x)
# ressam=fregre.sam.fr(yt~s(X4),data=ldatos,basis.y=b.y,basis.x=b.x)
# reskam=fregre.kam.fr(yt~X4,data=ldatos,par.metric=pmetric)


prlin=predict(reslin,ldatnew)
prsam=predict.fregre.sam.fr(ressam,ldatnew)
prsam2=predict.fregre.sam.fr(ressam2,ldatnew)
prkam=predict(reskam,ldatnew)

rr[i,1]=1-sum(norm.fdata(ldatnew$yt-prlin)^2)/rest
rr[i,2]=1-sum(norm.fdata(ldatnew$yt-prsam)^2)/rest
rr[i,3]=1-sum(norm.fdata(ldatnew$yt-prkam)^2)/rest

##### Refund
Y1D=yt$data;X1D=X1$data;X2D=X2$data;X3D=X3$data;X4D=X4$data

respff1=pffr(Y1D~sff(X4D,xind=tj,splinepars=list(bs="ps",m=c(2,2,2),k=7))
			    +sff(X2D,xind=tj,splinepars=list(bs="ps",m=c(2,2,2),k=7))
          #       +sff(X3D,xind=tj,splinepars=list(bs="ps",m=c(2,2,2),k=7))
				,yind=tj,bs.yindex=list(bs="ps",k=7,m=c(2,1)))

#respff1.fit=fdata(matrix(respff1$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)

respff1pc=pffr(Y1D~ffpc(X4D,xind=tj,npc.max=npc4,splinepars=list(k=11))
                  +ffpc(X2D,xind=tj,npc.max=npc2,splinepars=list(k=11))
#                  +ffpc(X3D,xind=tj,npc.max=npc3,splinepars=list(k=11))
					,yind=tj)
respff1pc.fit=fdata(matrix(respff1pc$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)

#prff1=fdata(predict(respff1,list(X4D=X4new$data,X2D=X2new$data,X3D=X3new$data)),argvals=tj)
prff1=fdata(predict(respff1,list(X4D=X4new$data,X2D=X2new$data)),argvals=tj)
# prff1=fdata(predict(respff1,list(X4D=X4new$data)),argvals=tj)

#prff1.pc=fdata(predict(respff1pc,list(X4D=X4new$data,X2D=X2new$data,X3D=X3new$data)),argvals=tj)
prff1.pc=fdata(predict(respff1pc,list(X4D=X4new$data,X2D=X2new$data)),argvals=tj)
# prff1.pc=fdata(predict(respff1pc,list(X4D=X4new$data)),argvals=tj)



rr[i,4]=1-sum(norm.fdata(ldatnew$yt-prff1.pc)^2)/rest
rr[i,5]=1-sum(norm.fdata(ldatnew$yt-prff1)^2)/rest

##### FRegSigCom
#XX=list(X4=X4$data,X2=X2$data,X3=X3$data)
#XXnew=list(X4=X4new$data,X2=X2new$data,X3=X3new$data)
#t.x.list=list(X4=tj,X2=tj,X3=tj)
XX=list(X4=X4$data,X2=X2$data)
XXnew=list(X4=X4new$data,X2=X2new$data)
 t.x.list=list(X4=tj,X2=tj)
# XX=list(X4=X4$data)
# XXnew=list(X4=X4new$data)
# t.x.list=list(X4=tj)

##### cv.nonlinear
resSC.NL1=cv.nonlinear(XX, Y1D, t.x.list, tj,s.n.basis=15,x.n.basis=11,t.n.basis=15,
			upper.comp=max(npc1,npc2,npc3,npc4))

resSC.NL1.fit=fdata(pred.nonlinear(resSC.NL1, XX),tj)
prSCNL1=fdata(pred.nonlinear(resSC.NL1,XXnew),tj)
#----- ---- ---- ----     cv.sigcom
resSC.Lin1=cv.sigcom(XX, Y1D, t.x.list, tj,s.n.basis=15,t.n.basis=15,upper.comp=max(npc1,npc2,npc3,npc4))
resSC.Lin1.fit=fdata(pred.sigcom(resSC.Lin1, XX),tj)
prSCLin1=fdata(pred.sigcom(resSC.Lin1,XXnew),tj)
#
rr[i,6]=1-sum(norm.fdata(ldatnew$yt-prSCLin1)^2)/rest
rr[i,7]=1-sum(norm.fdata(ldatnew$yt-prSCNL1)^2)/rest
print(rr[i,])
}

#' Boxplot of $R^2_p$ for the models (including the covariates in order of
#' importance) for predicting log(NCBR + 1).
boxplot(rr)
#' Median of $R^2_p$ for the models (including the covariates in order of
#' importance) for predicting log(NCBR + 1).
round(apply(rr,2,median),2)
# 2 variables 
round(apply(rr,2,median),3)
# FLMFR FSAMFR FKAMFR    PFR   FAMM    LSC   DISC 
# 0.502  0.608  0.636  0.512  0.406  0.486  0.636 
# 0.502  0.608   0.636 0.512  0.510  0.486  0.636
# # 1
# > apply(rr,2,median)
# FLMFR    FSAMFR    FKAMFR       PFR      FAMM 
# 0.5507157 0.6215387 0.6068991 0.5436763 0.5811326 
# LSC      DISC 
# 0.5441019 0.5919039 