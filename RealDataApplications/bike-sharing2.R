library(fda.usc.devel)
library(refund)
library(FRegSigCom)
data=read.csv("hour.csv")
l=which(data$weekday==6)

datan=data[l,]
tab=table(datan$dteday,datan$hr)
dayf=unique(rownames(which(tab==0,arr.ind=TRUE)))
l=which(data$weekday==6 & !(data$dteday %in% dayf))
datan=data[l,]
ddays=unique(datan$dteday)
lesp=which(ddays %in% c("2011-09-17","2012-01-21","2012-06-09"))
 tmin=-8; tmax=39
 atmin=-16; atmax=50
 fhum=100
 fws=67
tj=seq(0.5,23.5)
y0=fdata(matrix(datan$casual,ncol=24,byrow=TRUE),tj,c(0,24),
		names=list(main="NBCR",xlab="Hour",ylab="NBCR"))
y=fdata(matrix(log(datan$casual+1),ncol=24,byrow=TRUE),tj,c(0,24),
		names=list(main="log(NBCR)",xlab="Hour",ylab="log(NBCR)"))
x1=fdata(matrix(tmin+datan$temp*(tmax-tmin),ncol=24,byrow=TRUE),tj,c(0,24),
		names=list(main="Temperature",xlab="Hour",ylab="ºC"))
x2=fdata(matrix(fhum*datan$hum,ncol=24,byrow=TRUE),tj,c(0,24),
		names=list(main="Humidity",xlab="Hour",ylab=""))
x3=fdata(matrix(fws*datan$windspeed,ncol=24,byrow=TRUE),tj,c(0,24),
		names=list(main="Wind Speed",xlab="Hour",ylab=""))
x4=fdata(matrix(atmin+datan$atemp*(atmax-atmin),ncol=24,byrow=TRUE),tj,c(0,24),
		names=list(main="Feeling Temperature",xlab="Hour",ylab="ºC"))
ldatm=ldata(y=y,x1=x1,x2=x2,x3=x3,x4=x4)

pdf(file="bikesharing.pdf",width=10.67,height=6)
m=cbind(c(1,3,5),c(2,4,5))
layout(m)
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


nrep=20
rr1=matrix(NA,nrep,7)
rr2=matrix(NA,nrep,7)
rr3=matrix(NA,nrep,7)
colnames(rr1)=c("FLMFR","FSAMFR","FKAMFR","PFR","FAMM","LSC","DISC")
colnames(rr2)=c("FLMFR","FSAMFR","FKAMFR","PFR","FAMM","LSC","DISC")
colnames(rr3)=c("FLMFR","FSAMFR","FKAMFR","PFR","FAMM","LSC","DISC")


# one covariate X4
set.seed(20030101)
for (i in 1:nrep){
index=sample(1:102,82)
print(paste0("Repetición:",i,"/",nrep))
yt=y[index];npcy=6
X1=x1[index];npc1=2
X2=x2[index];npc2=4
X3=x3[index];npc3=6
X4=x4[index];npc4=2

ytnew=y[-index]
X1new=x1[-index]
X2new=x2[-index]
X3new=x3[-index]
X4new=x4[-index]

ldatos=ldata(X1=X1,X2=X2,X3=X3,X4=X4,yt=yt)
ldatnew=ldata(X1=X1new,X2=X2new,X3=X3new,X4=X4new,yt=ytnew)

rest=sum(norm.fdata(ldatnew$yt-func.mean(ldatos$yt))^2)

####### fda.usc.devel
b.x=list(X1=create.pc.basis(X1,1:npc1),X2=create.pc.basis(X2,1:npc2),
X3=create.pc.basis(X3,1:npc3),X4=create.pc.basis(X4,1:npc4))
b.y=create.pc.basis(yt,1:npcy)
#
pmetric=list(df=data.frame(idx=1:nrow(yt)),yt=list(metric=metric.lp,lp=2), 
X1=list(metric=metric.lp,lp=2),X2=list(metric=metric.lp,lp=2),
X3=list(metric=metric.lp,lp=2),X4=list(metric=metric.lp,lp=2))
p.np=list(X1=list(Ker=AKer.norm),X2=list(Ker=AKer.norm),
		  X3=list(Ker=AKer.norm),X4=list(Ker=AKer.norm))

reslin=fregre.mlm.fr(yt~X4,data=ldatos,basis.y=b.y,basis.x=b.x)
ressam=fregre.sam.fr(yt~s(X4),data=ldatos,basis.y=b.y,basis.x=b.x)
reskam=fregre.kam.fr(yt~X4,data=ldatos,par.metric=pmetric)


prlin=predict(reslin,ldatnew)
prsam=predict(ressam,ldatnew)
prkam=predict(reskam,ldatnew)

rr1[i,1]=1-sum(norm.fdata(ldatnew$yt-prlin)^2)/rest
rr1[i,2]=1-sum(norm.fdata(ldatnew$yt-prsam)^2)/rest
rr1[i,3]=1-sum(norm.fdata(ldatnew$yt-prkam)^2)/rest

##### Refund
Y1D=yt$data;X1D=X1$data;X2D=X2$data;X3D=X3$data;X4D=X4$data

respff1=pffr(Y1D~sff(X4D,xind=tj,splinepars=list(bs="ps",m=c(2,2,2),k=c(9,11,11)))
				,yind=tj,bs.yindex=list(bs="ps",k=7,m=c(2,1)))

respff1.fit=fdata(matrix(respff1$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)

respff1pc=pffr(Y1D~ffpc(X4D,xind=tj,npc.max=npc4,splinepars=list(k=11))
					,yind=tj)
respff1pc.fit=fdata(matrix(respff1pc$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)

prff1=fdata(predict(respff1,list(X4D=X4new$data)),argvals=tj)

prff1.pc=fdata(predict(respff1pc,list(X4D=X4new$data)),argvals=tj)

rr1[i,4]=1-sum(norm.fdata(ldatnew$yt-prff1.pc)^2)/rest
rr1[i,5]=1-sum(norm.fdata(ldatnew$yt-prff1)^2)/rest

########## FRegSigCom
XX=list(X4=X4$data)
XXnew=list(X4=X4new$data)
t.x.list=list(X4=tj)

#----- ---- ---- ----     cv.nonlinear
resSC.NL1=cv.nonlinear(XX, Y1D, t.x.list, tj,s.n.basis=15,x.n.basis=11,t.n.basis=15,
			upper.comp=max(npc4,4))
resSC.NL1.fit=fdata(pred.nonlinear(resSC.NL1, XX),tj)
prSCNL1=fdata(pred.nonlinear(resSC.NL1,XXnew),tj)
#----- ---- ---- ----     cv.sigcom
resSC.Lin1=cv.sigcom(XX, Y1D, t.x.list, tj,s.n.basis=15,t.n.basis=15,upper.comp=max(npc4,4))
resSC.Lin1.fit=fdata(pred.sigcom(resSC.Lin1, XX),tj)
prSCLin1=fdata(pred.sigcom(resSC.Lin1,XXnew),tj)
#
rr1[i,6]=1-sum(norm.fdata(ldatnew$yt-prSCLin1)^2)/rest
rr1[i,7]=1-sum(norm.fdata(ldatnew$yt-prSCNL1)^2)/rest

}
boxplot(rr1)
apply(rr1,2,median)

# Two covariates X4 y X2
set.seed(20030101)
for (i in 1:nrep){
index=sample(1:102,82)
print(paste0("Repetición:",i,"/",nrep))
yt=y[index];npcy=6
X1=x1[index];npc1=2
X2=x2[index];npc2=4
X3=x3[index];npc3=6
X4=x4[index];npc4=2

ytnew=y[-index]
X1new=x1[-index]
X2new=x2[-index]
X3new=x3[-index]
X4new=x4[-index]

ldatos=ldata(X1=X1,X2=X2,X3=X3,X4=X4,yt=yt)
ldatnew=ldata(X1=X1new,X2=X2new,X3=X3new,X4=X4new,yt=ytnew)

rest=sum(norm.fdata(ldatnew$yt-func.mean(ldatos$yt))^2)

####### fda.usc.devel
b.x=list(X1=create.pc.basis(X1,1:npc1),X2=create.pc.basis(X2,1:npc2),
X3=create.pc.basis(X3,1:npc3),X4=create.pc.basis(X4,1:npc4))
b.y=create.pc.basis(yt,1:npcy)
#
pmetric=list(df=data.frame(idx=1:nrow(yt)),yt=list(metric=metric.lp,lp=2), 
X1=list(metric=metric.lp,lp=2),X2=list(metric=metric.lp,lp=2),
X3=list(metric=metric.lp,lp=2),X4=list(metric=metric.lp,lp=2))
p.np=list(X1=list(Ker=AKer.norm),X2=list(Ker=AKer.norm),
		  X3=list(Ker=AKer.norm),X4=list(Ker=AKer.norm))

reslin=fregre.mlm.fr(yt~X4+X2,data=ldatos,basis.y=b.y,basis.x=b.x)
ressam=fregre.sam.fr(yt~s(X4)+s(X2,k=5),data=ldatos,basis.y=b.y,basis.x=b.x)
reskam=fregre.kam.fr(yt~X4+X2,data=ldatos,par.metric=pmetric)


prlin=predict(reslin,ldatnew)
prsam=predict(ressam,ldatnew)
prkam=predict(reskam,ldatnew)

rr2[i,1]=1-sum(norm.fdata(ldatnew$yt-prlin)^2)/rest
rr2[i,2]=1-sum(norm.fdata(ldatnew$yt-prsam)^2)/rest
rr2[i,3]=1-sum(norm.fdata(ldatnew$yt-prkam)^2)/rest

##### Refund
Y1D=yt$data;X1D=X1$data;X2D=X2$data;X3D=X3$data;X4D=X4$data

respff1=pffr(Y1D~sff(X4D,xind=tj,splinepars=list(bs="ps",m=c(2,2,2),k=c(9,11,11)))
			    +sff(X2D,xind=tj,splinepars=list(bs="ps",m=c(2,2,2),k=c(9,9,9)))
				,yind=tj,bs.yindex=list(bs="ps",k=7,m=c(2,1)))

respff1.fit=fdata(matrix(respff1$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)

respff1pc=pffr(Y1D~ffpc(X4D,xind=tj,npc.max=npc4,splinepars=list(k=11))
                  +ffpc(X2D,xind=tj,npc.max=npc2,splinepars=list(k=11))
					,yind=tj)
respff1pc.fit=fdata(matrix(respff1pc$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)

prff1=fdata(predict(respff1,list(X4D=X4new$data,X2D=X2new$data)),argvals=tj)

prff1.pc=fdata(predict(respff1pc,list(X4D=X4new$data,X2D=X2new$data)),argvals=tj)

rr2[i,4]=1-sum(norm.fdata(ldatnew$yt-prff1.pc)^2)/rest
rr2[i,5]=1-sum(norm.fdata(ldatnew$yt-prff1)^2)/rest

########## FRegSigCom
XX=list(X4=X4$data,X2=X2$data)
XXnew=list(X4=X4new$data,X2=X2new$data)
t.x.list=list(X4=tj,X2=tj)

#----- ---- ---- ----     cv.nonlinear
resSC.NL1=cv.nonlinear(XX, Y1D, t.x.list, tj,s.n.basis=15,x.n.basis=11,t.n.basis=15,
			upper.comp=max(npc4,npc2,4))
resSC.NL1.fit=fdata(pred.nonlinear(resSC.NL1, XX),tj)
prSCNL1=fdata(pred.nonlinear(resSC.NL1,XXnew),tj)
#----- ---- ---- ----     cv.sigcom
resSC.Lin1=cv.sigcom(XX, Y1D, t.x.list, tj,s.n.basis=15,t.n.basis=15,upper.comp=max(npc4,npc2,4))
resSC.Lin1.fit=fdata(pred.sigcom(resSC.Lin1, XX),tj)
prSCLin1=fdata(pred.sigcom(resSC.Lin1,XXnew),tj)
#
rr2[i,6]=1-sum(norm.fdata(ldatnew$yt-prSCLin1)^2)/rest
rr2[i,7]=1-sum(norm.fdata(ldatnew$yt-prSCNL1)^2)/rest

}
boxplot(rr2)
apply(rr2,2,median)

# Three covariates X4, X2, X3
set.seed(20030101)
for (i in 1:nrep){
index=sample(1:102,82)
print(paste0("Repetición:",i,"/",nrep))
yt=y[index];npcy=6
X1=x1[index];npc1=2
X2=x2[index];npc2=4
X3=x3[index];npc3=6
X4=x4[index];npc4=2

ytnew=y[-index]
X1new=x1[-index]
X2new=x2[-index]
X3new=x3[-index]
X4new=x4[-index]

ldatos=ldata(X1=X1,X2=X2,X3=X3,X4=X4,yt=yt)
ldatnew=ldata(X1=X1new,X2=X2new,X3=X3new,X4=X4new,yt=ytnew)

rest=sum(norm.fdata(ldatnew$yt-func.mean(ldatos$yt))^2)

####### fda.usc.devel
b.x=list(X1=create.pc.basis(X1,1:npc1),X2=create.pc.basis(X2,1:npc2),
X3=create.pc.basis(X3,1:npc3),X4=create.pc.basis(X4,1:npc4))
b.y=create.pc.basis(yt,1:npcy)
#
pmetric=list(df=data.frame(idx=1:nrow(yt)),yt=list(metric=metric.lp,lp=2), 
X1=list(metric=metric.lp,lp=2),X2=list(metric=metric.lp,lp=2),
X3=list(metric=metric.lp,lp=2),X4=list(metric=metric.lp,lp=2))
p.np=list(X1=list(Ker=AKer.norm),X2=list(Ker=AKer.norm),
		  X3=list(Ker=AKer.norm),X4=list(Ker=AKer.norm))

reslin=fregre.mlm.fr(yt~X4+X2+X3,data=ldatos,basis.y=b.y,basis.x=b.x)
ressam=fregre.sam.fr(yt~s(X4,k=5)+s(X2,k=5)+s(X3,k=5),data=ldatos,basis.y=b.y,basis.x=b.x)
reskam=fregre.kam.fr(yt~X4+X2+X3,data=ldatos,par.metric=pmetric)

prlin=predict(reslin,ldatnew)
prsam=predict(ressam,ldatnew)
prkam=predict(reskam,ldatnew)

rr3[i,1]=1-sum(norm.fdata(ldatnew$yt-prlin)^2)/rest
rr3[i,2]=1-sum(norm.fdata(ldatnew$yt-prsam)^2)/rest
rr3[i,3]=1-sum(norm.fdata(ldatnew$yt-prkam)^2)/rest

##### Refund
Y1D=yt$data;X1D=X1$data;X2D=X2$data;X3D=X3$data;X4D=X4$data

respff1=pffr(Y1D~sff(X4D,xind=tj,splinepars=list(bs="ps",m=c(2,2,2),k=c(9,11,11)))
			    +sff(X2D,xind=tj,splinepars=list(bs="ps",m=c(2,2,2),k=c(9,9,9)))
                 +sff(X3D,xind=tj,splinepars=list(bs="ps",m=c(2,2,2),k=c(5,5,5)))
				,yind=tj,bs.yindex=list(bs="ps",k=7,m=c(2,1)))

respff1.fit=fdata(matrix(respff1$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)

respff1pc=pffr(Y1D~ffpc(X4D,xind=tj,npc.max=npc4,splinepars=list(k=11))
                  +ffpc(X2D,xind=tj,npc.max=npc2,splinepars=list(k=11))
                  +ffpc(X3D,xind=tj,npc.max=npc3,splinepars=list(k=11))
					,yind=tj)
respff1pc.fit=fdata(matrix(respff1pc$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj)

prff1=fdata(predict(respff1,list(X4D=X4new$data,X2D=X2new$data,X3D=X3new$data)),argvals=tj)

prff1.pc=fdata(predict(respff1pc,list(X4D=X4new$data,X2D=X2new$data,X3D=X3new$data)),argvals=tj)

rr3[i,4]=1-sum(norm.fdata(ldatnew$yt-prff1.pc)^2)/rest
rr3[i,5]=1-sum(norm.fdata(ldatnew$yt-prff1)^2)/rest

########## FRegSigCom
XX=list(X4=X4$data,X2=X2$data,X3=X3$data)
XXnew=list(X4=X4new$data,X2=X2new$data,X3=X3new$data)
t.x.list=list(X4=tj,X2=tj,X3=tj)

#----- ---- ---- ----     cv.nonlinear
resSC.NL1=cv.nonlinear(XX, Y1D, t.x.list, tj,s.n.basis=15,x.n.basis=11,t.n.basis=15,
			upper.comp=max(npc4,npc2,npc3,4))
resSC.NL1.fit=fdata(pred.nonlinear(resSC.NL1, XX),tj)
prSCNL1=fdata(pred.nonlinear(resSC.NL1,XXnew),tj)
#----- ---- ---- ----     cv.sigcom
resSC.Lin1=cv.sigcom(XX, Y1D, t.x.list, tj,s.n.basis=15,t.n.basis=15,upper.comp=max(npc4,npc2,npc3,4))
resSC.Lin1.fit=fdata(pred.sigcom(resSC.Lin1, XX),tj)
prSCLin1=fdata(pred.sigcom(resSC.Lin1,XXnew),tj)
#
rr3[i,6]=1-sum(norm.fdata(ldatnew$yt-prSCLin1)^2)/rest
rr3[i,7]=1-sum(norm.fdata(ldatnew$yt-prSCNL1)^2)/rest

}
boxplot(rr3)
apply(rr3,2,median)

