#---------------------------
library(fda.usc.devel)
library(FRegSigCom)
library(refund)
#---------------------------

#=====================         load data air      ===============      
data("AirQuality")
names(AirQuality)
 
listfinal <- AirQuality
dim(listfinal$df)
al=lapply(listfinal[-1],function(v) which(fda.usc.devel:::is.na.fdata(v)))
al=unique(rapply(al,drop))
listfinal=subset(listfinal,-al)
dim(listfinal$df)

par(mfrow=c(2,3))
plot(listfinal$NMHC,col="grey")
lines(func.mean(listfinal$NMHC))
plot(listfinal$CO,col="grey")
lines(func.mean(listfinal$CO))
plot(listfinal$NOx,col="grey")
lines(func.mean(listfinal$NMHC))
plot(listfinal$O3,col="grey")
lines(func.mean(listfinal$O3))
plot(AirQuality$RH,col="grey")
lines(func.mean(listfinal$RH))
plot(listfinal$C6H6,col="grey")
lines(func.mean(listfinal$C6H6))

mdcor=diag(length(listfinal)-1)
colnames(mdcor)=names(listfinal)[-1]
rownames(mdcor)=colnames(mdcor)
for (i in 2:(length(listfinal)-1)){
 for (j in (i+1):length(listfinal)){
 mdcor[i-1,j-1]=dcor.xy(listfinal[[i]],listfinal[[j]])$estimate
 mdcor[j-1,i-1]=mdcor[i-1,j-1]
 }
}
#names(listfinal)
round(mdcor,3)
irow <- c(3,2,4,6,8)
icol <- c(1,3,2,4,6)
round(mdcor[irow,icol],3)
################################################################################

computeFAMM <- TRUE #Time consuming

nrep=100
mtrain=matrix(NA,ncol=7,nrow=nrep)
mpred=matrix(NA,ncol=7,nrow=nrep)
colnames(mpred)=c("FLMFR","FSAMFR","FKAMFR","PFFR","FAMM","LSIG","DISC")
colnames(mtrain)=colnames(mpred)
n <- nrow(listfinal$df)
for (i in 1:nrep){
  # i <- 1
  ltr=sample(1:n,floor(0.75*n))
  
  lgtr=1:n %in% ltr
  
  ltrain=subset(listfinal,lgtr)
  lpred=subset(listfinal,!lgtr)
  mutr=func.mean(ltrain$C6H6)
  MSEtr=mean(norm.fdata(ltrain$C6H6-mutr)^2)
  MSEpr=mean(norm.fdata(lpred$C6H6-mutr)^2)
  print(paste0("Rep.:",i," MSEtr:",round(MSEtr,3), " MSEpr:",round(MSEpr,3)))
  
  b.y=create.pc.basis(ltrain$C6H6)
  b.x=list(CO=create.pc.basis(ltrain$CO,1:5),NMHC=create.pc.basis(ltrain$NMHC,1:5),
  		 NOx=create.pc.basis(ltrain$NOx,1:5),NO2=create.pc.basis(ltrain$NO2,1:4),
  		 O3=create.pc.basis(ltrain$O3,1:5),Temp=create.pc.basis(ltrain$Temp,1:2),RH=create.pc.basis(ltrain$RH,1:3))
  			 
  pmetric=list(C6H6=list(metric=metric.lp,lp=2), CO=list(metric=metric.lp),
  	NMHC=list(metric=metric.lp,lp=2),NOx=list(metric=metric.lp),	
      NO2=list(metric=metric.lp),O3=list(metric=metric.lp),Temp=list(metric=metric.lp),RH=list(metric=metric.lp))
  lprob=c(seq(0.005,.1,len=51),.2,0.3,.4,0.5,.6,.7)	
  par.np=list(CO=list(h=c(h.default(ltrain$CO,prob=lprob),Inf)),
  	NMHC=list(h=c(h.default(ltrain$NMHC,prob=lprob),Inf)),
  	NOx=list(h=c(h.default(ltrain$NOx,prob=lprob),Inf)),
  	NO2=list(h=c(h.default(ltrain$NO2,prob=lprob),Inf)), 
  	O3=list(h=c(h.default(ltrain$O3,prob=lprob),Inf)),
  	Temp=list(h=c(h.default(ltrain$Temp,prob=lprob),Inf)),
  	RH=list(h=c(h.default(ltrain$RH,prob=lprob),Inf)))
  	
  #mod.lm.fr=fregre.mlm.fr(C6H6~NMHC,data=ltrain,basis.y=b.y,basis.x=b.x)
  #mod.sam.fr=fregre.sam.fr(C6H6~s(NMHC),data=ltrain,basis.y=b.y,basis.x=b.x)
  #mod.kam.fr=fregre.kam.fr(C6H6~NMHC,data=ltrain,par.metric=pmetric,par.np=par.np,control=list(trace=TRUE))
  
  #for (i in 3:length(ltrain)){print(paste0(names(ltrain)[i],":",round(dcor.xy(mod.lm.fr$residuals,listfinal[[i]])$estimate,3)))}
  #for (i in 3:length(ltrain)){print(paste0(names(listfinal)[i],":",round(dcor.xy(mod.sam.fr$residuals,listfinal[[i]])$estimate,3)))}
  #for (i in 3:length(ltrain)){print(paste0(names(listfinal)[i],":",round(dcor.xy(mod.kam.fr$residuals,listfinal[[i]])$estimate,3)))}
  
  #mod.lm.fr=fregre.mlm.fr(C6H6~NMHC+O3+NOx,data=listfinal,basis.y=b.y,basis.x=b.x)
  #mod.sam.fr=fregre.sam.fr(C6H6~s(NMHC)+s(O3)+s(NOx),data=listfinal,basis.y=b.y,basis.x=b.x)
  #mod.kam.fr=fregre.kam.fr(C6H6~NMHC+O3+NOx,data=listfinal,par.metric=pmetric,par.np=par.np,control=list(trace=TRUE))
  		 
  #for (i in 3:length(listfinal)){print(paste0(names(listfinal)[i],":",round(dcor.xy(mod.lm.fr$residuals,listfinal[[i]])$estimate,3)))}
  #for (i in 3:length(listfinal)){print(paste0(names(listfinal)[i],":",round(dcor.xy(mod.sam.fr$residuals,listfinal[[i]])$estimate,3)))}
  #for (i in 3:length(listfinal)){print(paste0(names(listfinal)[i],":",round(dcor.xy(mod.kam.fr$residuals,listfinal[[i]])$estimate,3)))}
  
  mod.lm.fr=fregre.mlm.fr(C6H6~NMHC,data=ltrain,basis.y=b.y,basis.x=b.x)
  #mod.lm.fr=fregre.mlm.fr(C6H6~NMHC+O3+CO,data=ltrain,basis.y=b.y,basis.x=b.x)
  #mod.lm.fr=fregre.mlm.fr(C6H6~NMHC+O3+CO+NO2+NOx,data=ltrain,basis.y=b.y,basis.x=b.x)
  
  mod.sam.fr=fregre.sam.fr(C6H6~s(NMHC),data=ltrain,basis.y=b.y,basis.x=b.x)
  #mod.sam.fr=fregre.sam.fr(C6H6~s(NMHC)+s(O3)+s(CO),data=ltrain,basis.y=b.y,basis.x=b.x)
  #mod.sam.fr=fregre.sam.fr(C6H6~s(NMHC)+s(O3)+s(CO)+s(NO2)+s(NOx),data=ltrain,basis.y=b.y,basis.x=b.x)
  
  mod.kam.fr=fregre.kam.fr(C6H6~NMHC,data=ltrain,par.metric=pmetric,par.np=par.np,control=list(trace=TRUE))
  #mod.kam.fr=fregre.kam.fr(C6H6~NMHC+O3+CO,data=ltrain,par.metric=pmetric,par.np=par.np,control=list(trace=TRUE))
  #mod.kam.fr=fregre.kam.fr(C6H6~NMHC+O3+CO+NO2+NOx,data=ltrain,par.metric=pmetric,par.np=par.np,control=list(trace=TRUE))
  		 
  #for (i in 3:length(ltrain)){print(paste0(names(ltrain)[i],":",round(dcor.xy(mod.lm.fr$residuals,ltrain[[i]])$estimate,3)))}
  #for (i in 3:length(listfinal)){print(paste0(names(ltrain)[i],":",round(dcor.xy(mod.sam.fr$residuals,ltrain[[i]])$estimate,3)))}
  #for (i in 3:length(listfinal)){print(paste0(names(ltrain)[i],":",round(dcor.xy(mod.kam.fr$residuals,ltrain[[i]])$estimate,3)))}
  
  predlm=predict(mod.lm.fr,lpred)
  predsam=predict(mod.sam.fr,lpred)
  predkam=predict(mod.kam.fr,lpred)
  
  mtrain[i,1]=1-mean(norm.fdata(ltrain$C6H6-mod.lm.fr$fitted.values)^2)/MSEtr
  mtrain[i,2]=1-mean(norm.fdata(ltrain$C6H6-mod.sam.fr$fitted.values)^2)/MSEtr
  mtrain[i,3]=1-mean(norm.fdata(ltrain$C6H6-mod.kam.fr$fitted.values)^2)/MSEtr
  
  mpred[i,1]=1-mean(norm.fdata(lpred$C6H6-predlm)^2)/MSEpr
  mpred[i,2]=1-mean(norm.fdata(lpred$C6H6-predsam)^2)/MSEpr
  mpred[i,3]=1-mean(norm.fdata(lpred$C6H6-predkam)^2)/MSEpr
    
  
  tj=0:23+.5
  s.n.basis=11;t.n.basis=11;x.n.basis=11
  nk=11;nkk=c(11,11)
  
  YD=ltrain$C6H6$data
  X1=ltrain$NMHC$data;X2=ltrain$O3$data;X3=ltrain$CO$data
  X4=ltrain$NO2$data;X5=ltrain$NOx$data
  
  mod.pffr=pffr(YD~ffpc(X1,xind=tj,npc.max=5)
  #			+ffpc(X2,xind=tj,npc.max=5)+ffpc(X3,xind=tj,npc.max=5)
  #			+ffpc(X4,xind=tj,npc.max=4)+ffpc(X5,xind=tj,npc.max=5)
  			,yind=tj)
  
  if (computeFAMM){
    mod.pffnl=pffr(YD~ sff(X1,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
    #				  +sff(X2,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
    #				  +sff(X3,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
    #				  +sff(X4,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
    #				  +sff(X5,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
    ,        yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
    			  )
    predsff=fdata(predict(mod.pffnl,list(X1=lpred$NMHC$data,X2=lpred$O3$data,
    		X3=lpred$CO$data,X4=lpred$NO2$data,X5=lpred$NOx$data)),argvals=tj,c(0,24))
    mtrain[i,5]=1-mean(norm.fdata(ltrain$C6H6-mod.sam.fr$fitted.values)^2)/MSEtr
    mpred[i,5]=1-mean(norm.fdata(lpred$C6H6-predsff)^2)/MSEpr
  }  
  
  predpffr=fdata(predict(mod.pffr,list(X1=lpred$NMHC$data,X2=lpred$O3$data,
  		X3=lpred$CO$data,X4=lpred$NO2$data,X5=lpred$NOx$data)),argvals=tj,c(0,24))
  
  fitt=fdata(matrix(mod.pffr$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj,c(0,24))
  mtrain[i,4]=1-mean(norm.fdata(ltrain$C6H6-fitt)^2)/MSEtr
  
  mpred[i,4]=1-mean(norm.fdata(lpred$C6H6-predpffr)^2)/MSEpr
  
  
  XX=list(X1=X1
  #		  ,X2=X2,X3=X3
  #		  ,X4=X4,X5=X5
  		  )
  XXnew=list(X1=lpred$NMHC$data
  #			,X2=lpred$O3$data,X3=lpred$CO$data
  #			,X4=lpred$NO2$data,X5=lpred$NOx$data
  			)
  		  
  t.x.list=list(X1=tj
  #			 ,X2=tj,X3=tj
  #			 ,X4=tj,X5=tj
  			 )
  
  mod.SCLin=cv.sigcom(XX,YD,t.x.list,tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
  mod.SCNL=cv.nonlinear(XX,YD,t.x.list,tj,s.n.basis=s.n.basis,
  				 x.n.basis=x.n.basis,t.n.basis=t.n.basis)
  			
  predSCLin=fdata(pred.sigcom(mod.SCLin,XXnew),tj,c(0,24))
  predSCNL=fdata(pred.nonlinear(mod.SCNL,XXnew),tj,c(0,24))  
  fitt3=fdata(pred.sigcom(mod.SCLin,XX),tj,c(0,24))
  fitt4=fdata(pred.nonlinear(mod.SCNL,XX),tj,c(0,24))
  mtrain[i,6]=1-mean(norm.fdata(ltrain$C6H6-fitt3)^2)/MSEtr
  mtrain[i,7]=1-mean(norm.fdata(ltrain$C6H6-fitt4)^2)/MSEtr
  
  mpred[i,6]=1-mean(norm.fdata(lpred$C6H6-predSCLin)^2)/MSEpr
  mpred[i,7]=1-mean(norm.fdata(lpred$C6H6-predSCNL)^2)/MSEpr
  print(i)
  print(mpred[i,])
}


# M1=X1
# M2=X1,X2,X3
# M3=X1,X2;X3,X4,X5
round(colMeans(mtrain),3)
round(colMeans(mpred),3)
#FLMFR FSAMFR FKAMFR   PFFR   FAMM   LSIG   DISC  PC-MOF
# 0.889  0.900  0.847  0.855  0.816  0.923  0.919 

# Model FLMFR FSAMFR FKAMFR PFR FAMM LSC DISC #PAPPER
# M1 0.889 0.900 0.848 0.855 0.815 0.923 0.920
# M2 0.892 0.901 0.858 0.864 0.817 0.926 0.912
# M3 0.891 0.900 0.859 0.862 0.815 0.924 0.899


