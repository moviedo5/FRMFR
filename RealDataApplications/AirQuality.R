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
al <- lapply(listfinal[-1],function(v) which(fda.usc.devel:::is.na.fdata(v)))
al <- unique(rapply(al,drop))
listfinal <- subset(listfinal,-al)
dim(listfinal$df)

par(mfrow=c(2,3))
plot(listfinal$NMHC,col="grey")
lines(func.mean(listfinal$NMHC))
plot(listfinal$CO,col="grey")
lines(func.mean(listfinal$CO))
plot(listfinal$NOx,col="grey")
lines(func.mean(listfinal$NOx))
plot(listfinal$O3,col="grey")
lines(func.mean(listfinal$O3))
plot(AirQuality$RH,col="grey")
lines(func.mean(listfinal$RH))
plot(listfinal$C6H6,col="grey")
lines(func.mean(listfinal$C6H6))
par(mfrow=c(1,1))

mdcor <- diag(length(listfinal)-1)
colnames(mdcor) <- names(listfinal)[-1]
rownames(mdcor) <- colnames(mdcor)

for (i in 2:(length(listfinal)-1)){
 for (j in (i+1):length(listfinal)){
 mdcor[i-1,j-1] <- dcor.xy(listfinal[[i]],listfinal[[j]])$estimate
 mdcor[j-1,i-1] <- mdcor[i-1,j-1]
 }
}
#names(listfinal)
round(mdcor,3)
irow <- c(3,2,4,6,8)
icol <- c(1,3,2,4,6)
round(mdcor[irow,icol],3)
################################################################################
# Steps to select M2 and M3
#mod.lm.fr <- fregre.mlm.fr(C6H6~NMHC,data=ltrain,basis.y=b.y,basis.x=b.x)
#mod.sam.f <- fregre.sam.fr(C6H6~s(NMHC),data=ltrain,basis.y=b.y,basis.x=b.x)
#mod.kam.fr <- fregre.kam.fr(C6H6~NMHC,data=ltrain,par.metric=pmetric,par.np=par.np,control=list(trace=tr))

#for (i in 3:length(ltrain)){print(paste0(names(ltrain)[i],":",round(dcor.xy(mod.lm.fr$residuals,listfinal[[i]])$estimate,3)))}
#for (i in 3:length(ltrain)){print(paste0(names(listfinal)[i],":",round(dcor.xy(mod.sam.fr$residuals,listfinal[[i]])$estimate,3)))}
#for (i in 3:length(ltrain)){print(paste0(names(listfinal)[i],":",round(dcor.xy(mod.kam.fr$residuals,listfinal[[i]])$estimate,3)))}

#mod.lm.fr <- fregre.mlm.fr(C6H6~NMHC+O3+NOx,data=listfinal,basis.y=b.y,basis.x=b.x)
#mod.sam.fr <- fregre.sam.fr(C6H6~s(NMHC)+s(O3)+s(NOx),data=listfinal,basis.y=b.y,basis.x=b.x)
#mod.kam.fr <- fregre.kam.fr(C6H6~NMHC+O3+NOx,data=listfinal,par.metric=pmetric,par.np=par.np,control=list(trace=tr))

#for (i in 3:length(listfinal)){print(paste0(names(listfinal)[i],":",round(dcor.xy(mod.lm.fr$residuals,listfinal[[i]])$estimate,3)))}
#for (i in 3:length(listfinal)){print(paste0(names(listfinal)[i],":",round(dcor.xy(mod.sam.fr$residuals,listfinal[[i]])$estimate,3)))}
#for (i in 3:length(listfinal)){print(paste0(names(listfinal)[i],":",round(dcor.xy(mod.kam.fr$residuals,listfinal[[i]])$estimate,3)))}

# for (i in 3:length(ltrain)){print(paste0(names(ltrain)[i],":",round(dcor.xy(mod.lm.fr$residuals,ltrain[[i]])$estimate,3)))}
# for (i in 3:length(listfinal)){print(paste0(names(ltrain)[i],":",round(dcor.xy(mod.sam.fr$residuals,ltrain[[i]])$estimate,3)))}
# for (i in 3:length(listfinal)){print(paste0(names(ltrain)[i],":",round(dcor.xy(mod.kam.fr$residuals,ltrain[[i]])$estimate,3)))}
################################################################################
tj <- 0:23 + .5  # Discretized points
rj <- c(0,24) # Range of discretized points

# Refund parameters
s.n.basis <- 11
t.n.basis <- 11
x.n.basis <- 11
nk <- 11
nkk <- c(11,11,11) # Model has more coefficients than data FAMM

# fregre.kam.fr parameter
tr <- FALSE

nrep <- 100
mtrain1 <- matrix(NA,ncol=7,nrow=nrep)
colnames(mtrain1) <- c("FLMFR","FSAMFR","FKAMFR","PFFR","FAMM","LSIG","DISC")
mtrain3 <- mtrain2 <- mpred3 <- mpred2 <- mpred1 <- mtrain1

compModel <- rep(TRUE,len=7)
names(compModel) <- colnames(mtrain1)
compModel["FAMM"] <- FALSE # Time consuming
compModel <- !compModel  # Only FAMM

compModel[1:7]<- FALSE
compModel[4]<- TRUE

compModel

n <- nrow(listfinal$df)

set.seed(20030101)
for (i in 1:nrep){
  # i <- 1
  ltr <- sample(1:n,floor(0.75*n))
  
  lgtr <- 1:n %in% ltr
  
  ltrain <- subset(listfinal,lgtr)
  lpred <- subset(listfinal,!lgtr)
  
  ytrain <- ltrain$C6H6
  ypred <- lpred$C6H6
  
  mutr <- func.mean(ltrain$C6H6)
  MSEtr <- mean(norm.fdata(ltrain$C6H6-mutr)^2)
  MSEpr <- mean(norm.fdata(lpred$C6H6-mutr)^2)
  print(paste0("Rep.:",i," MSEtr:",round(MSEtr,3), " MSEpr:",round(MSEpr,3)))
  
  if (any(compModel[c("FLMFR","FSAMFR")])){
    b.y <- create.pc.basis(ltrain$C6H6)
    b.x <- list(CO=create.pc.basis(ltrain$CO,1:5),
                NMHC=create.pc.basis(ltrain$NMHC,1:5),
          	    NOx=create.pc.basis(ltrain$NOx,1:5),
          		  NO2=create.pc.basis(ltrain$NO2,1:4),
          		  O3=create.pc.basis(ltrain$O3,1:5),
          		  Temp=create.pc.basis(ltrain$Temp,1:2),
          		  RH=create.pc.basis(ltrain$RH,1:3)
          		 )
  }
  
  # Model 1 fda.usc
  form <- C6H6~NMHC
  if (compModel["FLMFR"]){
    mod.lm.fr <- fregre.mlm.fr(form, data=ltrain, basis.y=b.y, basis.x=b.x)
    predlm <- predict(mod.lm.fr,lpred)
    mtrain1[i,1] <- 1 - mean(norm.fdata(ytrain - mod.lm.fr$fitted.values)^2)/MSEtr
    mpred1[i,1] <- 1 - mean(norm.fdata(ypred - predlm)^2)/MSEpr
  }
  
  if (compModel["FSAMFR"]){
      forms <- C6H6~s(NMHC)
      mod.sam.fr <- fregre.sam.fr(forms, data=ltrain,basis.y=b.y,basis.x=b.x)
      predsam <- predict(mod.sam.fr,lpred)
      mtrain1[i,2] <- 1 - mean(norm.fdata(ytrain - mod.sam.fr$fitted.values)^2)/MSEtr
      mpred1[i,2] <- 1 - mean(norm.fdata(ypred - predsam)^2)/MSEpr
  }
  
  if (compModel["FKAMFR"]){
    
    pmetric <- list(C6H6=list(metric=metric.lp,lp=2), CO=list(metric=metric.lp),
                    NMHC=list(metric=metric.lp,lp=2),NOx=list(metric=metric.lp),	
                    NO2=list(metric=metric.lp),O3=list(metric=metric.lp),
                    Temp=list(metric=metric.lp),RH=list(metric=metric.lp))
    lprob <- c(seq(0.005,.1,len=51),.2,0.3,.4,0.5,.6,.7)	
    par.np <- list(CO=list(h=c(h.default(ltrain$CO,prob=lprob),Inf)),
                   NMHC=list(h=c(h.default(ltrain$NMHC,prob=lprob),Inf)),
                   NOx=list(h=c(h.default(ltrain$NOx,prob=lprob),Inf)),
                   NO2=list(h=c(h.default(ltrain$NO2,prob=lprob),Inf)), 
                   O3=list(h=c(h.default(ltrain$O3,prob=lprob),Inf)),
                   Temp=list(h=c(h.default(ltrain$Temp,prob=lprob),Inf)),
                   RH=list(h=c(h.default(ltrain$RH,prob=lprob),Inf)))
    
    mod.kam.fr <- fregre.kam.fr(form, data=ltrain,
                               par.metric=pmetric,par.np=par.np,
                               control=list(trace=tr))
    predkam <- predict(mod.kam.fr,lpred)
    mtrain1[i,3] <- 1 - mean(norm.fdata(ytrain - mod.kam.fr$fitted.values)^2)/MSEtr
    mpred1[i,3] <- 1 - mean(norm.fdata(ypred - predkam)^2)/MSEpr
  }

  # Model 2 fda.usc
  form <- C6H6~NMHC+O3+CO
  if (compModel["FLMFR"]){
    mod.lm.fr <- fregre.mlm.fr(form, data=ltrain, basis.y=b.y, basis.x=b.x)
    predlm <- predict(mod.lm.fr,lpred)
    mtrain2[i,1] <- 1 - mean(norm.fdata(ytrain - mod.lm.fr$fitted.values)^2)/MSEtr
    mpred2[i,1] <- 1 - mean(norm.fdata(ypred - predlm)^2)/MSEpr
  }
  
  if (compModel["FSAMFR"]){
    forms <- C6H6~s(NMHC)+s(O3)+s(CO)
    mod.sam.fr <- fregre.sam.fr(forms, data=ltrain,basis.y=b.y,basis.x=b.x)
    predsam <- predict(mod.sam.fr,lpred)
    mtrain2[i,2] <- 1 - mean(norm.fdata(ytrain - mod.sam.fr$fitted.values)^2)/MSEtr
    mpred2[i,2] <- 1 - mean(norm.fdata(ypred - predsam)^2)/MSEpr
  }
  
  if (compModel["FKAMFR"]){
    mod.kam.fr <- fregre.kam.fr(form, data=ltrain,
                                par.metric=pmetric,par.np=par.np,
                                control=list(trace=tr))
    predkam <- predict(mod.kam.fr,lpred)
    mtrain2[i,3] <- 1 - mean(norm.fdata(ytrain - mod.kam.fr$fitted.values)^2)/MSEtr
    mpred2[i,3] <- 1 - mean(norm.fdata(ypred - predkam)^2)/MSEpr
  }
  
  # Model 3 fda.usc
  form <- C6H6~NMHC+O3+CO+NO2+NOx
  if (compModel["FLMFR"]){
    mod.lm.fr <- fregre.mlm.fr(form, data=ltrain, basis.y=b.y, basis.x=b.x)
    predlm <- predict(mod.lm.fr,lpred)
    mtrain3[i,1] <- 1 - mean(norm.fdata(ytrain - mod.lm.fr$fitted.values)^2)/MSEtr
    mpred3[i,1] <- 1 - mean(norm.fdata(ypred - predlm)^2)/MSEpr
  }
  
  if (compModel["FSAMFR"]){
    forms <- C6H6~s(NMHC)+s(O3)+s(CO)+s(NO2)+s(NOx)
    mod.sam.fr <- fregre.sam.fr(forms, data=ltrain,basis.y=b.y,basis.x=b.x)
    predsam <- predict(mod.sam.fr,lpred)
    mtrain3[i,2] <- 1 - mean(norm.fdata(ytrain - mod.sam.fr$fitted.values)^2)/MSEtr
    mpred3[i,2] <- 1 - mean(norm.fdata(ypred - predsam)^2)/MSEpr
  }
  
  if (compModel["FKAMFR"]){
    mod.kam.fr <- fregre.kam.fr(form, data=ltrain,
                                par.metric=pmetric,par.np=par.np,
                                control=list(trace=tr))
    predkam <- predict(mod.kam.fr,lpred)
    mtrain3[i,3] <- 1 - mean(norm.fdata(ytrain - mod.kam.fr$fitted.values)^2)/MSEtr
    mpred3[i,3] <- 1 - mean(norm.fdata(ypred - predkam)^2)/MSEpr
  }
  
  # refund
  if (any(compModel[c("PFFR","FAMM","LSIG","DISC")])){
    YD <- ltrain$C6H6$data
    X1 <- ltrain$NMHC$data
    X2 <- ltrain$O3$data
    X3 <- ltrain$CO$data
    X4 <- ltrain$NO2$data
    X5 <- ltrain$NOx$data
  }  
  
  # Model 1 refund
  if (compModel["PFFR"]){
    mod.pffr <- pffr(YD~ffpc(X1,xind=tj,npc.max=5),yind=tj)
    predsff <- fdata(predict(mod.pffr,list(X1=lpred$NMHC$data,X2=lpred$O3$data,
                     X3=lpred$CO$data,X4=lpred$NO2$data,X5=lpred$NOx$data)),
                     argvals=tj,rj)
    predpffr <- fdata(predict(mod.pffr,list(X1=lpred$NMHC$data)),argvals=tj,rj)
   
    fitt <- fdata(matrix(mod.pffr$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj,rj)
    mtrain1[i,4] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
    mpred1[i,4] <- 1 - mean(norm.fdata(ypred - predpffr)^2)/MSEpr
   
    # Model 2 refund
    mod.pffr <- pffr(YD~ffpc(X1,xind=tj,npc.max=5)                    
    			            +ffpc(X2,xind=tj,npc.max=5)+ffpc(X3,xind=tj,npc.max=5)
    			            ,yind=tj)
    predpffr <- fdata(predict(mod.pffr,list(X1=lpred$NMHC$data,X2=lpred$O3$data,
                                              X3=lpred$CO$data)),argvals=tj,rj)
    
    fitt <- fdata(matrix(mod.pffr$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj,rj)
    mtrain2[i,4] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
    mpred2[i,4] <- 1 - mean(norm.fdata(ypred - predpffr)^2)/MSEpr
    
    # Model 3 refund  
    mod.pffr <- pffr(YD~ffpc(X1,xind=tj,npc.max=5)                    
                      +ffpc(X2,xind=tj,npc.max=5)+ffpc(X3,xind=tj,npc.max=5)
                			+ffpc(X4,xind=tj,npc.max=4)+ffpc(X5,xind=tj,npc.max=5)
                			,yind=tj)
    predpffr <- fdata(predict(mod.pffr,list(X1=lpred$NMHC$data,X2=lpred$O3$data,
                                            X3=lpred$CO$data,X4=lpred$NO2$data,
                                            X5=lpred$NOx$data)),argvals=tj,rj)
    fitt <- fdata(matrix(mod.pffr$fitted.values,ncol=length(tj),byrow=TRUE),argvals=tj,rj)
    mtrain3[i,4] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
    mpred3[i,4] <- 1 - mean(norm.fdata(ypred - predpffr)^2)/MSEpr
  }
 
  if (compModel["FAMM"]){
    # Model 1 refund  
    mod.pffnl <- pffr(YD ~ sff(X1,xind=tj,splinepars = list(bs = "ps",
                                                            m = c(2, 2, 2), 
                                                            k=nkk))
                       ,yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
                       )
    fitt <- fdata(matrix(mod.pffnl$fitted.values,
                         ncol=length(tj),byrow=TRUE),argvals=tj,rj)
    
    predsff <- fdata(predict(mod.pffnl,list(X1=lpred$NMHC$data)),
                     argvals=tj,rj)
    mtrain1[i,5] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
    mpred1[i,5] <- 1 - mean(norm.fdata(ypred - predsff)^2)/MSEpr
    
    # Model 2 refund  
     mod.pffnl <- pffr(YD ~ sff(X1,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
                           + sff(X2,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
                           + sff(X3,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
                        , yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
                        )
     fitt <- fdata(matrix(mod.pffnl$fitted.values,
                          ncol=length(tj),byrow=TRUE),argvals=tj,rj)
     predsff <- fdata(predict(mod.pffnl,list(X1=lpred$NMHC$data,X2=lpred$O3$data,
                                             X3=lpred$CO$data)),argvals=tj,rj)
     mtrain2[i,5] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
     mpred2[i,5] <- 1 - mean(norm.fdata(ypred - predsff)^2)/MSEpr
     
    # # Model 3 refund  
     mod.pffnl <- pffr(YD ~ sff(X1,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
                         	+ sff(X2,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
                           + sff(X3,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
                           + sff(X4,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
                           + sff(X5,xind=tj,splinepars = list(bs = "ps",m = c(2, 2, 2), k=nkk))
                        ,yind=tj,bs.yindex=list(bs="ps",k=nk,m=c(2,1))
                        )
     fitt <- fdata(matrix(mod.pffnl$fitted.values,
    #                      ncol=length(tj),byrow=TRUE),argvals=tj,rj)
     predsff <- fdata(predict(mod.pffnl,list(X1=lpred$NMHC$data,X2=lpred$O3$data,
         X3=lpred$CO$data,X4=lpred$NO2$data,X5=lpred$NOx$data)),argvals=tj,rj)
     mtrain3[i,5] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
     mpred3[i,5] <- 1 - mean(norm.fdata(ypred - predsff)^2)/MSEpr
  }  
  
  if (compModel["LSIG"]){
    # Model 1 SigCom
    XX <- list(X1=X1)
    XXnew <- list(X1=lpred$NMHC$data)
    t.x.list <- list(X1=tj)
    mod.SCLin <- cv.sigcom(XX,YD,t.x.list,tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
    predSCLin <- fdata(pred.sigcom(mod.SCLin,XXnew),tj,rj)
    fitt <- fdata(pred.sigcom(mod.SCLin,XX),tj,rj)
    mtrain1[i,6] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
    mpred1[i,6] <- 1 - mean(norm.fdata(ypred - predSCLin)^2)/MSEpr
    
    # Model 2 SigCom
    XX <- list(X1=X1,X2=X2,X3=X3)
    XXnew <- list(X1=lpred$NMHC$data, X2=lpred$O3$data, X3=lpred$CO$data)
    		  
    t.x.list <- list(X1=tj, X2=tj, X3=tj)
    mod.SCLin <- cv.sigcom(XX,YD,t.x.list,tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
    predSCLin <- fdata(pred.sigcom(mod.SCLin,XXnew),tj,rj)
    fitt <- fdata(pred.sigcom(mod.SCLin,XX),tj,rj)
    mtrain2[i,6] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
    mpred2[i,6] <- 1 - mean(norm.fdata(ypred - predSCLin)^2)/MSEpr
    
    # Model 3 SigCom
    XX <- list(X1=X1, X2=X2, X3=X3, X4=X4, X5=X5)
    XXnew <- list(X1=lpred$NMHC$data, X2=lpred$O3$data, X3=lpred$CO$data,
                  X4=lpred$NO2$data, X5=lpred$NOx$data)
    
    t.x.list <- list(X1=tj, X2=tj, X3=tj, X4=tj, X5=tj)
    mod.SCLin <- cv.sigcom(XX,YD,t.x.list,tj,s.n.basis=s.n.basis,t.n.basis=t.n.basis)
    predSCLin <- fdata(pred.sigcom(mod.SCLin,XXnew),tj,rj)
    fitt <- fdata(pred.sigcom(mod.SCLin,XX),tj,rj)
    mtrain3[i,6] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
    mpred3[i,6] <- 1 - mean(norm.fdata(ypred - predSCLin)^2)/MSEpr
  }
  
  if (compModel["DISC"]){
    # Model 1 SigCom
    XX <- list(X1=X1)
    XXnew <- list(X1=lpred$NMHC$data)
    t.x.list <- list(X1=tj)
    mod.SCNL <- cv.nonlinear(XX,YD,t.x.list,tj,s.n.basis=s.n.basis,
                             x.n.basis=x.n.basis,t.n.basis=t.n.basis)
    predSCNL <- fdata(pred.nonlinear(mod.SCNL,XXnew),tj,rj)  
    fitt <- fdata(pred.nonlinear(mod.SCNL,XX),tj,rj)
    mtrain1[i,7] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
    mpred1[i,7] <- 1 - mean(norm.fdata(ypred - predSCNL)^2)/MSEpr
    
    # Model 2 SigCom
    XX <- list(X1=X1,X2=X2,X3=X3)
    XXnew <- list(X1=lpred$NMHC$data, X2=lpred$O3$data, X3=lpred$CO$data)
    
    t.x.list <- list(X1=tj, X2=tj, X3=tj)
    mod.SCNL <- cv.nonlinear(XX,YD,t.x.list,tj,s.n.basis=s.n.basis,
                             x.n.basis=x.n.basis,t.n.basis=t.n.basis)
    
    predSCNL <- fdata(pred.nonlinear(mod.SCNL,XXnew),tj,rj)  
    fitt <- fdata(pred.nonlinear(mod.SCNL,XX),tj,rj)
    mtrain2[i,7] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
    mpred2[i,7] <- 1 - mean(norm.fdata(ypred - predSCNL)^2)/MSEpr
    
    # Model 3 SigCom
    XX <- list(X1=X1, X2=X2, X3=X3, X4=X4, X5=X5)
    XXnew <- list(X1=lpred$NMHC$data, X2=lpred$O3$data, X3=lpred$CO$data,
                  X4=lpred$NO2$data, X5=lpred$NOx$data)
    
    t.x.list <- list(X1=tj, X2=tj, X3=tj, X4=tj, X5=tj)
    mod.SCNL <- cv.nonlinear(XX,YD,t.x.list,tj,s.n.basis=s.n.basis,
                             x.n.basis=x.n.basis,t.n.basis=t.n.basis)
    mod.SCNL <- cv.nonlinear(XX,YD,t.x.list,tj,s.n.basis=s.n.basis,
                             x.n.basis=x.n.basis,t.n.basis=t.n.basis)
    predSCNL <- fdata(pred.nonlinear(mod.SCNL,XXnew),tj,rj)  
    fitt <- fdata(pred.nonlinear(mod.SCNL,XX),tj,rj)
    mtrain3[i,7] <- 1 - mean(norm.fdata(ytrain - fitt)^2)/MSEtr
    mpred3[i,7] <- 1 - mean(norm.fdata(ypred - predSCNL)^2)/MSEpr
  }
  print(i)
  print(round(mpred1[i,],2))
  print(round(mpred2[i,],2))
  print(round(mpred3[i,],2))
}


# M1=X1
# M2=X1,X2,X3
# M3=X1,X2;X3,X4,X5
round(colMeans(mtrain1),3)
round(colMeans(mtrain2),3)
round(colMeans(mtrain3),3)

# Table 13
round(colMeans(mpred1),3)
round(colMeans(mpred2),3)
round(colMeans(mpred3),3)

# Table 13 (paper)
# Model FLMFR FSAMFR FKAMFR PFR FAMM LSC DISC #PAPER
# M1 0.889 0.900 0.848 0.855 0.815 0.923 0.920
# M2 0.892 0.901 0.858 0.864 0.817 0.926 0.912
# M3 0.891 0.900 0.859 0.862 0.815 0.924 0.899


