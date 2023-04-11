library(fda.usc.devel)
library(FRegSigCom)
library(refund)

best.estpred <- function(Y, 
                         Ynew, 
                         listYhat, 
                         listYpred){
  nam <- names(listYhat)
  Mat <- matrix(NA,nrow=nrow(Y),ncol=length(listYhat))
  Mpred <- matrix(NA,nrow=nrow(Ynew),ncol=length(listYpred))
  for (i in 1:length(listYhat)){
    Mat[,i] <- drop(norm.fdata(Y-listYhat[[i]])^2)
    Mpred[,i] <- drop(norm.fdata(Ynew-listYpred[[i]])^2)
  }
  colnames(Mat) <- nam
  colnames(Mpred) <- nam
  nullrss <- mean(norm.fdata(Y-func.mean(Y))^2)
  #nullpred <- mean(norm.fdata(Ynew-func.mean(Ynew))^2)
  nullpred <- mean(norm.fdata(Ynew-func.mean(Y))^2)
  R2 <- round(1-colMeans(Mat)/nullrss,3)
  names(R2) <- nam
  R2pr <- round(1-colMeans(Mpred)/nullpred,3)
  names(R2pr) <- nam
  return(list(est.R2=R2, pr.R2=R2pr, Err.est=Mat, Err.pred=Mpred))
}

best.estpred2 <- function(Y, 
                          Ynew, 
                          listYhat, 
                          listYpred){
  nam <- names(listYhat)
  Mat <- matrix(NA,nrow = nrow(Y),ncol=length(listYhat))
  Mpred <- matrix(NA,nrow=nrow(Ynew),ncol=length(listYpred))
  for (i in 1:length(listYhat)){
   Mat[,i] <- drop(norm.fdata(Y-listYhat[[i]])^2)
   Mpred[,i] <- drop(norm.fdata(Ynew-listYpred[[i]])^2)
  }
  colnames(Mat) <- nam
  colnames(Mpred) <- nam
  #nullrss <- mean(norm.fdata(Y-func.mean(Y))^2)
  nullrss <- drop(norm.fdata(Y-func.mean(Y)))
  #nullpred <- mean(norm.fdata(Ynew-func.mean(Ynew))^2)
  nullpred <- mean(norm.fdata(Ynew-func.mean(Y))^2)
  nullpred <- drop(norm.fdata(Ynew-func.mean(Y))^2)
  R2 <- round(1-colMeans(Mat/nullrss),3)
  names(R2) <- nam
  R2pr <- round(1-colMeans(Mpred/nullpred),3)
  names(R2pr) <- nam
  return(list(est.R2=R2 ,pr.R2=R2pr, Err.est=Mat, Err.pred=Mpred))
}

cteeps <- function(M,
                   error,
                   rho=0.9){
  (1-rho)*sum(norm.fdata(M-func.mean(M))^2)/(rho*sum(norm.fdata(error)^2))
  }

imodelfrfr <- function(imodel=1, 
                       ivar=2, 
                       N=200, 
                       Ntest=100, 
                       nT=71, 
                       nS=51,
                       rho=0.9){ 
  tj <- seq(0,1,len = nT)
  si <- seq(0,1,len = nS)
  delta.t <- diff(range(tj))/nT
  
  X1 <- rproc2fdata(N,si,sigma="OU")
  nc1 <- 4
  X1new <- rproc2fdata(Ntest,si,sigma="OU")
  X2 <- rproc2fdata(N,si,sigma="vexponential",par.list=list(scale=.5,theta=.7))
  nc2 <- 5
  X2new <- rproc2fdata(Ntest,si,sigma="vexponential",
                       par.list=list(scale=.5,theta=.7))
  
  TNLP  <-  function(X){X+exp(0.5-X^2/2)}
  TNL  <-  function(X){exp(0.5-X^2/2)}
  TNL1  <-  function(X){1+X^2/5+cos(2*pi*(X/2-1))*sin(2*pi*(X/2-1))}
  TNL2  <-  function(X){(X/3-1)*sin(2*pi*(X/3-1))}
  TNL3  <-  function(X){(X-1)*sin(2*pi*(X-1))}
  
  ML <- function(X,
                 beta,
                 tj){
    fdata(X$data%*%beta,tj)
    }
  MNL <- function(X,
                  TNL,
                  beta,
                  tj){
    fdata(TNL(X)$data%*%beta,tj)
    }
  
  epsilon  <-  rproc2fdata(N,tj,sigma="vexponential",
                           par.list=list(scale=0.05,theta=0.3))
  epsnew  <- rproc2fdata(Ntest,tj,sigma="vexponential",
                         par.list=list(scale=0.05,theta=0.3))
  
  #beta1 <- outer(si,tj,function(u,v){(exp(-(-0.5-(2*u-1)-(2*u-1)^3))+0.8* exp((1.8*v+0.4))/2.25)})
  #beta2 <- outer(si,tj,function(u,v){(-5*cos(2*pi*(2*u-1)*(1*v-1))*(1*v)^2)})
  
  # Also for easy setting we can consider this parameters 
  #beta1 <- outer(si,tj,function(u,v){(u+1)*cos(2*pi*v)})
  #beta2 <- outer(si,tj,function(u,v){6*sqrt(v*(u))*sin(4*pi*v)})
  etiqmodel <- c("E.L","H.L","ENL","HNL")
  print(paste0("Generating ",etiqmodel[imodel]))
  # 
  switch(imodel
   ,{  # -------EASY.L
  #beta1 <- outer(si,tj,function(u,v){(u+1)*cos(2*pi*v)})
  beta2 <- outer(si,tj,function(u,v){-(u*v+1)*cos(2*pi*sqrt(u*v))})
  beta1 <- outer(si,tj,function(u,v){6*sqrt(v*u)*sin(4*pi*v)})
  if (ivar==2) {
  MX <- ML(X1,beta1,tj)+ML(X2,beta2,tj)
  MXw <- ML(X1new,beta1,tj)+ML(X2new,beta2,tj)
   } else {
  MX <- ML(X1,beta1,tj)
  MXw <- ML(X1new,beta1,tj)
  }
   cte <- cteeps(delta.t*MX,epsilon,rho=rho)
   Y <- delta.t*MX+sqrt(cte)*epsilon; ncy=6
   Ynew <- delta.t*MXw+sqrt(cte)*epsnew}
  
   ,{ # -------HARD.L
  #beta1 <- outer(si,tj,function(u,v){7/2*(ifelse(v>= 0 & v< 0.35,
  #                  exp(u^3)* 1.6* ((v-0.1)/0.25)^3,
  #                     ifelse(v>= 0.35 & v< 0.8,
  #                            1.8*((v-0.575)/0.23)^5,
  #                         -cos(pi*u/2)*8*((v-1)/0.40)^2)))})
  beta1 <- outer(si,tj,function(u,v){7/2*(ifelse(v>= 0 & v< 0.333,
                    exp(u^3)* 1.6* ((v-0.1)/0.25)^3,
                       ifelse(v>= 0.333 & v< 0.75,
                              1.8*((v-0.6)/0.25)^5,
                           -cos(pi*u/2)*8*((v-1)/0.5)^2)))})
  
  #beta2 <- outer(si,tj,function(u,v){ifelse(v>= 0 & v< 0.5,
  #                                      -8*cos(2*pi*(2*u-1)*(2*v-1))*v^2,(2-3*v)^2)})
  beta2 <- outer(si,tj,function(u,v){(ifelse(v>= 0 & v< 0.5,
                              -5*cos(2*pi*(2*u-1)*(2*v-1))*(2*v)^2,(2-3*v)^2))})
  
  if (ivar==2) {
  MX <- ML(X1,beta1,tj)+ML(X2,beta2,tj)
  MXw <- ML(X1new,beta1,tj)+ML(X2new,beta2,tj)
   } else {
  MX <- ML(X1,beta1,tj)
  MXw <- ML(X1new,beta1,tj)
  }
  cte <- cteeps(delta.t*MX,epsilon,rho=rho)
  Y <- delta.t*MX+sqrt(cte)*epsilon; ncy <- 4
  Ynew <- delta.t*MXw+sqrt(cte)*epsnew}
  
  ,{ # -------EASY.NL
  beta2 <- outer(si, tj, function(u,v){-(u*v+1)*cos(2*pi*sqrt(u*v))})
  beta1 <- outer(si, tj, function(u,v){6*sqrt(v*u)*sin(4*pi*v)})
  if (ivar==2) {
   MX <- MNL(X1,TNL,beta1,tj)+MNL(X2,TNL1,beta2,tj)
   MXw <- MNL(X1new,TNL,beta1,tj)+MNL(X2new,TNL1,beta2,tj)
   } else {
  MX <- MNL(X1,TNL,beta1,tj)
  MXw <- MNL(X1new,TNL,beta1,tj)
  }
   cte <- cteeps(delta.t*MX,epsilon,rho=rho)
   Y <- delta.t*MX+sqrt(cte)*epsilon; ncy <- 3
   Ynew <- delta.t*MXw+sqrt(cte)*epsnew}
  
  ,{			# -------HARD.NL
  # beta1 <- outer(si,tj,function(u,v){(exp(-(-0.5-(2*u-1)-(2*u-1)^3))+ifelse(v>= 0 & v< 0.45,
                                                          # 0.8* exp(4*v+0.4)/2.25,
                                                                # ifelse(v>= 0.45 & v< 0.8,
                                                                      # 2.2*((2*v-1.25)/0.35)^5,
                                                                      # 3*((2*v-2)/0.28)^3+1)))})
  beta1 <- outer(si,tj,function(u,v){
    7/2*(ifelse(v>= 0 & v< 0.333,
    exp(u^3)* 1.6* ((v-0.1)/0.25)^3,
              ifelse(v>= 0.333 & v< 0.75,
                  1.8*((v-0.6)/0.25)^5,
                  -cos(pi*u/2)*8*((v-1)/0.5)^2)))
    }
    )
  #beta2 <- outer(si,tj,function(u,v){ifelse(v>= 0 & v< 0.5,
  #                                      -8*cos(2*pi*(2*u-1)*(2*v-1))*v^2,(2-3*v)^2)})
  beta2 <- outer(si,tj,
                 function(u,v){(ifelse(v>= 0 & v< 0.5,
                                -5*cos(2*pi*(2*u-1)*(2*v-1))*(2*v)^2,(2-3*v)^2))})
  if (ivar==2) {
    MX <- MNL(X1,TNL,beta1,tj) + MNL(X2,TNL1,beta2,tj)
    MXw <- MNL(X1new,TNL,beta1,tj) + MNL(X2new,TNL1,beta2,tj)
   } else {
    MX <- MNL(X1,TNL,beta1,tj)
    MXw <- MNL(X1new,TNL,beta1,tj)
  }									  
   cte <- cteeps(delta.t*MX,epsilon,rho=rho)
   Y <- delta.t*MX+sqrt(cte)*epsilon; ncy <- 4
   Ynew <- delta.t*MXw+sqrt(cte)*epsnew}
  ) # End switch
  ldatos <- ldata(X1=X1,X2=X2,Y=Y)
  ldatnew <- ldata(X1=X1new,X2=X2new,Y=Ynew)
  npc <- c(nc1,nc2,ncy)
  return(list(ldatos=ldatos,ldatnew=ldatnew,npc=npc))
}


ejecmodels=function(imodel,
                    N,
                    ivar=2,
                    Ntest=100,
                    nT=71,
                    nS=51,
                    rho=.9){
  
  dgen <- imodelfrfr(imodel=imodel,ivar=ivar,
                     N=N,Ntest=Ntest,nT=nT,
                     nS=nS,rho=rho)
  ###############
  #ldatos <- datagen$ldatos
  #ldatnew <- datagen$ldatnew
  nc1 <- dgen$npc[1]
  nc2 <- dgen$npc[2]
  ncy <- dgen$npc[3]
  # par(mfrow=c(1,3))
  # plot(ldatos$X1,main="X1(s)",xlab="s")
  # plot(ldatos$X2,main="X2(s)",xlab="s")
  # plot(ldatos$Y,main="Response Y(t)",xlab="s")
  
  
  YD <- dgen$ldatos$Y$data
  X1D=dgen$ldatos$X1$data
  X2D=dgen$ldatos$X2$data
  XX <- list(X1=dgen$ldatos$X1$data,X2=dgen$ldatos$X2$data)
  XXnew <- list(X1=dgen$ldatnew$X1$data,X2=dgen$ldatnew$X2$data)
  b.x <- list(X1=create.pc.basis(dgen$ldatos$X1,1:nc1),
              X2=create.pc.basis(dgen$ldatos$X2,1:nc2)) #,X3=create.pc.basis(X3,1:nc3))
  b.y <- create.pc.basis(dgen$ldatos$Y,1:ncy)
  
  # Model Estimation
  # fda.usc pkg
  pmetric <- list(df=data.frame(idx=1:nrow(dgen$ldatos$Y)),
                  Y=list(metric=metric.lp,lp=2),
                  X1=list(metric=metric.lp,lp=2),
                  X2=list(metric=metric.lp,lp=2))
                  #,X3=list(metric=metric.lp,lp=2))
  p.np <- list(X1=list(Ker=AKer.norm),
               X2=list(Ker=AKer.norm)) 
               #,X3=list(Ker=AKer.norm))
  
  n <- nrow(dgen$ldatos$X1)
  kpar <- floor(n/(2*mean(dgen$npc[1:2])))
  if (kpar>9) kpar <- -1
  
  reslin <- fregre.mlm.fr(Y~X1+X2,data=dgen$ldatos,basis.y=b.y,basis.x=b.x)
  ressam <- fregre.sam.fr(Y~s(X1,k=kpar)+s(X2,k=kpar),
                          data=dgen$ldatos,basis.y=b.y,basis.x=b.x)
  reskam <- fregre.kam.fr(Y~X1+X2,data=dgen$ldatos,par.metric=pmetric)
  
  prlin <- predict(reslin,dgen$ldatnew)
  prsam <- predict(ressam,dgen$ldatnew)
  prkam <- predict(reskam,dgen$ldatnew)
  
  
  # Refund pkg
  #respff.NL <- pffr(YD~sff(X1D,xind=ldatos$X1$argvals)
                  # +sff(X2D,xind=ldatos$X2$argvals),yind=ldatos$Y$argvals,
  #								bs.yindex=list(bs="ps",k=41,m=c(2,1)),bs.int=list(bs="ps",k=21,m=c(2,1)))
  respff.NL <- pffr(YD~sff(X1D,xind=dgen$ldatos$X1$argvals,
                    splinepars=list(bs="ps",k=c(11,5,5),m=c(2,2,2)))
  				          +sff(X2D,xind=dgen$ldatos$X2$argvals,
  				          splinepars=list(bs="ps",k=c(11,5,5),m=c(2,2,2))),
  				          yind=dgen$ldatos$Y$argvals,
  				          bs.yindex=list(bs="ps",k=5,m=c(2,1))
  				          )
  
  respffNL.fit <- fdata(matrix(respff.NL$fitted.values,
                               ncol=length(dgen$ldatos$Y$argvals),byrow=TRUE),
                        argvals=dgen$ldatos$Y$argvals)
  respffpc <- pffr(YD~ffpc(X1D,xind=dgen$ldatos$X1$argvals,npc.max=nc1)+
                    ffpc(X2D,xind=dgen$ldatos$X2$argvals,npc.max=nc2),
  		            	yind=dgen$ldatos$Y$argvals)
  respffpc.fit <- fdata(matrix(respffpc$fitted.values,
                               ncol=length(dgen$ldatos$Y$argvals),
                               byrow=TRUE),
                        argvals=dgen$ldatos$Y$argvals)
  
  prff.NL <- fdata(predict(respff.NL,
                           list(X1D=dgen$ldatnew$X1$data,X2D=dgen$ldatnew$X2$data)),
                   argvals=dgen$ldatos$Y$argvals)
  prff.pc <- fdata(predict(respffpc,
                           list(X1D=dgen$ldatnew$X1$data,X2D=dgen$ldatnew$X2$data)),
                         argvals=dgen$ldatos$Y$argvals)
  
  # FRegSigCom pkg
  t.x.list <- list(X1=dgen$ldatos$X1$argvals,X2=dgen$ldatos$X2$argvals)
  #----- ---- ---- ----     cv.sigcom
  resSC.Lin <- cv.sigcom(XX, YD, t.x.list, 
                         dgen$ldatos$Y$argvals,s.n.basis=15,t.n.basis=19)
  resSC.Lin.fit <- fdata(pred.sigcom(resSC.Lin, XX),dgen$ldatos$Y$argvals)
  prSCLin <- fdata(pred.sigcom(resSC.Lin,XXnew),dgen$ldatos$Y$argvals)
  
  #----- ---- ---- ----     cv.nonlinear
  resSC.NL <- cv.nonlinear(XX, YD, t.x.list, dgen$ldatos$Y$argvals,
                           x.n.basis=11,s.n.basis=15,t.n.basis=19)
  resSC.NL.fit <- fdata(pred.nonlinear(resSC.NL, XX),dgen$ldatos$Y$argvals)
  prSCNL <- fdata(pred.nonlinear(resSC.NL,XXnew),dgen$ldatos$Y$argvals)
  #
  lest <- list("Lin"=reslin$fitted.values, "SAM"=ressam$fitted.values,
               "KAM"=reskam$fitted.values, "PFF.PC"=respffpc.fit,
               "PFF.NL"=respffNL.fit,"SC.Lin"=resSC.Lin.fit,
               "SC.NL"=resSC.NL.fit)
  lpred <- list("Lin"=prlin, "SAM"=prsam,"KAM"=prkam,
                "PFF.PC"=prff.pc,"PFF.NL"=prff.NL,"SC.Lin"=prSCLin,
                "SC.NL"=prSCNL)
  resul1 <- best.estpred(dgen$ldatos$Y,dgen$ldatnew$Y,lest,lpred)
  #resul1$est.R2
  #resul1$pr.R2
  RR <- data.frame(method=c("Est","Pred"),imodel=rep(imodel,2),
                   ivar=rep(ivar,2),N=rep(N,2),
  	rho=rep(rho,2),rbind(resul1$est.R2,resul1$pr.R2))
  return(RR)
}

##############
# par(mfrow=c(2,2))
# for (i in 1:4){
# plot(ldatnew$Y[i],lwd=2,col="gray50",main=paste0("Pred. Curve:",i))
# lines(prlin[i],col="red")
# lines(prsam[i],col="red",lty=2)
# lines(prkam[i],col="red",lty=3)
# lines(prff.pc[i],col="blue")
# lines(prff.NL[i],col="blue",lty=2)
# lines(prSCLin[i],col="green")
# lines(prSCNL[i],col="green",lty=2)
# }

cesga <- FALSE
if (cesga) {
  library(doSNOW)
  cl <- getMPIcluster()
  registerDoSNOW(cl)
  } else{
  library(doParallel)
  ncores <- detectCores()
  cl <- makeCluster(ncores,type="PSOCK")
  registerDoParallel(cl)
}




# Numerical study
replic <- 100

time.ini <- Sys.time()
imod <- 1:4
ivar <- c(1,2)
N <- c(200)
rho <- c(.90,0.80,.65)
fname <- "res2_200.RData"
#rho <- .80
ejec <- expand.grid(imod=imod,ivar=ivar,N=N,rho=rho)

restotal <- NULL
#for (j in 1:nrow(ejec)){
for (ni in 1:replic){
  info <- Sys.info()[c("nodename", "machine")]
  print(paste(Sys.time(),"-Nodo:", info[1], "CPU:", info[2],"Replica:",ni))
  # resultado=foreach(ni=1:replic,.combine=rbind,.packages=c("fda.usc.devel","refund","FRegSigCom")) %dopar% {
  resultado=foreach(j=1:nrow(ejec),
                    .combine=rbind,
                    .packages=c("fda.usc.devel","refund","FRegSigCom")) %dopar% {
                        ejecmodels(imodel=ejec[j,"imod"],
                                   N=ejec[j,"N"],
                                   ivar=ejec[j,"ivar"],
                                   Ntest=100,
                                   nT=71,
                                   nS=51,
                                   rho=ejec[j,"rho"])
                    }
  j<-1
  res <- ejecmodels(imodel=ejec[j,"imod"],
                     N=ejec[j,"N"],
                     ivar=ejec[j,"ivar"],
                     Ntest=100,
                     nT=71,
                     nS=51,
                     rho=ejec[j,"rho"])
  #nam=paste0("Ej-",j,"-",ejec[j,"sigma"],".csv")
  #write.csv(resultado,file=nam)
  restotal=rbind(restotal,resultado)
}
stopCluster(cl)
rm(cl)
save(restotal,file=fname)
tabresul <- aggregate(restotal[,-c(1:5)],
                   by=list(rho=restotal$rho,method=restotal$method,
                           imodel=restotal$imodel,ivar=restotal$ivar,N=restotal$N)
                   ,mean)

time.fin <- Sys.time()
time.tot <- time.fin - time.ini 
print(time.tot)

registerDoSEQ()

# unregister_dopar <- function() {
#  env <- foreach:::.foreachGlobals
#  rm(list=ls(name=env), pos=env)
#}

#    cl <- makeCluster(2)
#    registerDoParallel(cl)
#    on.exit(stopCluster(cl))