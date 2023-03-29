library(fda.usc)
data(tecator)
ab=tecator$absorp.fdata
ab2=fdata.deriv(ab,2)
# ab <-ab2
ifat=factor(ifelse(tecator$y$Fat<8,1,ifelse(tecator$y$Fat>16,3,2)),label=c("F0","F8","F16"))

J <- 10
Lpred<- matrix(NA,J,6)
for (j in 1:J){
  set.seed(j)
  ll=sample(1:215,165)
  ltec=ldata(df=data.frame(ifat=ifat[ll]),ab=ab[ll],ab2=ab2[ll])
  ltecnew=ldata(df=data.frame(ifat=ifat[-ll]),ab=ab[-ll],ab2=ab2[-ll])

  
  # btec=list(ab=create.pc.basis(ltec$ab,1:4),ab2=create.pc.basis(ltec$ab2,1:4))
  bsp <- create.fdata.basis(ltec$ab,1:7)
  btec=list(ab=bsp,ab2=bsp)
  
  wei=fda.usc:::weights4class(ltec$df$ifat,type="inverse")
  wei=wei/sum(wei)
  ww=tapply(fda.usc:::weights4class(ltec$df$ifat,type="inverse"),ltec$df$ifat,mean)
  ww=ww/sum(ww)
  #print(tt<-table(ifat))
  #wei=1/tt[ifat]
  #wei=wei/sum(wei)
  #ww=(1/tt)/sum(1/tt)
  #cat("Frequencies:",tt,"\n")
  cat("Weights per group:",round(ww,3),"\n")
  
  ctrl2<- list(trace = FALSE, draw = FALSE,gray.scale=FALSE,fine=51)
  colores=c("green","orange","red")
  
  par(mfrow=c(1,2))
  plot(ab,lty=1,col=colores[ifat])
  legend("topleft",legend=c("F[0,8]","F[8,16]","F[16-]"),col=colores,lwd=2)
  plot(ab2,lty=1,col=colores[ifat],main="2nd derivative")
  legend("topleft",legend=c("F[0,8]","F[8,16]","F[16-]"),col=colores,lwd=2)
  
  type1 <- "1vsall"
  type2 <- "majority"
  type2 <- "1vsall"
  
  tec3=classif.DD(ltec$df$ifat,ltec$ab,depth="mode",classif="glm",control=ctrl2)
  tec.glm1=classif.glm(ifat~ab,ltec,basis.x=btec,type=type1)
  tec.glm2=classif.glm(ifat~ab,ltec,basis.x=btec,type=type2)
  tec.glmw=classif.glm(ifat~ab,ltec,basis.x=btec,type=type2,weights=wei)
  
  B = 10
  tec.ada=classif.adaboost(ifat~ab,ltec,classif="classif.glm",# weights = wei, 
                           B=B,par.classif=list(basis.x=btec,type=type2))
  B = 10
  N <- length(ltec$df$ifat) # *10
  N <- 1000
  boot = "global"
  #boot ="local"
  
  tec.boot=classif.bootstrap(ifat~ab,ltec,boot= boot
             , weights = wei
             , classif="classif.glm"
             , par.boot=list(B=B,N=N)
             , par.classif=list(basis.x=btec,type=type2))
  #summary(tec.boot)
  #tec.boot$prob.classification=fda.usc:::classifKgroups(ltec$df$ifat,tec.boot$prob.group,levels(ltec$df$ifat))$prob1
  
  M=rbind(tec3$prob.classification,tec.glm1$prob.classification,tec.glm2$prob.classification,tec.glmw$prob.classification,tec.ada$prob.classification,tec.boot$prob.classification)%*%ww
  #M=rbind(tec3$prob.classification,tec.glm1$prob.classification,tec.glm2$prob.classification,tec.glmw$prob.classification,tec.ada$prob.classification)%*%ww
  M=cbind(rbind(tec3$prob.classification,tec.glm1$prob.classification,tec.glm2$prob.classification,tec.glmw$prob.classification,tec.ada$prob.classification,tec.boot$prob.classification),M)
  #M=cbind(rbind(tec3$prob.classification,tec.glm1$prob.classification,tec.glm2$prob.classification,tec.glmw$prob.classification,tec.ada$prob.classification),M)
  rownames(M)=c("DD","GLM/OVA","GLM/OVO","GLM/OVO-W","Adaboost","Bootstrap")
  #rownames(M)=c("DD","GLM/OVA","GLM/OVO","GLM/OVO-W","Adaboost")
  colnames(M)=c("F0","F8","F16","Weig. Prob. Class.")
  # print(round(M,4))
  
  Mpred=matrix(0,nrow=nrow(M),ncol=ncol(M))
  pr.DD=predict(tec3,ltecnew$ab)
  #pr.DD=predict(tec.DD,ltecnew$ab)
  tt.DD=table(ltecnew$df$ifat,pr.DD)
  Mpred[1,1:3]=diag(tt.DD)/apply(tt.DD,1,sum)
  Mpred[1,4]=sum(Mpred[1,1:3]*ww)
  colnames(Mpred)=colnames(M)
  rownames(Mpred)=rownames(M)
  pr.glm1=predict(tec.glm1,ltecnew)
  tt.glm1=table(ltecnew$df$ifat,pr.glm1)
  Mpred[2,1:3]=diag(tt.glm1)/apply(tt.glm1,1,sum)
  Mpred[2,4]=sum(Mpred[2,1:3]*ww)
  pr.glm2=predict(tec.glm2,ltecnew)
  tt.glm2=table(ltecnew$df$ifat,pr.glm2)
  Mpred[3,1:3]=diag(tt.glm2)/apply(tt.glm2,1,sum)
  Mpred[3,4]=sum(Mpred[3,1:3]*ww)
  pr.glmw=predict(tec.glmw,ltecnew)
  tt.glmw=table(ltecnew$df$ifat,pr.glmw)
  Mpred[4,1:3]=diag(tt.glmw)/apply(tt.glmw,1,sum)
  Mpred[4,4]=sum(Mpred[4,1:3]*ww)
  pr.ada=predict(tec.ada,ltecnew)
  tt.ada=table(ltecnew$df$ifat,pr.ada)
  Mpred[5,1:3]=diag(tt.ada)/apply(tt.ada,1,sum)
  Mpred[5,4]=sum(Mpred[5,1:3]*ww)
  pr.boot=predict(tec.boot,ltecnew)
  tt.boot=table(ltecnew$df$ifat,pr.boot)
  Mpred[6,1:3]=diag(tt.boot)/apply(tt.boot,1,sum)
  Mpred[6,4]=sum(Mpred[6,1:3]*ww)
  Mpred
  Lpred[j,]<-Mpred[,4]
  cat(j)
}
colnames(Lpred)<-rownames(Mpred)
c(J,B,boot,type2)
colMeans(Lpred)


########################################
# Resultados:  ab
> c(J,B,boot,type2)
[1] "10"     "10"     "local"  "1vsall"
> colMeans(Lpred)
DD   GLM/OVA   GLM/OVO GLM/OVO-W  Adaboost Bootstrap 
0.4289931 0.7228660 0.7228660 0.7269221 0.7018647 0.7168631
[1] "10"       "10"       "local"    "majority"
> colMeans(Lpred)
DD   GLM/OVA   GLM/OVO GLM/OVO-W  Adaboost Bootstrap 
0.4289931 0.7228660 0.7322840 0.7298964 0.5440975 0.7269309 
[1] "10"       "10"       "global"   "majority"
> colMeans(Lpred)
DD   GLM/OVA   GLM/OVO GLM/OVO-W  Adaboost Bootstrap 
0.4289931 0.7228660 0.7322840 0.7298964 0.5440975 0.7270647 
[1] "10"     "10"     "global" "1vsall"
> colMeans(Lpred)
DD   GLM/OVA   GLM/OVO GLM/OVO-W  Adaboost Bootstrap 
0.4289931 0.7228660 0.7228660 0.7269221 0.7018647 0.7277917 



# Resultados: ab2
1] "50"       "20"       "local"    "majority"
> colMeans(Lpred)
DD   GLM/OVA   GLM/OVO GLM/OVO-W  Adaboost Bootstrap 
0.8395000 0.8835373 0.8899585 0.8958948 0.7577300 0.8917994 

> c(J,B,boot,type2)
[1] "50"       "20"       "local"    "majority"
> colMeans(Lpred)
DD   GLM/OVA   GLM/OVO GLM/OVO-W  Adaboost Bootstrap 
0.8395000 0.8835373 0.8899585 0.8958948 0.7577300 0.8917994 

#[1] "50"       "20"       "global"   "majority"
DD   GLM/OVA   GLM/OVO GLM/OVO-W  Adaboost Bootstrap 
0.8395000 0.8835373 0.8899585 0.8958948 0.7577300 0.8713656 

# local  B3     DD   GLM/OVA   GLM/OVO GLM/OVO-W  Adaboost Bootstrap 
# 0.8429040 0.8796117 0.9037658 0.9150508 0.7054710 0.9203397 
# global DD   GLM/OVA   GLM/OVO GLM/OVO-W  Adaboost Bootstrap 
# 0.8429040 0.8796117 0.9037658 0.9150508 0.7054710 0.8993394 
# data(phoneme)
# mlearn <- phoneme[["learn"]]
# glearn <- phoneme[["classlearn"]]
# dataf <- data.frame(glearn)
# dat = list("df" = dataf, "x" = mlearn)
# dat.boot <- convex.bootstrap(dat,"glearn")
# dim(dat$df)
# dim(dat.boot$df) 
# par(mfrow=c(1,2))
# plot(dat$x,col=dat$df$glearn)
# plot(dat.boot$x,col=dat.boot$df$glearn)
# 



# cat2alpha <- fda.usc:::cat2alpha
# predict.classif.bootstrap <- fda.usc:::predict.classif.bootstrap
# subset.fdata <- fda.usc:::subset.fdata
