##################################################### 
# cat2meas<-fda.usc:::cat2meas
# weights.inverse<-fda.usc:::weights.inverse
# #tab2waccuracy<-fda.usc:::tab2waccuracy
# tab2accuracy<-fda.usc:::tab2accuracy
# prob2classif<-fda.usc:::prob2classif
# weights4class<-fda.usc:::weights4class
# convex2group<-fda.usc:::convex2group
##################################################### 
library(fda.usc)
subset.ldata <- fda.usc:::subset.ldata
{
 # BBDD1
 data(phoneme)
 mlearn<-phoneme[["learn"]]
 glearn<-phoneme[["classlearn"]]
 mtest<-phoneme[["test"]]
 gtest<-phoneme[["classtest"]]
 dataf<-data.frame(glearn)
 ij<-c(41:70,101:200,241:250)
# ij<-c(1:140,191:210)
 dat=list("df"=dataf[ij,,drop=F],"x"=mlearn[ij,])
 newdat=list("df"=dataf[,,drop=F],"x"=mtest[,])
 #newdat<-dat
 table(dat$df$glearn)
 table(newdat$df$glearn)
 classif<-"classif.glm"

 ytest<-newdat$df$glearn
 
 formula <- glearn~x

 res0<-classif.glm(formula,data=dat)
 pred0<-predict(res0,newdat)

 res1<-classif.glm(formula,data=dat,weights="inverse")
 pred1<-predict(res1,newdat)
 cat2meas(ytest,pred0);cat2meas(ytest,pred1)
 
 
 B<-50
 res2<-classif.adaboost(formula,data=dat,classif=classif,B=B,
                        coeflearn ="Breiman",   weights = "equal")
 summary(res2)
 pred2<-predict.classif.adaboost(res2,newdat)
 res3<-classif.adaboost(formula,data=dat,classif=classif,B=B,
                        coeflearn ="Breiman",   weights = "inverse")
 summary(res3)
 pred3<-predict.classif.adaboost(res3,newdat)
 cat2meas(ytest,pred0);cat2meas(ytest,pred1)
 cat2meas(ytest,pred2);cat2meas(ytest,pred3)
 
 res2$alpha.boost;res3$alpha.boost
 plot(res2$alpha.boost,col=4,type="l")
 lines(res3$alpha.boost,col=2)
 }
 
##################################################### 
{
   # BBDD2
   data(tecator)    
cutpoint <- 13
#cutpoint <- 15
#cutpoint <- 18
cutpoint <- 11
tecator$y$class <-factor(ifelse(tecator$y$Fat<cutpoint,0,1))
table(tecator$y$class )
# todos los conjuntos de datos tendrian esta estructura:
# df con la variable class
# x  con la variable funcional
x<-tecator[[1]]
# x1  derivada
x1<-fdata.deriv(tecator[[1]])
x1<-rproc2fdata(nrow(x),1:50)
#x2<-rproc2fdata(nrow(x),x$argvals,mu=func.mean(x)$data[1,])
x2<-fdata.deriv(tecator[[1]],2)
plot(x2)
plot(x1)
data<- list("df"=tecator$y,x=x,x1=x1,x2=x2)
subset.ldata <- fda.usc:::subset.ldata
set.seed(1:4)
itrain<-sample(1:215,replace=F)<=129#165 # 129
ytrain <- data$df[,"class"]
dat<-subset.ldata(data,itrain)
newdat<-subset.ldata(data,!itrain)
########################################
formula0<-class ~ x
formula1<-class ~ x+x1
formula2<-class ~ x + x2


formula<-formula1
ytest<-newdat$df$class
#formula<- formula(class~s(x,k=3)+s(x1,k=3))
classif="classif.glm";

res0<-classif.glm(formula,data=dat,weights="equal")
pred0<-predict(res0,newdat)

res1<-classif.glm(formula,data=dat,weights="inverse")
pred1<-predict(res1,newdat)

B<-50
res2<-classif.adaboost(formula,data=dat,classif=classif,B=B,
                       coeflearn ="Breiman",   weights = "equal")
summary(res2)
pred2<-predict.classif.adaboost(res2,newdat)
res3<-classif.adaboost(formula,data=dat,classif=classif,B=B,
                       coeflearn ="Breiman",   weights = "inverse")
summary(res3)
pred3<-predict.classif.adaboost(res3,newdat)
cat2meas(ytest,pred0);cat2meas(ytest,pred1)
cat2meas(ytest,pred2);cat2meas(ytest,pred3)

res2$alpha.boost;res3$alpha.boost
plot(res2$alpha.boost,col=4,type="l")
lines(res3$alpha.boost,col=2)
}

{
  # BBDD3
  library(fdatasets)
  data(gun)    
  names(gun)
  dat<- list("df"=data.frame(class=gun$classlearn),x=gun$learn)
  newdat<- list("df"=data.frame(class=gun$classtest),x=gun$test)
  
  set.seed(1:4)
  ########################################
  formula0<-class ~ x
  formula<-formula0
  ytest<-newdat$df$class
  #formula<- formula(class~s(x))
  #classif="classif.nnet" # el ada con equal enpeora, con inverse mejora
  classif="classif.glm"
  par.classif<-list(formula=formula,data=dat,weights="equal")
  set.seed(1:4)
  res0<-do.call(classif,par.classif)
  summary(res0)
  par.classif$weights="inverse"
  set.seed(1:4)
  res1<-do.call(classif,par.classif)
  pred0<-predict(res0,newdat)
  pred1<-predict(res1,newdat)
  
  set.seed(1:4)
  B<-50
  res2<-classif.adaboost(formula,data=dat,classif=classif,B=B, weights = "equal")
  summary(res2)
  pred2<-predict.classif.adaboost(res2,newdat)
  set.seed(1:4)
  res3<-classif.adaboost(formula,data=dat,classif=classif,B=B,weights = "inverse")
  summary(res3)
  pred3<-predict.classif.adaboost(res3,newdat)
  pred3<-predict.classif.adaboost(res3,newdat)
  length(res3$list.fit)
  ts.plot(res3$alpha.boost)
  lines(res2$alpha.boost,col=2)
  
  cat2meas(ytest,pred0);cat2meas(ytest,pred1)
  cat2meas(ytest,pred2);cat2meas(ytest,pred3)
  table(ytest,pred2);table(ytest,pred3)

}