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
 #ij<-c(1:140,191:210)
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
 N <- 1000
 B <- 50
 nh1 <- 4
 nh2 <- 4
 
 nn1 <- 5
 nn2 <- 5

 
 set.seed(1:4)
 res21 <- classif.bootstrap(formula,data=dat, weights = "equal",  
                            boot = "global",classif = "classif.glm",
                            par.boot = list(B = B, N = N, Nhull = nh1))
 set.seed(1:4)
 res22 <- classif.bootstrap(formula,data=dat, weights = "inverse", 
                            boot = "global",classif = "classif.glm",
                            par.boot = list(B = B, N = N, Nhull = nh1))
 
 set.seed(1:4)
 res23 <- classif.bootstrap(formula,data=dat, weights = "equal", 
                            boot = "local",classif = "classif.glm",
                            par.boot = list(B = B, N = N, Nhull = nh2, Nnbh = nn2))
 set.seed(1:4)
 res24 <- classif.bootstrap(formula,data=dat, weights = "inverse", 
                            boot = "local",classif = "classif.glm",
                            par.boot = list(B = B, N = N, Nhull = nh2,  Nnbh = nn2))
 
 pred21 <- predict(res21,newdat)
 pred22 <- predict(res22,newdat)
 pred23 <- predict(res23,newdat)
 pred24 <- predict(res24,newdat)
 
 cat2meas(ytest,pred0);cat2meas(ytest,pred1)
 cat2meas(ytest,pred21);cat2meas(ytest,pred22);cat2meas(ytest,pred23);cat2meas(ytest,pred24)

 }
 
##################################################### 
{
   # BBDD2
data(tecator)    
cutpoint <- 13
#cutpoint <- 15
cutpoint <- 18
cutpoint <- 10

cutpoint <- 11; ntest <- 100 #COMBINACION BUENA PARA X
cutpoint <- 13; ntest <- 100 #COMBINACION BUENA PARA X
cutpoint <- 11; ntest <- 100 #COMBINACION BUENA PARA X
cutpoint <- 11; ntest <- 100 #COMBINACION BUENA PARA X
ntrain <- 215 - ntest
ntrain

tecator$y$class <-factor(ifelse(tecator$y$Fat<cutpoint,0,1))
table(tecator$y$class )
# todos los conjuntos de datos tendrian esta estructura:
# df con la variable class
# x  con la variable funcional
x<-tecator[[1]]
# x1  derivada
x1<-fdata.deriv(tecator[[1]])
x1<-rproc2fdata(nrow(x),1:ntest)
#x2<-rproc2fdata(nrow(x),x$argvals,mu=func.mean(x)$data[1,])
x2<-fdata.deriv(tecator[[1]],2)
plot(x2)
plot(x1)
data<- ldata("df"=tecator$y,x=x,x1=x1,x2=x2)

set.seed(1:4)
itrain<-sample(1:215,replace=F)<=ntrain#165 # 12
ytrain <- data$df[,"class"]

dat<-data[itrain,row=TRUE]
newdat<- data[!itrain,row=TRUE]

#dat<-subset.ldata(data,itrain)
#newdat<-subset.ldata(data,!itrain)
########################################
formula0<-class ~ x
formula1<-class ~ x1
formula2<-class ~ x2

formula <- formula0
ytest <- newdat$df$class
#formula<- formula(class~s(x,k=3)+s(x1,k=3))
classif <- "classif.glm"

res0<-classif.glm(formula,data=dat,weights="equal")
pred0<-predict(res0,newdat)
wei <- fda.usc:::weights.inverse(dat$df$class)

res1<-classif.glm(formula,data=dat,weights=wei)#"inverse")
pred1<-predict(res1,newdat)

N <- 1000
B <- 50
nh1 <- 4
nh2 <- 4

nn1 <- 5
nn2 <- 5

set.seed(1:4)
res21 <- classif.bootstrap(formula,data=dat, weights = "equal",  
                         boot = "global",classif = "classif.glm",
                         par.boot = list(B = B, N = N, Nhull = nh1))
set.seed(1:4)
res22 <- classif.bootstrap(formula,data=dat, weights = wei,#"inverse", 
                           boot = "global",classif = "classif.glm",
                           par.boot = list(B = B, N = N, Nhull = nh1))

set.seed(1:4)
res23 <- classif.bootstrap(formula,data=dat, weights = "equal", 
                           boot = "local",classif = "classif.glm",
                           par.boot = list(B = B, N = N, Nhull = nh2, Nnbh = nn2))
set.seed(1:4)
res24 <- classif.bootstrap(formula,data=dat, weights = wei,#"inverse", 
                           boot = "local",classif = "classif.glm",
                           par.boot = list(B = B, N = N, Nhull = nh2,  Nnbh = nn2))

pred21 <- predict(res21,newdat)
pred22 <- predict(res22,newdat)
pred23 <- predict(res23,newdat)
pred24 <- predict(res24,newdat)

cat2meas(ytest,pred0);cat2meas(ytest,pred1)
cat2meas(ytest,pred21);cat2meas(ytest,pred22);cat2meas(ytest,pred23);cat2meas(ytest,pred24)

}

{
  # BBDD3
  #devtools::install_github("moviedo5/fda.tsc")
  #devtools::install_github("moviedo5/fdatasets",auth_token ="e691d0b0e8d73c44adbe3279c8c18ad523a60bef")
  #devtools::install_github("moviedo5/fdatasets",auth_token ="3960ebe2eb4ae21b9ae80520b1dff673db08f0eb")
  #library(fda.tsc)
  # devtools::install("C:/Users/moviedo/github/fdatasets")
  library(fdatasets)
  data(gun)    
  names(gun)
  dat<- list("df"=data.frame(class=gun$classlearn),x=gun$learn)
  newdat<- list("df"=data.frame(class=gun$classtest),x=gun$test)
  
  data(kalivas)
  names(kalivas$df)
  
  kalivas$x1 <- fdata.deriv(kalivas$x)
  
  class(kalivas)<-c("ldata",class(kalivas))
  n <- nrow(kalivas$df)
  
  ntrain <- 50 #floor(n*.5)
  set.seed(1:4)
  itrain<-sample(1:n,replace=F)<=ntrain#165 # 12
  dat<-kalivas[itrain,row=TRUE]
  newdat<- kalivas[!itrain,row=TRUE]
  
  set.seed(1:4)
  ########################################
  #formula<- formula(class~s(x))
  formula0<-class ~ x
  formula0<-class ~ x1
  formula<-formula0
  ytest<-newdat$df$class
  classif="classif.glm"
  #classif="classif.gsam";
  
  
  N <- 500
  B <- 25
  nh1 <- 4
  nh2 <- 10
  
  nn1 <- 10
  nn2 <- 5
  
  par.classif<-list(formula=formula,data=dat,weights="equal")
  set.seed(1:4)
  res0<-do.call(classif,par.classif)
  summary(res0)
  par.classif$weights="inverse"
  set.seed(1:4)
  res1<-do.call(classif,par.classif)
  pred0<-predict(res0,newdat)
  pred1<-predict(res1,newdat)
  
  #res.np <- classif.knn(dat$df$class,dat$x1,knn=1)
  #res.np <- classif.np(dat$df$class,dat$x1)
  #pred.np <- predict(res.np,newdat$x1)
  res.np <- classif.gkam(class~x1,data=dat)
  res.np
  pred.np <- predict(res.np,newdat)
  
  set.seed(1:4)
  res21 <- classif.bootstrap(formula,data=dat, 
                             boot = "global",classif = classif, weights = "equal",
                             par.boot = list(B = 50, N = 1000,  Nhull = nh1))
  set.seed(1:4)
  res22 <- classif.bootstrap(formula,data=dat, 
                             boot = "global",classif = classif, weights = "inverse",
                             par.boot = list(B = 50, N = 1000,  Nhull = nh1))
  set.seed(1:4)
  res23 <- classif.bootstrap(formula,data=dat, 
                             boot = "local",classif = classif, weights = "inverse",
                             par.boot = list(B = 50, N = 1000, Nhull = nh1, Nnbh = nn1))
  set.seed(1:4)
  res24 <- classif.bootstrap(formula,data=dat, 
                             boot = "local",classif = classif, weights = "inverse",
                             par.boot = list(B = 50, N = 1000,  Nhull = nh1,Nnbh = nn2))
  
  pred21<-predict(res21,newdat)
  pred22<-predict(res22,newdat)
  pred23<-predict(res23,newdat)
  pred24<-predict(res24,newdat)
  
  #args(classif.bootstrap)
  cat2meas(ytest,pred0);cat2meas(ytest,pred1);cat2meas(ytest,pred.np)
  cat2meas(ytest,pred21);cat2meas(ytest,pred22);cat2meas(ytest,pred23);cat2meas(ytest,pred24)
  round((1-cat2meas(ytest,pred0))*100,2)
  round((1-cat2meas(ytest,pred24))*100,2)
  }
