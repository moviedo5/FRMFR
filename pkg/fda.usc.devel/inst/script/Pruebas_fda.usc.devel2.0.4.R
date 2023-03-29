## Comparación con fregre.basis se usa el ejemplo de la función fRegress de fdas
library(fda.usc.devel)
# source("./R/fregre.basis.R")
# source("./R/fregre.basis.cv.R")
# source("./R/fregre.glm.R")
# source("./R/fregre.gsam.R")
# source("./R/fdata2model.R")
yobs <- log10(apply(CanadianWeather$dailyAv[,,"Precipitation.mm"], 2,sum))
smallbasis  <- create.fourier.basis(c(0, 365), 5)
smallbasisx  <- create.fourier.basis(c(0, 365), 15)
smallbasis  <- create.bspline.basis(c(0, 365), 5)
smallbasisx  <- create.bspline.basis(c(0, 365), 15)
tt <- day.5 
tt <- 1:365
tempfd <- smooth.basis(tt, CanadianWeather$dailyAv[,,"Temperature.C"], smallbasis)$fd
# precip.Temp1 <- fRegress(y ~ tempfd) # Bug en fda si la respuesta se llama "y"
precip.Temp1 <- fRegress(yobs ~ tempfd) 
yest0 <- precip.Temp1$yhatfdobj[,1]
names(yest0)<-names(yobs)
plot(yobs, yest0, type="p", pch="o")
sqrt(mean((yobs-yest0)^2))


temp <- fdata(t(CanadianWeather$dailyAv[,,"Temperature.C"]))
res0<-fda.usc:::fregre.basis(temp,yobs,basis.b=smallbasis)
#is.na.fdata <- fda.usc:::is.na.fdata
res1<-fregre.basis(temp,yobs,basis.b=smallbasis)#,basis.x=smallbasisx)
res01<-fregre.basis(temp,yobs,basis.b=smallbasis,lambda=0.00000000001) # arreglar penalización!!!!
res01<-fregre.basis(temp,yobs,basis.b=smallbasis,lambda=1e5) # arreglar penalización!!!!

res0
summary(res01)
summary(res01)
res01          
res01$coefficients
cbind(res1$coefficients,res01$coefficients)

#Minverse <- fda.usc:::Minverse
res11<-fregre.basis.cv(temp,yobs,basis.b=smallbasis$nbasis,
                       basis.x=smallbasis$nbasis,lambda = 0)
res11$fregre.basis
res12<-fregre.basis.cv(temp,yobs,basis.x=smallbasis$nbasis,lambda = 1e5)
res13<-fda.usc:::fregre.basis.cv(temp,yobs,basis.x=smallbasis$nbasis,lambda = 1e6)
#  plot the estimated regression function
plot(precip.Temp1$betaestlist[[2]])
lines(res0$beta.est,col="red")
lines(res1$beta.est,col=4)
lines(res01$beta.est,col=5)

lines(res11$fregre.basis$beta.est,col=3)
lines(res12$fregre.basis$beta.est,col=5,lty=2,lwd=3)
lines(res13$fregre.basis$beta.est,col=6,lwd=2)

sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res0$beta.est$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res1$beta.est$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res11$fregre.basis$beta.est$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res12$fregre.basis$beta.est$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res13$fregre.basis$beta.est$coefs)^2)
precip.Temp1$betaestlist[[1]]$fd$coefs - res11$fregre.basis$coefficients[1]
precip.Temp1$betaestlist[[1]]$fd$coefs - res12$fregre.basis$coefficients[1]

ldf <- list("df"=data.frame(yobs=yobs),"temp"=temp)
lbspb <- list("temp"=smallbasis)
lbspx <- list("temp"=smallbasisx)

#lbspb <- list("temp"=tempfd$basis);lbspx <- list("temp"=tempfd$basis)

names(ldf)
names(lbspb)
res2<-fda.usc:::fregre.lm(yobs~temp,data=ldf,basis.b=lbspb,basis.x=lbspb)
res3<-fregre.lm(yobs~temp,data=ldf,basis.b=lbspb)#,basis.x=lbspx)

plot(precip.Temp1$betaestlist[[2]])
lines(res2$beta.l[[1]],col="red")
lines(res3$beta.l[[1]],col=4)

sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res2$beta.l[[1]]$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res3$beta.l[[1]]$coefs)^2)
precip.Temp1$betaestlist[[1]]$fd$coefs - res2$coefficients[1]
precip.Temp1$betaestlist[[1]]$fd$coefs - res3$coefficients[1]
###########################################################
precip.Temp1 <- fRegress(yobs ~ tempfd) 
res2<-fda.usc:::fregre.glm(yobs~temp,data=ldf,basis.b=lbspb)#,basis.x=lbspx)
res3<-fregre.glm(yobs~temp,data=ldf,basis.b=lbspb)#,basis.x=lbspx)

precip.Temp1$betaestlist[[2]]$fd$coefs -res2$beta.l[[1]]$coefs
precip.Temp1$betaestlist[[2]]$fd$coefs - res3$beta.l[[1]]$coefs
# par(mfrow=c(1,1))
plot(precip.Temp1$betaestlist[[2]],lty=2)
lines(res2$beta.l[[1]],col="red")
lines(res3$beta.l[[1]],col=4)

sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res2$beta.l[[1]]$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res3$beta.l[[1]]$coefs)^2)
precip.Temp1$betaestlist[[1]]$fd$coefs - res2$coefficients[1]
precip.Temp1$betaestlist[[1]]$fd$coefs - res3$coefficients[1]
###########################################################
###########################################################

res2<-fda.usc:::fregre.gsam(yobs~temp,data=ldf,basis.b=lbspb)#,basis.x=lbspx)
res3<-fregre.gsam(yobs~temp,data=ldf,basis.b=lbspb)#,basis.x=lbspx)

sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res2$beta.l[[1]]$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res3$beta.l[[1]]$coefs)^2)
precip.Temp1$betaestlist[[1]]$fd$coefs - res2$coefficients[1]
precip.Temp1$betaestlist[[1]]$fd$coefs - res3$coefficients[1]


##########################################
res2<-fregre.lm(logprec~temp,data=ldf,basis.b=list("temp"=smallbasis))
res3<-fregre.glm(logprec~temp,data=ldf,basis.b=list("temp"=smallbasis))
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res$beta.est$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res1$beta.l[[1]]$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res2$beta.l[[1]]$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res3$beta.l[[1]]$coefs)^2)
cbind(precip.Temp1$betaestlist[[2]]$fd$coefs,res$beta.est$coefs,res1$beta.l[[1]]$coefs,res2$beta.l[[1]]$coefs)
mean(annualprec)
c(precip.Temp1$betaestlist$const$fd$coefs,precip.Temp1$betaestlist$tempfd$fd$coefs[,1])
res2$coefficients

















###########################################################################
# fregre.lm
library(fda.usc.devel)
data(tecator)
x <- tecator$absorp.fdata
y <- tecator$y$Fat
tt <- x$argvals
rtt <- x$rangeval
dataf <- as.data.frame(tecator$y)
nbasis.x <- 5
x.d1 <- fdata.deriv(x,nderiv=1,class.out='fdata',nbasis=nbasis.x)
x.d2 <- fdata.deriv(x.d1,nderiv=1,class.out='fdata',nbasis=nbasis.x)
ldat <- ldata("df"=dataf,"x"=x,"x.d2"=x.d2)

###########################################
# BSP
bsp5 <- create.bspline.basis(ldat$x$rangeval,nbasis=nbasis.x)
lbsp <- list("x"=bsp5,"x.d1"=bsp5,"x.d2"=bsp5)
f1 <- formula(Fat ~ x)
res0 <- fda.usc:::fregre.basis(x,y,basis.x=lbsp$x,basis.b=lbsp$x)
res1 <- fregre.basis(x,y,basis.x=lbsp$x,basis.b=lbsp$x)
res2 <-fda.usc:::fregre.lm(f1,ldat,basis.x=lbsp,basis.b=lbsp)
res3 <-fregre.lm(f1,ldat,basis.x=lbsp,basis.b=lbsp)
# res0;res1; res2; res3
par(mfrow=c(1,1))
plot(res0$beta.est,col=1,lwd=3)
lines(res1$beta.l$x,col=2,lwd=2,lty=2)
lines(res2$beta.l$x,col=3,lwd=2,lty=2)
lines(res3$beta.l$x,col=4,lwd=2,lty=2)


###########################################
# BSP penalizado

#' lambda <- list("x"=.5)
#' P <- list("x"=c(0,0,1))
#' res1 <- fregre.basis(x,y,basis.x=lbsp$x,basis.b=lbsp$x)
#' res2 <- fregre.basis(x,y,basis.x=lbsp$x,basis.b=lbsp$x,lambda=10,Lfdobj = vec2Lfd(c(0,1), rtt))
#' res3 <-fda.usc:::fregre.lm(f1,ldat,basis.x=lbsp,basis.b=lbsp,lambda=lambda)
#' res4 <-fregre.lm(f1,ldat,basis.x=lbsp,basis.b=lbsp,P=P,lambda=lambda)
#' par(mfrow=c(1,1))
#' plot(res1$beta.est,col=1,lwd=3,ylim=c(-25,25))
#' lines(res2$beta.est,col=2,lwd=3,lty=2)
#' lines(res3$beta.l$x,col=3,lwd=3,lty=3)
#' lines(res4$beta.l$x,col=4,lwd=2,lty=4)
#' ###########################################
#' # BSP penalizado 2 variables
#' f2 <- formula(Fat ~ x + x.d2)
#' res1 <-fda.usc:::fregre.lm(f2,ldat,basis.x=lbsp,basis.b=lbsp)
#' res2 <-fregre.lm(f2,ldat,basis.x=lbsp,basis.b=lbsp)
#' 
#' lambda <- list("x"=1e-5)
#' P <- list("x"=c(0,0,1))
#' res3 <-fregre.lm(f2,ldat,basis.x=lbsp,basis.b=lbsp,P=P,lambda=lambda)
#' lambda <- list("x.d2"=1e-10)
#' P <- list("x.d2"=c(0,0,1))
#' res4 <-fregre.lm(f2,ldat,basis.x=lbsp,basis.b=lbsp,P=P,lambda=lambda)
#' 
#' lambda <- list("x"=1e-5,"x.d2"=1e-10)
#' P <- list("x"=c(0,0,1),"x.d2"=c(0,0,1))
#' res5 <-fregre.lm(f2,ldat,basis.x=lbsp,basis.b=lbsp,P=P,lambda=lambda)
#' par(mfrow=c(1,2))
#' plot(res1$beta.l$x,col=1,lwd=3)
#' lines(res2$beta.l$x,col=2,lwd=2,lty=2)
#' lines(res3$beta.l$x,col=3,lwd=2,lty=3)
#' lines(res4$beta.l$x,col=4,lwd=2,lty=3)
#' lines(res5$beta.l$x,col=5,lwd=2,lty=3)
#' ##############################################################
#' # PC
#' lpc <- list("x"=create.pc.basis(ldat$x,1:3),
#'             "x.d1"=create.pc.basis(ldat$x.d1,1:3),
#'                         "x.d2"=create.pc.basis(ldat$x.d2,1:3))
#' res1 <-fregre.pc(x,y,basis.x=lpc$x)
#' res2 <-fda.usc:::fregre.lm(f1,ldat,basis.x=lpc)
#' res3 <-fregre.lm(f1,ldat,basis.x=lpc)
#' res1; res2; res3
#' par(mfrow=c(1,1)) 
#' plot(res1$beta.est,col=1,lwd=3)
#' lines(res2$beta.l$x,col=4,lwd=2,lty=2)
#' lines(res3$beta.l$x,col=3,lwd=2,lty=4)
#' ###########################################
#' # PC penalized
#' lambda <- list("x"=10)
#' P <- list("x"=c(1))
#' res1 <-fregre.pc(x,y,basis.x=lpc$x)
#' res2 <-fregre.pc(x,y,basis.x=lpc$x,lambda=lambda$x,P=P$x)
#' res3 <-fda.usc:::fregre.lm(f1,ldat,basis.x=lpc,lambda=lambda)
#' res4 <-fregre.lm(f1,ldat,basis.x=lpc,P=P,lambda=lambda)
#' par(mfrow=c(1,1))
#' plot(res1$beta.est,col=1,lwd=3,ylim=c(-4,4))
#' lines(res2$beta.est,col=2,lwd=3,lty=2)
#' lines(res3$beta.l$x,col=3,lwd=3,lty=3)
#' lines(res4$beta.l$x,col=4,lwd=2,lty=4)

#' ###########################################
#' # PC penalized 2 covariates
#' f2 <- formula(Fat ~ x + x.d2)
#' res1 <-fda.usc:::fregre.lm(f2,ldat,basis.x=lpc)
#' res2 <-fregre.lm(f2,ldat,basis.x=lpc)
#' 
#' lambda <- list("x"=.1)
#' P <- list("x"=c(1))
#' res3 <-fregre.lm(f2,ldat,basis.x=lpc,P=P,lambda=lambda)
#' 
#' lambda <- list("x.d2"=5)
#' P <- list("x.d2"=c(0,0,1))
#' res4 <-fregre.lm(f2,ldat,basis.x=lpc,P=P,lambda=lambda)
#' 
#' lambda <- list("x"=1,"x.d2"=5)
#' P <- list("x"=c(1),"x.d2"=c(0,0,1))
#' res5 <-fregre.lm(f2,ldat,basis.x=lpc,P=P,lambda=lambda)
#' par(mfrow=c(1,2))
#' plot(res1$beta.l$x,col=1,lwd=3)
#' lines(res2$beta.l$x,col=2,lwd=2,lty=2)
#' lines(res3$beta.l$x,col=3,lwd=2,lty=3)
#' lines(res4$beta.l$x,col=4,lwd=2,lty=3)
#' lines(res5$beta.l$x,col=5,lwd=2,lty=3)
#' 
#' 
#' 
#' 
#' 
#' fregre.lm.fr
#' ##################################################  
# # YBSP - XBSP
# res11 <-  fda.usc.devel:::fregre.lm.fr(as.formula(y~x+x2),mdat)
# res22 <-  fregre.lm.fr(as.formula(y~x2),mdat)
# res11
# res22
# plot(res11$beta.l$x)
# plot(res22$beta.l$x2)
# 
# plot(res22$beta2d$x2,type="persp")
# 
# 
# 
# 
# bsp7 <- create.bspline.basis(mdat$y$rangeval,7)
# bsp5 <- create.bspline.basis(mdat$y$rangeval,5)
# 
# lbsp <- list("x"=bsp7,"x2"=bsp5)
# res22 <-  fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=bsp7,basis.x=lbsp)
# res22
# 
# plot(res22$beta2d$x)
# plot(res22$beta2d$x,type="persp")
# plot(res22$beta2d$x2,type="persp",phi=45,theta=60)
# 
# 
# ##################################################  
# # YPC - XPC
# pcy <- create.pc.basis(mdat$y,1:4)
# pcx <- create.pc.basis(mdat$x,1:5)
# pcx2 <- create.pc.basis(mdat$x2,1:3)
# lpc <- list("x"=pcx,"x2"=pcx2)
# res11 <-  fda.usc.devel:::fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=pcy,basis.x=lpc)
# res22 <-  fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=pcy,basis.x=lpc)
# 
# par(mfrow=c(1,2))
# plot(res11$beta.l$x)
# plot(res11$beta.l$x2)
# 
# # Falta hacer el objeto beta2d!!!!
# res11$beta.l$x2$data
# 
# res22 <-  fregre.lm.fr(as.formula(y~x),mdat,basis.y=pcy,basis.x=lpc)
# 
# 
# 
##################################################  
# # YPC - XBSP
# pcy <- create.pc.basis(mdat$y,1:3)
# res11 <-  fda.usc.devel:::fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=pcy)
# res11
# plot(res11$beta.l$x)
# 
# res22 <-  fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=pcy)
# bsp7 <- create.bspline.basis(mdat$y$rangeval,7)
# bsp5 <- create.bspline.basis(mdat$y$rangeval,5)
# 
# lbsp <- list("x"=bsp7,"x2"=bsp5)
# res22 <-  fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=bsp7,basis.x=lbsp)
# res22
# 
# plot(res22$beta2d$x)
# plot(res22$beta2d$x,type="persp")
# plot(res22$beta2d$x2,type="persp",phi=45,theta=60)



# ##################
# predict.fregre.lm.fr# ejemplo MLM
# x0 <- tecator$absorp.fdata$data
# y0 <- as.matrix(tecator$y)
# df0 <- list("y0"=y0,"x0"=x0)
# #df0 <- data.frame(y0=y0,x0=x0)
# names(df0)
# 
# res <- lm(as.formula(y0~x0),data=df0)
# res$coefficients
# res$residuals
# predict(res,df0)
# 
# #################################
# x0 <- tecator$absorp.fdata
# mdf <- mfdata(y0=fdata.deriv(x0),x0=x0)
# names(mdf)
# names(mdf)<-c("y0","x0")
# names(mdf)
# rr <- fregre.lm.fr(as.formula(y0~x0),data=mdf)
# names(rr)
# rr$residuals
#res0$basis.y.class
# pred1 <- predict.fregre.lm.fr(res0,mdat)
# class(pred1)
# dim(pred1)

# plot(pred1,col=1)
# lines(res0$fitted.values,col=2)

# plot(pred1-res0$fitted.values)















# fregre.gsam.vs
########################################################
# dist.list <- fda.usc.devel:::dist.list
# dcor.y <- fda.usc.devel:::dcor.y
# sp <- fda.usc.devel:::sp
# pvalue.anova <- fda.usc.devel:::pvalue.anova
########################################################
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat", include = covar)
# summary(res.gam2)
# names(ldat)
# bsp<-create.bspline.basis(c(0,1),7)
# basis<-list("x"=bsp,"x1"=bsp,"x2"=bsp)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis)
# summary(res.gam2)
# 
# basis2<-c("pc",12,TRUE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2,exclude="x")
# summary(res.gam2)
# 
# basis2<-c("pc",12,FALSE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# basis2<-c("pc",11,FALSE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# basis2<-c("pc",7,TRUE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat")
# 
# basis2<-c("bspline",13,FALSE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# 
# basis2<-c("bspline",7,TRUE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# basis2<-c("pls",5,FALSE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# basis2<-c("pls",4,TRUE)
# res.gam1 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# basis2<-c("pls",8,FALSE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# anova(res.gam1,res.gam2)

# dist.list <- fda.usc.devel:::dist.list
# dcor.y <-  fda.usc.devel:::dcor.y

########################################################
# dist.list <- fda.usc.devel:::dist.list
# dcor.y <- fda.usc.devel:::dcor.y
# sp <- fda.usc.devel:::sp
# pvalue.anova <- fda.usc.devel:::pvalue.anova

Alternativas a CQ con datos no normales e correlados--
  Liu, profundidad! HVAC panamá
  Control de la eficiencia energética en una tienda 
  Consumo, ilumanacion, climatización.--detecta sábados y 
  domingos con un horario diferente, también paradas y 
  averías del sistema.
  Curvas funcionales - fase 1 limpiar la muestra,
  2 calcular la profundidad
  3 bootstrap 2008 encontrar límite de control para detectar anomalias.
  Fase II, monitorizar curva, calcular profundidad respecto calibrado
  y obtengo gráfico de rango para los rangos.
  e intentar ver las causas asignables
  1) deteccion de anomialia aprendizaje máquina y boostrap LOCI
  2) métodos no parametrocos mandel interlaboratorio
  3) graficos de control para observaciones incorreladas.

datos temperatura para toda europa e***45 en pixeles usando
paquete eco***
en lugar de usar promedio usar evolución?
EUROMOMO MOMO

  