
data(aemet)
aemet$df$logprecsum<-rowSums(aemet$logprec$data)
aemet$df$logprecdia1<-aemet$logprec$data[,1]
aemet$df$logprecdia1<-aemet$logprec
ff0<-formula(logprecdia1~temp+altitude)
res0<-fregre.lm(ff0,data=aemet)

ff1<-formula(logprec~temp)
res1<-fregre.lm.fr(ff1,data=aemet)

names(aemet)
names(aemet$df)
# pred0<-predict(res0,aemet)
# pred1<-predict.fregre.lm.fr(res1,aemet)

ff2<-logprec~temp+altitude
res2<-fregre.lm.fr(ff2,data=aemet)
#pred2<-predict.fregre.lm.fr(res2,aemet2)
res3<-fregre.lm.fr(ff2,data=aemet,lambda=10000)
#pred3<-predict.fregre.lm.func(res3,aemet2)


ff2<-logprec~temp+altitude
# 
# lpc<-list("temp"=create.pc.basis(aemet$temp))
# ff1<-logprecdia1~temp+altitude
# res0<-fregre.lm(ff1,data=aemet,basis.x=lpc)
# res1<-fregre.lm.func(ff1,data=aemet,basis.x=lpc)
# pred0<-predict.fregre.lm(res0,aemet,basis.x=lpc)
# pred1<-predict.fregre.lm.func(res1,aemet,basis.x=lpc)

ff2<-logprec~temp+altitude
res2<-fregre.lm.fr(ff2,data=aemet,basis.x=lpc)
dim(res2$coefficients)

ff2<-lp~temp+altitude
aemet$lp<-fdata2fd(aemet$logprec)
(ff2)
names(aemet)
res2<-fregre.lm.fr(ff2,data=aemet)
dim(res2$coefficients)
res0<-lm(aemet$logprec$data[,1:3]~aemet$temp$data[,1:8])
res0
#pred2<-predict.fregre.lm.func(res2,aemet)
#res3<-fregre.lm.func(ff2,data=aemet,lambda=10000,basis.x=lpc)
#pred3<-predict.fregre.lm.func(res3,aemet2)        #no va Respuesta funcional con penalizacion

lpc<-list("temp"=create.pc.basis(aemet$temp),"wind.speed"=create.pc.basis(aemet$wind.speed,1:3))
res3<-fregre.lm.fr(logprec~temp+altitude+wind.speed,data=aemet,basis.x=lpc,lambda=1000)
res4<-fregre.lm.frc(logprec~temp+altitude+wind.speed,data=aemet,basis.x=lpc,lambda=c(100,1000))


lpc<-list("temp"=create.pc.basis(aemet$temp),"wind.speed"=create.pc.basis(aemet$wind.speed,1:3))
res3<-fregre.lm.fr(logprec~temp+altitude+wind.speed,data=aemet,basis.x=lpc,lambda=1000)
#ver penalizaciones (mlm no se puede, repetir npy veces!)
# penalizar las vnf
