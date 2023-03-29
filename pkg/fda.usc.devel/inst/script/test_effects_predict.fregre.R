# hacer testhat files
data(tecator)
ind<-1:129
x=tecator$absorp.fdata
x.d2<-fdata.deriv(x,nderiv=2)
tt<-x[["argvals"]]
dataf=as.data.frame(tecator$y)
nbasis.x=11;nbasis.b=7
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
basis.x=list("x.d2"=basis1)
basis.b=list("x.d2"=basis2)
ldat = ldata("df"=dataf[ind,],"x.d2"=x.d2[ind])
############## LM
res.lm=fregre.lm(Fat~Water+Protein+x.d2,data=ldata)
fake <- predict.fregre.lm(res.lm,type="terms")
head(fake)
fake <- predict.fregre.lm(res.lm,type="effects")
head(fake)
fake <- predict.fregre.lm(res.lm,ldata[1,row=T],type="terms")
head(fake)
fake <- predict.fregre.lm(res.lm,ldata[1,row=T],type="effects")
head(fake)
############## GLM
res.glm=fregre.glm(Fat~Water+x.d2,data=ldata,family=Gamma(link = "inverse"))
fake <- predict.fregre.glm(res.glm,type="terms")
head(fake)
fake <- predict.fregre.glm(res.glm,type="effects")
head(fake)
fake <- predict.fregre.glm(res.glm,ldata[1,row=T],type="terms")
head(fake)
fake <- predict.fregre.glm(res.glm,ldata[1,row=T],type="effects")
head(fake)
############## GSAM
res.gam=fregre.gsam(Fat~s(Water,k=4)+Protein+s(x.d2),data=ldata)
fake <- predict.fregre.gsam(res.gam,type="terms")
head(fake)
fake <- predict.fregre.gsam(res.gam,type="effects")
head(fake)
fake <- predict.fregre.gsam(res.gam,ldata[1,row=T],type="terms")
head(fake)
fake <- predict.fregre.gsam(res.gam,ldata[1,row=T],type="effects")
head(fake)
##############################