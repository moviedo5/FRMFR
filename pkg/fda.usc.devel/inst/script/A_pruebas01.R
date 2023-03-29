
##
## formula interface:  specify the model by a formula, the method
## fRegress.formula automatically sets up the regression coefficient functions,
## a constant function for the intercept, and a higher dimensional function
## for the inner product with temperature
##

precip.Temp1 <- fRegress(annualprec ~ tempfd)
#  the output is a list with class name fRegress, display names
names(precip.Temp1)
#  the vector of fits to the data is object  precip.Temp1$yfdPar,
#  but since the dependent variable is a vector, so is the fit
annualprec.fit1 <- precip.Temp1$yhatfdobj
#  plot the data and the fit
plot(annualprec.fit1, annualprec, type="p", pch="o")
lines(annualprec.fit1, annualprec.fit1, lty=2)
#  print root mean squared error
RMSE <- sqrt(mean((annualprec-annualprec.fit1)^2))
print(paste("RMSE =",RMSE))
#  plot the estimated regression function
plot(precip.Temp1$betaestlist[[2]])
#  This isn't helpful either, the coefficient function is too
#  complicated to interpret.
#  display the number of basis functions used:
print(precip.Temp1$betaestlist[[2]]$fd$basis$nbasis)
#  25 basis functions to fit 35 values, no wonder we over-fit the data

temp <- fdata(t(CanadianWeather$dailyAv[,,"Temperature.C"]))
res<-fregre.basis(temp,annualprec,basis.b=smallbasis)
ldf <- list("df"=data.frame(logprec=annualprec),"temp"=temp)
res1<-fda.usc:::fregre.lm(logprec~temp,data=ldf,basis.b=list("temp"=smallbasis))
res2<-fregre.lm(logprec~temp,data=ldf,basis.b=list("temp"=smallbasis))
res3<-fregre.glm(logprec~temp,data=ldf,basis.b=list("temp"=smallbasis))
plot(precip.Temp1$betaestlist[[2]])
lines(res$beta.est,col="red")
lines(res2$beta.l[[1]],col=4)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res$beta.est$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res1$beta.l[[1]]$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res2$beta.l[[1]]$coefs)^2)
sum((precip.Temp1$betaestlist[[2]]$fd$coefs-res3$beta.l[[1]]$coefs)^2)
cbind(precip.Temp1$betaestlist[[2]]$fd$coefs,res$beta.est$coefs,res1$beta.l[[1]]$coefs,res2$beta.l[[1]]$coefs)
mean(annualprec)
c(precip.Temp1$betaestlist$const$fd$coefs,precip.Temp1$betaestlist$tempfd$fd$coefs[,1])
res2$coefficients
res$coefficients


sqrt(mean((annualprec-res$fitted.values)^2))