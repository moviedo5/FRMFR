summary.fregre.func<-function(object,draw=TRUE,...){
    x<-object$x
    ttx=x$argvals
    y<-object$y
    tty=y$argvals
    n=nrow(x)
   if (object$call[[1]]=="fregre.pc") {
     cat(" *** Summary Functional Regression with Penalized Principal Components ***\n")
      object$lm$call<-object$call
      print(summary(object$lm))
#            cat("\n-R squared: ",object$r2)
#      cat("\n-Residual variance: ",
#            object$sr2,"on ",n-object$df," degrees of freedom\n")
       cat("-Lambda penalty: ",object$lambda)
       #     object$lm$call<-object$call
#     print(summary(object$lm))
     var.1<-apply(object$fdata.comp$x, 2, var)
     pr.x= var.1/sum(var.1)
 cat("\n-With",length(object$l),"Principal Components is explained ",round(sum(pr.x[object$l])*100
 ,2),"%\n of the variability of explicative variables. \n -Variability for each  principal components -PC- (%):\n")
    print(round(pr.x[object$l] * 100, 2))
    }
    if (object$call[[1]]=="fregre.pls") {
     cat(" *** Summary Functional Regression with Penalized Partial Least Squares ***\n")
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefs)
            cat("\n-R squared: ",object$r2)
#            cat("\n-Residual variance: ",object$sr2,"\n")
              cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")
       cat("-Lambda penalty: ",object$lambda)
      #     object$lm$call<-object$call
#     print(summary(object$lm))
#     var.1<-apply(object$fdata.comp$x, 2, var)
#     pr.x= var.1/sum(var.1)
# cat("\n-With",length(object$l),"Partial Least Squares is  explained ",round(sum(pr.x[object$l])*100
# ,2),"%\n of the variability of explicative variables. \n -Variability for each Partial Least Squares -PLS- (%):\n")
#    print(round(pr.x[object$l] * 100, 2))
    }
    if (object$call[[1]]=="fregre.basis.func") {
     cat(" *** Summary Functional Response Regression with representation in Basis *** \n")
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefs)
            cat("\n-R squared: ",object$r2)
            cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")

    }
    if (object$call[[1]]=="fregre.np") {
     cat(" *** Summary non-parametric Regression for Functional data *** \n\n")
     cat("-Call: ");    print(object$call)
     cat("\n-Bandwidth (h): ",object$h.opt)
    cat("\n-R squared: ",object$r2)
#    cat("\n-Residual variance: ",object$sr2,"\n")
cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")
    }
    if (!isfdata) {
    if (draw) {
      oldpar <- par()
      C<-match.call()
      lenC=length(C)
      j=1
      while (j<=lenC) {
        if (names(C)[j]=="ask") {
           ask=C[[j]]
           j=lenC +1             }
        else {      j=j+1
                    ask=FALSE             }
       }
       if (ask) {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          }
       else   par(mfrow=c(2,3))
 plot(object$fitted.values,y,xlab="Fitted values",main=paste("R-squared=",
     round(object$r2,2)))
 plot(object$fitted.values,object$residuals,ylab="Residuals",
    xlab="Fitted values",main="Residuals vs fitted.values")
    text(object$fitted.values[i.atypical],object$residuals[i.atypical],
    rownames(x)[i.atypical],cex=0.7)
    abline(h=mean(object$residuals),lwd=1,lty=2)
    abline(h=up,col=2,lwd=2,lty=2)
    abline(h=lo,col=2,lwd=2,lty=2)
#############
resid.sd=sqrt(abs(object$residuals/sd(object$residuals)))
main= "Scale-Location"
ylab23<-"Standardized residuals"
ylim <- c(0, max(resid.sd, na.rm = TRUE))
yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
plot(object$fitted.values,resid.sd, xlab = "Fitted values",
 ylab = yl, main = main,ylim = ylim)
 text(object$fitted.values[i.atypical],resid.sd[i.atypical],
 rownames(x)[i.atypical],cex=0.7)
 plot(diag(object$H),1:nrow(x),xlab="Leverage",ylab="Index.curves",
    main="Leverage")
text(diag(object$H)[i.influence],i.influence,
rownames(x)[i.influence],cex=0.7)
abline(v=times.influ*lim.influ,col=2,lwd=2,lty=2)
#  plot(density(object$residuals),main="Residuals")
    qqnorm(object$residuals,main="Residuals")
    boxplot(object$residuals,main="Residuals")
    par(mfrow=c(1,1))
    } }
#    cat("\n")
return(invisible(list("Influence"=influence,"i.influence"=i.influence,
"i.atypical"=i.atypical)))
}
#
#res1 <-  fregre.basis.func(tempfd, precfd)
#res2 <-  fregre.basis.func(tempfd, precfd,basis.s=smallbasisx,basis.t=smallbasisy)
#res1 <-  fregre.basis.func(tempfd, precfd)
# res1$coefficients
# names(res1)
#  [1] "call"          "alpha.est"     "beta.estbifd"  "fitted.values" "residuals"
#  [6] "x"             "y"             "lambda.s"      "lambda.t"      "Lfdobj.s"
# [11] "Lfdobj.t"      "Dmat"          "Cmat"          "R"             "S"
# [16] "HHCP"          "basis.tt"