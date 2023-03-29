#' @title Summarizes information from fregre.fr objects.
#' 
#' @description Summary function for  \code{\link{fregre.basis.fr}},
#' \code{\link{fregre.np.fr}} and \code{\link{fregre.np.cv.fr}} functions.
#' 
#' Shows:\cr \tabular{ll}{ \tab -Call.\cr \tab -R squared.\cr \tab -Residual
#' variance.\cr \tab -Index of possible atypical curves or possible}
# outliers.\cr \tab -Index of possible influence curves.\cr } If the
# \code{fregre.fr} object comes from the \code{\link{fregre.pc}} then shows:
# \tabular{ll}{ \tab -Variability of explicative variables explained by
# Principal Components.\cr \tab -Variability for each principal components
# -PC-.\cr }
#' 
#' If draw=TRUE plot: \cr \tabular{ll}{ \tab -y vs y fitted values.\cr \tab
#' -Residuals vs fitted values.\cr  -Functional residual boxplot.\cr }
#'  If \code{ask}=FALSE draw graphs in one window, by default. If \code{ask}=TRUE, draw each graph in a window,
#' waiting to confirm.
#' 
#' @aliases summary.fregre.fr  
#' @param object Estimated by functional response regression, \code{fregre.fr} object.
#' @param times.influ Limit for detect possible infuence curves.
# @param times.sigma Limit for detect possible oultiers or atypical curves.
#' @param draw =TRUE draw estimation and residuals graphics.
#' @param \dots Further arguments passed to or from other methods.
# @return 
# \itemize{
# \item {Influence}{ Vector of influence measures.} 
# \item {i.influence}{ Index of possible influence curves.} 
# \item {i.atypical}{ Index of possible atypical curves or possible outliers.}
# }
#' @author Manuel Febrero-Bande and Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
#' @seealso  \code{\link{fregre.basis.fr}}, \code{\link{fregre.np.fr}} and \code{\link{fregre.np.cv.fr}}.
#' @keywords print
# @examples
# \dontrun{
# }
#' 
#' @export 
summary.fregre.fr <- function(object,times.influ=3,#times.sigma=3,
                              draw=TRUE,...){
    if (is.fdata(object$y)){
#      print(class(object$y))
    n=nrow(object$y)
    ry <-object$y$rangeval
    t <- object$y$argvals
    y <- object$y
    yfit=object$fitted.values
    } else {
#      print(class(object$y))
    n=ncol(object$y$coefs)
    ry=object$y$basis$rangeval
    t=seq(ry[1],ry[2],len=max(object$y$basis$nbasis,101))
    y=fdata(object$y,argvals=t,rangeval=ry)
    yfit=fdata(object$fitted.values,argvals=t,rangeval=ry)      
    }
norm.e=norm.fdata(y-yfit)
norm.cen=norm.fdata(y-func.mean(y)) 
norm.ycfit=norm.fdata(yfit-func.mean(y))                         
rss.null = sum(norm.cen^2)
rss=sum(norm.e^2)
r2 <- 1-rss/rss.null
sr2 <- mean(norm.e)


    if (object$call[[1]] == "fregre.basis.fr") {
      cat(" *** Summary Functional Linear Regression with Functional Response (Basis) *** \n")
      cat("-Call: ");    print(object$call)
      cat("\n")
      print(object$coefs)
      ndf=sum(diag(object$H))
    }

    if (object$call[[1]] == "lm") {
    cat(" *** Summary Functional Linear Regression with Functional Response (Basis) *** \n")
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefficients)
            ndf=object$df
    }

    if (object$call[[1]]=="fregre.np") {
     cat(" *** Summary Functional Non-linear Regression with Functional Response*** \n")
     cat("-Call: ");    print(object$call)
     cat("\n-Bandwidth (h): ",object$h.opt)
     ndf=sum(diag(object$H))
      }

    if (object$call[[1]]=="fregre.np.cv.fr") {
#      dc <-  dcor.xy(yfit,y)
      cat(" *** Summary Functional Non-linear Regression with Functional Response *** \n")
      cat("-Call: ");    print(object$call)
      cat("\n-Bandwidth (h): ",object$h.opt)
      ndf=sum(diag(object$H))
#      cat("\n-Distance correlation between fitted values and response: ",dc$estimate)
          }

    cat("\n-1-sum(||Y-hat(Y)||^2)/sum(||Y-bar(Y)||^2): ",r2)
    cat("\n-Residual Mean Norm: ", sr2,"on ",n-ndf," degrees of freedom.\n")

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
      else   par(mfrow=c(2,2))
       #else   par(mfrow=c(2,3))
         
     i.influ <- diag(object$H) > (times.influ*ndf/n)
     
     if (sum(i.influ)>0) influ <- which(i.influ)else influ <- "No influence curves"
     
     cat("-Indices of influence curves:",influ,"\n")


      plot(object$residuals, ylab="Residuals",main=paste("R-squared=",round(r2,3)),col="gray50")
      lines((y-yfit)[i.influ],col=4,pch=21,lwd=1)
#      lines(object$residuals,col=1)

     plot(norm.ycfit,norm.cen,xlab="Centered Fitted values norm",ylab="Centered Response norm")
     points(norm.ycfit[i.influ],norm.cen[i.influ],col=4,pch=16)
     #points(norm.fit[ bb$outpoint],norm.y[ bb$outpoint],col=2)
     
     plot(norm.ycfit,norm.e,xlab="Centered Fitted values norm",ylab="Residuals norm")
     points(norm.ycfit[i.influ],norm.e[i.influ],col=4,pch=16)
     #points(norm.fit2[ bb$outpoint],norm.e[ bb$outpoint],col=2)
     
     plot(norm.cen,norm.e,xlab="Centered Response Norm",ylab="Residuals Norm")
     points(norm.cen[i.influ],norm.e[i.influ],col=4,pch=16)
     #points(norm.ycen[ bb$outpoint],norm.e[ bb$outpoint],col=2)
     
     #points(norm.ycen[out$outliers],norm.e[ out$outliers],col="pink",pch=16)
     
    par(oldpar)
        }
    res=list(rss=rss,rss.null=rss.null,ndf=ndf,R2=r2,RMN=sr2)
#    cat("\n")
#return(invisible(list("Influence"=influence,"i.influence"=i.influence,"i.atypical"=i.atypical)))
    return(invisible(res))
}


#' @export 
print.fregre.fr<-function (x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\n-Call: ", deparse(x$call), "\n", sep = "")
  if (length(coef(x))) {
    cat("\n-Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2,
                  quote = FALSE)
    if (x$call[[1]]=="fregre.lm")      print(x$beta.est[[2]])
  }
  else {
    if (x$call[[1]]=="fregre.np.fr" || x$call[[1]]=="fregre.np.cv.fr") {
      cat("\n-Bandwidth (h): ",x$h.opt)
    }
  }
  cat("\n-R squared: ",x$r2)
  cat("\n-Residual variance: ",x$sr2,"\n")
  invisible(x)
}


