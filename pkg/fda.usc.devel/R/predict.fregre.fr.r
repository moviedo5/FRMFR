#' Predict method for functional response model
#' 
#' @description Computes predictions for regression between functional explanatory variables
#' and functional response.
#' 
#' @aliases predict.fregre.fr
#' @param object \code{fregre.fr} object.
#' @param newdata New functional explanatory data of \code{fdata} class.
#' @param \dots Further arguments passed to or from other methods.
#' @return Return the predicted functional data.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as: \code{\link{fregre.basis.fr}}
#' @keywords regression
#' @examples 
#' \dontrun{ 
#' # CV prediction for CandianWeather data
#' rtt <- c(0, 365)
#' basiss  <- create.bspline.basis(rtt,7)
#' basist  <- create.bspline.basis(rtt,9)
#' nam <- dimnames(CanadianWeather$dailyAv)[[2]]
#' 
#' # fdata class (raw data)
#' tt <- 1:365
#' tempfdata <- fdata(t(CanadianWeather$dailyAv[,,1]),tt,rtt)
#' log10precfdata <- fdata(t(CanadianWeather$dailyAv[,,3]),tt,rtt)
#' rng <- range(log10precfdata) 
#' for (ind in 1:35){
#'  res1 <-  fregre.basis.fr(tempfdata[-ind], log10precfdata[-ind],
#'  basis.s=basiss,basis.t=basist)
#'  pred1 <- predict(res1,tempfdata[ind])
#'  plot( log10precfdata[ind],col=1,ylim=rng,main=nam[ind])
#'  lines(pred1,lty=2,col=2)
#'  Sys.sleep(1)
#' }
#' 
#' # fd class  (smooth data)
#' basis.alpha  <- create.constant.basis(rtt)
#' basisx  <- create.bspline.basis(rtt,65)
#' 
#' dayfd <- Data2fd(day.5,CanadianWeather$dailyAv,basisx) 
#' tempfd <- dayfd[,1]
#' log10precfd <- dayfd[,3]
#' for (ind in 1:35){
#'  res2 <- fregre.basis.fr(tempfd[-ind], log10precfd[-ind],
#'  basis.s=basiss,basis.t=basist)
#'  pred2 <- predict(res2,tempfd[ind])
#'  plot(log10precfd[ind],col=1,ylim=range(log10precfd$coef),main=nam[ind]) 
#'  lines(pred2,lty=2,col=2)
#'  Sys.sleep(.5)
#' }
#' }
#' 
#' @export
predict.fregre.fr <- function(object, newdata, ...){
  if (is.null(object)) stop("No fregre.fd object entered")
  if (missing(newdata)) return(object$fitted.values)
  output <- switch(as.character(object$call[[1]])
                     , fregre.basis.fr = predict.fregre.basis.fr(object, newdata)
                     , fregre.basis.fr.cv = predict.fregre.basis.fr(object, newdata)
                     , fregre.np.fr = predict.fregre.np.fr(object, newdata)
                     , fregre.np.cv.fr = predict.fregre.np.fr(object, newdata)
  )
  return(output)
}

#####################################
predict.fregre.basis.fr <- function (object, newdata){
  #if (object$call[[1]]=="fregre.basis.fr" || object$call[[1]]=="fregre.basis.fr.cv"){
    beta.est <- object$coefficients
    isfdx <- is.fd(newdata)
    if (isfdx) {
      xcoef <- newdata$coef
      ncurves <- ncol(xcoef)
    }
    else {
      xfdobj <- Data2fd(argvals =newdata$argvals, y = t(newdata$data)
                        , basisobj = object$basis.s)
      xcoef <- xfdobj$coef
      ncurves <- ncol(xcoef)
      if (any(newdata$argvals!=object$x$argvals)) stop("Incorrect argvals")
    }
    H = t(xcoef) %*% object$H 
    beta.xest = beta.est %*% t(H)
    beta.xfd   = fd(beta.xest, object$basis.t)
    if (isfdx) {
      #fitted.values  <- fd(coef=yhat, basisobj=object$y$basis, fdnames=object$y$fdnames)
      yhat = eval.fd(object$argvals.y,object$alpha.est) %*% 
                    matrix(1,1,ncurves) + eval.fd(object$argvals.y, beta.xfd)  
      fitted.values   <- smooth.basis(object$argvals.y, yhat, object$y$basis)$fd 
    }
    else {
      yhat = eval.fd(object$y$argvals,object$alpha.est) %*% matrix(1,1,ncurves) + eval.fd(object$y$argvals, beta.xfd)
      fitted.values <- fdata(t(yhat),newdata$argvals,newdata$rangeval,newdata$names)
    }
    return(fitted.values)
}

######################  Auxiliary function ############################

predict.fregre.np.fr <- function (object, newdata){
  y.mat <- object$y$data
  h <- object$h.opt
  x <- object$fdataobj
  n <- nrow(x)
  nn <- nrow(newdata)
  np <- ncol(x)
  if (is.null( rownames(newdata$data)))  
    rownames(newdata$data) <- 1:nn
  par.S <- object$par.S
  bs=as <- list() 
  Ker=object$Ker
  par.metric <- attr(object$mdist$x,"par.metric")
  par.metric[["fdata1"]] <- newdata
  par.metric[["fdata2"]] <- x
  a1 <- attr(object$mdist$x,"call")
  a2 <- attr(object$par.S,"call")
  nmdist <- do.call(a1,par.metric)
  par.S$tt <- nmdist
  par.S$cv=FALSE
  H <- do.call(a2,par.S)
  yp <- H %*% y.mat    
  return(fdata(yp,object$y$argvals,object$y$rangeval,object$y$names))       
}




