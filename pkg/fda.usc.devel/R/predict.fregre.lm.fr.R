#' @title Predict method for functional response model
#' 
#' @description 
#' Computes predictions for regression between functional (and non functional)
#' explanatory variables and functional response. 
#' \itemize{ 
#' \item \code{predict.fregre.lm.fr}, Predict method for functional linear model of
#' \code{\link{fregre.lm.fr}} fits object using basis or principal component
#' representation.
#' \item \code{predict.fregre.mlm.fr}, Predict method for functional linear model of
#' \code{\link{fregre.mlm.fr}} fits object using basis or principal component
#' representation.
#' \item \code{predict.fregre.sam.fr}, Predict method for functional additive model of
#' \code{\link{fregre.sam.fr}} fits object using basis or principal component
#' representation
#' }
#' 
#' \code{predict.fregre.lm.fr} uses the model fitting function \code{\link{lm}} properties.\cr
#' If using functional data derived, is recommended to use a number of bases
#' to represent beta  lower than the number of bases used to represent the
#' functional data. \cr
#' 
#' The first item in the \code{data} list of \code{newx} argument is called
#' \emph{"df"} and is a data frame with the response and non functional
#' explanatory variables, as \code{\link{lm}}. 
#' Functional variables (\code{fdata} and \code{fd} class)
#' are introduced in the following items in the \code{data} list of \code{newx}
#' argument.
#' 
#' @aliases predict.fregre.lm.fr
#' @param object \code{fregre.fr} object.
#' @param newx An optional data list in which to look for variables with which
#' to predict. If omitted, the fitted values are used. List of new explanatory
#' data.
#' @param type a character vector, Type of prediction: (\code{response}, \code{terms} for model terms or \code{effects}  for model terms where
#' the partial effects are summarized for each functional variable.
#' @param se.fit =TRUE (not default) standard error estimates are returned for
#' each prediction.
#' @param scale Scale parameter for std.err. calculation.
#' @param df Degrees of freedom for scale.
#' @param interval Type of interval calculation.
#' @param level Tolerance/confidence level.
#' @param pred.var the variance(s) for future observations to be assumed for
#' prediction intervals. See \code{link{predict.lm}} for more details.
#' @param weights variance weights for prediction. This can be a numeric vector
#' or a one-sided model formula. In the latter case, it is interpreted as an
#' expression evaluated in newdata
#' @param \dots Further arguments passed to or from other methods.
#' @return Return the predicted values and optionally:
#' \itemize{
#' \item {predict.lm.fr}{ produces a vector of predictions
#' or a matrix of predictions and bounds with column names fit, lwr, and upr if
#' interval is set. If se.fit is TRUE, a list with the following components is
#' returned: fit vector or matrix as above.} 
#' \item {se.fit}{ standard error of predicted means.} 
#' \item {residual.scale}{ residual standard deviations.}
#' \item {df}{ degrees of freedom for residual.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' 
#' @seealso See Also as: \code{\link{fregre.lm.fr}}. 
#' 
#' @references Febrero-Bande, M., Oviedo de la Fuente, M. (2012).
#' \emph{Statistical Computing in Functional Data Analysis: The R Package
#' fda.usc.} Journal of Statistical Software, 51(4), 1-28.
#' \url{https://www.jstatsoft.org/v51/i04/}
#' 
#' @keywords regression
#' @examples
#' \dontrun{
#' library(fda.usc.devel)
#' data(tecator)
#' absorp=tecator$absorp.fdata
#' ind=1:129
#' x=fdata.deriv(absorp[ind,])
#' y=fdata.deriv(x,1)
#' x2 <- fdata.deriv(y)

#' mdat<-mfdata("y"=y,"x"=x,"x2"=x2)
#' res0 <-  fregre.lm.fr(as.formula(y~x+x2),mdat)
#' pred0 <- predict.fregre.lm.fr(res0,mdat)
#' plot(res0$fitted.values-pred0)
#' 
#' res1 <-  fregre.mlm.fr(as.formula(y~x+x2),mdat)
#' pred1 <- predict.fregre.mlm.fr(res1,mdat)
#' plot(res0$fitted.values-pred0)

#' by<-create.bspline.basis(mdat$y$rangeval,5)
#' res0 <-  fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=by)
#' pred0 <- predict.fregre.lm.fr(res0,mdat)
#' plot(res0$fitted.values-pred0)
#' 
#' res1 <-  fregre.mlm.fr(as.formula(y~x+x2),mdat,basis.y=by)
#' pred1 <- predict.fregre.mlm.fr(res1,mdat)
#' plot(res1$fitted.values-pred0)
#' plot(pred1-pred0)

#' by<-create.pc.basis(mdat$y,1:5)
#' res0 <-  fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=by)
#' pred0 <- predict.fregre.lm.fr(res0,mdat)
#' plot(res0$fitted.values-pred0)
#' 
#' res1 <-  fregre.mlm.fr(as.formula(y~x+x2),mdat,basis.y=by)
#' pred1 <- predict.fregre.mlm.fr(res1,mdat)
#' plot(res1$fitted.values-pred0)
#' plot(pred1-pred0)
#' 
#' 
#' }                                                                                                              
#' @rdname predict.fregre.lm.fr
#' @export 
predict.fregre.lm.fr <- function(object, newx = NULL,
                                  type = "response",...){
  if (is.null(object)) stop("No fregre.lm.fr object entered")
  # print(1)
  # object <-res
  # newx <- ldata
  # type <- "response"
  if (is.null(newx)) {
    if (type == "effects"){
      # fake  = predict.lm(object, type = "terms", ...) 
      # yp <- effect.gam(object,fake)
    } else{
      # yp  = predict.gam(object, type = type, ...)
      yp <- object$fitted.values
    }
    return(yp)
  } else {
    # print(2)
    
    data=newx
    basis.x=object$basis.x
    basis.y=object$basis.y
    formula=object$formula.ini
    tf <- terms.formula(formula)
    terms <- attr(tf, "term.labels")
    if (length(terms)==0) return(rep(object$coefficient,length=nrow(newx[[1]])) ) 
    nt <- length(terms)
    if (attr(tf, "response") > 0) {
      response <- as.character(attr(tf, "variables")[2])
    }
    vtab<-rownames(attr(tf,"factors"))
    name.coef <- NULL
    vfunc <- object$vfunc
    vnf <- object$vnf
    nnf<-length(vnf)
    # print(3)
    if (!is.null(vnf)) {
      first=FALSE
      XX=data.frame(data[["df"]][,c(vnf)])
      names(XX)=vnf
    } else {  
      first=TRUE
      XX <- NULL
    }
    lenfunc<-length(vfunc)
    bsp <- TRUE
    if (object$basis.y.class== "fdata.comp")
    bsp <- FALSE
    
    raw <- object$raw
    # print(4)
    if (lenfunc>0) {
      k=1
      mean.list=vs.list=JJ=list()
      for (i in 1:lenfunc) {
        # print("i")
        # print(class(newx[[vfunc[i]]]))
        if(class(newx[[vfunc[i]]])[1]=="fdata"){
          #fdataobj<-data[[vfunc[i]]]
          fdat<-data[[vfunc[i]]]; 
          dat<-fdat$data
          tt<-fdat[["argvals"]]
          rtt<-fdat[["rangeval"]]
          if (nrow(dat)==1) rwn<-NULL         else rwn<-rownames(dat)
          # print(vfunc[i])
          # print(dim(fdat))
          xaux <- fdata2basis(fdat,basis.x[[vfunc[i]]])#,method ="inprod")
          name.coef[[vfunc[i]]] <- colnames(xaux$coefs) <- paste(vfunc[i],".",colnames(xaux$coefs),sep="")
          Z <- xaux$coefs
         # print(dim(Z))
          if (first) {
            XX=Z
            first=FALSE
          }   else {
            XX = cbind(XX,Z)} }
      }
    }
    # print(5)
    # print(head(XX))
    # print(head(object$XX$x2))
    yfdata <- object$data[[response]]
    tty <- yfdata[["argvals"]]
    rtty <- yfdata[["rangeval"]]
    namy <- yfdata[["names"]]
    # print(rtty)
    # print("rtty")
    npy<-NCOL(yfdata)
    if (!is.data.frame(XX)) 
      XX=data.frame(XX)
    # if (is.null(object$basis.y)){
    #   print("raw")
    #   yp<-object$fitted.values
    #   yp$data<-matrix(NA,nrow=nrow(XX),ncol=ncol(object$fitted.values))
    #   print(dim(yp))
    #   for (i in 1:npy)  {
    #     # print(i)
    #    object$result$coefficients <- object$coefficients[i,]
    #    yp$data[,i]=predict.lm(object=object$result,newdata=XX,type=type,...)
    #    # yp = object$a.est * rep(1, len = nn) + Z %*%  object$b.est
    #   }  
    #   return(yp)
    # }
    # print(6)
      npy <- NCOL(object$coefficients)
      yp <- matrix(NA,NROW(XX),npy)
      # print(npy)
      # print(dim(yp))
      zcoef <- object$coefficients
      for (i in 1:npy)  {
        object$coefficients <- zcoef[,i,drop=T]
#print(object$coefficients)
#print(head(XX))
        # data.frame(object$XX$x2)
        yp[,i]=predict.lm(object=object,newdata=XX,type=type,...)
      }  
      # print(7)
    if (object$raw) {
      yp <- fdata(yp,tty,rtty,namy)
    }    else{
      if (bsp) {
        yp <- fd(t(yp),basis.y)
        yp <- fdata(yp,tty,rtty,namy)
      } else  {
      # print("PC");        print(dim(yp))
        #yp <- fdata(yp,tty,rtty,namy)
#         yp <- gridfdata(yp,object$basis.y$basis)
        yp <- gridfdata(yp,object$basis.y$basis,basis.y$mean)
         #yp <- fdata(yp,tty,rtty,namy)
#         yp$data<-sweep(yp$data,2,matrix(basis.y$mean$data,ncol=1),"+")              
        # if (basis.y.class=="fdata.comp"){
        # print(dim(yp))
        # print(dim( basis.y$basis$data))
      #  yhatfdata <- (yp) %*% basis.y$basis$data
       #  yhatfdata<-sweep(yhatfdata,2,matrix(basis.y$mean$data,ncol=1),"+")      
        # yp <- fdata(yhatfdata,tty,rtty,namy)
        # }
        
        
      }
      }
    }                          
    return(yp)
}      
