#' @name Descriptive
#' @title Descriptive measures for functional data.
#' 
#' @description Central and dispersion measures for functional data.
#' 
#' @aliases func.mean.formula Descriptive func.mean func.trim.FM func.trim.mode
#' func.trim.RP func.trim.RT func.trim.RPD func.med.FM func.med.mode
#' func.med.RP func.med.RT func.med.RPD func.var func.trimvar.FM
#' func.trimvar.mode func.trimvar.RP func.trimvar.RT func.trimvar.RPD
#' @param formula a formula, such as y ~ group, where y is a fdata object to be
#' split into groups according to the grouping variable group (usually a
#' factor).
#' @param data List that containing the variables in the formula. The item
#' called \emph{"df"} is a data frame with the grouping variable. The item
#' called \emph{"y"} is a fdata object.
#' @param drop logical indicating if levels that do not occur should be dropped
#' (if f is a factor or a list).
#' @param fdataobj \code{\link{fdata}} class object.
#' @param x \code{\link{fdata}} or \code{\link{ldata}}  class object.
#' @param \dots Further arguments passed to or from other methods.  If the
#' argument \code{p} is passed, it used \code{\link{metric.lp}} function, by
#' default \code{p=2}.\cr If the argument \code{trim} (alpha of the trimming)
#' is passed, it used \code{\link{metric.lp}} function.\cr If the argument
#' \code{deriv} (number of derivatives to use) is passed. This parameter is
#' used in \code{\link{depth.RPD}} function, by default it uses \code{deriv
#' =(0,1)}.
#' @param niter Number of iterations.
#' @param eps epsilon.
#' @param metric metric function.
#' @return \code{\link{func.mean.formula}} The value returned from split is a
#' list of fdata containing the mean curves\cr for the groups. The components
#' of the list are named by the levels of f (after converting to a factor, or
#' if already a factor and drop = TRUE, dropping unused levels).\cr
#' 
#' \tabular{ll}{ \tab \code{\link{func.mean}} gives mean curve. \cr \tab
#' \code{\link{func.var}} gives variance curve. \cr \tab
#' \code{\link{func.trim.FM}} Returns the average from the \code{(1-trim)}\%
#' deepest curves following FM criteria. \cr \tab \code{\link{func.trim.mode}}
#' Returns the average from the \code{(1-trim)}\% deepest curves following mode
#' criteria. \cr \tab \code{\link{func.trim.RP}} Returns the average from the
#' \code{(1-trim)}\% deepest curves following RP criteria. \cr \tab
#' \code{\link{func.trim.RT}} Returns the average from the \code{(1-trim)}\%
#' deepest curves following RT criteria. \cr \tab \code{\link{func.trim.RPD}}
#' Returns the average from the \code{(1-trim)}\% deepest curves following RPD
#' criteria. \cr \tab \code{\link{func.med.FM}} Returns the deepest curve
#' following FM criteria. \cr \tab \code{\link{func.med.mode}} Returns the
#' deepest curve following mode criteria. \cr \tab \code{\link{func.med.RP}}
#' Returns the deepest curve following RP criteria. \cr \tab
#' \code{\link{func.med.RPD}} Returns the deepest curve following RPD criteria.
#' \cr \tab \code{\link{func.trimvar.FM}} Returns the marginal variance from
#' the deepest curves followinng FM criteria. \cr \tab
#' \code{\link{func.trimvar.mode}} Returns the marginal variance from the
#' deepest curves followinng mode criteria. \cr \tab
#' \code{\link{func.trimvar.RP}} Returns the marginal variance from the deepest
#' curves followinng RP criteria. \cr 
#' \tab \code{\link{func.trimvar.RT}}
#' Returns the marginal variance from the deepest curves followinng RT
#' criteria. \cr 
#' \tab \code{\link{func.trimvar.RPD}} Returns the marginal
#' variance from the deepest curves followinng RPD criteria. \cr 
#' \tab \code{\link{func.med.L1}}
#' Returns the curve that minimizes the Mean Absolute Error. 
#' \eqn{min_{a\in{S}_n}\sum_{i=1}^n d(X_i,a)}}
#'  
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @references Febrero-Bande, M., Oviedo de la Fuente, M. (2012).
#' \emph{Statistical Computing in Functional Data Analysis: The R Package
#' fda.usc.} Journal of Statistical Software, 51(4), 1-28.
#' \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords descriptive
#' @examples
#' \dontrun{
#' #' # Example with Montreal Daily Temperature (fda-package)
#' fdataobj<-fdata(MontrealTemp)
#' 
#' # Measures of central tendency by group
#' fac<-factor(c(rep(1,len=17),rep(2,len=17)))
#' ldat <- ldata("df"=data.frame(fac),"fdataobj"=fdataobj)
#' a1 <- func.mean.formula(fdataobj~fac,data=ldat)
#' plot(a1)
#' 
#' # Measures of central tendency
#' a1 <- func.mean(fdataobj)
#' a2 <- func.trim.FM(fdataobj)
#' a3 <- func.trim.mode(fdataobj)
#' a4 <- func.trim.RP(fdataobj)
#' a5 <- func.med.FM(fdataobj)
#' a6 <- func.med.mode(fdataobj)
#' a7 <- func.med.RP(fdataobj)
#' a8 <- func.med.L1(fdataobj)
#' 
#' par(mfrow = c(1,2))
#' plot(c(a1,a2,a3,a4), ylim = c(-26,29),
#'      main="Central tendency: trimmed mean")
#' plot(c(a1,a6,a7,a8), ylim = c(-26,29),
#'      main = "Central tendency: median")
#' 
#' ## Measures of dispersion
#' b1 <- func.var(fdataobj)
#' b2 <- func.trimvar.FM(fdataobj)
#' b3 <- func.trimvar.FM(fdataobj,trim=0.1)
#' b4 <- func.trimvar.mode(fdataobj)
#' b5 <- func.trimvar.mode(fdataobj,p=1)
#' b6 <- func.trimvar.RP(fdataobj)
#' b7 <- func.trimvar.RPD(fdataobj)
#' b8 <- func.trimvar.RPD(fdataobj)
#' dev.new()
#' par(mfrow=c(1,2))
#' plot(c(b1,b2,b3,b4),ylim=c(0,79),main="Measures of dispersion I")
#' plot(c(b1,b6,b7,b8),ylim=c(0,79),main="Measures of dispersion II")
#' }
#' 
#' @rdname Descriptive
#' @export
func.mean<-function (x) {
  #if (!is.fdata(x)) 
  #  x <- fdata(x)
#print("func.mean")  
  if (is.fdata(x)) {
     cnames <- colnames(x[["data"]])
     x[["data"]] <- matrix(colMeans(x[["data"]], 
                                   na.rm = TRUE), nrow = 1)
    
    if (!is.null(cnames)) colnames(x[["data"]]) <- cnames
    x$names$main <- "mean"
    xnew<-x
  } 
  if (is.data.frame(x)) {
    #p <- ncol(x)
    xnew <- x[1,,drop=F]
    nam <-names(x)
    clases <- sapply(x,class)
    iclases <- (clases %in% c("numeric","integer","double"))
    #x <- as.data.frame(t(colMeans(x[,iclases,drop=F],na.rm = TRUE)))
    x <- colMeans(---------------------x[,iclases,drop=F],na.rm = TRUE)
    names(x) <- nam[iclases]
    #attributes(x)$df.class <- iclases
    xnew[!iclases] <- NA
    xnew[iclases] <- x
    rownames(xnew)<-"mean"
  }
  if (is.matrix(x)) {
    nam <-colnames(x)
    x <- colMeans(x,na.rm = TRUE)
    xnew[!iclases] <- NA
    xnew[iclases] <- x
  }
  xnew
}



# func.mean <- function(fdataobj){
#   if (!is.fdata(fdataobj)) fdataobj<-fdata(fdataobj)
#   fdataobj[["data"]] <- matrix(colMeans(fdataobj[["data"]],na.rm=TRUE),nrow=1)
#   fdataobj$names$main<-"mean"
#   fdataobj
# }

#' @rdname Descriptive
#' @export
func.var<-function(fdataobj){
  if (!is.fdata(fdataobj)) fdataobj<-fdata(fdataobj)
  n<-dim(fdataobj)[1]
  fdataobj[["data"]]<-(n-1)*apply(fdataobj[["data"]],2,var)/n
  fdataobj[["data"]]<-matrix(fdataobj[["data"]],nrow=1)
  fdataobj$names$main<-"var"
  fdataobj
}

#' @rdname Descriptive
#' @export
func.trim.FM=function(fdataobj,...){depth.FM(fdataobj,...)$mtrim}

#' @rdname Descriptive
#' @export
func.trim.mode=function(fdataobj,...){depth.mode(fdataobj,...)$mtrim}

#' @rdname Descriptive
#' @export
func.trim.RP=function(fdataobj,...){depth.RP(fdataobj,...)$mtrim} 

#' @rdname Descriptive
#' @export
func.trim.RT=function(fdataobj,...){depth.RT(fdataobj,...)$mtrim} 
#' @rdname Descriptive
#' @export
func.trim.RPD=function(fdataobj,...){depth.RPD(fdataobj,...)$mtrim}
#' @rdname Descriptive
#' @export
func.med.FM=function(fdataobj,...){depth.FM(fdataobj,...)$median} 
#' @rdname Descriptive
#' @export
func.med.mode=function(fdataobj,...){depth.mode(fdataobj,...)$median}
#' @rdname Descriptive
#' @export
func.med.RP=function(fdataobj,...){ depth.RP(fdataobj,...)$median}
#' @rdname Descriptive
#' @export
func.med.RT=function(fdataobj,...){ depth.RT(fdataobj,...)$median}
#' @rdname Descriptive
#' @export
func.med.RPD=function(fdataobj,...){ depth.RPD(fdataobj,...)$median}
#' @rdname Descriptive
#' @export
func.trimvar.FM=function(fdataobj,...){
  lista=depth.FM(fdataobj,...)$ltrim
  func.var(fdataobj[lista,])
}
#' @rdname Descriptive
#' @export
func.trimvar.mode=function(fdataobj,...){
  lista=depth.mode(fdataobj,...)$ltrim
  func.var(fdataobj[lista,])
  }
#' @rdname Descriptive
#' @export
func.trimvar.RP=function(fdataobj,...){
 lista=depth.RP(fdataobj,...)$ltrim
 func.var(fdataobj[lista,])}
#' @rdname Descriptive
#' @export
func.trimvar.RPD=function(fdataobj,...){
 lista=depth.RPD(fdataobj,...)$ltrim
 func.var(fdataobj[lista,])}
#' @rdname Descriptive
#' @export
func.trim.RT=function(fdataobj,...){depth.RT(fdataobj,...)$mtrim}
#' @rdname Descriptive
#' @export
func.med.RT=function(fdataobj,...){ depth.RT(fdataobj,...)$median}
#' @rdname Descriptive
#' @export
func.trimvar.RT=function(fdataobj,...){
 lista=depth.RT(fdataobj,...)$ltrim
 func.var(fdataobj[lista,])}

#  for multivarate data
# func.trim.SD=function(fdataobj,...){depth.SD(fdataobj,...)$mtrim}
# func.trim.PD=function(fdataobj,...){depth.PD(fdataobj,...)$mtrim}
# func.trim.HD=function(fdataobj,...){depth.HD(fdataobj,...)$mtrim} 
# func.trim.MhD=function(fdataobj,...){depth.MhD(fdataobj,...)$mtrim}

# func.med.SD=function(fdataobj,...){depth.SD(fdataobj,...)$median} 
# func.med.PD=function(fdataobj,...){depth.PD(fdataobj,...)$median}
# func.med.HD=function(fdataobj,...){ depth.HD(fdataobj,...)$median}
# func.med.MhD=function(fdataobj,...){ depth.MhD(fdataobj,...)$median}

#' @rdname Descriptive
#' @export
func.med.L1 <- function(fdataobj, 
                        eps=0.001, 
                        niter=100, 
                       # draw=FALSE, 
                       metric=metric.lp,...){
  
  if (!is.fdata(fdataobj)) { stop("fdataobj is not a fdata object")}
    
  # if (draw) plot(fdataobj, col= "gray40")
  #y0 <- func.mean(fdataobj)
  y0 <- (fdataobj[1] + fdataobj[2])/2
  etol <- eps * norm.fdata(y0,metric=metric,...)
  conv <- FALSE; iter <- 0
  while(!conv){
    norm <- norm.fdata(fdataobj - y0, metric=metric,...)
    wi <- 1/norm
    wi <- wi/sum(wi,na.rm=TRUE)
    lna <- is.na(wi)
    if (any(lna)) {eta=1} else eta=0
    Rmu <- gridfdata(coef=matrix(wi[!lna],nrow=1),fdataobj[!lna]-y0)
    gammamu <- min(1,eta/norm.fdata(Rmu))
    y1 <- (1-eta)*gridfdata(coef=matrix(wi[!lna],nrow=1),fdataobj[!lna])+eta*y0
    # if (draw) {
    #   plot(fdataobj, col= "gray40")
    #   #	lines(y0, col= "white", lwd = 2)
    #   lines(y1, col= 2, lwd = 3)}
    if (norm.fdata(y0 - y1, metric=metric,...)>=etol & (iter < niter))
    {
      y0 <- y1
      iter <- iter+1
    } else {conv <- TRUE} 
  }
  #print(iter)
  attr(y1,"niter") <- iter
  return(y1)
}


