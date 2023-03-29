#' Cross-validation functional regression with functional response using kernel
#' estimation.
#' 
#' @description Computes functional regression between a functional explanatory variable and
#' a functional response using an asymmetric kernel estimation by cross-validation method.
#' 
#' @details The non-parametric functional regression model can be written as follows
#' \deqn{ y_i =r(X_i) + \epsilon_i } where the unknown smooth real function
#' \eqn{r} is estimated using kernel estimation by means of
#' \deqn{\hat{r}(X)=\frac{\sum_{i=1}^{n}{K(h^{-1}d(X,X_{i}))y_{i}}}{\sum_{i=1}^{n}{K(h^{-1}d(X,X_{i}))}}}
#' where \eqn{K} is an kernel function (see \code{Ker} argument), \code{h} is
#' the smoothing parameter and \eqn{d} is a metric or a semi-metric (see
#' \code{metric} argument).
#' 
#' The function estimates the value of smoothing parameter (also called
#' bandwidth) \code{h} through Generalized Cross-validation \code{GCV}
#' criterium, using \code{GCV.S.fr} or \code{CV.S.fr} internal functions.
#' 
#' It computes the distance between curves using the
#' \code{\link{metric.lp}}, although any other semimetric could be used (see
#' \code{\link{semimetric.basis}} or \code{\link{semimetric.NPFDA}} functions).
#' Different asymmetric kernels can be used, see
#' \code{\link{Kernel.asymmetric}}.\cr
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Scalar response with length \code{n}.
#' @param h Bandwidth, \code{h>0}. Default argument values are provided as the
#' sequence of length 25 from 2.5\%--quantile to 25\%--quantile of the distance
#' between \code{fdataobj} curves, see \code{\link{h.default}}.
#' @param Ker Type of asymmetric kernel used, by default asymmetric normal
#' kernel.
#' @param metric List with components \code{x},\code{y} containing the metric 
#' to be applied to covariate and response, respectively.
#' @param par.metric List of components \code{x},\code{y} containing optional 
#' parameters for each metric in \code{metric}.
#' @param type.CV Type of cross-validation. By default generalized
#' cross-validation \code{\link{GCV.S}} method.
#' @param type.S Type of smothing matrix \code{S}. By default \code{S} is
#' calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).
#' @param par.CV List of parameters for \code{type.CV}: \code{trim}, the alpha
#' of the trimming\cr and \code{draw=TRUE}.
#' @param par.S List of parameters for \code{type.S}: \code{w}, the weights.
#' @return Return:
#' \itemize{
#' \item \code{call}{ The matched call.} 
#' \item \code{residuals}{ \code{y} minus \code{fitted values}.} 
#' \item \code{fitted.values}{ Estimated scalar response.} 
#' \item \code{df}{ The residual degrees of freedom.} 
#' \item \code{r2}{ Coefficient of determination.} 
#' \item \code{sr2}{ Residual variance.} 
#' \item \code{H}{ Hat matrix.} 
#' \item \code{y}{ Response.} 
#' \item \code{fdataobj}{ Functional explanatory data.}
#' \item \code{mdist}{ Distance matrix between \code{x} and \code{newx}.} 
#' \item \code{Ker}{ Asymmetric kernel used.} 
#' \item \code{gcv}{ CV or GCV values.} 
#' \item \code{h.opt}{ smoothing parameter or bandwidth that minimizes CV or GCV method.} 
#' \item \code{h}{ Vector of smoothing parameter or bandwidth.} 
#' \item \code{fit.CV}{ List with the fitted values and residuals estimated by CV, without the same curve.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.np.fr}},
#' \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}} .\cr
#' Alternative method: \code{\link{fregre.basis.cv}} and
#' \code{\link{fregre.np.cv}}.
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York.
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples 
#' \dontrun{
#' data(tecator)
#' absorp=tecator$absorp.fdata
#' ind=1:129
#' x=absorp[ind,]
#' y=fdata.deriv(x,2)
#' Ker=AKer.tri
#' res.np=fregre.np.cv.fr(x,y,Ker=Ker,crit="Shibata")
#' summary(res.np)
#' res.np2=fregre.np.cv(x,y,type.CV=GCV.S,crit="Shibata")
#' summary(res.np2)
#' 
#' ## Example with other semimetrics (not run)
#' res.pca1=fregre.np.cv(x,y,Ker=Ker,metric=list(x=metric.lp,y=semimetric.pca),
#' par.metric=list(x=list(lp=1),y=list(q=3))
#' summary(res.pca1)
#' }
#' 
#' @export
fregre.np.cv.fr<-function (fdataobj, y, h = NULL, Ker = AKer.norm, metric = list(x=metric.lp,y=metric.lp), 
    type.CV = CV.S.fr, type.S = S.NW, par.CV = list(trim = 0.02), par.metric=list(x=list(), y=list()),
    par.S = list(w = 1)) 
{
    if (is.function(type.CV)){
        tcv <- deparse(substitute(type.CV))
    } else tcv <- type.CV
    if (is.function(type.S)){ 
        ty <- deparse(substitute(type.S))
    } else ty <- type.S
    if (!is.fdata(fdataobj)) fdataobj <- fdata(fdataobj)
    if (!is.fdata(y)) stop("y is not an fdata object")
    nasx <- is.na.fdata(fdataobj)
    nasy <- is.na.fdata(y)
    if (is.null(names(y$data))) names(y$data) <- seq_len(length(y))
    if (any(nasx) | any(nasy)) {
        bb <- !nasx & !nasy
        if (ops.fda.usc()$warning) 
            warning(sum(!bb), " covariate or response curves with NA are omited\n")
        fdataobj$data <- fdataobj$data[bb, ]
        y$data <- y$data[bb,]
    }
    x <- fdataobj[["data"]]
#    tt <- fdataobj[["argvals"]]
#    rtt <- fdataobj[["rangeval"]]
    C <- match.call()
    m <- match(c("x", "y", "h", "Ker", "metric", "par.metric","type.CV", "type.S", "par.CV", "par.S"), names(C), 0L)
    n <- nrow(x)
    ny <- nrow(y)
    np <- ncol(x)
    if (n != ny) stop("Size of covariate does not coincide with response")
    if (is.null(rownames(x))) rownames(x) <- seq_len(n)
    if (is.null(colnames(x))) colnames(x) <- seq_len(np)
    tty <- y$argvals
    rtty <- y$rangeval
#    types <- FALSE
    if (is.matrix(metric$x)){ 
        mdistx <- metric$x
		metric$x <- get(attributes(mdistx)$call)
		par.metric$x <- attributes(mdistx)$par.metric
    } else mdistx <- do.call(metric$x, c(list(fdata1=fdataobj, fdata2=fdataobj), par.metric$x))
	if (is.matrix(metric$y)){
	mdisty <- metric$y
	metric$y <- get(attributes(mdisty)$call)
	par.metric$y <- attributes(mdisty)$par.metric
	} else {
	mdisty <- do.call(metric$y, c(list(fdata1=y, fdata2=y), par.metric$y))
	}
#    ke <- deparse(substitute(Ker))
#    if (!is.function(Ker)) Ker <- get(Ker)
    if (is.character(Ker)) {
            nker <- function(u,mik=Ker){0.5*get(mik)(u)}
            Ker <- get(Ker)
            } else {
            nker <- function(u,mik=Ker){0.5*mik(u)}
            }    
    attr(par.S, "call") <- ty
    if (is.null(h)){ 
#        h = h.default(fdataobj,metric=metric$x)
		  nker <- get(paste0("Ker.",unlist(strsplit(deparse(substitute(Ker)),"[.]"))[2]))
#		  h = do.call(h.default,c(list(fdataobj=fdataobj,metric=mdistx,prob=c(0.02,.3),Ker=nker),par.metric$x))
		  h <- do.call(h.default, c(list(fdataobj=fdataobj, metric=mdistx, prob=c(seq(0.025,.2,len=41),.25,.3,.4,.5,.6), Ker=nker), par.metric$x))
          h <- c(h, Inf)
    } else {
    if (any(h <= 0)) 
            stop("Error: Invalid range for h")
    }
    lenh <- length(h)
    cv <- gcv1 <- gcv <- cv.error <- array(Inf, dim = c(lenh))
    par.S2 <- par.S
    if (is.null(par.S2$h)) par.S$h <- h
    if (is.null(par.S$Ker)) par.S$Ker <- Ker
#    y.est.cv <- y.est <- matrix(NA, nrow = ny , ncol = ncol(y$data))
    par.S$tt <- mdistx
    par.CV$metric <- metric$y
	par.CV$par.metric<- par.metric$y
    for (i in seq_len(length(h))) {
        par.S$h <- h[i]
        par.S$cv <- TRUE
#        H.cv <- do.call(ty, par.S)
        par.S$cv <- FALSE
        H = do.call(ty, par.S)
        par.CV$S <-H
        gcv[i] <- do.call(tcv,c(list(y=y),par.CV))
    }
  if (all(!is.finite(gcv))) stop(paste0("All GCV values are infinite. h too short:",paste0(round(range(h),3),collapse="-")," or reconsider the value for trim:",par.CV$trim))
    l <- which.min(gcv)
    h.opt <- h[l]
    if ((h.opt == min(h) || h.opt==max(h)) & ops.fda.usc()$warning) 
        warning(" Warning: h.opt is at the extreme of the range provided\n   provided, range(h)=",range(h), "\n")
    par.S$tt <- mdistx
    par.S$h <- h.opt
    par.S$cv <- FALSE
    H <- do.call(ty, par.S)
    yp <- H %*% y$data
    par.S$cv <- TRUE
    Hcv <- do.call(ty, par.S)
    ypcv <- Hcv %*% y$data
    df <- fdata.trace(H)
    names(gcv) <- h
    names(cv) <- h
    yp <- fdata(yp, tty, rtty)
    rownames(yp$data) <- rownames(y$data)
    e <- y - yp
    ypcv <- fdata(ypcv, tty, rtty)
    rownames(ypcv$data) <- rownames(y$data)
    ecv <- y - ypcv
    ycen <- fdata.cen(y)$Xcen
    norm.e <- drop(do.call("norm.fdata",c(list(fdataobj=e, metric = metric$y), par.metric$y)))^2
    norm.ycen<-do.call("norm.fdata",c(list(fdataobj=ycen,metric=metric$y), par.metric=par.metric$y))^2
    sr2 <- sum(norm.e)/(n - df)
    r2 <- 1 - sum(norm.e)/sum(norm.ycen)
	yp2 <- ecv$data %*% t(ecv$data)
        out <- list(call = C, fitted.values = yp, H = H, residuals = e, 
        df = df, r2 = r2, sr2 = sr2, var.y = yp2, y = y, 
        fdataobj = fdataobj, mdist = list(x=mdistx,y=mdisty), Ker = Ker, metric = metric,par.metric=par.metric, 
        type.S = type.S, par.S = par.S, gcv = gcv, h.opt = h.opt, 
        h = h, fit.CV = list(fitted.values = ypcv, residuals = ecv))
    class(out) <- "fregre.fr"
    return(out)
}

