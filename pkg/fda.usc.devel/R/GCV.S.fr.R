#' The generalized correlated cross-validation (GCCV) score
#' 
#' @description Compute the  generalized correlated cross-validation (GCV) score. 
#' @param y fdata class object.
#' @param S Smoothing matrix, see \code{\link{S.NW}}, \code{\link{S.LLR}} or
#' @param criteria The penalizing function. By default \emph{"Rice"} criteria. 
#'  Possible values are \emph{"GCCV1"}, \emph{"GCCV2"}, \emph{"GCCV3"}, \emph{"GCV"}. 
#' @param W Matrix of weights.
#' @param trim The alpha of the trimming.
#' @param andf Degree of freedom
#' @param draw =TRUE, draw the curves, the sample median and trimmed mean.
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param par.metric Further arguments passed to metric.
#' 
#' @return { Returns GCV score calculated for input parameters.  }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{fregre.np.cv.fr}} \cr Alternative method:
#' \code{\link{CV.S.fr}}
#' @keywords utilities
#' @export
GCV.S.fr=function(y, S, criteria = "GCV", W = NULL, trim = 0, andf=0,draw = FALSE,metric = metric.lp,par.metric=list(lp=2)) 
{
    if (!is.fdata(y)) stop("The function GCV.S.fr is for fdata objects")
    n = ncol(S)
    nn <- nrow(y)
	if (nn!=n) stop("Dimensions of S and y does not coincide")
   if (is.null(W)) W <- diag(nn)
	 if (criteria=="CV") {out=CV.S.fr(y=y,S=S,W=W,trim=trim,andf=andf,draw=draw,metric=metric,par.metric=par.metric)
	 } else {
    tab = list("GCV", "AIC", "FPE", "Shibata", "Rice")
    type.i = pmatch(criteria, tab)
  		 y.est = S %*% y$data
       y.est <- fdata(y.est, y$argvals, y$rangeval, y$names)
       e <- y - y.est
	   ee=drop(do.call("norm.fdata",c(list(fdataobj=e,metric=metric),par.metric)))^2
	   ee[is.nan(ee)]=Inf
        if (trim > 0) {
            e.trunc = quantile(ee, probs = (1 - trim), type = 4)
            ind <- ee <= e.trunc
            if (draw) plot(y, col = (2 - ind))
            l <- which(abs(ee) <= e.trunc)
            res = mean(ee[ind], na.rm = TRUE)
        }else {l=1:n;res = mean(ee, na.rm = TRUE)}
    mdf <- sum(diag(S)[l])
    u=(mdf+andf)/length(l)
    if (is.na(type.i)) {
        if (u > 0.5) 
            vv = Inf
        else vv = 1/(1 - 2 * u)
    }else {
        vv <- switch(type.i, `1` = if (type.i == 1) vv = (1 - u)^(-2),
             `2` = if (type.i == 2) vv = exp(2 * u),
             `3` = if (type.i == 3) vv = (1 + u)/(1 - u), 
            `4` = if (type.i == 4) vv = (1 + 2 * u), 
            `5` = if (type.i == 5) {
                if (u > 0.5) vv = Inf else vv = 1/(1 - 2 * u)
            })
    }
    out <- res * vv/length(l)
    attr(out, "mdf") <- mdf
	 attr(out, "criteria") <- criteria
	 }
    return(out)
}
