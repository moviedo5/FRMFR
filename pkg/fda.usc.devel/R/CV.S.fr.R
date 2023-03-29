#' The cross-validation (CV) score
#' 
#' @description Compute the leave-one-out cross-validation score.
#' 
#' @param y Matrix of set cases with dimension (\code{n} x \code{m}), where
#' \code{n} is the number of curves and \code{m} are the points observed in
#' each curve.
#' @param S Smoothing matrix, see \code{\link{S.NW}}, \code{\link{S.LLR}} or
#' \eqn{S.KNN}.
#' @param W Matrix of weights.
#' @param trim The alpha of the trimming.
#' @param draw =TRUE, draw the curves, the sample median and trimmed mean.
#' @param andf Degree of freedom
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param par.metric Further arguments passed to metric.
##' @return { Returns CV score calculated for input parameters.  }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{fregre.np.cv.fr}} \cr Alternative method:
#' \code{\link{GCV.S.fr}}
#' @keywords utilities
#' @export CV.S.fr
CV.S.fr=function (y, S, W = NULL, trim = 0, andf=0,draw = FALSE, metric = metric.lp,par.metric=list()) 
{
    n = ncol(S)
    if (!is.fdata(y)) stop("Object y is not fdata object")
	 nn <- nrow(y)
    if (is.null(W))  W <- diag(nn)
		S.cv=S
		diag(S.cv)=0
		ss=apply(S.cv,1,sum)
		S.cv=sweep(S.cv,1,ss,"/")
        y.cv = S.cv %*% y$data
	     y.cv[is.nan(y.cv)]=Inf
        y.cv <- fdata(y.cv, y$argvals, y$rangeval, y$names)
		e= y-y.cv
        ee <- drop(do.call(norm.fdata,c(list(fdataobj=e,metric = metric), par.metric)))^2
		  ee[is.nan(ee)]=Inf
        if (trim > 0) {
            e.trunc = quantile(ee, probs =  (1 - trim), type = 4)
            ind <- ee <= e.trunc
			   if (draw) plot(y, col = (2 - ind))
            mdf=sum(diag(S[ind,ind]))
			   if ((mdf+andf)/sum(ind)>0.66) {res=Inf} else {res = sum(ee[ind])/(sum(ind)-mdf)} 
			}
        else {
				mdf=sum(diag(S))
				if ((mdf+andf)/n>0.66) {res=Inf} else {res = sum(ee)/(n-mdf)}
				}
    if (is.nan(res)) res = Inf
	 attr(res, "mdf") <- mdf
	 attr(res, "crit") <- "CV"
    return(res)
}
