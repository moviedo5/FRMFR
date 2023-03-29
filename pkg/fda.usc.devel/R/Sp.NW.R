#' @name Sp.NW
#' @title Smoothing matrix for multivariate data
#' 
#' @description Provides the smoothing matrix \code{S} for multivariate data using a Gaussian Kernel
#' 
#' @details Options: 
#' \itemize{
#'  \item Nadaraya-Watson kernel estimator (S.NW) with bandwidth parameter \code{h}. 
#  \item Local Linear Smoothing (S.LLR) with bandwidth parameter \code{h}.
#  \item K nearest neighbors estimator (S.KNN) with parameter \code{knn}.
#  \item Polynomial Local Regression Estimator (S.LCR) with parameter of polynomial \code{p} and of kernel \code{Ker}.
#  \item Local Cubic Regression Estimator (S.LPR) with kernel \code{Ker}.
#'  }
#' @aliases Sp.NW
#' @param x \code{matrix}, Multivariate dataset
#' @param h Smoothing parameter or bandwidth. If \code{numeric} a fixed bandidth applied to the variance and covariance matrix of X.  If \code{matrix} is the bandwidth matrix.
#' @param w Optional case weights.
#' @param cv If \code{TRUE}, cross-validation is done (see \code{\link{CV.S}} and \code{\link{GCV.S}}).
#' @return Return the smoothing matrix \code{S} by Nadaraya-Watson estimator using Gaussian kernel.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
#' @seealso See Also \code{\link{S.np}} for functional data object and  \code{\link{fregre.np}} for  non-parametric regression wit scalar response.
#' @references
#' Scott, D. W. (2015). \emph{Multivariate density estimation: theory, practice, and visualization.} John Wiley \& Sons.
#' @keywords smooth
#' @examples
#' \dontrun{
#' data(tecator)
#' x <- tecator$absorp.fdata$data[,c(1,25,50,75,100)]
#' H1 <- Sp.NW(x, cv = T)
#' H2 <- Sp.NW(x, cv = F)
#' y <- tecator$y$Fat
#' yest1 <- H1 %*% y
#' yest2 <- H2 %*% y
#' plot(y,yest1,pch=19)
#' points(y,yest2,col=2,pch=19)
#' }
#' @rdname Sp.NW
#' @export 
Sp.NW <- function (x, h = .1, w = NULL, cv = FALSE) {
  k <- Kerd.mat(x,h)
  n <- NROW(x)
  if (cv) 
    diag(k) = 0
  if (is.null(w)) 
    w <- rep(1, len = n)
  k1 <- sweep(k, 2, w, FUN = "*")
  rw <- rowSums(k1, na.rm = TRUE)
  rw[rw == 0] <- 1e-28
  S = k1/rw
  return(S)
}

######################
mvdnorm=function(x, mu = rep(0,ncol(x)), Sig = diag(ncol(x))){
  if (is.vector(x)) x <- matrix(x,nrow=1)
  p <- ncol(x)
  mmu <- matrix(mu,ncol=p,nrow=nrow(x),byrow=TRUE)
  cte <- sqrt(2*pi)^p * sqrt(det(Sig))
  delta <- diag((x-mmu) %*% solve(Sig) %*% t(x-mmu))
  res <- exp(-0.5*delta)/cte
}
######################

Kerd.mat <- function(x, h = .1){
  if (is.matrix(h)) Sig <- h
  else  Sig <- h * var(x)
  n <- nrow(x)
  mv <- matrix(NA,n,n)
  for (i in 1:n){
    for (j in i:n){
      mv[j,i] <- mv[i,j]<-  mvdnorm(x[i,,drop=F],mu=x[j,,drop=F],Sig=Sig)
    }
  }
  return(mv)
}



