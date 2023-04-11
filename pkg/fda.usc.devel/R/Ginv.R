#' Generalized inverse
#' 
#' @description This function computes the Matrix inverse.
#' @param X a square numeric or complex matrix. 
#' @param tol the tolerance for detecting linear dependencies in the columns of \code{a}. 
#' The default is `.Machine$double.eps.`
#' @details By default it uses the \code{\link{solve}} function. 
#' If system is computationally singular, the matrix inverse is computed by
#'  \code{\link{svd}} function (using the effective rank).
#' @return {  Return the inverse of \code{a} }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{solve}}, \code{\link{svd}} and \code{\link{complex}} 
#' @keywords utilities
#' @export


Minverse <- function (X, tol = sqrt(.Machine$double.eps)) 
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Minv <-try(solve(X),silent=TRUE)  
  if (!is(Minv,"try-error")) {   
    Minv
  }
  else{ 
  
  Xsvd <- svd(X)
  if (is.complex(X)) 
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
#  warning("System is computationally singular:\n
#          -Rank of X is: ",dim(X)[2],
#"The inverse of X is computed by svd procedure, \n
#         -Effective Rank of X is: ",sum(Positive))  
  
  warning("System is computationally singular (rank  ",dim(X)[2],")\n
          The  matrix inverse is computed by svd (effective rank ",sum(Positive),")")    
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                         t(Xsvd$u[, Positive, drop = FALSE]))
  }   
}
