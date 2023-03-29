#' @rdname predict.fregre.fr
#' @export 
predict.fregre.kam.fr <- function (object, newdata = NULL,...)
{
#  tf <- terms.formula(object$formula)
#  if (attr(tf, "intercept") == 0) intercept = FALSE     else intercept = TRUE
  namesx <- names(object$result)
  nvars <- length(namesx)
  nr <- nrow(newdata[[namesx[1]]])
  X <- list()
  for (i in 1:nvars) {
      rownames(newdata[[namesx[i]]]$data) <-1:nr
      X[[namesx[i]]] <- predict(object$result[[namesx[i]]],newdata[[namesx[i]]],...)
      #aux <- predict.fregre.np.fr(object$result[[namesx[i]]],newdata[[namesx[i]]])
    }
#    if (intercept) X[[nvars + intercept]] = gridfdata(rep(1,nr),object$alpha)
    X[[nvars+1]]=gridfdata(rep(1,nr),object$alpha)
    yhat=Reduce('+',X)
#    yhat <- X / nvars
    return(yhat)
}
