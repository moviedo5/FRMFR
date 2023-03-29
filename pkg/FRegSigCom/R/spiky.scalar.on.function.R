#' @import fda
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom utils data
#' @useDynLib FRegSigCom
#' @name FRegSigCom
#########################################


cv.folds <- function(n,nfolds=5)
  ## Randomly split the n samples into folds
  ## Returns a list of nfolds lists of indices, each corresponding to a fold
{
  return(split(sample(n),rep(1:nfolds,length=n)))
}


#########################################

######################################
#' @export
cv.sof.spike=function(X, Y, t.x, K.cv=10, upper.level=10)
{
  if (!is.list(X))
  {stop("Error: X must be a list!")}
  if (!is.list(t.x))
  {stop("Error: t.x must be a list!")}

  n.curves=length(X)
  t.x.list=list()
  max_time=0
  scale_fac=rep(1, n.curves)
  for(i in 1:n.curves)
  {
    t.x.list[[i]]=seq(0,1,length.out=ncol(X[[i]]))
    max_time=max(max_time, ncol(X[[i]]))
  }
  nsample=length(Y)
  all.folds <- cv.folds(nsample, K.cv)
  lower_bound_set=c(100, 1, 1e-2, 1e-4, 1e-6)
  db_order=8
  lambda.set=c(1)
  alpha.set=c(0.1, 0.5, 1, 2)
  x.params=list()
  x.params[[1]]=t.x.list
  x.params[[2]]=alpha.set
  x.params[[3]]=lambda.set
  x.params[[4]]=upper.level
  x.params[[5]]=wave.vals.8
  upper=2*db_order-1
  x.params[[6]]=upper
  nbasis_levels=rep(0,x.params[[4]]+1)
  nbasis_levels[1]=2*upper
  for(i in 1:x.params[[4]])
  {
    nbasis_levels[i+1]=nbasis_levels[i]+2^i-1+upper
    if(nbasis_levels[i+1]>max_time)
    {
      break
    }
  }
  max.level=i
  x.params[[4]]=max.level
  nbasis_levels=nbasis_levels[1:(i+1)];
  x.params[[7]]=nbasis_levels
  x.params[[8]]=lower_bound_set

  c_fit_cv=c_cv_spiky_scalar_on_function(t.x.list, X, Y, x.params, all.folds)
  coef=c_fit_cv$coef


  return(list( mu=c_fit_cv$mu,  coef=coef, opt.level=c_fit_cv$opt_level, opt.shift=c_fit_cv$opt_shift, opt.lower.bound=c_fit_cv$opt_lower_bound, opt.alpha=c_fit_cv$opt_alpha,
              opt.lambda=c_fit_cv$opt_lambda,  fit.obj=c_fit_cv, x.params=x.params))
}


#########################################

######################################
#' @export
pred.sof.spike=function(fit.cv, X.test)
{

  coef.list=fit.cv$coef
  Y.pred= fit.cv$mu
  for(i in 1:length(X.test))
  {
    Y.pred=Y.pred+X.test[[i]]%*%coef.list[[i]]/(ncol(X.test[[i]]))
  }
  return(Y.pred)
}

