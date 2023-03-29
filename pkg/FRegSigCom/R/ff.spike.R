#' @import fda
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom utils data
#' @useDynLib FRegSigCom
#' @name FRegSigCom


#########################################


cv.folds.mat <- function(n,nfolds=5)
  ## Randomly split the n samples into folds
  ## Returns a list of nfolds lists of indices, each corresponding to a fold
{
  d=split(sample(n),rep(1:nfolds,length=n))
  out=matrix(0, n, nfolds)
  for(i in 1:nfolds)
  {
    out[d[[i]],i]=1
  }
  return(out)
}


#########################################

######################################
#' @export
cv.fof.spike=function(X, Y, t.x, t.y, K.cv=5, upper.comp=10, thresh=0.0001)
{
  if(!is.list(X))
  {stop("Error!!: X must be a list!")}
  if (sum(sapply(1:length(X),function(k){!is.matrix(X[[k]])})))
  {stop("Error!!: X must be a list and all its components must be matrix!")
  }
  if(!is.list(t.x))
  {stop("Error!!: t.x must be a list!")}
  if (length(X)!=length(t.x))
  {stop("Error!!: both X and t.x must be lists and they have the same numbers of  components!")
  }
  dim.1=sapply(1:length(X),function(k){dim(X[[k]])[1]})
  if((length(unique(dim.1))!=1))
  {stop("Error!!: all components of X must be matrix and have the same numbers of rows!")
  }
  if((dim(X[[1]])[1]!=dim(Y)[1]))
  {stop("Error!!: the number of observations of X (that is, the number of rows of each component of X) must be equal to the number of observations of Y (that is, the number of rows of Y)!")
  }
  if(sum(sapply(1:length(X), function(k){dim(X[[k]])[2]!=length(t.x[[k]])}))!=0)
  {stop("Error!!: The number of columns of each component of X must be equal to the length of the corresponsing component of  t.x!")
  }
  if(dim(Y)[2]!=length(t.y))
  {stop("Error!!: the number of columns of Y must be equal to the length of the vector t.y of the observation points!")
  }

  n.curves=length(X)
  t.x.list=list()
  max_time=0

  for(i in 1:n.curves)
  {
    t.x.list[[i]]=seq(0,1,length.out=ncol(X[[i]]))
    max_time=max(max_time, ncol(X[[i]]))
  }
  nsample=dim(Y)[1]
  all.folds <- cv.folds.mat(nsample, K.cv)
  db_order=8
  alpha.set.x=c(0.1, 0.5, 1, 2)
  x.params=list()
  x.params[[1]]=t.x.list
  x.params[[2]]=alpha.set.x
  tau_x=c(0.1, 1,10)
  x.params[[3]]=tau_x
  x.params[[4]]=10
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
  lambda.set.x=c(1e-9, 1e-6, 1e-3, 1)
  x.params[[8]]=lambda.set.x
  y.smooth.params = list()
  alpha.set.y=c(0.1, 0.5, 1, 2)
  y.params=list()
  y.params[[1]]=seq(0,1,length.out=ncol(Y))
  y.params[[2]]=alpha.set.y
  tau_y=c(0.1, 1,10)

  y.params[[3]]=tau_y
  y.params[[4]]=10
  y.params[[5]]=wave.vals.8
  upper=2*db_order-1
  y.params[[6]]=upper
  nbasis_levels_y=rep(0,y.params[[4]]+1)
  nbasis_levels_y[1]=2*upper
  for(i in 1:y.params[[4]])
  {
    nbasis_levels_y[i+1]=nbasis_levels_y[i]+2^i-1+upper
    if(nbasis_levels_y[i+1]>max_time)
    {
      break
    }
  }
  max.level.y=i
  y.params[[4]]=max.level.y
  nbasis_levels_y=nbasis_levels_y[1:(i+1)];
  y.params[[7]]=nbasis_levels_y
  lambda.set.y=c(1e-9, 1e-6, 1e-3, 1)
  y.params[[8]]=lambda.set.y

  c_fit_cv=c_cv_ff_spike(t.x.list, X, Y, x.params, y.params, all.folds, upper.comp, thresh)


  return(list(mu=c_fit_cv$mu, Beta=c_fit_cv$Beta_list, c_fit_cv=c_fit_cv))
}


#########################################

######################################
#' @export
pred.fof.spike=function(fit.cv, X.test)
{
  n.test=nrow(X.test[[1]])
  Beta=fit.cv$Beta
  Y.pred= rep(1, n.test)%*%t(fit.cv$mu)
  for(i in 1:length(X.test))
  {
    Y.pred=Y.pred+X.test[[i]]%*%Beta[[i]]/(ncol(X.test[[i]]))
  }
  return(Y.pred)
}

