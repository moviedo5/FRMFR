#' @import fda
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom stats approx
#' @useDynLib FRegSigCom
#' @name FRegSigCom

#########################################

cv.folds <- function(n,nfolds=5)
  ## Randomly split the n samples into folds
  ## Returns a list of nfolds lists of indices, each corresponding to a fold
{
  #return(split(sample(n),rep(1:nfolds,length=n)))
  return(split(1:n,rep(1:nfolds,length=n)))
}


#########################################

#' @export
cv.msof = function(X, Y, t.x.list, nbasis = 50, K.cv = 5, upper.comp = 10, thresh = 0.001)
{
  if (!is.list(X))
  {
    stop("Error!!: X must be a list!")
  }
  if (sum(sapply(1:length(X), function(k) {
    !is.matrix(X[[k]])
  })))
  {
    stop("Error!!: X must be a list and all its components must be matrix!")
  }
  if (!is.list(t.x.list))
  {
    stop("Error!!: t.x.list must be a list!")
  }
  if (length(X) != length(t.x.list))
  {
    stop("Error!!: both X and t.x.list must be lists and they have the same numbers of  components!")
  }
  dim.1 = sapply(1:length(X), function(k) {
    dim(X[[k]])[1]
  })
  if ((length(unique(dim.1)) != 1))
  {
    stop("Error!!: all components of X must be matrix and have the same numbers of rows!")
  }
  if ((dim(X[[1]])[1] != dim(Y)[1]))
  {
    stop(
      "Error!!: the number of observations of X (that is, the number of rows of each component of X) must be equal to the number of observations of Y (that is, the number of rows of Y)!"
    )
  }
  if (sum(sapply(1:length(X), function(k) {
    dim(X[[k]])[2] != length(t.x.list[[k]])
  })) != 0)
  {
    stop(
      "Error!!: The number of columns of each component of X must be equal to the length of the corresponsing component of  t.x.list!"
    )
  }
  if(is.vector(Y))
  {
    if((dim(X[[1]])[1]!=length(Y)))
    {
      stop("Error!!: the number of observations of X (that is, the number of rows of each component of X) must be equal to the number of observations of Y (that is, the number of rows of Y)!")
    }
  }else
  {
    if(is.matrix(Y))
    {
      if(dim(X[[1]])[1]!=dim(Y)[1])
      {
        stop("Error!!: the number of observations of X (that is, the number of rows of each component of X) must be equal to the number of observations of Y (that is, the number of rows of Y)!")
      }
    }else
    {
      stop("Error!!: Y must be a vector or matrix!")
    }
  }

  all.folds <- cv.folds(dim(Y)[1], K.cv)
  is_Y_vector=0
  if(is.vector(Y))
  {
    is_Y_vector=1
    Y=cbind(Y,Y)
    upper.comp=1
  }else
  {
    if(dim(Y)[2]==1)
    {
      is_Y_vector=1
      Y=cbind(Y,Y)
      upper.comp=1
    }else
    {
      upper.comp=min(upper.comp, ncol(Y))
    }
  }

  n.curves=length(X)
  n.sample=dim(Y)[1]
  x.smooth.params=list()

  x.smooth.params[[1]]=n.curves
  lambda.set=c(1e-9, 1e-6,  1e-3,  1,1e3)
  x.smooth.params[[2]]=lambda.set
  x.smooth.params[[3]]=n.sample
  basis.obj.x =create.bspline.basis(c(0,1), nbasis)

  tmp.1=sapply(1:length(t.x.list), function(i){length(t.x.list[[i]])})
  x.smooth.params[[4]]= lapply(1:length(tmp.1), function(i){t(eval.basis(seq(0,1,length=tmp.1[i]),    basis.obj.x))})

  tmp=list()
  tmp[[1]]=getbasispenalty(basis.obj.x,0)
  tmp[[2]]=getbasispenalty(basis.obj.x,2)
  x.smooth.params[[5]]= tmp

  tau.set=c(1e-7,1e-5, 1e-3, 1e-1, 1e1) #interaction code
  x.smooth.params[[6]]= tau.set
  x.smooth.params[[7]]=t.x.list
  x.smooth.params[[8]]=basis.obj.x
  fit.cv=C_cv_mof(t.x.list,   X, Y, x.smooth.params,   all.folds, upper.comp, thresh)


  return(list(fitted_model=fit.cv, is_Y_vector=is_Y_vector,  Y=Y, x.smooth.params=x.smooth.params))
}


######################################
#' @export
pred.msof <- function(fit.obj,  X.test)
{
  if(!is.list(X.test))
  {stop("Error!!: X.test must be a list!")}

  fit.cv=fit.obj$fitted_model
  x.smooth.params=fit.obj$x.smooth.params
  t.x.list=x.smooth.params[[7]]
  n.curves=length(t.x.list)
  if (length(X.test)!=n.curves)
  {stop("Error!!: X.test must be lists with the same number of components as the X.train!")
  }
  if (sum(sapply(1:length(X.test),function(k){!is.matrix(X.test[[k]])})))
  {stop("Error!!: X.test must be a list and all its components must be matrix!")
  }
  if(sum(sapply(1:length(X.test), function(k){ncol(X.test[[k]])!=length(t.x.list[[k]])}))!=0)
  {stop("Error!!: The number of columns of each component of X.test must be equal to the length of the corresponsing component of X.train!")
  }

  B.vals.list=x.smooth.params[[4]]



  mu= fit.cv$mu.expan.coef

  Y.pred=rep(1, nrow(X.test[[1]]))%*%t(mu)
  beta.expan.coef <- fit.cv$beta.expan.coef  #a list of Ls*Lt matrices
  beta.coefs <- list()
  for(i in 1:n.curves)
  {
    beta.coefs[[i]] <- crossprod(B.vals.list[[i]], beta.expan.coef[[i]])
    Y.pred=Y.pred+X.test[[i]]%*%beta.coefs[[i]]/ncol(X.test[[i]])
  }
  if(fit.obj$is_Y_vector==1)
  {
    Y.pred=Y.pred[,1]
  }
  return(Y.pred)
}

 
