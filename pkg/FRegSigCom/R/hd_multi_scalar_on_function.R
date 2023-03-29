#' @import fda
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom stats cov
#' @useDynLib FRegSigCom
#' @name FRegSigCom


#################################################################
#' @export
cv.msof.hd <- function(X, Y, t.x.list, n.basis=25,  K.cv=5, upper.comp=10, thresh=0.02)
{
  if(!is.list(X))
  {stop("Error!!: X must be a list!")}
  if (sum(sapply(1:length(X),function(k){!is.matrix(X[[k]])})))
  {stop("Error!!: X must be a list and all its components must be matrix!")
  }
  if(!is.list(t.x.list))
  {
    stop("Error!!: t.x.list must be a list!")}
  if (length(X)!=length(t.x.list))
  {
    stop("Error!!: both X and t.x.list must be lists and they have the same numbers of  components!")
  }
  dim.1=sapply(1:length(X),function(k){dim(X[[k]])[1]})
  if((length(unique(dim.1))!=1))
  {
    stop("Error!!: all components of X must be matrix and have the same numbers of rows!")
  }
  if(is.vector(Y))
  {
    if((dim(X[[1]])[1]!=length(Y)))
    {stop("Error!!: the number of observations of X (that is, the number of rows of each component of X) must be equal to the number of observations of Y (that is, the number of rows of Y)!")
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

  if(sum(sapply(1:length(X), function(k){dim(X[[k]])[2]!=length(t.x.list[[k]])}))!=0)
  {stop("Error!!: The number of columns of each component of X must be equal to the length of the corresponsing component of t.x.list!")
  }

  all.folds <- split(sample(nrow(X[[1]])),rep(1:K.cv,length=nrow(X[[1]])))

  tol=1e-11
  tol1=1e-11
  x.params=list()
  x.basis<- create.bspline.basis(c(0,1), n.basis, 4)

  tmp=lapply(1:length(t.x.list), function(k){seq(0,1,length.out=length(t.x.list[[k]]))})
  t.x.list=tmp
  x.params[[1]]=t.x.list
  eta.set=c(1e-6, 1e-3, 1)
  x.params[[2]]=eta.set
  J.a <- getbasispenalty(x.basis,0)
  J2.a <- getbasispenalty(x.basis, 2)
  nvarX=length(X)
  x.params[[3]]=J.a
  x.params[[4]]=J2.a
  x.params[[5]]=x.basis$nbasis
  x.params[[6]]=nvarX
  x.params[[7]]=cbind(c(0.1, 1,10, 100), c(0.1,0.2,0.3, 0.4))

  nbasis=x.basis$nbasis
  temp.fun=function(k)
  {
    return(eval.basis(t.x.list[[k]], x.basis)/length(t.x.list[[k]]))
  }
  wb<-lapply(1:nvarX, function(k){temp.fun(k)})
  if(is.null(all.folds))
  {
    all.folds <- cv.folds(dim(X[[1]])[1], K.cv)
  }
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
      upper.comp=min(upper.comp, dim(Y)[2])
    }
  }

  Y.scale=sqrt(sum(diag(cov(Y)))/dim(Y)[2])
  #Y.scale=sqrt(max(diag(cov(Y))))

    Y=Y/Y.scale

  fit_c=cv_hd_msof(X, Y, wb, x.params, all.folds, upper.comp, thresh)

  return(list(opt.K=fit_c$opt.K, opt.tau=fit_c$opt.tau, opt.lambda=fit_c$opt.lambda,
              opt.eta=fit_c$opt.eta, min.error=fit_c$min.error, errors=fit_c$errors, mu.x=fit_c$mu.x,Beta=fit_c$Beta,
              is_Y_vector=is_Y_vector, X.scale=fit_c$X.scale, Y.scale=Y.scale,
               XbTransInv=fit_c$XbTransInv, normTransInv=fit_c$normTransInv, nbasis=nbasis, eta.set=eta.set,
              wb=wb,nvarX=nvarX, nabsis=fit_c$nabsis, Y=Y))
}



##############################################################

#' @export
pred.msof.hd <- function(fit.cv, X.test)
{
  X.scale=fit.cv$X.scale
  XbTransInv=fit.cv$XbTransInv
  Beta=(fit.cv$Beta)[,1:(fit.cv$opt.K)]
  if(is.vector(Beta))
  {
      Beta=matrix(Beta, ncol=1)
  }
  normTransInv=fit.cv$normTransInv

  nvarX=length(X.test)
  if(is.vector(X.test[[1]]))
  {
    tmp=lapply(1:nvarX, function(k){rbind(X.test[[k]], X.test[[k]])})
    X.test=tmp
    nobs=2
  }else
  {
    nobs=dim(X.test[[1]])[1]
  }
  nbasis=fit.cv$nbasis
  wb=fit.cv$wb

  Y.pred=get_pred_msof(X.test, X.scale, XbTransInv, Beta, normTransInv, fit.cv$mu.x, fit.cv$wb, fit.cv$Y, nobs, nbasis, nvarX)
  if(is.vector(X.test[[1]]))
  {
    Y.pred=as.vector(Y.pred[1,])
  }
  Y.pred=Y.pred*fit.cv$Y.scale
  if(fit.cv$is_Y_vector==1)
  {
    Y.pred=Y.pred[,1]
  }

  return(Y.pred)
}

