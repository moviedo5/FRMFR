#' @import fda
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @useDynLib FRegSigCom
#' @importFrom stats sd
#' @name FRegSigCom



cv.folds <- function(n,nfolds=5)
  ## Randomly split the n samples into folds
  ## Returns a list of nfolds lists of indices, each corresponding to a fold
{
  #return(split(sample(n),rep(1:nfolds,length=n)))
  return(split(1:n,rep(1:nfolds,length=n)))
}


#########################################

###################################################################


#######################################################################
#' @export
cv.sigcom=function(X, Y, t.x, t.y, Z=NULL, s.n.basis=50, t.n.basis=50, K.cv=5, upper.comp=20, thresh=0.001, basis.type.x="Bspline", basis.type.y="Bspline")
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
  
  all.folds <- cv.folds(dim(Y)[1], K.cv)
  
  X.org=X
  Y.org=Y
  tmp=scale(Y.org, scale=FALSE)
  Y.scale=max(abs(tmp))
  Y=Y.org/Y.scale
  
  X.scale=list()
  for(i in 1:length(X))
  {
    tmp=scale(X.org[[i]], scale=FALSE)
    X.scale[[i]]=max(abs(X[[i]]))
    X[[i]]=X.org[[i]]/X.scale[[i]]
  }
  
  n.curves=length(X)
  n.sample=dim(Y)[1]
  x.smooth.params=list()
  
  x.smooth.params[[1]]=n.curves
  lambda.set=c(1e-8,1e-6, 1e-4, 1e-2, 1, 1e2) #tuning parameters in interaction code
  x.smooth.params[[2]]=lambda.set
  x.smooth.params[[3]]=n.sample
  if(basis.type.x=="Bspline"){
    basis.obj.x =create.bspline.basis(c(0,1), s.n.basis)
  }else{
    if(basis.type.x=="Fourier"){
      if(s.n.basis%%2==0)
      {
        s.n.basis=s.n.basis+1
        print("In **create.fourier.basis(c(0, 1), s.n.basis)** s.n.basis must be an odd integer; since s.n.basis is even now, it will be increased by 1")
      }
      basis.obj.x =create.fourier.basis(c(0,1), s.n.basis)
    }
  }
  tmp.1=sapply(1:length(t.x), function(i){length(t.x[[i]])})
  x.smooth.params[[4]]= lapply(1:length(tmp.1), function(i){t(eval.basis(seq(0,1,length=tmp.1[i]),    basis.obj.x))})
  
  tmp=list()
  tmp[[1]]=getbasispenalty(basis.obj.x,0)
  tmp[[2]]=getbasispenalty(basis.obj.x,2)
  x.smooth.params[[5]]= tmp
  
  tau.set=c(1e-3,1e-1,1e1, 1e3) #interaction code
  x.smooth.params[[6]]= tau.set
  x.smooth.params[[7]]=t.x
  x.smooth.params[[8]]=basis.obj.x
  
  
  y.smooth.params = list()
  if(basis.type.y=="Bspline"){
    y.smooth.params[[1]] = create.bspline.basis(c(0, 1), t.n.basis)
  }else{
    y.smooth.params[[1]] = create.fourier.basis(c(0, 1), t.n.basis)
  }
  y.smooth.params[[2]] = c(1e-11, 1e-9,1e-7, 1e-5,1e-3,1e-1,1e1, 1e3)
  y.smooth.params[[3]] =  length(t.y)
  y.smooth.params[[4]] = t(eval.basis(seq(0, 1, length = length(t.y)), y.smooth.params[[1]]))
  tmp = getbasispenalty(y.smooth.params[[1]], 2)
  y.smooth.params[[5]] =  (tmp + t(tmp)) / 2
  B.vals = y.smooth.params[[4]]
  K.w = y.smooth.params[[5]]
  y.weights.aver = 1 / y.smooth.params[[3]]
  
  B.vals.weig = B.vals * y.weights.aver
  y.penalty.inv = list()
  kappa.set = y.smooth.params[[2]]
  tmp=list()
  tmp[[1]]=B.vals.weig %*% t(B.vals)
  tmp[[2]]=K.w *y.weights.aver
  
  y.smooth.params[[6]] = B.vals.weig
  y.smooth.params[[7]] = tmp
  y.smooth.params[[8]]=t.y
  
  if(is.null(Z))
  {
    is.null.Z=TRUE
    fit.cv=C_cv_ff_sig(t.x,   X, Y, x.smooth.params, y.smooth.params, all.folds, upper.comp, thresh)
    fit.cv[[1]]=Y.scale*fit.cv[[1]]
    for(i in 1:length(X))
    {
      fit.cv[[2]][[i]]=Y.scale*fit.cv[[2]][[i]]/X.scale[[i]]
    }
  }else
  {
    is.null.Z=FALSE
    if(is.vector(Z))
    {
      if(length(Z)!=n.sample)
      {
        stop("Error!!: when Z is a vector, the size of Z must be equal to the sample size!")
      }else
      {
        Z=matrix(Z,ncol=1)
      }
    }else
    {
      if(is.matrix(Z))
      {
        if(nrow(Z)!=n.sample)
        {
          stop("Error!!: when Z is a matrix, the number of rows of Z must be equal to the sample size!")
        }
      }else
      {
        stop("Error!!: Z must be a matrix or a vector!")
      }
    }
    delta.set=c(1e-8, 1e-5, 1e-2, 1e1) #interaction code
    x.smooth.params[[9]]=delta.set
    Z.org=Z
    Z.scale=apply(Z,2,sd)
    for(i in 1:ncol(Z))
    {
      Z[,i]=Z[,i]/Z.scale[i]
    }
    fit.cv=C_cv_ff_sig_with_scalar(t.x, X, Y, Z, x.smooth.params, y.smooth.params, all.folds, upper.comp, thresh)
    fit.cv[[1]]=Y.scale*fit.cv[[1]]
    for(i in 1:length(X))
    {
      fit.cv[[2]][[i]]=Y.scale*fit.cv[[2]][[i]]/X.scale[[i]]
    }
    for(i in 1:ncol(Z))
    {
      fit.cv[[3]][,i]=Y.scale*fit.cv[[3]][,i]/Z.scale[i]
    }
  }
  
  
  return(list(fitted_model=fit.cv, is.null.Z=is.null.Z, Z=Z, Y=Y, x.smooth.params=x.smooth.params, y.smooth.params=y.smooth.params))
}


######################################
#' @export
pred.sigcom <- function(fit.obj,  X.test, t.y.test=NULL, Z.test=NULL)
{
  if(!is.list(X.test))
  {stop("Error!!: X.test must be a list!")}
  
  fit.cv=fit.obj$fitted_model
  x.smooth.params=fit.obj$x.smooth.params
  t.x=x.smooth.params[[7]]
  n.curves=length(t.x)
  if (length(X.test)!=n.curves)
  {stop("Error!!: X.test must be lists with the same number of components as the X.train!")
  }
  if (sum(sapply(1:length(X.test),function(k){!is.matrix(X.test[[k]])})))
  {stop("Error!!: X.test must be a list and all its components must be matrix!")
  }
  if(sum(sapply(1:length(X.test), function(k){ncol(X.test[[k]])!=length(t.x[[k]])}))!=0)
  {stop("Error!!: The number of columns of each component of X.test must be equal to the length of the corresponsing component of X.train!")
  }
  y.smooth.params=fit.obj$y.smooth.params
  t.y <- y.smooth.params[[8]]
  B.vals.list=x.smooth.params[[4]]
  
  if(is.null(t.y.test))
  {
    B.vals.y=y.smooth.params[[4]]
  }else
  {
    if(!is.vector(t.y.test))
    {
      stop("Error!!: t.y.test must be a vector!")
    }else
    {
      if((min(t.y.test)<min(t.y))|(max(t.y.test)>max(t.y)))
      {
        stop("Error!!: the range of t.y.test must be within the range of the time points ofresponse curve in training data!")
      }else
      {
        tmp=seq(0,1, length.out=length(t.y))
        t.tmp=approx(t.y, tmp, xout=t.y.test)$y
        B.vals.y=t(eval.basis(t.tmp, y.smooth.params[[1]]))
      }
    }
  }
  is.null.Z=fit.obj$is.null.Z
  Z=fit.obj$Z
  if(!is.null.Z)
  {
    if(is.null(Z.test))
    {
      stop("Error!!: do not find Z.test, the test data for the scalar predictors!")
    }else
    {
      if(!is.matrix(Z.test))
      {
        if(is.vector(Z.test))
        {
          Z.test=matrix(Z.test,ncol=1)
        }else
        {
          stop("Error!!: Z.test must be a matrix or a vector!")
        }
      }
      if(ncol(Z)!=ncol(Z.test))
      {
        stop("Error!!: the number of variables in Z.test must be the same as that of Z!")
      }else
      {
        if(nrow(X.test[[1]])!=nrow(Z.test))
        {
          stop("Error!!: the number of rows of Z.test must be the same as that of each matrix in X.test!")
        }
      }
    }
  }
  
  mu.expan.coef= fit.cv$mu.expan.coef
  mu= crossprod(B.vals.y, mu.expan.coef)
  Y.pred=rep(1, nrow(X.test[[1]]))%*%t(mu)
  beta.expan.coef <- fit.cv$beta.expan.coef  #a list of Ls*Lt matrices
  beta.coefs <- list()
  for(i in 1:n.curves)
  {
    tmp <- tcrossprod(beta.expan.coef[[i]], t(B.vals.y))  #Ls*T'
    beta.coefs[[i]] <- crossprod(B.vals.list[[i]], tmp)
    Y.pred=Y.pred+X.test[[i]]%*%beta.coefs[[i]]/ncol(X.test[[i]])
  }
  if(!is.null.Z)
  {
    gamma.expan.coef= fit.cv$gamma.expan.coef
    gamma= t(crossprod(B.vals.y, gamma.expan.coef))
    Y.pred=Y.pred+Z.test%*%gamma
  }
  return(Y.pred)
}


######################################
#' @export
getcoef.sigcom <- function(fit.obj, t.x.coef=NULL, t.y.coef=NULL)
{
  fit.cv=fit.obj$fitted_model
  x.smooth.params=fit.obj$x.smooth.params
  t.x=x.smooth.params[[7]]
  n.curves=length(t.x)
  y.smooth.params=fit.obj$y.smooth.params
  t.y <- y.smooth.params[[8]]
  
  
  if(is.null(t.x.coef))
  {
    B.vals.list=x.smooth.params[[4]]
  }else
  {
    B.vals.list=list()
    if(!is.list(t.x.coef))
    {
      stop("Error!!: t.x.coef must be a list!")
    }else
    {
      if(length(t.x.coef)!=n.curves)
      {
        stop("Error!!: the number of components ofthe list t.x.coef must equal to the number of predictor curves!")
      }
      for(i in 1:n.curves)
      {
        if((min(t.x.coef[[i]])<min(t.x[[i]]))|(max(t.x.coef[[i]])>max(t.x[[i]])))
        {
          stop("Error!!: the range each component of t.x.coef must be within the range each component of t.x!")
        }else
        {
          tmp=seq(0,1, length.out=length(t.x.coef[[i]]))
          t.tmp=approx(t.x[[i]], tmp, xout=t.x.coef[[i]])$y
          B.vals.list[[i]]=t(eval.basis(t.tmp, x.smooth.params[[8]]))
        }
      }
    }
  }
  
  if(is.null(t.y.coef))
  {
    B.vals.y=y.smooth.params[[4]]
  }else
  {
    if(!is.vector(t.y.coef))
    {
      stop("Error!!: t.y.coef must be a vector!")
    }else
    {
      if((min(t.y.coef)<min(t.y))|(max(t.y.coef)>max(t.y)))
      {
        stop("Error!!: the range of t.y.coef must be within the range of the time points ofresponse curve in training data!")
      }else
      {
        tmp=seq(0,1, length.out=length(t.y))
        t.tmp=approx(t.y, tmp, xout=t.y.coef)$y
        B.vals.y=t(eval.basis(t.tmp, y.smooth.params[[1]]))
      }
    }
  }
  
  is.null.Z=fit.obj$is.null.Z
  
  
  mu.expan.coef= fit.cv$mu.expan.coef
  mu= as.vector(crossprod(B.vals.y, mu.expan.coef))
  
  beta.expan.coef <- fit.cv$beta.expan.coef  #a list of Ls*Lt matrices
  beta.coefs <- list()
  for(i in 1:n.curves)
  {
    tmp <- tcrossprod(beta.expan.coef[[i]], t(B.vals.y))  #Ls*T'
    beta.coefs[[i]] <- crossprod(B.vals.list[[i]], tmp)/(max(t.x[[i]])-min(t.x[[i]]))
  }
  gamma=NULL
  if(!is.null.Z)
  {
    gamma.expan.coef= fit.cv$gamma.expan.coef
    gamma= t(crossprod(B.vals.y, gamma.expan.coef))
  }
  return(list(mu=mu, beta=beta.coefs, gamma=gamma))
}



