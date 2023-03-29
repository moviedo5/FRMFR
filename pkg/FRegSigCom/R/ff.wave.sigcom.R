#' @import fda
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm
#' @useDynLib FRegSigCom
#' @name FRegSigCom

################################################################

determineMaxKcomp <- function(X, Y, params.set, upp.comp, thresh)
{
  repeat.solution=1
  tol=1e-12
  tol.1=1e-8
  p=dim(params.set)[1]
  nobs=dim(X)[1]
  nvar <- dim(X)[2]
  nvarY <- dim(Y)[2]
  max.m=min(upp.comp,nvarY)
  K.comp <- rep(max.m,p)
  x <- X
  y <- Y
  nobs <- dim(x)[1]
  muy <- drop(rep(1,nobs) %*% y) / nobs
  z.y <- scale(y,muy,scale=FALSE)
  mux <- drop(rep(1,nobs) %*% x) / nobs
  z.x <- scale(x,mux,scale=FALSE)
  
  Lambda.mtx <- cbind(diag(nobs),-z.x)
  
  
  Psi.0=t(find_orth_basis(Lambda.mtx))
  
  sigma.xy <- crossprod(z.x, z.y)/sqrt(nobs)
  sigma.yx <- t(sigma.xy)
  sigma.xx <- crossprod(z.x, z.x)
  
  obj.vals <- matrix(0,p,max.m)
  for(j in 1:p)
  {
    tau=params.set[j,1]
    lambda=params.set[j,2]
    Beta=NULL
    Psi=Psi.0
    for(ncomp in 1:max.m)
    {
      if(!is.null(Beta))
      {
        temp <- rbind(matrix(0,nobs,1),sigma.xx %*% beta)
        tt <- temp-Psi%*%(t(Psi)%*%temp)
        Psi <- cbind(Psi, tt/sqrt(sum(tt^2)))
      }
      qq=rep(0,repeat.solution)
      mqq=0
      for(jj in 1:repeat.solution)
      {
        gamma=rnorm(nobs+nvar,0,1)
        gamma.old=rep(0, nobs+nvar)
        count=0
        value=sum((sigma.yx%*%gamma[-(1:nobs)])^2)
        value.old=0
        while((value>value.old)&(min(t(gamma-gamma.old)%*%(gamma-gamma.old), t(gamma+gamma.old) %*% (gamma+gamma.old)) > sum(gamma^2)*tol.1))
        { if(count>0){value.old=value}
          count=count+1
          temp=sigma.xy %*% (sigma.yx %*%gamma[-(1:nobs)])
          t0=-crossprod(Psi, gamma)
          temp.1=temp/sqrt(sum(temp^2))
          ret=max_H(c(rep(0,nobs),temp.1), Psi, t0, nvar, tau, lambda)
          gamma.old=gamma
          gamma=ret$x
          value=sum((sigma.yx%*%gamma[-(1:nobs)])^2)
        }
        qq[jj]=value
        if(qq[jj]>mqq)
        {
          mqq=qq[jj]
          beta=gamma[-(1:nobs)]
        }
      }
      obj.vals[j,ncomp] <- mqq
      Beta=cbind(Beta,beta)
      if(mqq/sum(obj.vals[j,])<thresh){
        K.comp[j] <- ncomp
        break;
      }
    }#end of for(ncomp in 1:(nvarY-1))
  }
  
  return(list( obj.vals=obj.vals, max.K=K.comp))
}
################################################################



getAllComponents <- function(x, y, K.comp, tau, lambda, tol=10^(-12), repeat.solution=1, tol.1=10^{-8})
{
  nvar <- dim(x)[2]
  nobs <- dim(x)[1]
  muy <- drop(rep(1,nobs) %*% y) / nobs
  z.y <- scale(y,muy,scale=FALSE)
  mux <- drop(rep(1,nobs) %*% x) / nobs
  z.x <- scale(x,mux,scale=FALSE)
  Lambda.mtx <- cbind(diag(nobs),-z.x)
  Psi.0=t(find_orth_basis(Lambda.mtx))
  
  sigma.xy <- t(z.x) %*% z.y/sqrt(nobs)
  sigma.yx <- t(sigma.xy)
  sigma.xx <- t(z.x) %*% z.x
  
  Beta=NULL
  Psi=Psi.0
  for(ncomp in 1:K.comp)
  {
    if(!is.null(Beta))
    {
      temp <- rbind(matrix(0,nobs,1),sigma.xx %*% beta)
      tt <- temp-Psi%*%(t(Psi)%*%temp)
      Psi <- cbind(Psi, tt/sqrt(sum(tt^2)))
    }
    
    qq=rep(0,repeat.solution)
    mqq=0
    for(jj in 1:repeat.solution)
    {
      pt <- proc.time()
      gamma=rnorm(nobs+nvar,0,1)
      gamma.old=rep(0, nobs+nvar)
      count=0
      value=sum((sigma.yx%*%gamma[-(1:nobs)])^2)
      value.old=0
      
      while((value>value.old)&(min(t(gamma-gamma.old)%*%(gamma-gamma.old), t(gamma+gamma.old) %*% (gamma+gamma.old)) > sum(gamma^2)*tol.1))
      { 
        if(count>0)
          value.old=value
        count=count+1
        temp=sigma.xy %*% (sigma.yx %*%gamma[-(1:nobs)])
        t0=-t(Psi)%*%gamma
        temp.1=temp/sqrt(sum(temp^2))
        ret=max_H(c(rep(0,nobs),temp.1), Psi, t0, nvar, tau, lambda)
        gamma.old=gamma
        gamma=ret$x
        value=sum((sigma.yx%*%gamma[-(1:nobs)])^2)
      }
      qq[jj]=value
      if(qq[jj]>mqq)
      {mqq=qq[jj]
      beta=gamma[-(1:nobs)]
      }
    }
    Beta=cbind(Beta,beta)
  }
  return(list(z.x=z.x, z.y=z.y, mux=mux, muy=muy, Beta=Beta))
}



################################################################
getErrors.smooth <- function(t.y, Y, Y.test, T, T.test, J,K, smooth.params)
{
  
  q=length(smooth.params)
  error= matrix(0,q, dim(T)[2])
  t.mtx=cbind(rep(1,dim(T)[1])/sqrt(dim(T)[1]),T)
  t.test.mtx=cbind(rep(1,dim(T.test)[1])/sqrt(dim(T)[1]),T.test)
  c.w.0= J%*%t(Y)%*%t.mtx
  for(i in 1:q)
  {
    lambda=smooth.params[i]
    c.w= solve(J%*%t(J)+lambda*K)%*%c.w.0
    for(ncomp in 1:dim(T)[2])
    {
      Y.pred=t.test.mtx[,1:(1+ncomp)]%*%t(c.w)[1:(1+ncomp),]%*%J
      
      error[i,ncomp]=sum((Y.pred-Y.test)^2)
    }
  }
  return(error)
}

################################################################


cv.folds=function (n, folds = 10)
{
  split(sample(1:n), rep(1:folds, length = n))
}
################################################################

######################################
#' @export

cv.fof.wv <- function(X, Y, t.y, K.cv=5, upp.comp=10, thresh=0.01)
{
  if(!is.matrix(X))
  {stop("Error!!: X must be a matrix!")}
  if(!is.matrix(Y))
  {stop("Error!!: Y must be a matrix!")}
  if(!is.vector(t.y))
  {stop("Error!!: t.y must be a vector!")}
  
  if((dim(X)[1]!=dim(Y)[1]))
  {stop("Error!!: the number of observations of X (that is, the number of rows of each component of X) must be equal to the number of observations of Y (that is, the number of rows of Y)!")
  }

  if(dim(Y)[2]!=length(t.y))
  {stop("Error!!: the number of columns of Y must be equal to the length of the vector t.y of the observation points!")
  }
  
  repeat.solution=1
  tol=10^(-12)
  tol.1=10^{-8}
  params.set <- cbind(c(1e-4,1e-3, 0.001,0.01, 0.1  ,1, 10), 
                      c(1e-4,1e-3, 0.001,0.01, 0.1, 0.2, 0.3))  
  t.y=seq(0, 1, length =ncol(Y))
  y.smooth.params=list()
  y.smooth.params[[1]] = ncol(Y)
  y.smooth.params[[2]]=create.bspline.basis(c(0,1), ncol(Y))
  y.smooth.params[[3]]= c(1e-12, 1e-10,1e-8,1e-6,1e-4,1e-2) 
  
  W.basis <- y.smooth.params[[2]]
  J=t(eval.basis(t.y,W.basis))
  K=getbasispenalty(W.basis, 2)
  p=dim(params.set)[1]
  
  scale.x=max(abs(X))
  X=X/scale.x
  t0 <- proc.time()[3]
  max.K <- determineMaxKcomp(X, Y, params.set, upp.comp, thresh)$max.K
  
  print("The maximum components for all parameters")
  print(max.K)
  errors=list()
  for(j in 1:p)
    errors[[j]]= rep(0,max.K[j])
  nvar <- dim(X)[2]
  all.folds <- cv.folds(dim(Y)[1],K.cv)
  
  for(i in seq(K.cv))
  {   print(c("The CV fold ", i))
    omit <- all.folds[[i]]
    x <- X[-omit,]
    y <- Y[-omit,]
    x.valid <- X[omit,]
    y.valid <- Y[omit,]
    nobs <- dim(x)[1]
    muy <- drop(rep(1,nobs) %*% y) / nobs
    z.y <- scale(y,muy,scale=FALSE)
    z.y.valid <- scale(y.valid,muy,scale=FALSE)
    mux <- drop(rep(1,nobs) %*% x) / nobs
    z.x <- scale(x,mux,scale=FALSE)
    z.x.valid <- scale(x.valid,mux,scale=FALSE)
    
    Lambda.mtx <- cbind(diag(nobs),-z.x)
    Psi.0=t(find_orth_basis(Lambda.mtx))
    sigma.xy <- crossprod(z.x, z.y)/sqrt(nobs)
    sigma.yx <- t(sigma.xy)
    sigma.xx <- crossprod(z.x, z.x)
    
    for(j in 1:p)
    {
      tau=params.set[j,1]
      lambda=params.set[j,2]
      K.comp <- max.K[j]
      Beta=NULL
      Psi=Psi.0
      for(ncomp in 1:K.comp)
      {
        if(!is.null(Beta))
        {
          temp <- rbind(matrix(0,nobs,1),sigma.xx %*% beta)
          tt <- temp-Psi%*%(t(Psi)%*%temp)
          Psi <- cbind(Psi, tt/sqrt(sum(tt^2)))
        }
        
        qq=rep(0,repeat.solution)
        mqq=0
        for(jj in 1:repeat.solution)
        {
          gamma=rnorm(nobs+nvar,0,1)
          gamma.old=rep(0, nobs+nvar)
          count=0
          value=sum((sigma.yx%*%gamma[-(1:nobs)])^2)
          value.old=0
          while((value>value.old)&(min(t(gamma-gamma.old)%*%(gamma-gamma.old), t(gamma+gamma.old) %*% (gamma+gamma.old)) > sum(gamma^2)*tol.1))
          { if(count>0){value.old=value}
            count=count+1
            temp=sigma.xy %*% (sigma.yx %*%gamma[-(1:nobs)])
            t0=-t(Psi)%*%gamma
            temp.1=temp/sqrt(sum(temp^2))
            
            ret=max_H(c(rep(0,nobs),temp.1), Psi, t0, nvar, tau, lambda)
            gamma.old=gamma
            gamma=ret$x
            value=sum((sigma.yx%*%gamma[-(1:nobs)])^2)
          }
          qq[jj]=value
          if(qq[jj]>mqq)
          {
            mqq=qq[jj]
            beta=gamma[-(1:nobs)]
          }
        }
        Beta=cbind(Beta,beta)
      }#end of for(ncomp in 1:(K.comp-1))
      
      Beta <- matrix(Beta, nrow=nvar, ncol=K.comp)
      T=z.x%*%Beta
      tmp=sqrt(diag(t(T)%*%T))
      Beta=scale(Beta,center=FALSE,scale=tmp)
      T=z.x%*%Beta
      T.valid=z.x.valid%*%Beta
      errors[[j]]=errors[[j]]+getErrors.smooth(t.y, y, y.valid,  T, T.valid,  J, K, y.smooth.params[[3]])
    }
  }
  
  min.error=rep(0,p)
  for(j in 1:p)
    min.error[j]=min(errors[[j]],na.rm = TRUE)
  temp=which.min(min.error)
  min.tau1 <- params.set[temp,1]
  min.mu1 <- params.set[temp,2]
  opt.K1=which.min(errors[[temp]])
  min.error=rep(0,p)
  for(j in 1:p)
    min.error[j]=min(errors[[j]],na.rm = TRUE)
  temp=(which.min(min.error))[1]
  
  min.tau1 <- params.set[temp,1]
  min.mu1 <- params.set[temp,2]
  error=errors[[temp]]
  a=which(error==min(error),arr.ind=TRUE)
  opt.K1=a[1,2]
  opt.smooth=y.smooth.params[[3]][a[1,1]]
  
  return(list(min.error=min(min.error,na.rm = TRUE), scale.x=scale.x,  X=X,Y=Y, params.set=params.set, error1=errors, opt.K=opt.K1, opt.smooth=opt.smooth, min.tau=min.tau1, min.lambda=min.mu1, J=J, K=K))
  
}

################################################################



######################################
#' @export
pred.fof.wv<- function(cv.obj, X.test)
{
  J=cv.obj$J
  K=cv.obj$K
  X.test=X.test/cv.obj$scale.x
  alltrain10.ret1 <- getAllComponents(cv.obj$X, cv.obj$Y, cv.obj$opt.K, cv.obj$min.tau, cv.obj$min.lambda)
  z.x.test <- scale(X.test, alltrain10.ret1$mux, scale=FALSE)
  z.x=alltrain10.ret1$z.x
  Y=cv.obj$Y
  Beta=alltrain10.ret1$Beta

  T=z.x%*%Beta

  tmp=sqrt(diag(t(T)%*%T))
  Beta=scale(Beta,center=FALSE,scale=tmp)
  T=z.x%*%Beta
  lambda= cv.obj$opt.smooth
  T.test=z.x.test%*%Beta
  t.mtx=cbind(rep(1,dim(T)[1])/sqrt(dim(T)[1]),T)
  t.test.mtx=cbind(rep(1,dim(T.test)[1])/sqrt(dim(T)[1]),T.test)
  c.w.0= J%*%t(Y)%*%t.mtx
  c.w= solve(J%*%t(J)+lambda*K)%*%c.w.0
  Y.pred=t.test.mtx%*%t(c.w)%*%J

  return(Y.pred)
}


