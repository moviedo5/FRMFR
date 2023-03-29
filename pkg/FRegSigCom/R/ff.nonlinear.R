#' @import fda
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @useDynLib FRegSigCom
#' @name FRegSigCom
#########################################


cv.folds <- function(n,nfolds=5)
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

eval.tensor.spline.basis=function(t.x, x, bspline.x.obj, bspline.s.obj)
{

  lenghtX=length(t.x)
  x=as.vector(t(x))

  A=t(eval.basis(t.x,bspline.s.obj))
  B=t(eval.basis(x,bspline.x.obj))
  tmp=sapply(1:(length(x)/lenghtX), function(k){c_prod(A, t(B[,(k-1)*lenghtX+1:lenghtX]))/lenghtX})

  return(tmp)
}


#########################################

#########################################
#' @export
cv.nonlinear=function(X, Y, t.x.list, t.y, s.n.basis = 40, x.n.basis=40, t.n.basis = 40, K.cv=5, upper.comp=10, thresh=0.01)
{

  if(!is.list(X))
  {stop("Error!!: X must be a list!")}
  if (sum(sapply(1:length(X),function(k){!is.matrix(X[[k]])})))
  {stop("Error!!: X must be a list and all its components must be matrix!")
  }
  if(!is.list(t.x.list))
  {stop("Error!!: t.x.list must be a list!")}
  if (length(X)!=length(t.x.list))
  {stop("Error!!: both X and t.x.list must be lists and they have the same numbers of  components!")
  }
  dim.1=sapply(1:length(X),function(k){dim(X[[k]])[1]})
  if((length(unique(dim.1))!=1))
  {stop("Error!!: all components of X must be matrix and have the same numbers of rows!")
  }
  if((dim(X[[1]])[1]!=dim(Y)[1]))
  {stop("Error!!: the number of observations of X (that is, the number of rows of each component of X) must be equal to the number of observations of Y (that is, the number of rows of Y)!")
  }
  if(sum(sapply(1:length(X), function(k){dim(X[[k]])[2]!=length(t.x.list[[k]])}))!=0)
  {stop("Error!!: The number of columns of each component of X must be equal to the length of the corresponsing component of  t.x.list!")
  }
  if(dim(Y)[2]!=length(t.y))
  {stop("Error!!: the number of columns of Y must be equal to the length of the vector t.y of the observation points!")
  }


  n.curves=length(X)
  n.sample=dim(Y)[1]
  t.x.list.org=t.x.list
  t.x.list=lapply(1:n.curves, function(k){seq(0,1,length.out=ncol(X[[k]]))})
  shift.x.list=list()
  scale.x.list=list()

  for(i in 1:n.curves)
  {
    shift.x.list[[i]]=min(X[[i]])
    scale.x.list[[i]]=max(X[[i]])-min(X[[i]])
    X[[i]]=(X[[i]]-shift.x.list[[i]])/scale.x.list[[i]]
  }


  bspline.s.obj=create.bspline.basis(c(0,1), s.n.basis, 4)
  bspline.x.obj=create.bspline.basis(c(0, 1), x.n.basis, 4)


  K.x0=getbasispenalty(bspline.x.obj, 0)
  K.x1=getbasispenalty(bspline.x.obj, 1)
  K.x2=getbasispenalty(bspline.x.obj, 2)

  K.s0=getbasispenalty(bspline.s.obj, 0)
  K.s1=getbasispenalty(bspline.s.obj, 1)
  K.s2=getbasispenalty(bspline.s.obj, 2)

  J0=K.x0%x%K.s0
  J0=(J0+t(J0))/2

  J2=K.x2%x%K.s0+K.x1%x%K.s1+K.x0%x%K.s2
  J2=(J2+t(J2))/2

  x.params=list()
  x.params[[1]]=n.curves
  x.params[[2]]=J0
  x.params[[3]]=J2
  x.params[[4]]=nrow(J0)
  lambda.set=c(1e-10,1e-8, 1e-6,1e-4,1e-2,1, 1e2, 1e4)
  x.params[[5]]=lambda.set
  tau.set=c(1e-4,1e-1, 100)

  x.params[[6]]=tau.set
   
  d=split(sample(n.sample),rep(1:K.cv,length=n.sample))
  all.folds=matrix(0, n.sample, K.cv)
  for(i in 1:K.cv)
  {
    all.folds[d[[i]],i]=1
  }
  
  G=NULL
  G.mean.list=list()
  for(i in 1:n.curves)
  {
    tmp=t(eval.tensor.spline.basis(t.x.list[[i]], X[[i]], bspline.x.obj, bspline.s.obj))
    G.mean.list[[i]]=apply(tmp,2,mean)
    tmp=scale(tmp, scale=FALSE)/sqrt(n.sample)
    G=rbind(G, t(tmp))
  }

  y.params = list()
  y.params[[1]] = create.bspline.basis(c(0, 1), t.n.basis)

  y.params[[2]] = c(1e-10,1e-8,1e-6,1e-4,1e-2,1,100)
  y.params[[3]] =  length(t.y)
  y.params[[4]] = t(eval.basis(seq(0, 1, length.out = length(t.y)), y.params[[1]]))
  tmp = getbasispenalty(y.params[[1]], 2)
  y.params[[5]] =  (tmp + t(tmp)) / 2
  B.vals = y.params[[4]]
  K.w = y.params[[5]]
  y.weights.aver = 1 / y.params[[3]]

  B.vals.weig = B.vals * y.weights.aver
  y.penalty.inv = list()
  kappa.set = y.params[[2]]
  tmp=list()
  tmp[[1]]=B.vals.weig %*% t(B.vals)
  tmp[[2]]=K.w *y.weights.aver

  y.params[[6]] = B.vals.weig
  y.params[[7]] = tmp

  fit.1=C_cv_nonlinear_ff(G, Y, x.params, y.params, all.folds, upper.comp, thresh)

  return(list(opt.K=fit.1$opt_K, opt.lambda=fit.1$opt_lambda, opt.tau=fit.1$opt_tau, opt.kappa=fit.1$opt_kappa,
              opt.T=fit.1$opt_T, opt.z=fit.1$opt_z, bspline.x.obj=bspline.x.obj, bspline.s.obj=bspline.s.obj,
              shift.x.list=shift.x.list, scale.x.list=scale.x.list, y.params=y.params, Y=Y, t.x.list=t.x.list,
              t.x.list.org=t.x.list.org, G.mean.list=G.mean.list))
}


######################################
#' @export
pred.nonlinear <- function(fit.cv, X.test, t.y.test=NULL){

  t.x.list=fit.cv$t.x.list
  y.params=fit.cv$y.params
  n.curves=length(X.test)
  bspline.x.obj=fit.cv$bspline.x.obj
  bspline.s.obj=fit.cv$bspline.s.obj

  shift.x.list=fit.cv$shift.x.list
  scale.x.list=fit.cv$scale.x.list
  for(k in 1:n.curves)
  {

    X.test[[k]]=(X.test[[k]]-shift.x.list[[k]])/scale.x.list[[k]]
    X.test[[k]]=((X.test[[k]]+1)-abs(X.test[[k]]-1))/2
    X.test[[k]]=((X.test[[k]])+abs(X.test[[k]]))/2
  }
  Y=fit.cv$Y
  n.sample=dim(Y)[1]

  opt.K=fit.cv$opt.K
  opt.lambda=fit.cv$opt.lambda
  kappa=fit.cv$opt.kappa

  t.y=seq(0,1,length.out=ncol(Y))
  y.basis=y.params[[1]]


  J.w= t(eval.basis(t.y,y.basis))
  K.w=getbasispenalty(y.basis, 2)
  y.int.weights=diag(rep(1/length(t.y),length(t.y)))
  y.weights.aver=mean(diag(y.int.weights))

  T=fit.cv$opt.T
  z=fit.cv$opt.z
  ncol=opt.K
  G.mean.list=fit.cv$G.mean.list

  G.test.list=list()
  for(k in 1:n.curves)
  {
    tmp=t(eval.tensor.spline.basis(t.x.list[[k]], X.test[[k]],  bspline.x.obj, bspline.s.obj))
    G.test.list[[k]]= scale(tmp, center=G.mean.list[[k]], scale=FALSE)/sqrt(n.sample)
  }
  G.test=do.call(cbind,G.test.list)
  T.test=c_prod(G.test, z)

  T=as.matrix(T)
  T.test=as.matrix(T.test)
  tmp.1 <- sqrt(as.numeric(apply(T^2,2,sum)))
  T <- scale(T, center=FALSE, scale=tmp.1)
  T.test <- scale(T.test, center=FALSE, scale=tmp.1)


  t.mtx <- cbind(1/sqrt(dim(T)[1]),T)

  # print(t(t.mtx)%*%t.mtx)
  t.test.mtx <- cbind(1/sqrt(dim(T)[1]),T.test)

  coef.w.0= J.w%*%y.int.weights%*%t(Y)%*%t.mtx
  coef.w <- solve(J.w%*%y.int.weights%*%t(J.w)+kappa*K.w*y.weights.aver)%*%coef.w.0
  if(is.null(t.y.test)){
    Y.pred <- t.test.mtx%*%t(coef.w)%*%J.w
  }else{
    Y.pred <- t.test.mtx%*%t(coef.w)%*%t(eval.basis(t.y.test,y.basis))
  }
  return(Y.pred=Y.pred)
}


##########################################
#' @export
getcoef.nonlinear=function(fit.cv, n.x.grid=50)
{

  t.x.list.org=fit.cv$t.x.list.org
  t.x.list=fit.cv$t.x.list
  n.curves=length(t.x.list.org)
  range.t.x.list=lapply(1:n.curves, function(k){max(t.x.list.org[[k]])-min(t.x.list.org[[k]])})
  X.grid=list()
  for(k in 1:n.curves)
  {
    a=fit.cv$scale.x.list[[k]]
    b=fit.cv$shift.x.list[[k]]
    X.grid[[k]]=seq(b,a+b, length.out=n.x.grid)
  }
  t.y=seq(0,1,length.out=ncol(fit.cv$Y))

  tmp.list=lapply(1:n.curves, function(k){matrix(0, 1,length(t.x.list[[k]]))})
  Y.pred.0=pred.nonlinear(fit.cv, tmp.list)
  mu=as.vector(Y.pred.0)

  F=list()
  X.new=lapply(1:n.curves, function(k){A=diag(length(t.x.list[[k]])); X.grid[[k]]%x%A})
  for(k in 1:n.curves)
  {
    X.0=lapply(1:n.curves, function(j){matrix(0, dim(X.new[[k]])[1], dim(X.new[[j]])[2])})
    X.0[[k]]=X.new[[k]]
    Y.pred=pred.nonlinear(fit.cv, X.0)
    Y.pred=t(sapply(1:dim(Y.pred)[1],function(k){Y.pred[k,]}))
    Y.pred=t(sapply(1:dim(Y.pred)[1],function(k){Y.pred[k,]-Y.pred.0}))
    F[[k]]=array(Y.pred, c(length(t.x.list[[k]]),length(X.grid[[k]]),length(t.y)))
    F[[k]]=aperm(F[[k]], c(2,1,3))*length(t.x.list[[k]])/(range.t.x.list[[k]])
  }
  return(list(mu=mu, F=F, X.grid=X.grid, t.x.list=t.x.list.org, t.y=t.y))
}




