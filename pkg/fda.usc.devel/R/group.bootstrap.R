#' @title Bootstrap for Functional Data by factor.
#'
#' @description  Convex and Local Convex, smooth and naive  bootstrap for functional data (fdata and ldata class) by group label
#' 
#' @param x \code{list} that containing the variables in the model. The first item in the
#'  \code{data} list is called \emph{"df"} and is a data frame with the response and non 
#'  functional explanatory variables, as  \code{\link{glm}}.  Functional covariates of class 
#' @param response \code{character}, name of the factor containing the group labels.
# @param boot \code{character},  type of bootstrap. Values are "smooth" and "convex".
#' @param nb number of total resamples (if vector of length 1), number of resample for each class level  (if vector of length nlevels).
#' @param smo The smoothing parameter for the bootstrap samples as a proportion of the sample variance matrix.
#' @param Nhull \code{character}, size of the convex hull for generating a new sample. If \code{Nhull=0} a naive bootstrap is done.
#' @param Nnbh \code{character}, size of the neighborhood for selecting a local convex hull. Only needed for local convex procedures. \code{Nnbh}\eqn{>=}\code{Nhull}
#' @param weights Vector of weights. Only used when \code{nb} is a vector of length 1, then, the number of resamples per class is weighted. 
#' @param metric metric to be applied to \code{x}. Only needed for local convex bootstrap, when \code{Nnbh>0}.
#' 
#' @return A \code{ldata} class object.  This \code{list}  contains:
#'  A \code{data.frame} object (called \emph{"df"})  with the bootstrap resamples for response and multivariate variables.
#'  The following elements of the list are \code{fdata} class objects with bootstrap resamples for the functional variables.
#'
# @references
#' 
#' @aliases group.bootstrap
#' 
#' @details {
#' Different types of options:
#' \itemize{
#' \item Smooth bootstrap,  \code{boot="smooth"} and \code{smo>0} (\code{Nhull} and \code{Nnbh} are not used).
#' \item Naive  bootstrap,  \code{boot="smooth"} and \code{smo=0} (\code{Nhull} and \code{Nnbh} are not used).
#' \item Global convex bootstrap,  \code{boot="convex"} and \code{Nnbh=NULL}   (\code{smo} is not used).
#' Select additional \eqn{X_J} elements of the same group and generate a random convex combination 
#' \item Local convex bootstrap,  \code{boot="convex"} and \code{Nnbh>0}   (\code{smo} is not used).  The same as previous item but now the additional elements are selected from the neighborhood of Xi in group Gi
#' }
#' }
#' 
# @seealso See also  \code{\link{classif.bootstrap}}
#' @noRd  
#' @author
#' Manuel Febrero-Bande, Manuel Oviedo de la Fuente 
#' @examples
#' \dontrun{
#' data(phoneme)
#' mlearn <- phoneme[["learn"]]
#' glearn <- phoneme[["classlearn"]]
#' dataf <- data.frame(glearn)
#' dat = list("df" = dataf, "x" = mlearn)
#' dat.globalconvex <- group.bootstrap(dat,"glearn", boot = "convex")
#' dat.localconvex <- group.bootstrap(dat,"glearn", boot = "convex",Nnbh=10)
#' dat.smooth <- group.bootstrap(dat,"glearn", boot = "smooth", smo = 0.05)
#' par(mfrow=c(2,2))
#' plot(dat$x,col=dat$df$glearn)
#' plot(dat.globalconvex$x,col=dat.globalconvex$df$glearn)
#' plot(dat.localconvex$x,col=dat.localconvex$df$glearn)
#' plot(dat.smooth$x,col=dat.smooth$df$glearn)
#' # number of resamples by class level
#' dat.naive <- group.bootstrap(dat,"glearn", boot = "smooth", smo = 0, nb = c(5,5,25,5,5))
#' plot(dat.naive,var.name="glearn")
#' }
# @export group.bootstrap
group.bootstrap <-function(x, response #,boot = "convex"
                           , nb = NULL
                           , smo = 0, Nhull = NULL, Nnbh = NULL
                           , weights, metric = NULL){
  #,expand=1
#print(" entra convex bootstrap")  
  # print(boot)
  ndatos <- NROW(x[[1]])
  y <- x$df[,response]
  iresponse <- which(names(x$df)==response)
  nlev <- nlevels(y)
  lev  <- levels(y)
  if (nlev == length(nb)){
    #print(nlev);    print(nb)
    Bgroup <- nb
    nb <- sum(Bgroup)
    group.boot <- rep(lev, times = Bgroup)
  }else{
  if (missing(weights))   {
    # warning("missing weights in convex bootstrap")
    weights <- weights4class(y)
  }
  if (is.null(nb)) 
    nb = length(y)
  if (sum(weights)!=1) 
    weights<- weights/sum(weights)
    group.boot <- sample(y, size = nb, replace = TRUE, prob = weights )
    Bgroup <- table(group.boot)
  }
  tgroup <- table(y)
  xboost<-NULL
  idf<-(names(x)=="df")
  lenl<-length(x)#-1
  #if (is.null(par.local$metric))
#    par.local$metric <- metric.ldata(fbb,method="euclidean") Calcular D al pprio y despues seleccion DD<-D[ij,ij]
  if (is.null(metric))
    is.null.metric<-TRUE  else is.null.metric-FALSE
  for (i in 1:nlev){
#cat("convex 2 group - level ",i,"*************************** \n")
    ij <- y == lev[i]
    faa<- subset.ldata(x,ij)
    fbb <- faa
    if (NCOL(faa$df)==1) {
       fbb<-faa[-idf]
       lenl <- length(fbb)
       idf2 <- FALSE
    }    else    {
      fbb$df <- faa$df[,-iresponse,drop=F]
      idf2 <-TRUE
    } 
    # print("smo, nNhull y N")
    #print(smo);    print(Nhull);    print(Nnbh)
    if ( is.null(Nhull) & is.null(Nnbh)){   
# print(" smo or anive RESAMP")    
      fx <- boot.smooth.ldata(x = fbb, nb = Bgroup[i], smo =  smo)
    } else{
    if (!is.null(Nhull) & is.null(Nnbh)){                 
#print(" GLOBAL CONVEX RESAMPLE")
          fx<-boot.convex.global(fbb, Nhull=Nhull, nb=Bgroup[i])
    }      
    if ( !is.null(Nnbh)){         
#print(" LOCAL CONVEX RESAMPLE")
          #if (is.null(Nnbh)) 
           # Nnbh<- min(tgroup[i],Nhull*3)
         if (is.null(Nhull))            Nhull = min(4, floor(Nnbh/2))
         if (is.null.metric)
            metric <- metric.ldata(fbb, method = "euclidean")
          #fx<-do.call("boot.convex.local",par.local)
          fx = boot.convex.local(fbb, Nhull=Nhull, nb=Bgroup[i], Nnbh=Nnbh, metric=metric)
        }
    }
    if ( i > 1){
      for (j in 1:lenl){
        if (is.fdata(xboost[[j]])) {
          xboost[[j]]$data <- rbind(xboost[[j]]$data,fx[[j]]$data)
        }
        else xboost[[j]] <- rbind(xboost[[j]],fx[[j]])
      }
    }    else xboost <- fx

  }
  
  df <- data.frame(sort(group.boot))
  names(df)<-response
  if (idf2){
    idf2<-which(names(xboost)=="df")
    df<-cbind(df,xboost[[idf2]]) #mejor ordear los datos
    #output <- c(list("df"=df),xboost[-idf2])
    output <- ldata(df,mfdata=xboost[-idf2])
  }
  else{
    #output <- c(list("df"=df),xboost)
    output <- ldata(df,mfdata=xboost)
  }
  return(output)
}
##################################################################
# dat.naive <- group.bootstrap(dat,"glearn", boot = "smooth", smo = 0, nb = c(2,2,2,2,6))
# plot(dat.naive,var.name="glearn",lwd=5)
# dat.naive <- group.bootstrap(dat,"glearn", boot = "smooth", smo = 0, nb = 10)
##################################################################
# # funcion para borrar
# bconvex=function(x,J=4,B=5000){
#   n=nrow(x)
#   xm=colMeans(x)        #* 
#   xc=sweep(x,2,xm,"-")  #Nuevo quitar la media 
#   w=matrix(rnorm(B*(J+1)),nrow=B,ncol=J+1)
#   w=w^2
#   sumsq=apply(w,1,sum)
#   w=sweep(w,1,sumsq,"/")
#   w=sweep(w,1,sqrt(apply(w^2,1,sum)),"/") # Ajuste para varianza similar
#   ix=matrix(sample(1:n,size=(J+1)*B,replace=TRUE),ncol=J+1)
#   xx=matrix(0,ncol=ncol(x),nrow=B)
#   for (j in 1:(J+1)){
#     xx=xx+sweep(xc[ix[,j],],1,w[,j], "*") #Con datos centrados
#   }
#   xx=sweep(xx,2,xm,"+")  #Reponer la media #*
#   return(xx)
# }
##################################################################

##################################################################
boot.convex.global=function(x, Nhull = 4, nb = 1000){
#  print("global")
  islist <- FALSE
  idf<-NULL
  names.x <- names(x)
  if (is.list(x)) {
    isldata<-any(names.x=="df")
    laux<-lx<-x
    if (isldata){
      idf<-which(names.x=="df")
      lxaux<-lx
    } else  {
      x<-lx[[1]]
      laux<-lx
    }
    islist<-TRUE
  } else {
    if (is.data.frame(x))   laux<-lx<-list("df"=x)
    else     laux<-lx<-list("x"=x)
    isldata<- FALSE
    islist <- FALSE
  }
  if (is.null(idf))    lx <- mfdata.cen(lx)
  else  lx <- ldata.cen(lx)
  meanX <-lx$meanX
  lx <- lx$Xcen
  name.lx<-names(lx)
  nl <-length(lx)
  isfdata<-lapply(lx,is.fdata)
  n = NROW(lx[[1]])
  if (n<Nhull)  
    stop("A small sample size")#stop("The sample size is too small")
  J1 = Nhull+1
  w = matrix(rnorm(nb*J1),nrow=nb,ncol=J1)
  w <- w^2
  sumsq <- rowSums(w)
  w <- sweep(w,1,sumsq,"/") #* expand
  w <- sweep(w,1,sqrt(rowSums(w^2)),"/") # Ajuste para varianza similar
  ix=matrix(sample(1:n,size=J1*nb,replace=TRUE),ncol=J1)
  lxx <- lx
  for (i in 1:nl){
    if ( isfdata[[i]])       x<-lx[[i]]$data    
    else       x<-lx[[i]]
    
    p <- NCOL(lx[[i]])
    xx = matrix(0,ncol=p,nrow=nb)
    for (j in 1:J1){
      xx=xx+sweep(x[ix[,j],,drop=F],1,w[,j], "*")
    }
    if ( isfdata[[i]]) {
         lxx[[name.lx[i]]]$data <- sweep(xx,2,meanX[[i]]$data,"+")
    } else    {
      lxx[["df"]]= sweep(xx,2,as.matrix(meanX$df),"+")  
    }
  }
  #if (islist)  names(lxx)<-names.x
  #else lxx <- lxx[[1]]
  # print("sale boot.convex.global")
  return(lxx)
}
##################################################################
boot.convex.local=function(x, Nhull = 4,nb = 1000, Nnbh = NULL, metric =NULL){
  J1 <-Nhull+1
  w=matrix(rnorm(nb*(J1)),nrow=nb,ncol=J1)
  w=w^2
  sumsq=apply(w,1,sum)
  w=sweep(w,1,sumsq,"/")
  w=sweep(w,1,sqrt(apply(w^2,1,sum)),"/")
  
  clase <- class(x)[1]
  #isfdata<-is.fdata(x);  #isldata<-is.ldata(x)  #ismatrix<-is.matrix(x)
  isldata<-isfdata<-ismatrix<-isdf<-FALSE
  xaux  <- x
  switch(clase,
         list={
           D=metric.mfdata(x, method="euclidean")
           n=nrow(x[[1]])
           isldata <- TRUE}, ###
         ldata={
           D=metric.ldata(x, method="euclidean")
           n=nrow(x[[1]])
           isldata <- TRUE},
         matrix={ 
           D=metric.dist(x)
           n=nrow(x)
           ismatrix<-TRUE},
         data.frame={ 
           D=metric.dist(x)
           n=nrow(x)
           x <- as.matrix(x)
           isdf<-TRUE},
         fdata={
           D=metric.lp(x)
           n=nrow(x)
           isfdata<-TRUE}         )
  qu=apply(D,1,rank)
  ix=sample(1:n,size=nb,replace=TRUE)
  if (isldata){
    lenl<-length(x)
    clases<-sapply(x,class)
    lx <- x
    #if (isfdata)  { x<-x$data}
    xx <- list()
    for (j in 1:lenl)
      xx[[j]]=matrix(0,ncol=ncol(x[[j]]),nrow=nb)
    
    for (i in 1:nb){
      nix = sample(which(qu[ix[i],]<=Nnbh),size=J1,replace=TRUE)
      # print(i)
      for (j in 1:lenl){
        p <- ncol(x[[j]])
        if (clases[j]=="fdata"){
            x0 <- x[[j]]$data  
        } else x0 <- as.matrix(x[[j]])
        xm = colMeans(x0[nix,,drop=F])
        xx[[j]][i,]=matrix(w[i,],nrow=1)%*%(x0[nix,]-matrix(xm,byrow=TRUE,ncol=p,nrow=J1))+xm
      }}
  }else{
    p<-ncol(x)
    xx=matrix(0,ncol=p,nrow=nb)
    if (isfdata)  { x<-x$data}
    for (i in 1:nb){
      nix=sample(which(qu[ix[i],]<=Nnbh),size=J1,replace=TRUE)
      xm=colMeans(x[nix,,drop=FALSE])
      xx[i,]=matrix(w[i,],nrow=1)%*%(x[nix,]-matrix(xm,byrow=TRUE,ncol=p,nrow=J1))+xm
    }
    if (isfdata)  { 
      xaux$data<-xx
      xx<-xaux
  }}
  if (isldata){
    xaux <- x
    for (j in 1:lenl){
      if (clases[j]=="fdata")
        xaux[[j]]$data <- xx[[j]]
      else xaux[[j]] <- xx[[j]]
      }
    xx <- xaux
  }
  return(xx)
}
##################################################################
boot.smooth.ldata <- function(x, nb = 1000, smo = 0){
#print("SMOOTH bootstrap")
  #print(smo)
  islist <- FALSE
  idf<-NULL
  names.x <- names(x)
  if (is.list(x)) {
    isldata<-any(names.x=="df")
    laux<-lx<-x
    if (isldata){
      idf<-which(names.x=="df")
      lxaux<-lx
    } else  {
      x<-lx[[1]]
      laux<-lx
    }
    islist<-TRUE
  } else {
    if (is.data.frame(x))   laux<-lx<-list("df"=x)
    else     laux<-lx<-list("x"=x)
    isldata<- FALSE
    islist <- FALSE
  }
  name.lx<-names(lx)
  nl <-length(lx)
  isfdata<-sapply(lx,is.fdata)
  n = NROW(lx[[1]])
  ix = sample(1:n,size=nb,replace=TRUE)
  lxx <- lx
  ncl<-sapply(lx,NCOL)
  nrl<-sapply(lx,NROW)
  nl <- length(lx)
  
  for (i in 1:nl){
    p <- ncl[i]
    n <- nrl[i]
    err <- 0
    if ( isfdata[i]){
      x<-lx[[i]]$data    
    if (smo>0) {
        err <- mvrnorm(n=nb,rep(0,p),var(x)*smo) 
        #err <- rproc2fdata(n=nb,lx[[i]]$argvals,Sigma=(var(x)*smo))$data
        }
      lxx[[name.lx[i]]]$data <- x[ix,,drop=F] + err
    }           
    else       {
      x<-lx[[i]]
      if (smo>0) {err <- mvrnorm(n=nb,rep(0,p),var(x)*smo)}
      lxx[["df"]]= x[ix,,drop=F] + err
    }
  }
  return(lxx)
}

##################################################################
boot.naive.ldata <- function(x,nb = 1000){
  #print("Naive bootstrap**************************************")
  islist <- FALSE
  idf<-NULL
  names.x <- names(x)
  if (is.list(x)) {
    isldata<-any(names.x=="df")
    laux<-lx<-x
    if (isldata){
      idf<-which(names.x=="df")
      lxaux<-lx
    } else  {
      x<-lx[[1]]
      laux<-lx
    }
    islist<-TRUE
  } else {
    if (is.data.frame(x))   laux<-lx<-list("df"=x)
    else     laux<-lx<-list("x"=x)
    isldata<- FALSE
    islist <- FALSE
  }
  name.lx<-names(lx)
  nl <-length(lx)
  isfdata<-lapply(lx,is.fdata)
  n = NROW(lx[[1]])
  ix=sample(1:n,size=nb,replace=TRUE)
  lxx <- lx
  for (i in 1:nl){
    if ( isfdata[[i]])       x<-lx[[i]]$data    
    else       x<-lx[[i]]
    p <- NCOL(lx[[i]])
    if ( isfdata[[i]]) {
      lxx[[name.lx[i]]]$data <- x[ix,,drop=F]
    } else    {
      lxx[["df"]]= x[ix,,drop=F]
    }
  }
  return(lxx)
}
##################################################################
boot.naive.raw<-function(x, nb = 1000,nr = nrow(x)){
  return(x[sample(nr,size=nb,replace=TRUE),])
}
##################################################################
# n: tamaño de la muestra a devolver
# Sigma: varianza si se usa smoothing
boot.smooth.df<-function(x, smo = 0.05, n = NROW(x), Sigma=var(x)){
  if (smo>0){
  return(x[sample(NROW(x),size=n,replace=TRUE),]+ mvrnorm(n,rep(0,NCOL(x)),Sigma*smo))}
  else {
    return(x[sample(NROW(x),size=n,replace=TRUE),])
  }
}
##################################################################
boot.smooth.fdata<-function(x, smo = 0.05, n = NROW(x), Sigma=cov(x$data)){
  x$data <- boot.smooth.df(x$data,smo=smo,n=n,Sigma=Sigma)
  return(x)  #func.var(x)
#  return(x[sample(nr,size=nb,replace=TRUE),]+ rproc2fdata(n=nb,x$argvals,Sigma=var(x$data)*smo))
}


##################################################################
# smooth.boot=function(x,nb=1000,smo = 0.05){
#   #  print("Naive bootstrap")
#   islist <- FALSE
#   idf<-NULL
#   names.x <- names(x)
#   if (is.list(x)) {
#     isldata<-any(names.x=="df")
#     laux<-lx<-x
#     if (isldata){
#       idf<-which(names.x=="df")
#       lxaux<-lx
#     } else  {
#       x<-lx[[1]]
#       laux<-lx
#     }
#     islist<-TRUE
#   } else {
#     if (is.data.frame(x))   laux<-lx<-list("df"=x)
#     else     laux<-lx<-list("x"=x)
#     isldata<- FALSE
#     islist <- FALSE
#   }
#   name.lx<-names(lx)
#   nl <-length(lx)
#   isfdata<-lapply(lx,is.fdata)
#   n = NROW(lx[[1]])
#   ix=sample(1:n,size=nb,replace=TRUE)
#   lxx <- lx
#   for (i in 1:nl){
#     
#     if ( isfdata[[i]])       x<-lx[[i]]$data    
#     else       x<-lx[[i]]
#     n = NROW(lx[[i]])
#     p <- NCOL(lx[[i]])
#     #smt <- mvrnorm(n=nb,rep(0,p),var(x) * abs(smo))
#     smt <- mvrnorm(n=nb,rep(0,p),var(x[ix,]) * abs(smo))
#     print(dim(x))
#     print(dim(smt))
#     if ( isfdata[[i]]) {
#       lxx[[name.lx[i]]]$data <- x[ix,,drop=F]+smt
#     } else    {
#       lxx[["df"]]= x[ix,,drop=F]+smt
#     }
#   }
#   return(lxx)
# }
##################################################################
# 
# library(fda.usc)
# data(tecator)
# absorp<-tecator$absorp.fdata
# 
# 
# plot(absorp,col=1)
# set.seed(1)
# out.naive<-boot.naive(absorp,nb=20)
# lines(out.naive,col=2,lwd=2)
# 
# set.seed(1)
# out.smo <- boot.smooth.fdata(absorp,nb=20,smo=0.0)
# lines(out.smo,col=4,lwd=2)
# 
# set.seed(1)
# out.smo<-boot.smooth.prueba(list("absorp"=absorp),nb=20,smo=0.00)
# lines(out.smo$absorp,col=3,lwd=2)

# cambiar mvrnom por rproc2fdata

# set.seed(1)
# out.boot=fdata.bootstrap(absorp,statistic=func.trim.FM,nb=20,draw=TRUE)
# lines(out.boot$fdataobj.boot,col=2)
#############

###########################################
# # funcion para borrar
# bconvex=function(x,J=4,B=5000){
#   n=nrow(x)
#   xm=colMeans(x)        #* 
#   xc=sweep(x,2,xm,"-")  #Nuevo quitar la media 
#   w=matrix(rnorm(B*(J+1)),nrow=B,ncol=J+1)
#   w=w^2
#   sumsq=apply(w,1,sum)
#   w=sweep(w,1,sumsq,"/")
#   w=sweep(w,1,sqrt(apply(w^2,1,sum)),"/") # Ajuste para varianza similar
#   ix=matrix(sample(1:n,size=(J+1)*B,replace=TRUE),ncol=J+1)
#   xx=matrix(0,ncol=ncol(x),nrow=B)
#   for (j in 1:(J+1)){
#     xx=xx+sweep(xc[ix[,j],],1,w[,j], "*") #Con datos centrados
#   }
#   xx=sweep(xx,2,xm,"+")  #Reponer la media #*
#   return(xx)
# }
# ###########################################
# 
# naive.boot=function(x,nb=1000){
# #  print("Naive bootstrap")
#   islist <- FALSE
#   idf<-NULL
#   names.x <- names(x)
#   if (is.list(x)) {
#     isldata<-any(names.x=="df")
#     laux<-lx<-x
#     if (isldata){
#       idf<-which(names.x=="df")
#       lxaux<-lx
#     } else  {
#       x<-lx[[1]]
#       laux<-lx
#     }
#     islist<-TRUE
#   } else {
#     if (is.data.frame(x))   laux<-lx<-list("df"=x)
#     else     laux<-lx<-list("x"=x)
#     isldata<- FALSE
#     islist <- FALSE
#   }
#   name.lx<-names(lx)
#   nl <-length(lx)
#   isfdata<-lapply(lx,is.fdata)
#   n = NROW(lx[[1]])
#   ix=sample(1:n,size=nb,replace=TRUE)
#   lxx <- lx
#   for (i in 1:nl){
#     if ( isfdata[[i]])       x<-lx[[i]]$data    
#     else       x<-lx[[i]]
#     p <- NCOL(lx[[i]])
#     if ( isfdata[[i]]) {
#       lxx[[name.lx[i]]]$data <- x[ix,,drop=F]
#     } else    {
#       lxx[["df"]]= x[ix,,drop=F]
#     }
#   }
#   return(lxx)
# }
###############################################################################################

# Delete
# nb<-5
# nc<-3
# x0<-x$data
# smo=.01
# tt<-1:10/10
# x<-rproc2fdata(1000,tt,2*sin(tt*2*pi))
# x1<-boot.smooth.fdata(x, smo = 0.5, nb = 1000, nr = NROW(x), nc = NCOL(x))
# x2<-boot.naive.raw(x, nr = NROW(x))
# plot(x,col=1)
# lines(x2,col=4,lwd=2)
# lines(x1,col=2)
# lines(x,col=1)
# x1<-boot.smooth.ldata(list("x"=x), smo = 0.5, nb = 1000)
# class(x1)
# x2<-boot.convex.local(x, Nhull = 4,nb = 1000)
# dim(x2$data)