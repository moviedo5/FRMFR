#' Create Basis Set for Functional Data of fdata class
#' 
#' @description Compute basis for functional data.
#' 
#' @aliases create.fdata.basis create.pc.basis create.pls.basis
#' create.raw.fdata 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Vector of response (scalar).
#' @param l Vector of basis index.
#' @param maxl maximum number of basis
#' @param type.basis Type of basis (see create.basis function).
#' @param rangeval A vector of length 2 giving the lower and upper limits of
#' the range of permissible values for the function argument.
#' @param norm If \code{TRUE} the norm of eigenvectors \code{basis} is 1.
#' @param class.out =="fd" basisfd class, =="fdata" fdata class.
#' @param basis "fd" basis object.
#' @param lambda Amount of penalization. Default value is 0, i.e. no
#' penalization is used.
#' @param P If P is a vector: coefficients to define the penalty matrix object.
#' By default P=c(0,0,1) penalize the second derivative (curvature) or
#' acceleration.  If P is a matrix: the penalty matrix object.
#' @param \dots Further arguments passed to or from other methods.
#' @return 
#' \itemize{
#' \item \code{basis}{ basis} 
#' \item \code{x}{ if \code{TRUE} the value of the rotated data (the centred data multiplied by the basis matrix) is returned}
#' \item \code{mean}{ functional mean of \code{fdataobj}} 
#' \item \code{df}{ degree of freedom} 
#' \item \code{type}{ type of basis}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \link[fda]{create.basis} and \code{\link{fdata2pc}}.
#' @references Ramsay, James O. and Silverman, Bernard W. (2006),
#' \emph{Functional Data Analysis}, 2nd ed., Springer, New York.
#' 
#' N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). Penalized Partial Least
#' Squares with Applications to B-Spline Transformations and Functional Data.
#' Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69.
#' \doi{10.1016/j.chemolab.2008.06.009}
#' @keywords multivariate
#' @examples
#' \dontrun{
#' data(tecator)
#' basis.pc<-create.pc.basis(tecator$absorp.fdata,c(1,4,5))
#' plot(basis.pc$basis,col=1)
#' summary(basis.pc)
#' basis.pls<-create.pls.basis(tecator$absorp.fdata,y=tecator$y[,1],c(1,4,5))
#' summary(basis.pls)
#' plot(basis.pls$basis,col=2)
#' summary(basis.pls)
#' 
#' basis.fd<-create.fdata.basis(tecator$absorp.fdata,c(1,4,5),
#' type.basis="fourier")
#' plot(basis.pc$basis)
#' basis.fdata<-create.fdata.basis(tecator$absorp.fdata,c(1,4,5),
#' type.basis="fourier",class.out="fdata")
#' plot(basis.fd,col=2,lty=1)
#' lines(basis.fdata,col=3,lty=1)
#' }

#' @export
create.fdata.basis <-function(fdataobj, l = 1:5, maxl = max(l), type.basis = "bspline", 
          rangeval = fdataobj$rangeval, class.out = "fd") 
{
  aa1 <- paste("create.", type.basis, ".basis", sep = "")
  if (type.basis == "pc") {
    as <- list()
    as$fdataobj <- fdataobj
    as$l <- l
    basis = do.call(aa1, as)
  }
  if (type.basis %in% c("bspline", "fourier", "constant", "exponential", 
                        "polygonal", "power")) {
    as <- list()
    as$rangeval <- rangeval
    as$nbasis <- maxl
    basis = do.call(aa1, as)
    #basis$params <- diff(rangeval)
    basis$dropind <- setdiff(1:maxl, l)
    if (class.out == "fdata") {
      nam <- basis$names[intersect(1:maxl, l)]
      basis = fdata(t(eval.basis(fdataobj$argvals, basis)), 
                    fdataobj$argvals, fdataobj$rangeval)
      rownames(basis$data) <- nam
      basis$type <- type.basis
      basis$nbasis <- maxl
      basis$dropind <- as$dropind
    }
  }
  basis
}


# scores.basis <-function(fdataobj,l=1:5,maxl=max(l),type.basis="bspline", lambda ){
#       aa1 <- paste("create.",type.basis,".basis", sep = "")
#       if (type.basis=="pc"){
#         
#         as <- list() 
#         as$fdataobj <- fdataobj
#         as$l <- l
#         if (missing(lambda)) lambda = 0
#         as$lambda <- lambda
#         basis=do.call(aa1,as)
#       }
#       if (type.basis %in% c("bspline","fourier","constant","exponential"
#                             ,"polygonal","power")){
#         
# ################
# #        fdataobj<-tecator$absorp.fdata
# #        l=1:5
# #        type.basis="bspline"
# #        lambda = NULL
# ################        
#         # as <- list()
#         maxl=max(l)
#         if (missing(lambda)) lambda = NULL
#         #as$lambda <- lambda
#         #basis=do.call(aa1,as)
#         basis0 = fdata2fd(fdataobj, type.basis, nbasis = maxl, lambda = lambda)
#         basis<-list()
#         basis$basis <- basis0$basis
#         basis$basis.fdata=fdata(t(eval.basis(fdataobj$argvals,basis0$basis)),
#                           fdataobj$argvals,fdataobj$rangeval)
#         basis$coefs <- t(basis0$coefs)  
#         basis$fdataobj <- fdataobj
#         basis$l <- l
#         basis$lambda <- as$lambda
#         basis$type <- type.basis
#         basis$nbasis <-maxl
#         basis$dropind<-setdiff(1:maxl,l)
#         basis$rangeval <- fdataobj$rangeval
#         basis$mean <- func.mean(fdataobj)
#         
#         #x.fd = Data2fd(argvals = tt,
#          #              y = t(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]]$data), 
#           #             basisobj = object$basis.x[[vfunc[i]]]$basis, 
#            #            fdnames = fdnames)
#         
#         
#         #as$nbasis <-maxl
#         #as$dropind<-setdiff(1:maxl,l)
#         #as$rangeval <- fdataobj$rangeval
#         #nam<-basis$names[intersect(1:maxl,l)]
#         #rownames(basis$data)<-nam
#         #basis$type<-type.basis
#         #basis$nbasis<-maxl
#         #basis$dropind<-as$dropind
#       }
#   basis
# } 
#######################
     
#' @rdname create.fdata.basis
#' @export
create.pc.basis<-function(fdataobj,l=1:5,norm=TRUE,basis=NULL,
                          lambda=0,P = c(0, 0, 1),...){
 tt<-fdataobj$argvals
 rtt<-fdataobj$rangeval
 dropind=NULL
 if (lambda>0) pc <- fdata2pc(fdataobj,norm=norm,ncomp=max(l),lambda=lambda,P=P,...)
 else  pc <- fdata2pc(fdataobj,norm=norm,ncomp=max(l),...)
 
 vs<-pc$basis$data    
 lenl<-length(l) 
 pc.fdata<-pc$u[,l,drop=FALSE]%*%(diag(lenl)*pc$d[l])%*%vs[l,,drop=FALSE]
 pc.fdata<-sweep(pc.fdata,2,matrix(pc$mean$data,ncol=1),"+")
 basis.pc = pc$basis[l, ,drop=FALSE]
 rownames(basis.pc$data) <- paste("PC", l, sep = "")
 # basisobj <- pc
 fdnames<- colnames(pc$coefs[,l,drop=FALSE])
 if (is.null(basis)) {
   pc.fdata <- fdata(pc.fdata,tt,rtt,fdataobj$names)
   # out <- list(basis = basis.pc, coefs = pc$coefs, mean = pc$mean,
   #             fdataobj.pc=pc.fdata, fdataobj.cen = pc$fdataobj.cen,
   #             fdataobj = fdataobj,l = l,norm=norm,lambda = pc$lambda,
   #             lambda=lambda, P=P , type = "pc",call="fdata2pc",values=pc$d)
   # class(out) <- "fdata.comp"
   pc$fdata.est <- pc.fdata
   names(pc$d) <- paste("PC", 1:length(pc$d), sep = "")
   pc$l <- paste("PC", l, sep = "")
   pc$coefs <- pc$coefs[,l,drop=F]
   pc$lambda <- lambda
   pc$type <- "pc"
   pc$basis <- basis.pc
   pc$call <- match.call()
   return(pc)
   }
 else {
      fdobj<- Data2fd(argvals = tt, y = t(pc.fdata),basisobj = basis)
      out<-list()
      out$harmonics<-fdobj
      colnames(out$harmonics$coefs)<-rownames(fdataobj$data)
      out$values<-pc$newd^2
      #out$scores<-pc$coefs[,l,drop=FALSE]
      #rownames(out$scores)<-rownames(fdataobj$data)
      out$coefs<-pc$coefs[,l,drop=FALSE]
      rownames(out$coefs)<-rownames(fdataobj$data)
      out$varprop<-out$values[l]/sum(out$values)
      out$meanfd<- Data2fd(argvals = tt, y = pc$mean$data[1,],basisobj = basis)
      out$call <- match.call()
      class(out) <- "fdata2pc"
      return(out) 
      }
}



#' @rdname create.fdata.basis
#' @export
create.pls.basis<-function(fdataobj, y, l=1:5, norm=TRUE,
                           lambda=0, P = c(0, 0, 1),...){
if (lambda>0) pls<-fdata2pls(fdataobj,y,norm=norm,ncomp=max(l),lambda=lambda,P=P,...)
 else  pls<-fdata2pls(fdataobj,y,norm=norm,ncomp=max(l),...)
     basis=pls$basis[l,,drop=FALSE]
     rownames(basis$data)<-paste("PLS",l,sep="")
     fdata.est <- gridfdata(pls$coefs, pls$basis, pls$mean)
     #fdata.est <- fdata(pl$coefs %*% pl$basis$data,fdataobj$argvals,fdataobj$rangeval)
pls$call <- match.call()
pls$fdata.est <- fdata.est
pls$l <-l
pls$type <- "pls"
#out <- list(call="fdata2pls","basis"=basis,"coefs"=pls$coefs,"mean"=pls$mean,"df"=pls$df,
#"fdataobj.cen"=pls$fdataobj.cen,"fdataobj"=fdataobj,norm=norm,
#"l"=l,"type"="pls","y"=y,fdata.est=fdata.est)
#class(pls) <- "fdata.comp"
return(pls)
} 


#' @rdname create.fdata.basis
#' @export
create.raw.fdata=function (fdataobj, l = 1:nrow(fdataobj))
{
    return(list(basis =fdataobj[l,] , type = "raw"))
}

#########################
create.mfdata.basis <- function(mfdata, l = 1:5
                                , type.basis = "bspline"
                                , class.out = "fd") 
{
  nvar <- length(mfdata)
  nam <- names(mfdata)
  aa1 <- paste("create.", type.basis, ".basis", sep = "")
  basis.x <- NULL
  maxl <- max(l)
  if (type.basis == "pc") {
    for (i in 1:nvar) {
      as <- list()
      as$fdataobj <- mfdata[[nam[i]]]
      as$l <- l
      basis = do.call(aa1, as)
      basis.x[[nam[i]]] <- basis
  }}
  
  if (type.basis %in% c("bspline", "fourier", "constant", "exponential", 
                        "polygonal", "power")) {
    for (i in 1:nvar) {
    as <- list()
    fdataobj <- mfdata[[nam[i]]]
    rangeval = fdataobj$rangeval
    as$rangeval <- rangeval
    as$nbasis <- maxl
    as$dropind <- setdiff(1:maxl, l)
    basis = do.call(aa1, as)
    
    if (class.out == "fdata") {
      nam <- basis$names[intersect(1:maxl, l)]
      basis = fdata(t(eval.basis(fdataobj$argvals, basis)), 
                    fdataobj$argvals, fdataobj$rangeval)
      rownames(basis$data) <- nam
      basis$type <- type.basis
      basis$nbasis <- maxl
      basis$dropind <- as$dropind
    }
   basis.x[[nam[i]]] <- basis
   }
  }
  return(invisible(basis.x))
}

#########################
create.ldata.basis <- function(x, l = 1:5
                               , type.basis = "bspline"
                               , class.out = "fd") 
{
  clases <- sapply(x,class)
  ifdata <- which(clases == "fdata")
  basis.x <- create.mfdata.basis(x[ifdata], l = l, 
             type.basis = type.basis, class.out = "fd") 
  return(basis.x)
}  


# create.wavelets.basis <- function(
#     rangeval = NULL, 
#     nbasis = NULL, 
#     lev = NULL
#     ){
#   # Creación de una base de Haar
#   
#  if (is.null(lev)) lev <- 2
#  if (is.null(rangeval)) rangeval <- c(0,1)
#  else {
#    if (length(rangeval) == 1) {
#      if (rangeval <= 0) 
#        stop("'rangeval' a single value that is not positive, is ", 
#             rangeval)
#      rangeval = c(0, rangeval)
#    }
#    if (length(rangeval) > 2) {
#      if (!is.null(breaks)) 
#        stop("breaks can not be provided with length(rangeval) > 2;  ", 
#             " length(rangeval) = ", length(rangeval), 
#             " and length(breaks) = ", length(breaks))
#      breaks <- rangeval
#      rangeval <- range(breaks)
#    }
#    if (rangeval[1] >= rangeval[2]) 
#      stop("rangeval[1] must be less than rangeval[2];  instead ", 
#           "rangeval[1] = ", rangeval[1], c("==", 
#           ">")[diff(rangeval) < 0], " rangeval[2] = ", rangeval[2])
#  
#    }
#  if (is.null(nbasis)) nbasis <- 2^lev
#  else   {
#    if (!is.numeric(nbasis)) 
#      stop("nbasis must be numeric, is ", class(nbasis))
#    if ((lnb <- length(nbasis)) > 1) 
#      stop("nbasis must be a single positive integer;  ", 
#           "length(nbasis) = ", lnb, " > 1;  first 2 elements = ", 
#           nbasis[1], ", ", nbasis[2])
#    if ((nbasis%%2) > 0) 
#      stop("nbasis is not  a power of 2")
#   }
# 
#  t01 <- seq(rangeval[1], rangeval[2], len = nbasis + 1)
#  # función wavelet madre 
#   mother <- function(x){
#     # modificar para que funcionen en otro grid
#    ifelse(x >= 0 & x <= .5, 1, ifelse(x >= .5 & x <= 1, -1, 0))
#  } 
#  wavelet <- c(fdata(rep(1, nbasis+1), t01),
#               fdata(mother(t01), argvals = t01))
#  nam <- c("Cte","H0,0")
#  nk <- c(1,1)
#  for (k in 1:lev){ 
#    for(j in 0:(2^k-1)){ 
#      wavelet <- c(wavelet, 
#                   fdata( 2^(k/2) * mother(2^k * t01 - j), argvals = t01)
#                   )
#      nam <- c(nam, paste("H", k,",", j, sep =""))
#      nk <- c(nk,k)
#    }
#  }
#  rownames(wavelet$data) <- nam
#  wavelet # lista?
# }
 
# kk<-10
# b1 <- create.wavelets.basis(c(0,1),2^kk,2)
# plot(b1)  
# dim(b1)
# b1$type <- "wvl"
# tt <- b1$argvals
# x <- rproc2fdata(30,tt,sin(pi*tt),sigma=0.1)
# plot(x)
# a1 <- fdata2basis(x,b1,"inprod")
# a2 <- fdata2basis(x,b1,"grid")
# # A = t(inprod.fdata(x, b1))
# # B = inprod.fdata(b1)
# # coefs = t(solve(B) %*% A)
# xg1 <- xg2 <- x
# xg1$data <- a1$coefs %*% b1$data
# xg2$data <- a2$coefs %*% b1$data
# xg1 <- xg1+func.mean(x)
# plot(x,col=1)
# lines(xg1,col=2)
# xg2 <- xg2 + func.mean(x)
# lines(xg2,col=4)
# # plot(x-xg)

# comparar con la implementación 
# D:\Users\moviedo\OneDrive - Universidade da Coruña\FDA_documentos\Codigo R\Wavelets\R

# # basisobj <- basisfd()
# t01 <- seq(0, 1, len = 513)
# # función wavelet madre 
# mother <- function(x){
#   ifelse(x >= 0 & x <= .5, 1, ifelse(x >= .5 & x <= 1, -1, 0))
# } 
# wavelet <- c(fdata(rep(1, length(t01)), t01),
#              fdata(mother(t01), argvals = t01))
# nam <- c("Cte","H0,0")
# nk <- c(1,1)
# for (k in 1:2){ 
#   for(j in 0:(2^k-1)){ 
#     wavelet <- c(wavelet, fdata(2^(k/2) * mother(2^k*t01-j), argvals = t01))
#     nam <- c(nam,paste("H",k,",",j,sep = ""))
#     nk <- c(nk,k)
#   }
# }
# rownames(wavelet$data) <- nam
# dim(wavelet)
# plot(wavelet)
# round(inprod.fdata(wavelet), 3)
