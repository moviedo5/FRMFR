# funcion optimizada que solo llama 1 vez al lm()
# en lugar npy veces

# hay que revisar el codigo de las beta funcionales siguiendo 
# el beta2d del fregre.basis.fr

# Añadir names en fdata cuando es 2d!:
# beta2d[[vfunc[i]]] <-  fdata(dd,list(data[[vfunc[i]]]$argvals,tty),fdata2d=TRUE)

#' Fitting Functional Linear Models with functional response using basis representation.
#' 
#' @description Computes functional regression between functional (and non functional)
#' explanatory variables and functional response using basis representation.
#' 
#' @details This section is presented as an extension of the functional linear regression models with scalar response:
#' \code{\link{fregre.lm}}, \code{\link{fregre.pc}}, \code{\link{fregre.pls}} and
#' \code{\link{fregre.basis}}. Now, the functional response \eqn{Y} is estimated by
#' more than one functional covariate \eqn{X^j(t)} and also more than one non
#' functional covariate \eqn{Z^j}. The regression model is given by:
#' \deqn{E[Y|X,Z]=\alpha+\sum_{j=1}^{p}\beta_{j}Z^{j}+\sum_{k=1}^{q}\frac{1}{\sqrt{T_k}}\int_{T_k}{X^{k}(t)\beta_{k}(t)dt}
#' }{E[Y|X,Z]=\alpha+\sum_j \beta_j Z^j + \sum_k <X^k,\beta_k>}
#' 
#' where \eqn{Y(s)=\left[ Y^{1}(s_1),\cdots,Y^{p}(s_p)\right]}{Y(s)=[Y^1(s),...,Y^r(s)}, \eqn{Z=\left[ Z^1,\cdots,Z^p \right]}{Z=[Z^1,...,Z^p]} are the non
#' functional covariates, \eqn{X(t)=\left[ X^{1}(t_1),\cdots,X^{q}(t_q)
#' \right]}{X(t)=[X^1(t),...,X^q(t)]} are the functional ones and
#' \eqn{\epsilon} are random errors with mean zero , finite variance
#' \eqn{\sigma^2} and \eqn{E[X(t)\epsilon]=0}{E[X(t)\epsilon]=0}.  
#' 
#' \code{basis.y}, is a basis for represent the response. Parameter option:
#' \itemize{
#' \item Basis representation.  By default a basis of bsplines is used 
#' as in the \code{basis.x} parameter definition.\cr

#  Type of  basis objects :
# \code{\link{create.pc.basis}}, \code{\link{pca.fd}}
# \code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} or
# \code{\link{create.basis}}.

#' \item Raw representation. 
#' The response is estimated directly using the  \code{\link{fdata}} object.
#' }
#' 
#' \code{basis.x} is a list of basis for represent each functional covariate.
#' The basis object can be created by the function:
#' \code{\link{create.pc.basis}}, \code{\link{pca.fd}}
#' \code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} or
#' \code{\link{create.basis}}.\cr \cr
#' 
# The user can penalty the basis elements by: (i) \code{lambda} is a list of
# rough penalty values for the second derivative of each functional covariate,
# see \code{\link{fregre.basis}} for more details.\cr (ii) \code{rn} is a list
# of Ridge penalty value for each functional covariate, see
# \code{\link{fregre.pc}}, \code{\link{fregre.pls}} and
# \code{\link{P.penalty}} for more details.\cr Note: For the case of the
# Functional Principal Components basis two penalties are allowed (but not the
# two together). \cr
#' 
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' details of model specification are given under \code{Details}.
#' @param data List that containing the variables in the model. Allowed objects:
#' \itemize{
#' \item A list with \code{\link{fdata}} objects and the "df" object for non-functional covariates.
#' \item An object of class \code{\link{ldata}}. (equivalent to the previous one)
#' \item An object of class \code{\link{mfdata}} (only for \code{\link{fdata}} objects).
#' }
#' @param basis.y Basis for functional response variable.
#' @param basis.x List of basis for functional explanatory data estimation.
#' @param lambda List, indexed by the names of the functional covariates, which contains the Roughness penalty parameter.
#' @param P List, indexed by the names of the functional covariates, which contains the parameters for the creation of the penalty matrix.
# @param weights weights
#' @param \dots Further arguments passed to or from other methods.
#' @return Return \code{lm} object (MLM output) plus:
#' \itemize{
#' \item \code{sr2}{ Residual variance.}
# \item \code{Vp}{ Estimated covariance matrix for the parameters.} 
# \item \code{lambda}{ A roughness penalty.} 
#' \item \code{basis.y}{ Basis used for response.} 
#' \item \code{basis.x}{ Basis used for \code{fdata} or \code{fd} covariates.} 
# \item \code{beta.l}{ List of estimated beta parameter of functional covariates.}
#' \item \code{data}{ List that containing the variables in the model.}
#' \item \code{formula}{ formula.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as: \code{\link{predict.fregre.lm}} and
#' \code{\link{summary.lm}}.\cr Alternative method: \code{\link{fregre.glm}}.
#' @references Chiou, J. M., Muller, H. G. and Wang, J. L. (2004). Functional response models. Statistica Sinica, 675--693.
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' absorp <- tecator$absorp.fdata
#' x <- fdata.deriv(absorp)
#' y <- fdata.deriv(x,1)
#' x2 <- fdata.deriv(y)
#' mdat <- mfdata("y"=y,"x"=x,"x2"=x2)
#' # YBSP - XBSP
#' ff <- as.formula(y~x2)
#' res1 <-  fregre.lm.fr(ff,mdat)
#' plot(res1$beta.l$x2)
#' 
#' bsp7 <- create.bspline.basis(mdat$y$rangeval,7)
#' bsp5 <- create.bspline.basis(mdat$y$rangeval,5)
#' lbsp <- list("x"=bsp7,"x2"=bsp5)
#' res2 <-  fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=bsp7,basis.x=lbsp)
#' plot(res2$beta2d$x)
#' plot(res2$beta2d$x2)
#' plot(res2$beta2d$x2,type="persp",phi=30,theta=30)
#' # YPC - XPC
#' pcy <- create.pc.basis(mdat$y,1:4)
#' pcx <- create.pc.basis(mdat$x,1:5)
#' pcx2 <- create.pc.basis(mdat$x2,1:3)
#' lpc <- list("x"=pcx,"x2"=pcx2)
#' res3 <-  fregre.lm.fr(as.formula(y~x2),mdat,basis.y=bsp7,basis.x=lpc)
#' res3 <-  fregre.lm.fr(as.formula(y~x2),mdat,basis.y=pcy,basis.x=lpc)
#' res3 <-  fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=pcy,basis.x=lpc)
#' 
#' par(mfrow=c(1,2))
#' plot(res1$beta.l$x)
#' plot(res2$beta.l$x2)
#' 
#' # YPC - XBSP
#' pcy <- create.pc.basis(mdat$y,1:3)
#' res4 <-  fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=pcy)
#' bsp7 <- create.bspline.basis(mdat$y$rangeval,7)
#' bsp5 <- create.bspline.basis(mdat$y$rangeval,5)
#' lbsp <- list("x"=bsp7,"x2"=bsp5)
#' res5 <-  fregre.lm.fr(as.formula(y~x+x2),mdat,basis.y=bsp7,basis.x=lbsp)
#' plot(res5$beta2d$x,type="persp")
#' plot(res5$beta2d$x2,type="persp",phi=45,theta=60)
#' }
#' 
#' @aliases fregre.lm.fr
#' @export

fregre.lm.fr <- function(formula,data,basis.y=NULL
                          ,basis.x=NULL
                          ,lambda=NULL,P=NULL
                          ,...){
  # res0 <-  fregre.mlm.fr(as.formula(y~x),dat)
  # formula<-as.formula(y~x)
  # data<-mdat
  ##### lambda=0
  ##### P=c(1,0,0)
  # res2 <-  fregre.mlm.fr(as.formula(y~x),mdat,basis.y=basis.pc)
  # formula <- as.formula(y~x)
  # data <- mdat
  # basis.y=basis.pc
  # basis.x=NULL
  # lambda=NULL ;P=NULL
  if (is.null(basis.y)) raw=TRUE
  else raw=FALSE
  # print(raw)
  beta2d <- list() # se guardan las beta funcionales en 2D (superficies)
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  datanames<-names(data)
  if (any(datanames=="df")) 
    vnf=intersect(terms,names(data$df))  else vnf <- NULL
  vfunc=setdiff(terms,vnf)
  off<-attr(tf,"offset") 
  name.coef=nam=par.fregre=beta.l=list()
  ind.name<-NULL
  kterms <- 1
  isfdata <- FALSE
  XX <- list()
  lambdap <- XX2 <- NULL
  RR <- NULL
  penalization <- FALSE    
  # opciones para la respuesta (Bsp, pc o raw)
  if (any(datanames==response)) {
# print("respuesta funcional")
    isfdata<-TRUE#is.fdata(data[[response]])
    yfdata <- data[[response]]
    ydat <- yfdata$data       # si es de la clase fd??
    #       print(paste("Functional response:",response))
    #XX[[response]]<- yfdata$data  
    tty <- yfdata[["argvals"]]
    rtty <- yfdata[["rangeval"]]
    namy <- yfdata[["names"]]
    n <- nrow(ydat)
    npy <- NCOL(ydat)
    yy <- fdata.cen(data[[response]])
    #tty<-yy$argvals
    #ydat <- yy$data
    if (is.null(basis.y)) basis.y.class <- "raw"
    else basis.y.class <-class(basis.y)
    
# print(basis.y.class)    
    # if (is.null(basis.y)){
    #  basis.y<-create.fdata.basis(yfdata,l=1:7)
    #     print(class(basis.y))
    #} 
    ycoef <- ydat
    #basis.y.class <- class(basis.y)
    
    if (basis.y.class == "basisfd"){
      # print("entra if basisfd")
      #basis.y <- basis.x[[1]]
      #names(basis.x)
      fdnames=list("time"=tty,"reps"=rownames(ydat),"values"="values")
# print("peta basisfd1")
      yfd= Data2fd(argvals = tty, y = t(ydat),
                   basisobj = basis.y,fdnames=fdnames,lambda=0)
      yfd2= fdata2fd(yfdata,type.basis=basis.y$type,nbasis =basis.y$nbasis )
      par(mfrow=c(1,2))
      plot(yfd,col=1)
      lines(yfd2,col=2)
      plot(yfd-yfd2,col=2)
  # print("peta basisfd2")
      ycoef <- t(yfd$coefs)
      # print(dim(ycoef))
    }
    
    # basis.y <- create.pc.basis(yfdata,1:4)
    # class(basis.pc)
    if (basis.y.class=="fdata.comp"){
      #  # print("peta pc1")
      ycoef <- basis.y$coefs[, basis.y$l]
      #  # print("peta pc2")
    }
    
    if (basis.y.class=="character"){
      #   # print("peta raw1")
      ycoef <- ydat
      #    # print("peta raw2")
    }
    #XX[[response]]<- yfdata$data  
    XX[[response]]<- ycoef
  }
  else { 
    if (any(names(data[["df"]])==response)) {
      stop(paste0("No functional response,",response," object must be fdata class"))
      yy<-as.matrix(data[["df"]][response])
      XX[[response]]<-yy    
    }
    else stop("Response not found")      
  }
 # print("aaaaaaaaaaaaaaaaa")
  ####################  
 npy2 <- ncol(ycoef)
  nvnf<-length(vnf)
  nn<-1
  if (nvnf>0) {
    XX2<-XX[[vnf]]<-as.matrix(data$df[,vnf])
    nn<-nvnf+1
    lambdap<-rep(0,nn)
  }
  else lambdap<-0
  lvfunc<-length(vfunc)
  # penalty <- TRUE
  # if (is.null(P) | is.null(lambda)){
  #   P <- as.list(numeric(lvfunc))
  #   names(P)<- vfunc
  #   lambda <- P
  #   penalty <- FALSE
  # }
  
  #####if (length(lambda)==1) lambda<-rep(lambda,lvfunc)
  #lambda<-rep(0,lvfunc)
  # # print(paste("Functional covariate/s:",vfunc))
  if (lvfunc>0) { 
    mean.list<-vs.list<-list()
     # print(names(XX))
    # print("antes 2model")
    out <- fdata2model.fr(vfunc = vfunc, vnf=vnf, response=response, XX=XX
                          , data=data, basis.x = basis.x, #basis.y=basis.y, 
                          pf=pf, tf=tf, lambda=lambda, P=P)
      }
    XX <- out$XX
    # print(names(out))
    # print("despues 2model")
    # print(out$bsp1)
    bsp1 <- out$bsp1
    ## print(bsp1)
    name.coef <- out$name.coef
    name.coef2 <- unlist(name.coef)
    basis.x <- out$basis.x
    vs.list <- out$vs.list
    mean.list <- out$mean.list
  
  ####################################
  if (penalization){
    stop("not implemented yet")
    # z=list()
    # rownames(b.est)<-colnames(XX2)
    # z$coefficients<-b.est
    # #z$R<-R
    # z$rank<-df
    # z$df.residual<-n-df
    # if (nvnf>0) colnames(XX2)[2:(1+nvnf)]=vnf
    # coeff<-b.est[-1,,drop=F]
    # if (isfdata) {       z$residuals<-fdata(e,tty,rtty,namy)
    # z$fitted.values<-fdata(yp,tty,rtty,namy)
    # }    else {
    #   z$residuals<-e
    #   z$fitted.values<-yp
    # }                      
  }
  else { 
    #z=lm(formula=as.formula(formula),data=XX,...) 
    # # print(3)
    # # print(out$pf)
    formula <- as.formula(out$pf)
    # print(4)
    # print("y original vs ybasis")
    # print(npy)
    # print(npy2)
    zfitted <- zresid <- matrix(NA,n,npy2)
    pcoef <- zcoef<-NULL
    result<-list()
    #XX o ydatg o ycoef
    # Xmat <- data.frame(ydat[,1,drop=F])
    # print(head(ycoef))
    # print("XX")
    # print(head(XX))
    Xmat <- data.frame(ycoef[,1,drop=F])
    if (length(vnf)>0) {
      Xmat<-data.frame(Xmat,data$df[,vnf])
      names(Xmat)<-c(response,vnf)        
      #Xmat <- data$df[,vnf,drop=F]
    }  else       names(Xmat) <- response                   
    
    # # print("names(Xmat)1")
    # # print(names(Xmat))
    # # print(class(Xmat))
    
    # # print("names(XX)2")
    # print(names(XX))
    # print(class(XX))
    
    if (is.list(XX)){
      for (i in 2:length(XX)){
        Xmat<-cbind(Xmat,XX[[i]])
      }
    } else   Xmat <- (cbind(Xmat,XX))
    # print("names(Xmat)3")
    # print(colnames(Xmat))
    
    
    # # print("name.coef4")
    # # print(name.coef)
    
    # # print("name.coef5")
    # # print(name.coef2)
    
    # # print(response)
    # print(dim(ydat))
    # print(dim(ycoef))
#                                  Xmat[,response] <- ydat[,1] ########################
    #Xmat[,response] <- ycoef[,1]
     # print(formula)
     # print(dim(Xmat))
     # print(Xmat[1:3,1:3])
     z <- lm(formula=formula,data=Xmat)  
     # print("saleeeeeeeeeeeeeee")
     
    zcoef <- z$coefficients
    zfitted[,1] <-z$fitted.values
    #  # print(names(z))
    # print("bucle for npy")
  
    for (i in 2:npy2){      # Evita llamar npy2-1 veces al lm()
      #result[[i]] <- z
     # zcoef <- rbind(zcoef, qr.coef(z$qr,ydat[,i])) 
      zcoef <- rbind(zcoef, qr.coef(z$qr,ycoef[,i]))
      zfitted[,i] <-qr.fitted(z$qr,ycoef[,i], k = z$qr$rank)
      #cat(i,ycoef[1:3,i],"\n")
      # H[,i] <-  hatvalues(z)
    }
    #colnames(H) <- nam.y
    # # print(rownames(zcoef))# <- nam.y
    #  # print(bsp1)
    out <- list()
    ## print(raw)
    if (!isfdata){
      if (bsp1) {
        # print("bsp1")
        # # print(bsp1)
        # # print(dim(zfitted))
        # # print(basis.y)
        
        zfitted <- fd(t(zfitted),basis.y)
        out$fitted.values <- fdata(zfitted,tty,rtty,namy)
        # H <- fd(t(H),basis.y)
        # H <- fdata(H,tty,rtty,namy)
      } else {
        # print("no bsp1")
#        if (class(basis.y)=="fdata.comp"){
        if (inherits(basis.y,"fdata")) {
# print("fdata.comp basis.y")
          # print(zfitted[1:3,1:4])
          out$fitted.values <- gridfdata(zfitted,out$basis)
          # print("ffffffffffffffff")
      }
        #else  out$fitted.values <- fdata(zfitted,tty,rtty,namy)
      }}else  {
 # print("fdata")
        if (raw) out$fitted.values <- fdata(zfitted,tty,rtty,namy)
        else      {
          #if (class(basis.y)=="fdata.comp"){
          if (inherits(basis.y,"fdata")) {
# print("fdata.comp!! basis.y")
# print(dim(zfitted))
# print(dim(basis.y$basis))
            #out$fitted.values <- gridfdata(zfitted,basis.y$basis)# esta mal
            ##out$fitted.values <- gridfdata(zfitted,yfdata)# esta mal
            yhatfdata <- zfitted %*% basis.y$basis$data
            yhatfdata<-sweep(yhatfdata,2,matrix(basis.y$mean$data,ncol=1),"+")      
            out$fitted.values <- fdata(yhatfdata,tty,rtty,namy)
            # print("ffffffffffffffff")
          } else
            out$fitted.values <- fdata(fd(t(zfitted),basis.y),tty,rtty,namy)
        }
        # print("eeeeeeeeeeeeeeeettttttttttttttte")
        #  out$fitted.values <- fdata(fd(t(zfitted),basis.y),tty,rtty,namy)
        # print("eeeeeeeeeeeeeeeettttttttttttttte")
        # print(out$fitted.values$data[1:3,1:4])
        # print(dim(zfitted))
        # print(dim(out$fitted))
      }
  
    # if (is.vector(z$coefficients)) {
    #   z$coefficients<-matrix(z$coefficients,ncol=1)
    #   colnames(z$coefficients)<- paste(response,".",rownames(basis.y$basis$data),sep="")
    #   # # print("PC coefs")
    #   # # print(z$coefficients)
    #   # # print(name.coef)
    #   rownames(z$coefficients)<- c("(Intercept)",unlist(name.coef))
    #   #rownames(z$coefficients)<- colnames(XX)
    #   # # print(names(XX))
    #   
    # }
    # npy <- NCOL(z$coefficients)
    # z$call <- z$call[1:2]
    # df <- z$rank
    # colnames(z$coefficients) <- paste(response,".",colnames(z$coefficients),sep="")
    # rownames(z$coefficients)[-1] <- name.coef2
    # coeff <- z$coefficients
    # 
    # if (basis.y.class=="basisfd"){
    #   z$yfit<-z$fitted.values  
    #   yhatfd <- fd(t(z$fitted.values), basis.y)
    #   z$fitted.values <- fdata(yhatfd,tty,rtty,namy)
    #   z$residuals<-yfdata-z$fitted.values
    #   #z$fitted.values<-fdata(z$fitted.values,tty,rtty,namy)
    # }
    # if (basis.y.class=="fdata.comp"){
    #   z$yfit<-z$fitted.values  
    #   # yhatfdata <- z$fitted.values %*% (diag(length(basis.y$l))*basis.y$values[basis.y$l]) %*% basis.y$basis$data
    #   yhatfdata <- z$fitted.values %*% basis.y$basis$data
    #   # # print(dim(yhatfdata))
    #   # # print("fdata.comp")
    #   yhatfdata<-sweep(yhatfdata,2,matrix(basis.y$mean$data,ncol=1),"+")      
    #   z$fitted.values <- fdata(yhatfdata,tty,rtty,namy)
    #   z$residuals<-yfdata-z$fitted.values
    #   
    #   #z$residuals<-fdata(z$residuals,tty,rtty,namy)
    #   #z$fitted.values<-fdata(z$fitted.values,tty,rtty,namy)
    #   
    #   #pc.fdata<-pc$u[,l,drop=FALSE]%*%(diag(lenl)*pc$d[l])%*%vs[l,,drop=FALSE]
    #   #pc.fdata<-sweep(pc.fdata,2,matrix(pc$mean$data,ncol=1),"+")
    #   
    # }
    # if (basis.y.class=="character"){
    #   z$yfit<-z$fitted.values  
    #   z$fitted.values<-fdata(z$fitted.values,tty,rtty,namy)
    #   z$residuals<-fdata(z$residuals,tty,rtty,namy)
    # }
    
  }
  # print("bbbb")
  # print(isfdata)
  
  for (i in 1:length(vfunc)) { 
    #      if (penalization)      name.coef2<-name.coef[[vfunc[i]]]
    #      else name.coef2<-paste(vfunc[i],name.coef[[vfunc[i]]],sep="")
    name.coef2 <- name.coef[[vfunc[i]]]
    if (bsp1) {
      if (isfdata) {  
        # print("cccccccccccccccccccccccccccccccc1")
        # print(dim(zcoef))
        # print(class(zcoef))
        # print( colnames(zcoef))
        # print( rownames(zcoef))
        # print(vfunc)
        # print(name.coef)
        # print(name.coef[[vfunc[i]]])
        #print(zcoef[,name.coef[[vfunc[i]]]])
        # print("cccccccccccccccccccccccccccccccc2")
        beta.l[[vfunc[i]]]=fd(t(zcoef[,name.coef[[vfunc[i]]],drop=F]),
                              basis.x[[vfunc[i]]])
        # print("cccccccccccccccccccccccccccccccc3")
        #plot(beta.l[[vfunc[i]]])
        
      }
      else {beta.l[[vfunc[i]]] <- fd(z[["coefficients"]][name.coef2],
                                     basis.x[[vfunc[i]]])}                 
    }
    else{
      #if (class(data[[vfunc[i]]])[1] == "fdata") {
      if (inherits(data[[vfunc[i]]],"fdata")) {
        
        # print("bbb2")
        if (isfdata) {
          #  print(dim(zcoef))
          # print(name.coef)
          beta.est <- drop(as.numeric(zcoef[,name.coef2]))* vs.list[[vfunc[i]]]
        }
        else {
          beta.est <- z$coefficients[name.coef2] * vs.list[[vfunc[i]]]
        }
        beta.est$data <- apply(beta.est$data, 2, sum)
        beta.est$names$main <- "beta.est"
        beta.est$data <- matrix(as.numeric(beta.est$data), 
                                nrow = 1)
        beta.l[[vfunc[i]]] <- beta.est
        if (npy > 1) {
          for (j in 1:npy) {
            if (isfdata) 
              beta.est <- drop(as.numeric(zcoef[,name.coef2])) * vs.list[[vfunc[i]]]
            beta.est$data <- apply(beta.est$data, 2, 
                                   sum)
            beta.est$names$main <- "beta.est"
            beta.est$data <- matrix(as.numeric(beta.est$data), 
                                    nrow = 1)
            if (basis.x[[vfunc[1]]]$type == "pls") {
              if (basis.x[[vfunc[1]]]$norm) {
                sd.X <- sqrt(apply(data[[vfunc[i]]]$data, 
                                   2, var))
                beta.est$data <- beta.est$data/sd.X
              }
            }
            beta.l[[vfunc[i]]] <- c(beta.l[[vfunc[i]]], 
                                    beta.est)
          }
        }
      }
      else {
        if (length(name.coef2)<basis.x[[vfunc[i]]]$harmonics$basis$nbasis)
          basis.x[[vfunc[i]]]$harmonics$basis$dropind<-(length(name.coef2)+1):basis.x[[vfunc[i]]]$harmonics$basis$nbasis
        beta.l[[vfunc[i]]]<-fd(z$coefficients[name.coef2,],basis.x[[vfunc[i]]]$harmonics$basis)
      }
    }
    # print(class(basis.y))
    if (bsp1 & inherits(basis.y,"basisfd")) {  
   # print("biiiiifffffffff")
      cc <- beta.l[[vfunc[i]]]$coefs
      # print("biiiiifffffffff1")
      # print(dim(cc))
      # print((basis.x[[vfunc[i]]]$nbasis))
      # print((basis.y$nbasis))
      # print(dim(zcoef[,name.coef[[vfunc[i]]],drop=F]))
      # bi <- bifd(t(cc[,]),  basis.x[[vfunc[i]]], basis.y)
      # print("biiiiifffffffff2")
      ## dd <- eval.bifd(data[[vfunc[i]]]$argvals,data[[response]]$argvals,bi) # modificar grid.fdata para que haga lo mismo
      # print("biiiiifffffffff3")
      #beta2d[[vfunc[i]]] <-  fdata(dd,list(data[[vfunc[i]]]$argvals,data[[response]]$argvals),fdata2d=TRUE)
      
      
      bi <-bifd(cc,  basis.x[[vfunc[i]]], basis.y)
      dd <- eval.bifd(data[[vfunc[i]]]$argvals,tty,bi) # modificar grid.fdata para que haga lo mismo
      beta2d[[vfunc[i]]] <-  fdata(dd,list(data[[vfunc[i]]]$argvals,tty),fdata2d=TRUE)
    } else{
      beta2d <- NULL
    #   print("parte por completar")  **********************
      # eta<-object$basis.x$[[vfunc[i]]]$coefs
      # coef <-object$coefficients[,-1]
      # dd <-eta %*% t(coef) %*% t(phi)
      # beta2d[[vfunc[i]]] <-  fdata(#zcoef[,name.coef[[vfunc[i]]]]
      # (beta.l[[vfunc[i]]]$coefs)
      #    ,list(data[[vfunc[i]]]$argvals,data[[response]]$argvals),fdata2d=TRUE)
    }
  }
  # if (class(basis.y)!="basisfd"){  
    if (!inherits(basis.y,"basisfd")) {    
    # print(basis.y.class)
    if (basis.y.class=="raw")
      rownames(zcoef) <- paste(response,".",tty,sep="")
    else 
      rownames(zcoef) <- paste(response,".",rownames(basis.y$basis$data),sep="")
    
    #sume<-colSums(z$residuals$data^2)}
    #else {
    #print(dim(z$residuals))
    #yhatmat = eval.fd(y$argvals, alphafd) %*% matrix(1, 1, ncurves)+   eval.fd(y$argvals,  xbetafd)
    #z$residuals = eval.fd(tty, z$residuals)
    #z$residuals<-fdata(t(z$residuals),tty,rtty,namy)
    #z$residuals <- fd(z$coefficients[name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]])
    #sume<-colSums(z$residuals$data^2)
  }  else       rownames(zcoef) <- paste(response,".",basis.y$names,sep="")

  # print("eeeeeeeeeeeeeeeeeeeeeeee")
  #z$sr2=fdata(sume/(n-df),tty,rtty,namy)
  #z$r2<-fdata(1-sume/colSums(yy$Xcen$data^2),tty,rtty,namy)
  #}
  # else { #caso resp escalar
  #  z$sr2<-sum(z$residuals^2)/z$df.residual
  #  ###### z$Vp=z$sr2*S
  #}
  z$H <-hatvalues(z)
  z$beta.l <- beta.l
  z$formula <- pf
  z$mean <- mean.list
  
  z$formula.ini <- formula
  z$basis.y <- basis.y
  z$basis.x <- basis.x
  z$y<-data[[response]]
  #z$JJ<-JJ
  z$data <- z$data
  z$XX <- XX
  z$vs.list <- vs.list
  eY <- z$y-mean(z$y)
  nullrss <- sum(norm.fdata(eY)^2)
  z$fitted.values <- out$fitted.values 
  z$residuals <- z$y- out$fitted.values
  
  #rss <- sum(norm.fdata(z$residuals)^2)
  ndf <- NROW(z$coefficients)
  # z$sr2 <- rss/z$df.residual
  # z$r2  <-  1 - rss/nullrss
  z$response <- response
  # class(z)<-c(class(z),"fregre.fr")
  z$beta2d <- beta2d
  z$rn <- FALSE
  z$bsp<-bsp1
  z$raw<-raw
  z$vfunc <-vfunc
  z$vnf <- vnf
  z$data <- data
  
  z$coefficients <- t(zcoef)
  
  z$basis.y.class <- basis.y.class
  class(z)<-c("fregre.lm.fr",class(z))
  z
}
# 
# res0=fregre.lm.fr(y~x,ldat2,  basis.x=basis.x) #OK
# 
# traceback()
# 
# basis.y<-create.pc.basis(ldat2$y,1:3)
# res0=fregre.lm.fr(y~x,ldat2,  basis.x=basis.x,basis.y=basis.y)
# res0$coefficients
# 
# res0=fregre.sam.fr(y~x,ldat2,  basis.x=basis.x,basis.y=basis.y)
# res0$coefficients
# 
# res0=fregre.lm.fr(y~x,ldat2,  basis.x=basis.x)
# res0$coefficients
# 
# res00=fregre.mlm.fr(y~x,ldat2,  basis.x=basis.x,basis.y=basis.y)
# tty<-ldat2$y$rangeval
# bsp <-create.bspline.basis(tty,11)
# res00=fregre.mlm.fr(y~x,ldat2,basis.y=bsp)
# dim(res00$beta2d$x$data)
# dim(res00$coefficients)
# dim(res0$coefficients)
# 
# eta<-res0$basis.x$x$coefs
# dim(eta)
# coef <-res0$coefficients[,-1]
# res0$coefficients[1:3,]
# coef <-res0$coefficients[,-1]
# t(res00$coefficients[,1:3])
# dim(coef)
# dd <-eta %*% t(coef)
# phi <- basis.y$coefs
# dim(phi)
# dim(dd)
# beta<-dd %*% t(phi)
# dim(beta)
# par(mfrow=c(1,1))
# image(beta)
# persp(beta,phi=43,theta=26)
# 
# 
# dim(res0$beta.l[[1]]$data)
# t(eta)%*%B%*%phi
# 
# Ahora lo que ocurre es que eta es la base de componentes de X 
# y phi la de Y. Si Y fuese RAW tendrías directamente (creo)
# el estimador t(eta)%*%B que debería ser una matriz nSxnT

