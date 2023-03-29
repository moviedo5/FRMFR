#' Fitting Functional Generalized Spectral Additive Models with functinal response
#' 
#' @description Computes functional GAM model between functional covariate
#' \eqn{(X^1(t_1),\cdots,X^{q}(t_q))}{(X(t_1),...,X(t_q))} (and non functional
#' covariate \eqn{(Z^1,...,Z^p)}{(Z1,...,Zp)}) and functional response
#'  \eqn{Y(s)=\left[ Y^{1}(s_1),\cdots,Y^{p}(s_p)\right]}{Y(s)=[Y^1(s),...,Y^r(s)}.
#' 
#' This function is an extension of the functional generalized additve
#' regression model \code{\link{fregre.gsam}} where now the response is of functional nature.
#' 
# \deqn{E[Y|X,Z]) = \sum_{l=1}^{L}{f_{l}(\theta_l)} =g^{-1}(\alpha+\sum_{i=1}^{p}f_{i}(Z^{i})+\sum_{k=1}^{q}\sum_{j=1}^{k_q}{f_{j}^{k}(\xi_j^k)})}{E[Y|X,Z]=\sum_l\sum_{l=1}{f_j(\theta_l)}=g^{-1}(\alpha+\sum_i f_i(Z_{i})+\sum_k^q\sum_{j=1}^{k_q}{f_j^k(\xi_j^k)})} where
#' \deqn{E[Y|X,Z]) = \sum_{l=1}^{L}{f_{l}(\theta_l)}}{E[Y|X,Z]=\sum_l\sum_{l=1}{f_j(\theta_l)}} where
#' \eqn{\theta_l}{\xi_l}   is the coefficient of the basis function expansion of Y,  
#' \eqn{\xi_j^k}{\xi_j^k} is the coefficient of the basis function expansion of
#' \eqn{X^k}, (in PCA analysis \eqn{\xi_j^k}{\xi_j^k} is the score of the
#' \eqn{j}-functional PC of \eqn{X^k}) 
#' 
#'  \deqn{E[Y|X,Z]) = \alpha+\sum_{i=1}^{p}f_{i}(Z^{i})+\sum_{k=1}^{q}\sum_{j=1}^{k_q}{f_{j}^{k}(\xi_j^k)}}{E[Y|X,Z]=\alpha+\sum_i f_i(Z_{i})+\sum_k^q\sum_{j=1}^{k_q}{f_j^k(\xi_j^k)}}
#' 
#' 
#' Note that if \code{basis.y} is null, the coefficients of the basis are 
#' the raw values of the response observed in the grid (argvas of Y).
#' 
#' The smooth functions \eqn{f(\cdot)}{f(.)} can be added to the right hand
#' side of the formula to specify that the linear predictor depends on smooth
#' functions of predictors using smooth terms \code{\link{s}} and
#' \code{\link{te}} as in \code{\link{gam}} (or linear functionals of these as
#' \eqn{Z\beta} and \eqn{\big<X(t),\beta\big>}{< X(t),\beta(t) >} in
#' \code{\link{fregre.glm}}).
#' 
#' The first item in the \code{data} list is called \emph{"df"} and is a data
#' frame with the non functional explanatory variables, as
#' \code{\link{gam}}.\cr
#' 
#' Functional covariates and functional response of class \code{fdata} are introduced in
#' the following items in the \code{data} list.\cr
#' \code{basis.x} is a list of basis for represent each functional covariate. The basis object can be
#' created by the function: \code{\link{create.pc.basis}} or  \code{\link{create.fdata.basis}}.
#' 
#' \code{basis.y} is a  basis for represent the functional response.
#'  \itemize{
#' \item If it is \code{NULL} (by default) the response variable is not represented in any base. 
#' \item If it is fixed basis (\code{basisfd} class object) or data-driven basis (\code{fdata.comp}
#'  class for functional PCA and PLS basis)
#' }
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' details of model specification are given under \code{Details}.
# @param family a description of the error distribution and link function to
# be used in the model. This can be a character string naming a family
# function, a family function or the result of a call to a family function.
# (See \code{\link{family}} for details of family functions.)
#' @param data \code{list}, \code{ldata} or \code{mfdata} class object that containing the variables in the model.
#' @param weights weights
#' @param basis.y Basis for functional response. 
#' @param basis.x List of basis for functional explanatory data estimation.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @details A \code{\link{fregre.gsam}} model is adjusted for
#' each of the coefficients functions derived from the basis representation on a functional response.
#' 
#' @return Return \code{gam} object plus:
#' \itemize{
#' \item {fitted.values}{ \code{fdata} object with the tted model predictions of expected value for each datum.} 
#' \item {residuals}{ \code{fdata} object with the model residuals.} 
#' \item {coefficients}{\code{matrix}, with coefficients of each fitted
#'  model (by row) and parametric coefficients (by col).}
#' \item {edf.coef}{\code{matrix},  with estimated degrees of freedom 
#' for each model fitted model (by row) and parametric coefficients (by col).
#' the NA value is assigned to non-significant linear and smooth terms (p-value> 0.05).} 
#' \item{result}{\code{list}, with each fitted \code{\link{fregre.gsam}} fitted model
#' for each dimension of functional response basis representation.}
#' \item {basis.x}{ Basis used for \code{fdata} covariates.} 
#' \item {basis.y}{ Basis used for \code{fdata} response.} 
#' \item {data}{ List that containing the variables in the model.} 
#' \item {formula}{ formula.} 
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{predict.fregre.sam.fr}}.
#'  Alternative methods: \code{\link{fregre.lm.fr}}
#' and \code{\link{fregre.kam.fr}}.
#' @references Muller HG and Stadtmuller U. (2005). \emph{Generalized
#' functional linear models.} Ann. Statist.33 774-805.
#' 
#' Wood (2001) \emph{mgcv:GAMs and Generalized Ridge Regression for R}. R News
#' 1(2):20-25.
#' 
#' Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data
#' Analysis}, 2nd ed., Springer, New York.
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' x=tecator$absorp.fdata
#' x.d1<-fdata.deriv(x)
#' x.d2<-fdata.deriv(x.d1)
#' tt<-x[["argvals"]]
#' dataf=as.data.frame(tecator$y)
#' nbasis.x=19;nbasis.y=15
#' basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
#' basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.y)
#' basis.x=list("x"=basis1,"x.d1"=basis1)
#' basis.y=basis2
#' 
#' ldat=ldata("df"=dataf,"x"=x,"x.d1"=x.d1,"x.d2"=x.d2)
#' res=fregre.sam.fr(x.d2~+x+s(x.d1),ldat,
#' basis.x=basis.x,basis.y=basis.y)
#' image(res$edf.coef)
#' }
#' @export
fregre.sam.fr <- function (formula
                           #, family = gaussian()
                           , data = list(), weights = NULL, 
                            basis.y = NULL, basis.x = NULL, ...) 
{
  
  ########## borrar
  # formula <- as.formula(x.d2~s(x.d1))
  #  data <- ldata
  #  family=gaussian() # warning
  # basis.x=basis.x
  # basis.y=basis.y
  ########## borrar
  # if (family!=gaussian() ) stop("Only implemented for gaussian family")
  # print(1)
  family <- gaussian()
  nam.data <- names(data)
  nam.df <- names(data$df)
  nam.func <- setdiff(nam.data,"df")
# print(nam.data)  ;print(nam.df);print(nam.func)
  tf <- terms.formula(formula, specials = c("s", "te", "t2"))
  terms <- attr(tf, "term.labels")
  special <- attr(tf, "specials")
  nt <- length(terms)
  specials <- rep(NULL, nt)
  if (!is.null(special$s)) 
    specials[special$s - 1] <- "s"
  if (!is.null(special$te)) 
    specials[special$te - 1] <- "te"
  if (!is.null(special$t2)) 
    specials[special$t2 - 1] <- "t2"
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  }   else pf <-rf <- "~"
  vtab <- rownames(attr(tf, "fac"))
  gp <- interpret.gam(formula)
  len.smo <- length(gp$smooth.spec)
  nterms <- length(terms)
  speci <- NULL
  specials1 <- specials2 <- fnf1 <- fnf2 <- fnf <- bs.dim1 <- bs.dim2 <- vfunc2 <- vfunc <- vnf <- NULL
  func <- nf <- sm <- rep(0, nterms)
  names(func) <- names(nf) <- names(sm) <- terms
  ndata <- length(data) - 1
  # if (ndata > 0) {
  #   names.vfunc <- rep("", ndata)
  #   for (i in 1:ndata) names.vfunc[i] <- names(data)[i + 1]
  # }   else names.vfunc <- NULL
  nam.func <- setdiff(nam.data,response)
  
  covar <- gp$fake.names
  if (length(covar)==0) {
# print("aaaaaaaaaaaaaaaaaa")    
    z = gam(formula = gp$pf, data = data$df, family = family)
    z$data <- data
    z$formula.ini <-gp$pf
    z$formula <- gp$pf
    #z$nnf <- nnf
    class(z) <- c(class(z), "fregre.gsam")
    return(z)
  }

  vnf<-vfunc<-NULL
  for (i in 1:length(covar)){
    if (covar[i] %in% nam.df) vnf<-c(vnf,covar[i])
    if (covar[i] %in% nam.func) vfunc<-c(vfunc,covar[i])
  }
  nnf<-length(vnf)
  nfunc<-length(vfunc)
  bsp1 <-  TRUE
  bspy <-  TRUE
  raw <- FALSE

# cat("response ");print(response)
  # cat("scalar ");print(vnf)
  # cat("funcional ");print(vfunc)

#fnf2<-rep(0,nfunc)


  ######## manejando la respuesta  
  if (any(nam.data==response)) {
    isfdata<-is.fdata(data[[response]])
    if (isfdata) {
      #    print(paste("Functional response:",response))
      yfdata <-  data[[response]]
      
      if (is.null(basis.y)) {
# print("base y nula")        
        y <-yfdata$data
        bspy=0
        raw <- TRUE
        aux<-list("coefs"=NULL,"basis"=NULL)
        nam.y <- colnames(y)
        if (is.null(nam.y)) nam.y <-  yfdata$argvals
        
      }      else {
        aux <- fdata2basis(yfdata,basis.y) # si PCA que vaya centrada! meanY.list()
        y <- aux$coefs
        if (basis.y$type == "pc" | basis.y$type == "pls") 
          bspy <- 1
        else bspy <- 2
        nam.y <- colnames(y)
        #if (is.null(nam.y)) nam.y <-  yfdata$argvals
       # print(bsp1)
      }
      ndatos <- NROW(y)
      npy<-NCOL(y)
      
      #  XX[[response]] <- data[[response]][["data"]]
      tty <- yfdata[["argvals"]]
      rtty <- yfdata[["rangeval"]]
      namy <- yfdata[["names"]]
    }     else stop("Response must be of fdata class")
  }  else stop("Response not found in data object")     
  ########################################
  # print("len.smo")
  # print(len.smo)
  # print(specials)
#cat("names.vfunc");print(names.vfunc)
  if (len.smo == 0) {
    specials <- rep(NA, nterms)
    speci <- rep("0", nterms)
    gp$smooth.spec[1:nterms] <- NULL
    gp$smooth.spec[1:nterms]$term <- NULL
    fnf2<-rep(0,nterms)
    #    vnf <- terms
    fnf1 <- fnf <- nterms
  }  else {
    for (i in 1:nterms) if (!is.na(specials[i])) 
      speci <- c(speci, specials[i])
    for (i in 1:nterms) {
    #  cat("terms[i]")        ;      print(terms[i])        

      if (any(terms[i] == nam.df)) {
        #print(1)
        #        vnf <- c(vnf, terms[i])

        sm[i] <- nf[i] <- 1
        fnf1 <- c(fnf1, 0)
        bs.dim1 <- c(bs.dim1, 0)
        specials1 <- c(specials1, "0")
      }      else {
        #print(2)
        if (any(terms[i] == nam.func)) {
          # print(3)
          #          vfunc <- c(vfunc, terms[i])
          func[i] <- 1
          fnf2 <- c(fnf2, 0)
          bs.dim2 <- c(bs.dim2, 0)
          specials2 <- c(specials2, "0")
        }
      }
    }
    for (i in 1:len.smo) {
# cat("i");print(i)
# print(speci[i])      
      if (speci[i] != "s") {
        # print(4)
        if (any(gp$smooth.spec[[i]]$term == nam.df)) {
          #   print(5)
          #          vnf <- c(vnf, gp$smooth.spec[[i]]$margin[[1]]$term)
          bs.dim1 <- c(bs.dim1, gp$smooth.spec[[i]]$margin[[1]]$bs.dim)
          fnf <- c(fnf, 1)
          fnf1 <- c(fnf1, 1)
          fnf2 <- c(fnf2, 0)
          specials1 <- c(specials1, speci[i])
        }
        else {
          # print(6)
          if (any(gp$smooth.spec[[i]]$term == nam.func)) {
            # print(7)
            # vfunc <- c(vfunc, gp$smooth.spec[[i]]$margin[[1]]$term)
            bs.dim2 <- c(bs.dim2, gp$smooth.spec[[i]]$margin[[1]]$bs.dim)
            fnf <- c(fnf, 2)
            fnf2 <- c(fnf2, 1)
            #fnf2[i] <- 1
            specials2 <- c(specials2, speci[i])
          } else  fnf2 <- c(fnf2, 0)
        }
      }      else {
        # print(8)
        if (any(gp$smooth.spec[[i]]$term == nam.df)) {
          #           print(9)
          #   vnf <- c(vnf, gp$smooth.spec[[i]]$term)
          bs.dim1 <- c(bs.dim1, gp$smooth.spec[[i]]$bs.dim)
          fnf <- c(fnf, 1)
          fnf1 <- c(fnf1, 1)
          fnf2 <- c(fnf2, 1)
          specials1 <- c(specials1, speci[i])
        }        else {
          # print(10)
          if (any(gp$smooth.spec[[i]]$term == nam.func)) {
            #    vfunc <- c(vfunc, gp$smooth.spec[[i]]$term)
            bs.dim2 <- c(bs.dim2, gp$smooth.spec[[i]]$bs.dim)
            fnf <- c(fnf, 2)
            fnf2 <- c(fnf2, 1)
            specials2 <- c(specials2, speci[i])
          } else  fnf2 <- c(fnf2, 0)
        }
      }
    }
  }
  #nfunc <- sum(func)
  
  name.coef = nam = par.fregre = beta.l = list()
  kterms = 1
  if (nnf > 0) {
    XX = data[["df"]]
    if (attr(tf, "intercept") == 0) {
      pf <- paste(pf, -1, sep = "")
    }
    for (i in 1:nnf) {
      if (fnf1[i] == 1 & len.smo != 0) 
        sm1 <- TRUE
      else sm1 <- FALSE
      if (sm1) {
        pf <- paste(pf, "+", specials1[i], "(", vnf[i], 
                    ",k=", bs.dim1[i], ")", sep = "")
      }
      else pf <- paste(pf, "+", vnf[i], sep = "")
      kterms <- kterms + 1
    }
  }   else {
    XX <- NULL
    #   XX = data.frame(data[["df"]][, response])
    #   names(XX) = response
  }
  lenfunc <- length(vfunc) 
  ifunc <- lenfunc > 0
  mean.list = basis.list = list()
  if (ifunc) {
    k = 1

  for (i in 1:lenfunc ) {
      if (is(data[[vfunc[i]]], "fdata")) {
        tt <- data[[vfunc[i]]][["argvals"]]
        rtt <- data[[vfunc[i]]][["rangeval"]]
        fdat <- data[[vfunc[i]]]
        nms <- data[[vfunc[i]]]$names
        #dat <- data[[vfunc[i]]]$data
        if (is.null(basis.x[[vfunc[i]]])) 
          basis.x[[vfunc[i]]] <- 
          create.fdata.basis(fdat, l = 1:7) else
            if (basis.x[[vfunc[i]]]$type == "pc" | basis.x[[vfunc[i]]]$type == "pls") 
              bsp1 = FALSE
        
        xaux <- fdata2basis(data[[vfunc[i]]],basis.x[[vfunc[i]]])
        name.coef[[vfunc[i]]] <- colnames(xaux$coefs) <- paste(vfunc[i],".",colnames(xaux$coefs),sep="")
        Z <- xaux$coefs
        lencoef <- length(colnames(Z))
        XX = cbind(XX, Z)

        if (fnf2[i] == 1)    
          sm2 <- TRUE    else sm2 <- FALSE
        for (j in 1:lencoef) {
          if (sm2) {
            pf <- paste(pf, "+", specials2[i], "(", 
                        name.coef[[vfunc[i]]][j], ",k=", bs.dim2[i],")", sep = "")
          }     else pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
          kterms <- kterms + 1
        }       
        basis.list[[vfunc[i]]] <- xaux$basis
        # J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
        #   vs.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$basis
        if (bspy!=0) 
          mean.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$mean
        else {
          xcc <- fdata.cen(data[[vfunc[i]]])
          mean.list[[vfunc[i]]] = xcc[[2]]
        }
       }      else {
        stop("Please, enter functional covariate of fdata class object")
      }
    }
  }
  #par.fregre$formula=as.formula(pf)
  #par.fregre$data=XX
  formula=as.formula(pf)
  names(par.fregre)
  nx <- ncol(XX)
  Ymat<-y
  H <- zfitted <- zresid <- matrix(NA,ndatos,npy)
  pcoef<-zcoef<-NULL
  result<-list()
  Xmat<-data.frame(y[,1,drop=F])
  if (length(vnf)>0) {
    Xmat<-data.frame(Xmat,data$df[,vnf])
    names(Xmat)<-c(response,vnf)         
  }  else       names(Xmat) <- response                   
  # for (i in 2:length(XX)){
  #    XX1<-data.frame(XX1,XX[[i]])
  #  }
  Xmat <- (cbind(Xmat,XX))
  # print(class(zfitted))
  
  for (i in 1:npy){
    Xmat[,response] <- Ymat[,i]
    #        z=gam(formula=as.formula(par.fregre$formula),data=XX1,family=family,offset=rep(1,len=nrow(XX[[1]])))
    #if (missing(offset))   
    z=gam(formula=formula,data=Xmat,family=family)  
    # else   { descomentar ***** y missing(offset)
    #   off<-offset
    #   z=gam(formula=formula,data=Xmat,family=family,offset=off)
    # }
    result[[i]] <- z
    ss <- summary.sam(z)
    pcoef <- rbind(pcoef, ss)
    zcoef <- rbind(zcoef, z$coefficients)
    zfitted[,i] <-z$fitted.values
    H[,i] <- z$hat
  }
  rownames(pcoef) <- nam.y
  rownames(zcoef) <- nam.y
#  print(bsp1)
  out <- list()
  if (bspy==2) {
    zfitted <- fd(t(zfitted),basis.y)
    out$fitted.values <- fdata(zfitted,tty,rtty,namy)
    # Funciona pero se devuelve un promedio para que sea 
    # el mismo tipo de output que para el resto de modelos
    # H <- fd(t(H),basis.y)
    # H <- fdata(H,tty,rtty,namy)
  } else {
    if (inherits(basis.y,"fdata")){ #bspy==1
      # print(aux$mean$data)
      out$fitted.values <- gridfdata(zfitted,aux$basis,aux$mean)
      # print(dim(H))
      # Salen valores negativos
      # H <- gridfdata(H,aux$basis)

    }    else  { #bspy==0
# bspy=0 Raw basis.y
      if (bspy==0)      out$fitted.values <- fdata(zfitted,tty,rtty,namy)
       else out$fitted.values <- gridfdata(zfitted,aux$basis,aux$mean)
      # H <- fdata(H,tty,rtty,namy)
    }
  }
  # se devuelve un promedio
  H <-rowMeans(H)
  out$fitted.values$names$main <- "Fitted values"
  #z$fitted.values <- zfitted #fdata(zfitted,tty,rtty,namy)
  out$residuals <- yfdata - out$fitted.values 
  out$residuals$names$main <- "Residuals"  
  out$coefficients <- zcoef
  out$edf.coef <- pcoef
  out$result<-result
  out$formula <- pf
  out$mean <- mean.list
  out$formula.ini=formula
  out$basis.x=basis.x
  out$basis.y=basis.y
  out$data=data
  out$basis <- aux$basis
  #out$basis.list <- basis.list
  out$raw <- raw
  out$bsp <- bsp1
  out$bspy <- bspy
  out$vfunc <- vfunc;   out$nnf <- nnf;   out$vnf <- vnf
  out$H <- H
  class(out) <- c("fregre.sam.fr",class(out))
  # tabla con elementos significativos 
  # por fila los modelos y por columnas las mismas
  # componentes (edf si s p-valor por *)
  out
}


################## FIXED BASIS
# library(fda.usc.devel)
# data(tecator)
# x=tecator$absorp.fdata
# x.d1<-fdata.deriv(x)
# x.d2<-fdata.deriv(x.d1)
# tt<-x[["argvals"]]
# dataf=as.data.frame(tecator$y)
# nbasis.x=7;nbasis.y=5
# basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
# basis.x=list("x"=basis1,"x.d1"=basis1,"x.d2"=basis1)
# 
# nbasis.y <-25
# basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.y)
# basis.y=basis2
# 
# ldat=ldata("df"=dataf,"x"=x,"x.d1"=x.d1,"x.d2"=x.d2)
# res.bsp=fregre.sam.fr(x.d2~+x+s(x.d1),data=ldat,#family=gaussian(),
# basis.x=basis.x,basis.y=basis.y)
# plot(ldat$x.d2,col=1)
# lines(res.bsp$fitted.values,col=2)
# plot(res.bsp$residuals)
# 
# ################### RAW BASIS for response (basis.y=NULL)
# 
# res.raw=fregre.sam.fr(x.d2~+s(x.d1),data=ldata,#family=gaussian(),
#                    basis.x=basis.x,basis.y=NULL)
# par(mfrow=c(1,2))
# plot(ldata$x.d2,col=1)
# lines(res.raw$fitted.values,col=2)
# plot(res.raw$residuals)
# range(res.raw$residuals)
# range(res.raw$fitted.values)
# res.raw$mean
# 
# ################### PC BASIS for response (basis.y=NULL)
#  ypc <- create.pc.basis(ldata$x.d2,1:10)
#  pc0 <- create.pc.basis(ldata$x,1:1)
#  pc1 <- create.pc.basis(ldata$x.d1,1)
#  xpc <- list("x"=pc0,"x.d1"=pc1)
#  # fda.usc.devel:::create.mfdata.basis(ldata,1:3,type.basis="pc",class.out = "fdata")
#  res.pc=fregre.sam.fr(x.d2~+x+s(x.d1),data=ldata,#family=gaussian(),
#                     basis.x=xpc
#                     #,basis.y=ypc
#                     )
#  par(mfrow=c(1,2))
#  plot(ldata$x.d2,col=1)
#  lines(res.pc$fitted.values,col=2)
#  plot(res.pc$fitted.values-ldat$x.d2)
# lines(res.pc$residuals,col="blue")
# plot(Ycen,col=2)
# lines(res.pc$fitted.values-ldat$x.d2)
# 
# res.pc$rss.null
# res.pc$rss

############# scalar covariates
# res=fregre.sam.fr(x.d2~+s(Protein)+s(Fat)+s(x),ldata,family=gaussian(),
#                   basis.y=basis.y)


# hacer un foreach()?
# alias basis2fdata de gridfdata


# # lapply(res$result,summary)
# coefs <- res$coefficients
# coefs <- res$fitted
# yfdata <- ldata$x.d2

# if (!raw) 
# aa<-fd(t(res$fitted),basis.y)
# faa <- fdata(aa,yfdata$argvals,yfdata$rangeval)
# residuals <- yfdata - faa 
# res$coefficients

# plot(yfdata,col=1)
# lines(faa,col=4)

# inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
# dim(res$fitted)
#yy<-gridfdata(t(res$fitted),faa)
# lines(yy,col=3)

# coefs <- t(res$coefficients[2:10,])
# dim(coefs)
# basis2d2fdata(coefs,res$basis.x$x.d1,basis.y)

# cc <- res$fitted
# plot(gridfdata(t(cc),yfdata))

# bi <-bifd(cc,  basis.x[[vfunc[i]]], basis.y)
# dd <- eval.bifd(data[[vfunc[i]]]$argvals,data[[response]]$argvals,bi) # modificar grid.fdata para que haga lo mismo
# beta2d[[vfunc[i]]] <-  fdata(dd,list(data[[vfunc[i]]]$argvals,data[[response]]$argvals),fdata2d=TRUE)


# res1=fregre.sam.fr(x.d2~+s(x),data=mda,basis.y=basis.y)
# summary(res1$result[[1]])
# 
# res2=fregre.sam.fr(x.d2~+x,data=mda,basis.y=basis.y)
# summary(res2$result[[1]])
# 
# 
# mda=mfdata("x"=x,"x.d1"=x.d1,"x.d2"=x.d2)
# res=fregre.sam.fr(x.d2~+x+s(x.d1),data=mda,family=gaussian(),
#                     basis.x=basis.x,basis.y=basis.y)
# summary(res$result[[1]])
# 
# res=fregre.sam.fr(x.d2~s(Fat)+x+s(x.d1),data=ldata,basis.y=basis.y)
# summary(res$result[[1]])
# 
# res=fregre.sam.fr(x.d2~Fat+x+s(x.d1),data=ldata,basis.y=basis.y)
# summary(res$result[[1]])
# 
# res=fregre.sam.fr(x.d2~+Fat,data=ldata, basis.x=basis.x,basis.y=basis.y)
# summary(res$result[[1]])
# 
# 
# res=fregre.sam.fr(x.d2~+s(Fat),data=ldata, basis.x=basis.x,basis.y=basis.y)
# summary(res$result[[1]])
# traceback()


summary.sam <- function (object, alpha=0.05, ...) 
{
  ss <- summary(object)
  #aa <- c(ss$p.table[,4,drop=F],ss$s.table[,4,drop=F]) < alpha
  #ii<-ss$s.table[,4,drop=F] < alpha
  df <- as.numeric(ss$p.table[,4,drop=F]<alpha)
  df[df==0] <- NA
  edf <- ss$s.table[,1]
  edf[ss$s.table[,4]>alpha]<-NA
  bb <- c(df,edf)
  names(bb) <- c(rownames(ss$p.table),rownames(ss$s.table))
  round(bb,4)
}

# Ycen <- fdata.cen(yfdata)$Xcen
# Pasar código al summary(invisble=FALSE->devuelve texto,=TRUE devuelve una lista
# out$rss.null <- sum(norm.fdata(Ycen)^2)
# out$rss.null <- sum(norm.fdata(yfdata)^2)
# out$rss <-  1- (sum(norm.fdata(out$residuals)^2) /  out$rss.null )

