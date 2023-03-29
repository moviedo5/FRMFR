# pendiente
# h es lag para el delta de X (X(t_j)-X(t_{j-h}))
# paralelizar el for para que vaya + rapido
# modificar ayuda

#' Fitting Functional Concurrent (Generalized) Spectral Additive Models with functinal response
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
#' @param h \code{integer}, if \eqn{h=0} (by default) same discretization point for the functional
#'  response and the functional covariates is, for example \eqn{Y(t_j)=f(X(t_j))+\varepsilon}. 
#'  if \eqn{h>0}, the values observed in the ball are also used (j-h,j+h),
#'   for example \eqn{Y(t_j)=f(X(t_{j-1}),X(t_j),X(t_{j+1})+\varepsilon} with \eqn{h=1}.
#@param basis.y Basis for functional response. 
#@param basis.x List of basis for functional explanatory data estimation.
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
# \item {basis.x}{ Basis used for \code{fdata} covariates.} 
# \item {basis.y}{ Basis used for \code{fdata} response.} 
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
#' library(fda.usc.devel)
#' data(aemet)
#' class(aemet)<-c("ldata","list")
#' ltrain <- aemet[1:55,row=T]
#' ltest <- aemet[-(1:55),row=T]
#' ff <- logprec~+s(temp,k=3)+s(wind.speed,k=3)#+s(latitude,k=3)
#' ff <- logprec~+s(temp,k=3)+s(latitude,k=3)

#' names(ltrain)
#' res0=fregre.gsam.cl(ff,ltrain,h=0)
#' res1=fregre.gsam.cl(ff,ltrain,h=1)
#' res2=fregre.gsam.cl(ff,ltrain,h=2)
#' res3=fregre.gsam.cl(ff,ltrain,h=3)
#' summary(res0$result[[11]])
#' summary(res1$result[[11]])
#' summary(res2$result[[11]])
#' summary(res3$result[[11]])
#' 
#' newy<-ltest$logprec
#' y <- ltrain$logprec
#' pred0 <- predict.fregre.gsam.cl(res0,ltest)
#' pred1 <- predict.fregre.gsam.cl(res1,ltest)
#' pred2 <- predict.fregre.gsam.cl(res2,ltest)
#' pred3 <- predict.fregre.gsam.cl(res3,ltest)
#' rest=sum(norm.fdata(newy -func.mean(y))^2)
#' 1-sum(norm.fdata(newy -pred0)^2)/rest
#' 1-sum(norm.fdata(newy -pred1)^2)/rest
#' 1-sum(norm.fdata(newy -pred2)^2)/rest
#' 1-sum(norm.fdata(newy -pred3)^2)/rest
#' plot(newy -pred0,col=1)
#' lines(newy-pred1,col=2)
#' lines(newy-pred2,col=3)
#' lines(newy-pred3,col=4)
#' 
#' data(tecator)
#' x=tecator$absorp.fdata
#' x.d1<-fdata.deriv(x)
#' x.d2<-fdata.deriv(x.d1)
#' tt<-x[["argvals"]]
#' dataf=as.data.frame(tecator$y)
#' nbasis.x=19;nbasis.y=15
#' ldat=ldata("df"=dataf,"x"=x,"x.d1"=x.d1,"x.d2"=x.d2)
#' ltrain <- ldat[1:129,row=T]
#' ltest <- ldat[-(1:129),row=T]
#' res0=fregre.gsam.cl(x.d2~+s(x),ltrain,h=0)
#' res1=fregre.gsam.cl(x.d2~+s(x),ltrain,h=1)
#' res2=fregre.gsam.cl(x.d2~+s(x),ltrain,h=2)
#' res3=fregre.gsam.cl(x.d2~+s(x),ltrain,h=3)
#' newy<-ltest$x.d2
#' y <- ltrain$x.d2
#' pred0 <- predict.fregre.gsam.cl(res0,ltest)
#' pred1 <- predict.fregre.gsam.cl(res1,ltest)
#' pred2 <- predict.fregre.gsam.cl(res2,ltest)
#' pred3 <- predict.fregre.gsam.cl(res3,ltest)
#' rest=sum(norm.fdata(newy -func.mean(y))^2)
#' 1-sum(norm.fdata(newy -pred0)^2)/rest
#' 1-sum(norm.fdata(newy -pred1)^2)/rest
#' 1-sum(norm.fdata(newy -pred2)^2)/rest
#' 1-sum(norm.fdata(newy -pred3)^2)/rest
#' plot(newy -pred1,col=1)
#' lines(newy-pred2,col=2)
#' lines(newy-pred3,col=3)
#' }
#' @export
fregre.gsam.cl <- function (formula
                         #  , family = gaussian()
                           , data = list(), weights = NULL
                           , h = 1,...) 
{
  ########## borrar
  # formula <- as.formula(x.d2~s(x.d1))
  #  data <- ldat
  #family=gaussian() # warning
  # basis.x=basis.x
  # basis.y=basis.y
  ########## borrar
  # if (family!=gaussian() ) stop("Only implemented for gaussian family")
  # print(1)
  #family <- gaussian()
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
    z = gam(formula = gp$pf, data = data$df,...)
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
      
#      if (is.null(basis.y)) {
# print("base y nula")        
        y <- yfdata$data
        # bspy=0
        raw <- TRUE
        aux<-list("coefs" = NULL, "basis" = NULL)
        nam.y <- colnames(y)
        if (is.null(nam.y)) nam.y <-  yfdata$argvals
 #     }    
      # else {
      #   aux <- fdata2basis(yfdata,basis.y) # si PCA que vaya centrada! meanY.list()
      #   y <- aux$coefs
      #   if (basis.y$type == "pc" | basis.y$type == "pls") 
      #     bspy <- 1
      #   else bspy <- 2
      #   nam.y <- colnames(y)
      # }
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
          #fnf2 <- c(fnf2, 0)
          fnf2[i]<-1
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
          # fnf2 <- c(fnf2, 0)
          
          specials1 <- c(specials1, speci[i])
        }
        else {
          # print(6)
          if (any(gp$smooth.spec[[i]]$term == nam.func)) {
            # print(7)
            # vfunc <- c(vfunc, gp$smooth.spec[[i]]$margin[[1]]$term)
            bs.dim2 <- c(bs.dim2, gp$smooth.spec[[i]]$margin[[1]]$bs.dim)
            fnf <- c(fnf, 2)
            #fnf2 <- c(fnf2, 1)
            fnf2[i] <- 1
            specials2 <- c(specials2, speci[i])
          } #else  fnf2 <- c(fnf2, 0)
        }
      }      else {
        # print(8)
        if (any(gp$smooth.spec[[i]]$term == nam.df)) {
          #           print(9)
          #   vnf <- c(vnf, gp$smooth.spec[[i]]$term)
          bs.dim1 <- c(bs.dim1, gp$smooth.spec[[i]]$bs.dim)
          fnf <- c(fnf, 1)
          fnf1 <- c(fnf1, 1)
          #fnf2 <- c(fnf2, 1)
          specials1 <- c(specials1, speci[i])
        }        else {
          # print(10)
          if (any(gp$smooth.spec[[i]]$term == nam.func)) {
            #    vfunc <- c(vfunc, gp$smooth.spec[[i]]$term)
            bs.dim2 <- c(bs.dim2, gp$smooth.spec[[i]]$bs.dim)
            fnf <- c(fnf, 2)
            #fnf2 <- c(fnf2, 1)
            fnf2[i]<-1
            specials2 <- c(specials2, speci[i])
          } #else  fnf2 <- c(fnf2, 0)
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
  }
  
  Ymat<-y
  H <- zfitted <- zresid <- matrix(NA,ndatos,npy)
  pcoef<-zcoef<-NULL
  result<-list()
  
  Xmat0<-data.frame(y[,1,drop=F])
  names(Xmat0) <- response                   
  lenvnf <- length(vnf)
  if (lenvnf>0) {
    Zmat<-data$df[,vnf,drop=F]
    Xmat0 <- (cbind(Xmat0,Zmat))
  } 
  hh <- list()
  for (i in 1:npy){
    # i<-1
    Xmat <- Xmat0
    Xmat[,response] <- Ymat[,i,drop=F]
    pf0<-pf
   # cat("i",i,"h",h,"\n")
      if (h == 0 | i <= h)         { 
        #print("i<h o h=0")
        for (j in 1:lenfunc){
          nam.var <- paste(vfunc[j],".t",i,sep="",collpase="")
        if (specials2[j] != "0") {
            nam.var0 <- paste(specials2[j], 
                              "(",nam.var,",k=", bs.dim2[j],")", 
                              sep = "",collapse="+")
   
            }     else {
            nam.var0 <- paste(nam.var,sep = "",collapse="+")
            }
          # Xmat <- cbind(Xmat,data[[vfunc[j]]]$data[,hh])
          pf0 <- paste(pf0, "+", nam.var0,collapse="")
          Xmat[nam.var] <- data[[vfunc[j]]]$data[,i]
          #Xmat <- cbind(Xmat,data[[vfunc[j]]]$data[,h])
        }
      } else{
        # print("i>=h y h>0")
        # hh <-1
        pf0 <-pf
        # hh[[i]]<-i - h
        for (j in 1:lenfunc){
          #Xmat <- cbind(Xmat,data[[vfunc[j]]]$data[,hh])
          nam.var <- paste(vfunc[j],".t",i,sep="",collpase="")
          nam.var11 <- paste(vfunc[j],".t",i,"d",h,sep="",collpase="")
          if (specials2[j] != "0") {
             nam.var0 <- paste(specials2[j], 
                              "(",nam.var,",k=", bs.dim2[j],")", 
                              sep = "",collapse="+")
             nam.var1 <- paste(specials2[j], 
                               "(",nam.var11,",k=", bs.dim2[j],")", 
                               sep = "",collapse="+")
             
            }     else {
              nam.var0 <- paste(nam.var,sep = "",collapse="+")
              nam.var1 <- paste(nam.var11,sep = "",collapse="+")
            }
        pf0 <- paste(pf0, "+", nam.var0,"+",nam.var1,collapse="")
        Xmat[nam.var] <- data[[vfunc[j]]]$data[,i]
        Xmat[nam.var11] <- data[[vfunc[j]]]$data[,i] - data[[vfunc[j]]]$data[,i-h]
        }
      }
  # print(names(Xmat));  print(names(Xmat0))
  # print(pf0);  print(Xmat[1:3,])
    
   z=gam(formula=as.formula(pf0),data=Xmat,...)  
   result[[i]] <- z
   ss <- summary.sam(z)
   if (h == 0){
     pcoef <- rbind(pcoef, ss)
     zcoef <- rbind(zcoef, z$coefficients)
   } else{
     #print(ss);     print("ss");     print(pcoef)
     # con i =1 hh es mas corto que con i>h
     pcoef <- suppressWarnings(rbind(pcoef, ss))
     zcoef <- suppressWarnings(rbind(zcoef, z$coefficients))
   }
   zfitted[,i] <-z$fitted.values
   H[,i] <- z$hat
  }
  rownames(pcoef) <- nam.y
  rownames(zcoef) <- nam.y
#  print(bsp1)
  out <- list()
  out$fitted.values <- fdata(zfitted,tty,rtty,namy)
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
  #out$mean <- mean.list
  out$formula.ini=formula
  #out$basis.x=basis.x
  #out$basis.y=basis.y
  #out$basis <- aux$basis
  #out$basis.list <- basis.list
  out$data=data
  out$raw <- raw
  #out$bsp <- bsp1;  out$bspy <- bspy
  out$vfunc <- vfunc;   out$nnf <- nnf;   out$vnf <- vnf
  out$H <- H
  out$h <- h
  class(out) <- c("fregre.gsam.cl",class(out))
 out
}

# https://cran.r-project.org/web/packages/fcr/index.html
# https://cran.r-project.org/web/packages/fcr/vignettes/dynamic-prediction.html
# https://rdrr.io/cran/fcr/man/fcr.html
# comentar que hace cada una de las propuestas implementadas!!
