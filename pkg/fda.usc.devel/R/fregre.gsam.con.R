# Otra opcion, todos medidos en el mismo [a,b] pero las x's sparse
# paralelizar el for para que vaya + rapido
# modificar ayuda
# si h=0 modelo concurrente
# si h>0 se la matriz de disenyo se amplia (num.h veces por fila X un factor por columna)
# otra version es incluir un peso w en lugar del factor.

#' Fitting Functional Concurrent (Generalized) Spectral Additive Models 
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
# @param weights weights
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
#' @noRd 
#' @examples
#' \dontrun{
#' data(aemet)
#' class(aemet)<-c("ldata","list")
#' set.seed(0)
#' ii <- sample(nrow(aemet$df),60)
#' ltrain <- aemet[ii,row=T]
#' ltest <- aemet[-ii,row=T]
#' ff <- logprec~+s(temp,k=3)+s(wind.speed,k=3)#+s(latitude,k=3)
#' ff <- logprec~+s(temp,k=3)+s(latitude,k=3)
#' res0=fregre.gsam.con(ff,ltrain,h=0)
#' res1=fregre.gsam.con(ff,ltrain,h=1)
#' res2=fregre.gsam.con(ff,ltrain,h=2)
#' res3=fregre.gsam.con(ff,ltrain,h=3)
#' summary(res0$result[[151]])
#' summary(res2$result[[1]])
#' newy<-ltest$logprec
#' y <- ltrain$logprec
#' pred0 <- predict.fregre.gsam.con(res0,ltest)
#' pred1 <- predict.fregre.gsam.con(res1,ltest)
#' pred2 <- predict.fregre.gsam.con(res2,ltest)
#' pred3 <- predict.fregre.gsam.con(res3,ltest)
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
#' ldat=ldata("df"=dataf,"x"=x,"x.d1"=x.d1,"x.d2"=x.d2)
#' ltrain <- ldat[1:186,row=T]
#' ltest <- ldat[-(1:186),row=T]
#' res0=fregre.gsam.con(x.d2~+s(x),ltrain,h=0)
#' res1=fregre.gsam.con(x.d2~+s(x),ltrain,h=1)
#' res2=fregre.gsam.con(x.d2~+s(x),ltrain,h=2)
#' res3=fregre.gsam.con(x.d2~+s(x),ltrain,h=3)
#' res4=fregre.gsam.con(x.d2~+s(Water)+s(x),ltrain,h=0)
#' res5=fregre.gsam.con(x.d2~+s(Water)+s(x),ltrain,h=1)
#' res6=fregre.gsam.con(x.d2~+s(Water)+s(x),ltrain,h=2)
#' res7=fregre.gsam.con(x.d2~+s(Water)+s(x),ltrain,h=3) 
#' 
#' newy<-ltest$x.d2
#' y <- ltrain$x.d2
#' pred0 <- predict.fregre.gsam.con(res0,ltest)
#' pred1 <- predict.fregre.gsam.con(res1,ltest)
#' pred2 <- predict.fregre.gsam.con(res2,ltest)
#' pred3 <- predict.fregre.gsam.con(res3,ltest)
#' pred4 <- predict.fregre.gsam.con(res4,ltest)
#' pred5 <- predict.fregre.gsam.con(res5,ltest)
#' pred6 <- predict.fregre.gsam.con(res6,ltest)
#' pred7 <- predict.fregre.gsam.con(res7,ltest)
#' rest=sum(norm.fdata(newy -func.mean(y))^2)
#' 1-sum(norm.fdata(newy -pred0)^2)/rest; 1-sum(norm.fdata(newy -pred1)^2)/rest
#' 1-sum(norm.fdata(newy -pred2)^2)/rest; 1-sum(norm.fdata(newy -pred3)^2)/rest
#' 1-sum(norm.fdata(newy -pred4)^2)/rest; 1-sum(norm.fdata(newy -pred5)^2)/rest
#' 1-sum(norm.fdata(newy -pred6)^2)/rest; 1-sum(norm.fdata(newy -pred7)^2)/rest 
#' plot(newy -pred0,col=1)
#' lines(newy-pred1,col=2)
#' lines(newy-pred2,col=3)
#' }
#' 
# @export
fregre.gsam.con <- function (formula
                         #  , family = gaussian()
                           , data = list()#, weights = NULL
                           , h = 1,...) 
{
  ########## borrar
  # formula <- logprec~+s(temp,k=3)+s(latitude,k=3) #as.formula(x.d2~s(x.d1))
  #  data <- ltrain
  #family=gaussian() # warning
  # basis.x=basis.x
  # basis.y=basis.y
  ########## borrar
  # if (family!=gaussian() ) stop("Only implemented for gaussian family")
  # print(1)
  # family <- gaussian()
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
    # weights=NULL o ponderar despues?
    
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
  ######## manejando la respuesta  
  if (any(nam.data==response)) {
    isfdata<-is.fdata(data[[response]])
    if (isfdata) {
      #    print(paste("Functional response:",response))
      yfdata <-  data[[response]]
        y <- yfdata$data
        aux<-list("coefs" = NULL, "basis" = NULL)
        nam.y <- colnames(y)
        if (is.null(nam.y)) nam.y <-  yfdata$argvals
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
  if (ifunc)     k = 1
  Ymat<-y
  nh <- (1+h*2)
  H <- zfitted <- zresid <- matrix(NA,ndatos,npy)
  pcoef<-zcoef<-NULL
  result<-list()
  w <- NULL
  Xmat0<-data.frame(y[,1,drop=F])
  names(Xmat0) <- response                   
  lenvnf <- length(vnf)
  if (lenvnf>0) {
    Zmat<-data$df[,vnf,drop=F]
    Xmat0 <- (cbind(Xmat0,Zmat))
  } 
  hh <- list()
  Xmat <- Xmat0
  znf <-list()
 # print("nuevo"); print(h); print(pf)
  if (h>0) {
    Xmat0<-data.frame(y=rep(Xmat[,response],len=ndatos*nh))
    # w <- rep(c(1+0:(h),(h):1)/nh,each=ndatos)^2  # opcion con pesos
    names(Xmat0)<-response
    Xmat0[,"delta"] <- factor(rep(c(-h:(0),1:h),each=ndatos))
    pf <- paste(pf,"+ delta",collapse="")
    #print(pf)
  if (nnf > 0) {
    for (indnf in 1:lenvnf){
      Xmat0[[vnf[indnf]]]<-rep(Xmat[,vnf[indnf]],len=ndatos*nh)
    }
  }
  }
  pf
  #print("antes for npy")
  
  for (i in 1:npy){
    # i<-1
    #print("for npy")
    Xmat[,response] <- Ymat[,i,drop=F]
    pf0<-pf
   #print("antes h0")
      if (h == 0)         { 
        for (j in 1:lenfunc){
          nam.var <- paste(vfunc[j],".h",i,sep="",collpase="")
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
        
          Xmat0<-Xmat
        }
      } else{
        # hh <-1
        pf0 <-pf
        aux <- (i-h):(i+h)
        hh[[i]]<-aux[aux>0 & aux<=npy]
        # Xmat[,response] <- Ymat[,i,drop=F]
        nh <- length(hh[[i]])
        #print(nh)
        inn <- 1:length(c(Ymat[,hh[[i]]]))
        Xmat0[inn,response] <- c(Ymat[,hh[[i]]])
        Xmat["delta"] <- factor(rep(0,len=ndatos),levels=levels(Xmat0$delta))
        #print(table(Xmat$delta))
        #print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
        #print(dim(Xmat0));        print("antes del for")
        for (j in 1:lenfunc){
          #Xmat <- cbind(Xmat,data[[vfunc[j]]]$data[,hh])
          #nam.var <- paste(vfunc[j],".",i,"h",hh[[i]],sep="",collpase="")
          nam.var <- paste(vfunc[j],".h",i,sep="",collpase="")
          if (specials2[j] != "0") {
            nam.var0 <- paste(specials2[j], 
                              "(",nam.var,",k=", bs.dim2[j],")", 
                              sep = "",collapse="+")
          }     else {
            nam.var0 <- paste(nam.var,sep = "",collapse="+")
          }
        # Xmat <- cbind(Xmat,data[[vfunc[j]]]$data[,hh])
        #Xmat[nam.var] <- c(data[[vfunc[j]]]$data[,hh[[i]]])
          # print(hh[i]);          print(names(Xmat0))
        Xmat0[inn,nam.var]=c(data[[vfunc[j]]]$data[,hh[[i]]])
        Xmat[nam.var]=data[[vfunc[j]]]$data[,i] # para la prediccion
        pf0 <- paste(pf0, "+", nam.var0,collapse="")
        }  
      }
    # print("antes gam");    print(pf0);    print(names(Xmat0))
    # print(sapply(Xmat0,class))
   #z=gam(formula=as.formula(pf0),data=Xmat0, weights = w,...)
    z=gam(formula=as.formula(pf0),data=Xmat0,...)  
   # print("despues gam")
   result[[i]] <- z
   ss <- summary.sam(z)
   if (h == 0){
     pcoef <- rbind(pcoef, ss)
     zcoef <- rbind(zcoef, z$coefficients)
   } else{
     # print(ss);     print("ss");     print(pcoef)
     # con i =1 hh es mas corto que con i>h
     pcoef <- suppressWarnings(rbind(pcoef, ss))
     zcoef <- suppressWarnings(rbind(zcoef, z$coefficients))
   }
   # print("sale if")
   # print(dim(zfitted));   print(length(z$fitted.values));   print(nh)
   if (h>0){
   #ww <- matrix(rep(h,len=nh)/nh,ncol=1)
   # if (h==NCOL(tt))     ww <- matrix(c(1:(h+1),h:1)/nh,ncol=1)
     # print(head(Xmat))
    z$fitted.values <-  predict(z,Xmat)
    # print("bbbbbbbbbbbbbbbbbbbbbbbbbbbbb")
    tt<-matrix(z$hat,nrow=ndatos,ncol=nh)
    z$hat <-  (tt)%*% matrix(rep(h,len=nh)/nh,ncol=1)
   }
   zfitted[,i] <-z$fitted.values
   # zfitted[,i] <-z$fitted.values
   H[,i] <- z$hat
   
  }
  
  rownames(pcoef) <- nam.y
  rownames(zcoef) <- nam.y
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
  out$hh <- hh
  #out$w <-ww
  class(out) <- c("fregre.cgsam.fr",class(out))
 out
}

# https://cran.r-project.org/web/packages/fcr/index.html
# https://cran.r-project.org/web/packages/fcr/vignettes/dynamic-prediction.html
# https://rdrr.io/cran/fcr/man/fcr.html

# interpret.gam<-mgcv:::interpret.gam
# is.fdata<-fda.usc.devel:::is.fdata 
# summary.sam <- fda.usc.devel:::summary.sam

