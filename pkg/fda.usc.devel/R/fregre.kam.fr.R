#' Functional Kernel Additive regression with functional response.
#' 
#' @description Computes Functional Kernel Additive Regression between functional explanatory variables and
#' a functional response 
#' 
#' @details The non-parametric functional regression model can be written as follows
#' \deqn{ y_i =r_1(X^{(1)}_i)+...+ r_p(X^{(p)}_i)+\epsilon_i } where the unknown smooth real functions
#' \eqn{r_j} are estimated using kernel estimation and backfitting: 
#' \deqn{\hat{r}(X_0)=\hat{\alpha}+\hat{r}_1(X^{(1)}_0)+...+\hat{r}_p(X^{(p)}_0)} with
#' \deqn{\hat{r_j}(X^{(j)}_0)=\frac{\sum_{i=1}^{n}{K(h_{j}^{-1}d_j(X^{(j)}_0,X^{(j)}_{i}))(y_{i}-\hat{y}^{(-j)}_{i})}}{\sum_{i=1}^{n}{K(h_{j}^{-1}d_j(X^{(j)}_0,X^{(j)}_i))}}}
#' where \eqn{K} is a kernel function (see \code{Ker} argument), \code{h} is the vector of the smoothing parameter
#' and \eqn{d_j} is a metric or a semi-metric (see \code{metric} argument).
#' 
#' The function estimates the value of smoothing parameter for each component (also called
#' bandwidth) \code{h} through Generalized Cross-validation \code{GCV} criterium (see \code{GCV.S.fr} or \code{CV.S.fr}) and
#' computing the distances using a metric or a semimetric (see, for instance, \code{\link{metric.lp}}). 
#' Different asymmetric kernels can be used, see \code{\link{Kernel.asymmetric}}.\cr
#' 
#' @param formula an object of class "\link{formula}": a symbolic description of the model to be fitted (as a linear model).
#' @param data an \code{ldata} object (a list with the collection of functional covariates and response) 
#' @param weights an optional vector of weights to be used in the estimation process. 
#' @param par.metric list with components with the names in \code{formula} containing the details of each metric/semimetric to be applied 
#' to each covariate/response. The first component of that list called \code{metric} is the name of the function to be applied, the rest are the possible parameters for that metric
#' @param par.np List with components with the names of the covariates like in \code{formula} containing the needed components
#' for the smoothing task: \code{Ker} kernel, \code{type.S} type of smoothing, \code{par.S} parameters for smoothing, \code{h} bandwidths. 
#' Values by default are provided when \code{par.np[[]]} is NULL.  
#' @param control List of components \code{maxit}:maximum number of iterations, \code{epsilon}: relative change in objective functions, 
#' \code{trace}: TRUE/FALSE to show intermediate results and \code{inverse}: "solve"/"svd" way of computing the inverse. 
#' @return Return:
#' \itemize{
#' \item \code{result}{ List with the output of each component.} 
#' \item \code{residuals}{ \code{y} minus \code{fitted values}.} 
#' \item \code{fitted.values}{ Estimated scalar response.}
#' \item \code{H}{ The hat matrix.}
#' \item \code{effects}{ List with the contribution to \code{fitted.values} of each component.}
#' \item \code{alpha}{ Estimation of the intercept.}
#' \item \code{metric}{ List with the metric employed with each component.}
#' \item \code{par.metric}{ Options for the metric for each component.}
#' \item \code{RSS}{ Functional Residual Sum of Squares.}
#' \item \code{null.RSS} { Functional Residual Sum of Squares of the null model: \eqn{y-\bar{y}}}
#' \item \code{sr2}{ Residual variance.} 
#' \item \code{df}{ The residual degrees of freedom computed from \code{H}} 
#' \item \code{iter}{ Number of iterations consumed.} 
#' \item \code{weights}{ weights}. 
#' \item \code{eqrank}{ Degrees of freedom consumed by each component.} 
#' \item \code{converged}{ Final status of the iterations.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.np.cv.fr}},
#' \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}} .\cr
#' Alternative method: \code{fregre.lm.fr} and \code{\link{fregre.np.cv.fr}}.
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples 
#' \dontrun{
#' N=100
#' T=51
#' S=31
#' tj = seq(0,1,len=T)
#' si = seq(0,5,len=S)
#' X1 = rproc2fdata(N,tj,sigma="OU")
#' X2 = rproc2fdata(N,tj,sigma="vexponential")
#' X3 = rproc2fdata(N,tj,sigma="brownian")
#' beta1<-fdata(tj^2+3*tj^4,argvals=tj)
#' beta2<-fdata(4*(tj-0.5)^2,argvals=tj)
#' beta3<-fdata(exp(-tj^2/2),argvals=tj)
#' beta4<-outer(si,tj,function(si,tj){exp(-5*abs(tj-si/5))})
#' beta5<-outer(si,tj,function(si,tj){3*tj*si/5+tj^2*exp(-si/5)})
#' beta6<-outer(si,tj,function(si,tj){-3*abs(tj-si/5)+3*tj*si/5})
#' expX1 = exp(X1)
#' X3d=X3
#' X3d$data=sweep(sweep(X3d$data+0.001,2,beta2$data,"*"),2,beta3$data,"+")
#' FX1=fdata(t(beta4%*%t(expX1$data)),si)
#' FX2=fdata(t(beta5%*%t(X2$data)),si)
#' FX3=fdata(t(beta6%*%t(X3d$data)),si)

#' epsilon = rproc2fdata(N,si,sigma="vexponential",par.list=list(scale=.1,theta=2))
#' yt=FX1+FX2+FX3+epsilon
#' y1=FX1+epsilon
#' y2=FX2+epsilon
#' y3=FX3+epsilon

#' ldata=ldata(df=data.frame(idx=1:nrow(yt)),yt=yt,y1=y1,y2=y2,y3=y3,X1=X1,X2=X2,X3=X3)
#' pmetric=list(yt=list(metric=metric.lp,lp=1),y1=list(metric=metric.lp,lp=1),
#'              y2=list(metric=metric.lp,lp=1),y3=list(metric=metric.lp,lp=1),
#'	X1=list(metric=metric.lp),X2=list(metric=metric.lp),X3=list(metric=metric.lp))
#' modt=fregre.kam.fr(yt~X1+X2+X3,data=ldata,par.metric=pmetric,control=list(trace=TRUE))
#' modt0=fregre.kam.fr(yt~X1+X2+X3,data=ldata,control=list(trace=TRUE))
#' 
#' modt1=fregre.kam.fr(y1~X1+X2+X3,data=ldata,par.metric=pmetric,control=list(trace=TRUE))
#' modt2=fregre.kam.fr(y2~X1+X2+X3,data=ldata,par.metric=pmetric,control=list(trace=TRUE))
#' modt3=fregre.kam.fr(y3~X1+X2+X3,data=ldata,par.metric=pmetric,control=list(trace=TRUE))

#' }
#' 
#' @export
fregre.kam.fr=function (formula, data, weights = rep(1,nobs), par.metric = NULL, 
	par.np = NULL, control = list(maxit = 100,epsilon = 0.001, trace = FALSE, inverse = "solve")) 
{
# Only allowed  functional covariates
#	 w<-weights
    tf <- terms.formula(formula)
    terms <- attr(tf, "term.labels")
    nt <- length(terms)
    if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    }
    else pf <- rf <- "~"
	 y <- data[[response]]
    ynames <- rownames(y$data)
	 nobs <- nrow(y)
    vtab <- rownames(attr(tf, "factors")) 
    vnf <- intersect(terms, names(data$df)) # Non functional
    vnf2 <- intersect(vtab[-1], names(data$df)[-1]) # Non functional
    vfunc2 <- setdiff(terms, vnf) 
    vint <- setdiff(terms, vtab)
    vfunc <- setdiff(vfunc2, vint) # Functionals 
    vnf <- c(vnf2, vint) #Non functionals
	 if (length(vnf)>0) { 
	 print(paste("There are",length(vnf), " non functional variates."))
	 stop("Only allowed functional covariates by now")}
	 if (control$trace) cat(paste0("Response:",response, "Cov.:",paste0(vfunc,collapse="/"),"\n"))
    if (attr(tf, "intercept") == 0) intercept <- FALSE
    else intercept <- TRUE
    if (is.null(control$maxit)) control$maxit <- 100
    if (is.null(control$epsilon)) control$epsilon <- 0.001
    if (is.null(control$trace)) control$trace <- FALSE
    if (is.null(control$inverse)) control$inverse <- "solve"
    if (is.null(data$df)) xlist <- data else xlist<-data[-1]

    eps <- control$epsilon
    namesx <- vfunc
    nvars <- length(vfunc)
    conv <- FALSE
 
    EMPTY <- nvars == 0
    unless.null <- function(x, if.null) {
        if (is.null(x)) 
            if.null
        else x
    }
	 X <- vector("list",nvars+intercept)
	 metric <- pmetric<-vector("list",nvars+1)
	 names(metric) <- c(namesx,response)
	 names(pmetric) <- c(namesx,response)
	 result <- vector("list",nvars)
	 names(result) <- namesx
    eqrank <- c(rep(0, nvars), 1)
	names(X) <- ifelse(intercept,c(namesx, "Intercept"),namesx)
    names(eqrank) <- ifelse(intercept, c(namesx, "Intercept"),namesx) 
#    par.np2 <- par.np
    if (intercept) X[[nvars + intercept]] <- gridfdata(rep(1,nobs), mean(y))
	 for (i in 1:nvars){
	 X[[i]] <- fdata(matrix(0,nrow=nobs,ncol=ncol(y)), argvals=y$argvals)
	 }
#    if (control$trace) cat("----Computing the distance matrices ----\n")
#     pmetric.def<-list(metric=metric.lp, lp=2)
    if (is.null(par.metric)){
       par.metric<-vector("list",nvars)
       names(par.metric)=namesx
       for (i in seq_len(nvars)){
        #   par.metric[[i]]$metric<-metric.lp(xlist[[namesx[i]]],xlist[[namesx[i]]],lp=2)
           par.metric[[i]]$metric<-metric.lp
           par.metric[[i]]$lp<-2
           
       } 
    }    
    if (is.null(par.np)){
       par.np<-vector("list",nvars)
       names(par.np)=namesx
       for (i in seq_len(nvars)){
           par.np[[i]]<-list(Ker = AKer.norm, type.S = "S.NW",par.S = list(w = weights))
       } 
    }    

	 for (i in 1:nvars) {
        if (is.null(par.metric[[namesx[i]]])) {
            metric[[i]] <- metric.lp(xlist[[namesx[i]]], xlist[[namesx[i]]])
			pmetric[[i]] <- list(metric=metric.lp, lp=2)
        }
        else {
		if (is.matrix(par.metric[[namesx[i]]]$metric)){
		    metric[[i]] <- par.metric[[namesx[i]]]$metric
		    pmetric[[i]] <- c(list(metric=get(attributes(par.metric[[namesx[i]]]$metric)$call), 
                                                    attributes(par.metric[[namesx[i]]])$par.metric))
			   } else  {
			metric[[i]] <- do.call(par.metric[[namesx[i]]]$metric,c(list(fdata1=xlist[[namesx[i]]]),par.metric[[namesx[i]]][-1]))
			pmetric[[i]] <- par.metric[[namesx[i]]]
			   }
        }

        if (is.null(par.np[[namesx[i]]])) {
            par.np[[namesx[i]]] = list(Ker = AKer.norm, type.S = "S.NW",par.S = list(w = weights))
        } else {
        if (is.null(par.np[[namesx[i]]]$Ker)) par.np[[namesx[i]]]$Ker <- AKer.norm
        if (is.null(par.np[[namesx[i]]]$type.S)) par.np[[namesx[i]]]$type.S <- "S.NW"
        }
        if (is.null(par.np[[namesx[i]]]$h)) {
            if (is.character(par.np[[namesx[i]]]$Ker)){
                nker <- function(u,mik=par.np[[namesx[i]]]$Ker){0.5*get(mik)(u)}
            } else {
                nker <- function(u,mik=par.np[[namesx[i]]]$Ker){0.5*mik(u)} 
            }
#            iker <- 1
#            nker[[iker]] <- function(u){.5*par.np2[[namesx[i]]]$Ker(u)}
#         if (is.character(par.np[[namesx[i]]]$Ker)) {
#            iker <- 2
#            nker[[iker]] <- function(u,mik=par.np[[namesx[i]]]$Ker){0.5*get(mik)(u)}
#            } else if (is.function(par.np[[namesx[i]]]$Ker)){
#              iker <- 3
#              nker[[iker]] <- function(u,mik=par.np[[namesx[i]]]$Ker){0.5*mik(u)}
#            } else {
#              iker <- 4
#              nker[[iker]] <- function(u,mik=par.np2[[namesx[i]]]$Ker){0.5*mik(u)} 
#            }
#            nker <- nker[[iker]] 
#            if (is.function(par.np2[[namesx[i]]]$Ker)) nker=par.np2[[namesx[i]]] else {
#            nker=get(paste0("Ker.",unlist(strsplit(deparse(substitute(par.np2[[namesx[i]]]$Ker)),"[.]"))[2]))}
            par.np[[namesx[i]]]$h <- c(h.default(xlist[[namesx[i]]], 
                prob = c(seq(0.025,0.20,len=41),0.25,0.3,0.4,0.5,0.6,0.7,0.8), metric = metric[[namesx[i]]],
                type.S=par.np[[namesx[i]]]$type.S, Ker=nker), Inf)
        }
    }
	if (is.null(par.metric[[response]])) {
		par.metric[[response]]$metric <- metric.lp
		par.metric[[response]]$lp <- 2
		}
	if (is.matrix(par.metric[[response]]$metric)){
	 metric[[nvars+1]] <- get(attributes(par.metric[[response]]$metric)$call)
	 pmetric[[nvars+1]] <- attributes(par.metric[[response]]$metric)$par.metric
	 } else {
	 metric[[nvars+1]] <- par.metric[[response]]$metric
	 pmetric[[nvars+1]] <- par.metric[[response]][-1]
	 }
    conv <- FALSE
	yhat <- Reduce('+',X)
	resids <- y-yhat
	eY <- y-mean(y)
    nullrss <- sum(do.call("norm.fdata",c(list(fdataobj=eY,metric=metric[[nvars+1]]),pmetric[[nvars+1]]))^2*weights)
    rss <- sum(do.call("norm.fdata",c(list(fdataobj=resids,metric=metric[[nvars+1]]),pmetric[[nvars+1]]))^2*weights)
    if (control$trace) cat("Inicio RSS:", nullrss, "\n")	
    names(weights) <- ynames
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
	 cambio <-  rep(Inf,length(X))
	 names(cambio) <- ifelse(intercept, c(namesx,"Intercept"), namesx)
    for (iter in 1L:control$maxit) {
        rssold <- rss
        Xold <- X
#	     if (intercept) X[[nvars + intercept]] = gridfdata(rep(1,nobs),mean(y-Reduce('+',X[1:nvars]))) 
        if (control$trace) cat("#------------------------------------------------\n")
        if (control$trace) cat("Iter:", iter, "/", nobs, "RSS:", rss, "\n")

        for (i in 1:nvars) {
            off <- Reduce('+',X[-i]) 
            z <- y - off
            offdf <- sum(eqrank[-i])
            xfunc <- xlist[[namesx[i]]]
            h <- par.np[[namesx[i]]]$h
            if (control$trace) cat(namesx[i],"/Range h:", range(h), length(h), "\n")
            Ker <- par.np[[namesx[i]]]$Ker
            type.S <- par.np[[namesx[i]]]$type.S
            parS <- par.np[[namesx[i]]]$par.S
            parS$w <- weights
            if (is.function(type.S)) {ty <- deparse(substitute(type.S))} 
            else {ty <- type.S}
				W=diag(parS$w)
#			   cat(paste0("Entering fregre.np.cv.fr variable:",namesx[i],"\n")) 
                res = fregre.np.cv.fr(xfunc, z, h = h, type.CV = "GCV.S.fr", 
                Ker = Ker, type.S = ty, par.S = parS,metric=list(x=metric[[i]],y=metric[[nvars+1]]), 
                par.metric=list(x=pmetric[[i]],y=pmetric[[nvars+1]]),
				par.CV = list(andf = offdf, W = W,metric=metric[[nvars+1]]))
            if (control$trace) cat("Var:", namesx[[i]], " h.opt:", res$h.opt, " df:", res$df, "\n")
            eqrank[i] <- res$df
            X[[i]] <- res$fitted.values
            result[[i]] <- res
        }
        if (intercept) X[[nvars + intercept]] = gridfdata(rep(1,nobs),mean(y-Reduce('+',X[1:nvars])))
		  yhat=Reduce('+',X)
		  resids=y-yhat
        rss <- sum(do.call("norm.fdata",c(list(fdataobj=resids,metric=metric[[nvars+1]]),pmetric[[nvars+1]]))^2*weights)
        if (control$trace) {
				par(mfrow=c(1,nvars+1))
            plot(resids,main="Functional Residuals")
            for (i in 1:nvars) {
				if (intercept) {phat=yhat-X[[i]]-X[[nvars+intercept]]} else {phat=yhat-X[[i]]}
                plot(phat, col = i+1, ylab = paste0("yhat - ",namesx[i]), 
                  xlab = "argvals", main = paste(namesx[i],"EqPar:", round(eqrank[i], 1)))
                  abline(h = 0)
            }
        }

		  for (i in 1L:nvars){
		  eX=X[[i]]-Xold[[i]]
		  cambio[i]=sum(do.call("norm.fdata",c(list(fdataobj=eX,metric=pmetric[[i]]$metric),pmetric[[i]][-1]))^2*weights)
		  }
#		  eX=X[[nvars+1]]-Xold[[nvars+1]]
          if (intercept) eX=X[[nvars+1]]-Xold[[nvars+1]] else eX=fdata(rep(0,length(y$argvals)),y$argvals)
		  cambio[nvars+1]=sum(do.call("norm.fdata",c(list(fdataobj=eX,metric=metric[[nvars+1]]),pmetric[[nvars+1]]))^2*weights)
        if (control$trace) {
            cat("Shift Iter:", iter, "EqRank:", sum(eqrank)," RSS:",rss,"/", nobs, "\n")
            print(paste(round(cambio,3),collapse="/"))
        }
        if (any(cambio[1:nvars] > control$epsilon)) {
            conv <- FALSE
        }
        else {
            conv <- TRUE
            break
        }
        if (control$trace) 
            cat("RSS =", rss, "RSSold =",rssold, "Iterations -",iter, "\n")
        if (abs(rss - rssold)/(0.1 + abs(rss)) < control$epsilon | rss < control$epsilon) {
            conv <- TRUE
            break
        }
    }
    if (!conv) warning("kgam.fit: algorithm did not converge", call. = FALSE)
    rownames(resids$data) <- ynames
	rownames(yhat$data) <- ynames
#    for (i in 1L:nvars){
#		  cambio[i]=sum(do.call("norm.fdata",c(list(fdataobj=X[[i]]-Xold[[i]],metric=pmetric[[i]]$metric),pmetric[[i]][-1]))^2*weights)
#		  }

    names(result) <- namesx
    H <- kgam.H(result, control$inverse)
	ndf <- fdata.trace(H)
    sr2 <- rss/(nobs - ndf)
	r2 <- 1 - rss/nullrss
     if (intercept) alpha=func.mean(X[[nvars+intercept]]) else alpha=fdata(rep(0,length(y$argvals)),y$argvals)
    res <- list(result = result, residuals = resids, fitted.values = yhat, H=H,
        effects = X, alpha = alpha, metric=metric, par.metric=pmetric, 
        RSS = rss, null.RSS = nullrss, sr2 = sr2, r2=r2, df=ndf, iter = iter, weights = weights, 
		eqrank = eqrank, converged = conv
		,formula=formula)
    class(res) <- "fregre.kam.fr"
    res
}
