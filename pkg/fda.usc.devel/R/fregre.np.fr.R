#' @title Functional regression with functional response using non-parametric kernel
#' estimation
#' 
#' @description Computes functional regression between a functional explanatory variables and
#' functional response using kernel estimation. 
#' 
#' @details The non-parametric functional regression model can be written as follows \deqn{y_i =r(X_i)+\epsilon_i}{ y
#' = r(X) + \epsilon} where the unknown smooth real function \eqn{r} is
#' estimated using kernel estimation by means of
#' \deqn{\hat{r}(X)=\frac{\sum_{i=1}^{n}{K(h^{-1}d(X,X_{i}))y_{i}}}{\sum_{i=1}^{n}{K(h^{-1}d(X,X_{i}))}}}{\hat{r}(X)=(\sum_i
#' K(d(X,X_i))y_i/h) / (\sum_i K(d(X,X_i)/h)) i=1,...,n} where \eqn{K} is an
#' kernel function (see \code{Ker} argument), \code{h} is the smoothing
#' parameter and \eqn{d} is a metric or a semi-metric (see \code{metric}
#' argument).
#' 
#' The distance between curves is calculated using a metric (for instance, \code{\link{metric.lp}})
#' or a semimetric (see \code{\link{semimetric.basis}} or \code{\link{semimetric.NPFDA}} functions).
#' The kernel is applied to a metric or semi-metric that provides non-negative
#' values, so it is common to use asymmetric kernels. Different asymmetric
#' kernels can be used, see \code{\link{Kernel.asymmetric}}.\cr
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Functional response \code{\link{fdata}} class object.
#' @param h Bandwidth, \code{h>0}. Default argument values are provided as the
#' 5\%--quantile of the distance between \code{fdataobj} curves, see
#' \code{\link{h.default}}.
#' @param Ker Type of asymmetric kernel used, by default asymmetric normal
#' kernel.
#' @param metric List with components \code{x},\code{y} containing the metric 
#' to be applied to covariate and response, respectively.
#' @param par.metric List of components \code{x},\code{y} containing optional 
#' parameters for each metric in \code{metric}.
#' @param type.S Type of smothing matrix \code{S}. By default \code{S} is
#' calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).
#' @param par.S List of parameters for \code{type.S}: \code{w}, the weights.
#' @return Return:
#' \itemize{
#' \item {call}{ The matched call.} 
#' \item {fitted.values}{ Estimated scalar response.} 
#' \item {H}{ Hat matrix.} 
#' \item {residuals}{ \code{y} minus \code{fitted values}.} 
#' \item {df}{ The residual degrees of freedom.} 
#' \item {r2}{ Coefficient of determination.} 
#' \item {sr2}{ Residual variance.} 
#' \item {y}{ Response.} 
#' \item {fdataobj}{ Functional explanatory data.} 
#' \item {mdist}{ Distance matrix between \code{x} and \code{newx}.}
#' \item {Ker}{ Asymmetric kernel used.} 
#' \item {h.opt}{ Smoothing parameter or bandwidth.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.np.cv.fr}},
#' \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}} .\cr
#' Alternative method: \code{\link{fregre.basis.fr}}.
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York. \cr
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' ab=tecator[[1]]
#' ab2=fdata.deriv(ab,2)
#' ind=1:129
#' x=ab[ind]
#' y=ab2[ind]#' 
#' res.np=fregre.np.fr(x,y,Ker=AKer.epa)
#' summary(res.np)
#' }
#' 
#' @export
fregre.np.fr<-function(fdataobj,y,h=NULL,Ker=AKer.norm,metric=list(x=metric.lp,y=metric.lp),
par.metric=list(x=list(),y=list()), type.S=S.NW, par.S=list(w=1)){

if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
if (!is.fdata(y)) y=fdata(y)
nasx<-is.na.fdata(fdataobj)
nasy<-is.na.fdata(y)
if (is.null(names(y$data))) names(y$data)<-seq_len(nrow(y))
if (any(nasx) | any(nasy)) {
   bb<-!nasx & !nasy
   cat("Warning: ",sum(!bb)," curves with NA in covariate or response are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
                              
x<-fdataobj[["data"]]
#tt<-fdataobj[["argvals"]]
#rtt<-fdataobj[["rangeval"]]
C<-match.call()
m<-match(c("fdataobj","y","h","Ker","metric","par.metric","type.S","par.S"),names(C),0L)
n = nrow(x)
#np <- ncol(x)   
tty<-y$argvals
rtty<-y$rangeval
ny = nrow(y$data)
if (n != ny) stop("Different number of curves for covariate and response")
   

if (is.matrix(metric$x)) {
		mdistx<-metric$x
		metric$x=get(attributes(mdistx)$call)
		par.metric$x=attributes(mdistx)$par.metric
} else mdistx=do.call(metric$x,c(list(fdata1=fdataobj,fdata2=fdataobj),par.metric$x))
if (is.matrix(metric$y)) {
		mdisty<-metric$y
		metric$y=get(attributes(mdisty)$call)
		par.metric$y=attributes(mdisty)$par.metric
}
#ke<-deparse(substitute(Ker))
#if (!is.function(Ker)) Ker<-get(Ker)
ty<-deparse(substitute(type.S))
attr(par.S, "call") <- ty
#print(h)
if (is.null(h)) 
{
     if (is.character(Ker)) {
            nker <- function(u,mik=Ker){0.5*get(mik)(u)}
            Ker<-get(Ker)
#            ke=get(paste0("Ker.",unlist(strsplit(par.np[[namesx[i]]]$Ker,"[.]"))[2]))
            } else {
            nker <- function(u,mik=Ker){0.5*mik(u)}
            }
#    	nker=get(paste0("Ker.",unlist(strsplit(deparse(substitute(Ker)),"[.]"))[2]))
		h = do.call(h.default,c(list(fdataobj=fdataobj,metric=mdistx,prob=0.1,Ker=nker),par.metric$x))
}
    par.S$tt<-mdistx
    if (is.null(par.S$Ker))  par.S$Ker<-Ker
    if (is.null(par.S$h))  par.S$h<-h
    #if  (type.S=="S.KNN")  par.S$cv<-TRUE
    H=do.call(type.S,par.S)
    par.S$cv<-TRUE
    H.cv=do.call(type.S,par.S)
    df=trace.matrix(H)
    yp=H%*%y$data
    yp<-fdata(yp,tty,rtty,names=y$names)
    ypcv=H.cv %*% y$data
    rownames(yp$data)=rownames(y$data)
    e=y-yp
    ecv=y-ypcv
#    ecv2<-ecv$data %*% t(ecv$data)
    ycen=fdata.cen(y)$Xcen
    norm.e<-do.call("norm.fdata",c(list(fdataobj=e,metric=metric$y),par.metric$y))^2
    norm.ycen<-do.call("norm.fdata",c(list(fdataobj=ycen,metric=metric$y),par.metric$y))^2
    sr2=sum(norm.e)/(n-df)
 	 r2=1-sum(norm.e)/sum(norm.ycen)
    yp2=ecv$data%*%t(ecv$data)
out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df"=df,"r2"=r2,
"sr2"=sr2,"y"=y,"fdataobj"=fdataobj,"mdistx"=mdistx,"Ker"=Ker,var.y=yp2,
"metric"=metric,"par.metric"=par.metric,"type.S"=type.S,"par.S"=par.S,"h.opt"=h, 
fit.CV = list(fitted.values = ypcv, residuals = ecv))
     
class(out) <- "fregre.fr"
return(out)
}
