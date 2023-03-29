#' @title   Lag a multivariate functional data (ldata or mfdata class object)
#' 
#' @description Compute a lagged version of a ldata object, shifting the time base back by a given number of observations.
#' 
# @details  
#' 
#' @aliases ldata.lagp
#' 
#' @param x data (ldata or mfdata  class)	
#' @param p vector of integers with he number of lags (in units of observations)
#' @param var.name not used yet
#' 
#' @return 
#' returns a \code{ldata} with the lagged data marked with se suffixes "_p" at the end of its name.
#' 
#' @author Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es} 
#' 
#' @seealso \code{\link{ldata}} amd \code{\link{subset.ldata}}.
#' 
# @keywords htest multivariate nonparametric
#' @examples
#' \dontrun{ 
#' # no example yet
#' library(fda.usc.devel)
#' data(poblenou)
#' names(poblenou)
#' x <- ldata(poblenou$df,nox=poblenou$nox)
#' x0 <- ldata.lag(x,lag=0)
#' x5 <- ldata.lag(x,lag=5)
#' dim(x0[[1]])
#' dim(x5[[1]])
#' xx0 <- ldata.lagp(x,0)
#' xx1 <- ldata.lagp(x,1)
#' xx2 <- ldata.lagp(x,2)
#' xx12 <- ldata.lagp(x,c(1,2))
#' xx12 <- ldata.lagp(x,c(1,12))
#' names(xx0);names(xx1);names(xx2);names(xx12)
#' x$df[1:4,];xx0$df[1:4,];xx1$df[1:4,];xx2$df[1:4,];xx12$df[1:4,]
#' sapply(x,nrow);sapply(xx0,nrow);sapply(xx1,nrow);sapply(xx2,nrow);sapply(xx12,nrow)
#' # poner un indice que sera quen que no se repita !!!
#' }
#' 
#' @rdname ldata.lagp
#' @export  ldata.lagp
#####################################################
ldata.lagp <- function(x, p = 2,  var.name){
   # p<-2
  nvar <- length(x)
  names.x <- names(x)
  if (missing(var.name)){
    var.name <- setdiff(names(x),"df")
  }
  nvar <- length(var.name)
  np <- length(p)
  maxp <- max(p)
  xnew <- ldata.lag(x,lag=0, from=(1+maxp))
  names.xnew <- names(xnew)
  idf <- which(names.x != "df")
  names(xnew)[idf] <- paste0(names(xnew)[idf],"_p",0)
  if (maxp > 0){
    cls <- class(xnew)
    for (i in p){
    # i<-p
    xaux <- ldata.lag(x,lag=i, from=(maxp-i+1))
    idf <- which(names(xaux) != "df")
    names(xaux)[idf] <- paste0(names.xnew[idf],"_p",i)
    class(xnew)<-"list"
    names(xaux[["df"]]) <- paste0(names(xaux[["df"]]) ,"_p",i)
    xnew$df <- cbind(xnew$df,xaux$df)
    xnew <- c(xnew,xaux[idf])
    }
   class(xnew)<-cls
    }
  xnew
}
#####################################################
# @export  ldata.lag # no exported yet
ldata.lag <- function(x,lag= 0, from = 1, lag.var = NULL){
  n <- NROW(x[[1]])
  # comprovar que (lag+from) < n 
  ilag <- (from):(n-lag)
  length(ilag)
  xlag <- x[ilag,row=T]
}
#####################################################

#####################################################
# ldata.lagp <- function(x, p = 2, index = NULL, var.name){
#   # p<-2
#   # var.name solo para varriables funcionales
#   if (is.ldata(x)){
#     if  (class(x)[1]=="list") class(x)<-c("ldata","list")
#   } else stop("x is not a ldata class object")
#   
#   nvar <- length(x)
#   names.x <- names(x)
#   if (missing(var.name)){
#     var.name <- setdiff(names(x),"df")
#   }
#   nvar <- length(var.name)
#   
#   # if (is.null(index))  index <- (from):(n-lag)
#   # else if (length(index)>n) stop("wrong index values") #index <- index[(from):(n-lag)]
#   
#   n <- NROW(x[[1]])
#   np <- length(p)
#   maxp <- max(p)
#   lag <- 0
#   # if (is.null(index))  index <- (from):(n-lag)
#   # else if (length(index)>n) stop("wrong index values") #index <- index[(from):(n-lag)]
#   xnew <- ldata.lag(x,lag=0, from=(1+maxp),index= index)
#   
#   names.xnew <- names(xnew)
#   idf <- which(names.x != "df")
#   names(xnew)[idf] <- paste0(names(xnew)[idf],"_p",0)
#   if (maxp > 0){
#     cls <- class(xnew)
#     for (i in p){
#       # i<-p
#       xaux <- ldata.lag(x,lag = i, from=(maxp-i+1))
#       idf <- which(names(xaux) != "df")
#       names(xaux)[idf] <- paste0(names.xnew[idf],"_p",i)
#       class(xnew)<-"list"
#       names(xaux[["df"]]) <- paste0(names(xaux[["df"]]) ,"_p",i)
#       xnew$df <- cbind(xnew$df,xaux$df)
#       xnew <- c(xnew,xaux[idf])
#     }
#     class(xnew)<-cls
#   }
#   xnew
# }
# #####################################################
# # @export  ldata.lag # no exported yet
# ldata.lag <- function(x,lag= 0, from = 1, index = NULL){
#   if (is.ldata(x)){
#     if  (class(x)[1]=="list") class(x)<-c("ldata","list")
#   } else stop("x is not a ldata class object")
#   n <- NROW(x[[1]])
#   # comprovar que (lag+from) < n
#   if (is.null(index))  index <- (from):(n-lag)
#   else if (length(index)>n) stop("wrong index values") #index <- index[(from):(n-lag)]
#   # print(n);  print(lag)
#     print(index)
#   xlag <- x[index,row=T]
# }
# #####################################################
# Time Series Lag Plots
# Plot time series against lagged versions of themselves. Helps visualizing auto-dependence even when auto-correlations vanish.
# ldeaths
# lag(ldeaths, 12) # starts one year earlier
# 
# 
# lag.plot(nhtemp, 8, diag.col = "forest green")
# lag.plot(nhtemp, 5, main = "Average Temperatures in New Haven")
# ## ask defaults to TRUE when we have more than one page:
# lag.plot(nhtemp, 6, layout = c(2,1), asp = NA,
#          main = "New Haven Temperatures", col.main = "blue")
# 
# # que sería un lag plot funcional (f(t)-f(t-1)), norm(f(t)),norm(f(t-1)),
# plot(aemet$temp[1:2]-aemet$temp[2:3])
# plot(aemet$temp[1:2]/aemet$temp[2:3])
# plot(norm.fdata(aemet$temp[1:20]),norm.fdata(aemet$temp[2:21]))
# 
# ## Multivariate (but non-stationary! ...)
# lag.plot(freeny.x, lags = 3)
# 
# ## no lines for long series :
# lag.plot(sqrt(sunspots), set = c(1:4, 9:12), pch = ".", col = "gold")
# class(aemet)
# # aemet1 <- ldata.lag(aemet,lag= 0, from = 1, index = NULL)
# 
# aemet1 <- ldata.lag(aemet,lag= 0, from = 3, index = 20:77)
# head(aemet1$df);tail(aemet1$df)
# 
# 
# library(fda.usc.devel)
# data(aemet)
# names(aemet$df)
# aemet$df$index<-1:nrow(aemet$df)
# aemet1<-ldata.lag(aemet,lag= 0, from = 1, index = NULL)
# 
# 
# 
# lag(1:22,3)
# index<-c(1:5,10,15,21,22:)
# index1<-lag(index)
# 
# 
# 
# library(dplyr)
# lead(1:5)
# lead(1:5, default = 6)
# x <- 1:5
# tibble(behind = lag(x), x, ahead = lead(x))
# 
# 
# scrambled <- slice_sample(tibble(year = 2000:2005, value = (0:5) ^ 2), prop = 1)
# 
# # If you want to look more rows behind or ahead, use `n`
# lag(ts(1:5),  1)
# lag(1:5, -2)
# 
# 
# install.packages("lagged")
# library(lagged)
# v_lagged <- Lagged(0:6)                           # 1d object
# m_lagged <- Lagged(matrix(1:12, nrow = 4))        # 2d object
# a_lagged <- Lagged(array(1:24, dim = c(4,3,2)))  # 3d object
# 
# 
# 
# v_lagged[1:20]