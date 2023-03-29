\name{cv.sigcom}
\alias{cv.sigcom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for linear function-on-function regression
%%  ~~function to do ... ~~
}
\description{This function is used to perform cross-validation and build the final model using the signal compression approach for the following linear function-on-function regression model with or without scalar predictors:
\deqn{Y(t)= \mu(t)+Z^T\alpha(t)+\sum_{i=1}^p\int_{a_i}^{b_i}X_i(s)\beta_i(s,t)ds+\varepsilon(t)}{Y(t)= \mu(t)+Z'\alpha(t)+\int X_1(s)\beta_1(s,t)ds+...+\int X_p(s)\beta_p(s,t)ds+\epsilon(t),}
where \eqn{\mu(t)}{\mu(t)} is the intercept function, \eqn{Z}{Z} is a
multivariate predictor and \eqn{\alpha(t)}{\alpha(t)} is the vector of
corresponding coefficient functions. When \eqn{Z=}NULL, there is no
scalar predictor. The \eqn{\{X_i(s),1\le i\le p\}}{{X_i(s),1\le i\le p}}
are \eqn{p} functional predictors and \eqn{\{\beta_i(s,t),1\le i\le
  p\}}{{\beta_i(s,t),1\le i\le p}} are their corresponding coefficient
functions, where \eqn{p} is a positive integer. When \eqn{p} is large
(e.g., greater than 6), in addition to smoothness penalty, one may want to consider to impose sparsity penalty (see
\code{\link{cv.hd}}()). The \eqn{\epsilon(t)}{\epsilon(t)} is the noise function.

We require that all the sample curves of each functional predictor are observed in a common dense grid of time points, but the grid can be different for different predictors. All the sample curves of the functional response are observed in a common dense grid.

%This method estimates a special decomposition of the coefficient
%functions: \deqn{\beta_i(s,t)=\sum_{k=1}^\infty\psi_{ik}(s)w_k(t),
%1\le i\le p}{\beta_i(s,t)=\psi_{i1}(s)w_1(t)+...+\psi_{iK}(s)w_K(t),
%1\le i\le p.} We first estimate \eqn{\{\psi_{ik}(s), 1\le i\le
%p\}}{{\psi_{ik}(s), 1\le i\le p}} for each \eqn{k>0} by solving a
%generalized functional eigenvalue problem with penalty
%\deqn{\lambda\sum_{i=1}^p\{|| \psi_{ik}||^2+ \tau
%||\psi''_{ik}||^2\}}{\lambda{|| \psi_{1k}||^2+...+  ||\psi_{pk}^2||^2}
%+ \lambda\tau{|| \psi''_{1k}||^2+...+||\psi''_{pk}||^2},} and then estimate \eqn{w_{k}(t), k>0}{w_{k}(t), k>0} by penalized least square with smoothness tuning parameter \eqn{\kappa}.
}
\usage{
cv.sigcom(X, Y, t.x, t.y, Z = NULL, s.n.basis = 50, t.n.basis = 50,
           K.cv = 5, upper.comp = 20, thresh = 0.001,
           basis.type.x="Bspline", basis.type.y="Bspline")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of length \eqn{p}, the number of functional predictors. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of observation time points for  \eqn{X_i(s)}{X_i(s)}.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{the \eqn{n\times m}{n*m} data matrix for the functional response \eqn{Y(t)}{Y(t)}, where  \eqn{n}  is the sample size and  \eqn{m}  is the number of the observation time points for \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{Y} here~~
}
  \item{t.x}{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of obesrvation time points of the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}.
%%     ~~Describe \code{t.x} here~~
}
  \item{t.y}{the vector of obesrvation time points of the functional response \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{t.y} here~~
}
  \item{Z}{the \eqn{n\times q}{n*q} data matrix for multivariate scalar predictors, where  \eqn{n}  is the sample size and  \eqn{q}  is the number of the scalar predictors. Default is NULL, indicating no scalar predictors.
%%     ~~Describe \code{Z} here~~
}
  \item{s.n.basis}{the number of B-spline basis functions used for estimating the functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s}. Default is 50.
%%     ~~Describe \code{s.n.basis} here~~
}
  \item{t.n.basis}{the number of B-spline basis functions used for estimating the functions \eqn{w_{k}(t)}{w_{k}(t)'s}. Default is 50.
%%     ~~Describe \code{t.n.basis} here~~
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}
\item{upper.comp}{the upper bound for the maximum number of components to be calculated. Default is 20.
%%     ~~Describe \code{upper.comp} here~~
}
  \item{thresh}{a number between 0 and 1 used to determine the maximum number of components we need to calculate. The maximum number is between one and the "upp.comp" above. The optimal number of components will be chosen between 1 and this maximum number, together with other tuning parameters by cross-validation. A smaller thresh value leads to a larger maximum number of components and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.001.
%%     ~~Describe \code{thresh} here~~
}
   \item{basis.type.x}{the type of basis functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s} . Only "BSpline" (default) and "Fourier" are supported.
}
   \item{basis.type.y}{the type of basis functions \eqn{w_{k}(t)}{w_{k}(t)'s}. Only "BSpline" (default) and "Fourier" are supported.
}
}
\details{
We consider the model with \eqn{p=1} and \eqn{Z}=NULL. Details
for the model with scalar and multiple functional predictors can be
found in the reference and its supplementary material.
We consider a special decomposition of the coefficient function:
\deqn{\beta(s,t)=\sum_{k=1}^\infty\psi_{k}(s)w_k(t).}{\beta(s,t)=\psi_{1}(s)w_1(t)+...+\psi_{K}(s)w_K(t).}
Let \eqn{\{(x_1(s), y_1(t)),...,(x_n(s),y_n(t))\}} denote the samples. We first estimate \eqn{\psi_{k}(s)}{\psi_{k}(s)} for each \eqn{k>0} by solving a penalized generalized functional
eigenvalue problem
\deqn{max_{\psi} \frac{\int\int
    \psi(s)\hat{B}(s,s')\psi(s')dsds'}{\int\int \psi(s)\hat{\Sigma}(s,s')\psi(s')dsds'+P(\psi)}}
\deqn{  {\rm{s.t.}}\quad \int\int \psi(s) \hat{\Sigma}(s,s')\psi(s')dsds'=1}
\deqn{  {\rm{and}}\quad  \int\int \psi(s)
  \hat{\Sigma}(s,s')\psi_{k'}(s')dsds'=0 \quad {\rm{for}} \quad k'<k}
where \eqn{\hat{B}(s,s')=\sum_{\ell=1}^n\sum_{\ell'=1}^n
  \{x_{\ell}(s)-\bar{x}(s)\}\int\{y_{\ell}(t)-\bar{y}(t)\}\{y_{\ell'}(t)-\bar{y}(t)\}dt
  \{x_{\ell'}(s')-\bar{x}(s')\}/n^2}, \eqn{\hat{\Sigma}(s,s')=\sum_{\ell=1}^n
  \{x_{\ell}(s)-\bar{x}(s)\}\{x_{\ell}(s')-\bar{x}(s')\}/n},
and penalty
\deqn{P(\psi)=\lambda\{|| \psi_{k}||^2+ \tau ||\psi''_{k}||^2\}.}{\lambda{|| \psi_{k}||^2} + \lambda\tau{|| \psi''_{k}||^2}.}
Then we estimate \eqn{\{w_{k}(t), k>0\}}{{w_{k}(t), k>0}} by regressing
\eqn{\{y_{\ell}(t)\}_{\ell=1}^n} on \eqn{\{\hat{z}_{\ell
    1},... \hat{z}_{\ell K}\}_{\ell=1}^n} with penalty \eqn{\kappa \sum_{k=0}^K
  \|w''_k\|^2} tuned by the smoothness parameter
\eqn{\kappa}. Here \eqn{\hat{z}_{\ell k}= \int (x_{\ell}(s)-\bar{x}(s))\hat{\psi}_{k}(s)ds}.
}
\value{An object of the ``cv.sigcom'' class, which is used in the function \code{\link{pred.sigcom}} for prediction and \code{\link{getcoef.sigcom}} for extracting the estimated coefficient functions.
\item{t.x.list}{the input t.x.}
\item{t.y}{the input t.y.}
\item{Z}{the input Z}
\item{opt.index}{index of the optimal \eqn{lambda}.}
\item{opt.lambda}{optimal value for \eqn{lambda}, the parameter tuning the penalty on \eqn{{\psi_{ik}(s), 1\le i \le p}}.}
\item{opt.tau}{optimal value for \eqn{tau}, the smoothness parameter for \eqn{{\psi_{ik}(s), 1\le i \le p}}.}
\item{opt.kappa}{optimal value for \eqn{kappa}, the smoothness parameter for \eqn{w_k(t)}{w_k(t)}.}
\item{opt.K}{optimal number of components to be selected.}
\item{min.error}{minimum CV error.}
\item{errors}{list for CV errors. }
\item{x.smooth.params}{a list for internal use.}
\item{y.smooth.params}{a list for internal use.}
\item{fit.1}{a list for internal use.}
\item{is.null.Z}{a logic value indicating whether \eqn{Z} is null.}
}
\references{Ruiyan Luo and Xin Qi (2017) Function-on-Function Linear Regression by Signal Compression, Journal of the American Statistical Association.  112(518), 690-705. \url{ https://doi.org/10.1080/01621459.2016.1164053}
%% ~put references to the literature/web site here ~
}
\author{Ruiyan Luo and Xin Qi
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#########################################################################
# Example 1: linear function-on-function regresion with one predictor
#             curve and two scalar predictors
#########################################################################
ptm <- proc.time()
library(FRegSigCom)
library(refund)
data(DTI)
tmp=1*(DTI$sex=="male")
Z.0=cbind(DTI$case, tmp)

I=which(is.na(apply(DTI$cca,1,mean)))
Y=DTI$cca[-I,] # functional response
X=list(DTI$rcst[-I,-(1:12)]) #functional predictor
Z=Z.0[-I,] # scalar predictor

t.x <- list(seq(0,1,length=dim(X[[1]])[2]))
t.y <- seq(0,1,length=dim(Y)[2])
# randomly split all the observations into a training set with 200 observations
# and a test set.
train.id=sample(1:nrow(Y), 100)
X.train <- list(X[[1]][train.id,])
Y.train <- Y[train.id, ]
Z.train <- Z[train.id, ]
X.test <- list(X[[1]][-(train.id),])
Y.test <- Y[-(train.id), ]
Z.test <- Z[-(train.id), ]
fit.cv=cv.sigcom(X.train, Y.train, t.x, t.y, Z=Z.train)
Y.pred=pred.sigcom(fit.cv, X.test, Z.test=Z.test)
error<- mean((Y.pred-Y.test)^2)
print(c(" prediction error=", error))
coef.est.list=getcoef.sigcom(fit.cv)
mu.est=coef.est.list[[1]] # intercept curve
beta.est=coef.est.list[[2]] # coefficient functions of functional predictors
gamma.est=coef.est.list[[3]] # coefficient functions of scalar predictors
print(proc.time()-ptm)



#########################################################################
# Example 2: linear function-on-function regresion with four predictor curves
#########################################################################


ptm <- proc.time()
library(FRegSigCom)
data(ocean)
Y=ocean[[1]]
Y.train=Y[1:80,]
Y.test=Y[-(1:80),]
t.y=seq(0,1, length.out = ncol(Y))
X.list=list()
X.train.list=list()
X.test.list=list()
t.x.list=list()
for(i in 1:4)
{
  X.list[[i]]=ocean[[i+1]]
  X.train.list[[i]]=X.list[[i]][1:80,]
  X.test.list[[i]]=X.list[[i]][-(1:80),]
  t.x.list[[i]]=seq(0,1, length.out = ncol(X.list[[i]]))
}
fit.cv=cv.sigcom(X.train.list, Y.train, t.x.list, t.y)
Y.pred=pred.sigcom(fit.cv, X.test.list)
error<- mean((Y.pred-Y.test)^2)
print(c(" prediction error=", error))
coef.list=getcoef.sigcom(fit.cv)
mu.est=coef.list[[1]] # intercept curve
beta.est=coef.list[[2]] # coefficient functions of functional predictors
print(proc.time()-ptm)
}
