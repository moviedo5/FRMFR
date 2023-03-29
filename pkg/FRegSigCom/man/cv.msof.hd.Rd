\name{cv.msof.hd}
%hd_multi_scalar_on_function.R
\alias{cv.msof.hd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for high dimensional linear multivariate scalar-on-function regression
%%  ~~function to do ... ~~
}
\description{This function is used to perform cross-validation and build
  the final model using the signal compression approach for the
  following linear multivariate scalar-on-function regression model with
  a large number of functional predictors:
\deqn{\bold{Y}= \bold{\mu}+\sum_{i=1}^p\int_{a_i}^{b_i}X_i(s)\bold{\beta}_i(s)ds+\bold{\epsilon}}{Y= \mu+\int X_1(s)\beta_1(s)ds+...+\int X_p(s)\beta_p(s)ds+\epsilon,}
where \eqn{\bold{Y}} is an \eqn{m}-dimensional multivariate response variable, \eqn{\bold{\mu}}{\mu} is the \eqn{m}-dimensional intercept vector. The \eqn{\{X_i(s),1\le i\le p\}}{{X_i(s),1\le i\le p}}
are \eqn{p} functional predictors and \eqn{\{\bold{\beta}_i(s),1\le i\le
  p\}}{{\beta_i(s),1\le i\le p}} are their corresponding \eqn{m}-dimensional vector of coefficient
functions, where \eqn{p} is a positive integer. The \eqn{\bold{\epsilon}}{\epsilon} is the random noise vector.

We require that all the sample curves of each functional predictor are observed in a common dense grid of time points, but the grid can be different for different predictors. All the sample curves of the functional response are observed in a common dense grid.
}
\usage{
cv.msof.hd(X, Y, t.x.list, n.basis = 25, K.cv = 5, upper.comp = 10,
           thresh = 0.02)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of length \eqn{p}, the number of functional predictors. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of observation time points for  \eqn{X_i(s)}{X_i(s)}.
%%     ~~Describe \code{X} here~~
}
 \item{Y}{an \eqn{n\times q}{n*q} data matrix for the response \eqn{Y} or a \eqn{n}-dimensional vector if there is only one scalar response, where  \eqn{n}  is the sample size, and \eqn{q} is the number of scalar response variables.
%%     ~~Describe \code{Y} here~~
}
  \item{t.x.list}{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of observation time points of the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}.
%%     ~~Describe \code{t.x} here~~
}
  \item{n.basis}{the number of basis functions used for estimating the vector of functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s} (see the reference for details). Default is 50.
%%     ~~Describe \code{s.n.basis} here~~
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}
  \item{upper.comp}{the upper bound for the maximum number of components to be calculated. Default is 10.
%%     ~~Describe \code{upper.comp} here~~
}
  \item{thresh}{a number between 0 and 1 used to determine the maximum number of components we need to calculate. The maximum number is between one and the "upp.comp" above. The optimal number of components will be chosen between 1 and this maximum number, together with other tuning parameters by cross-validation. A smaller thresh value leads to a larger maximum number of components and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.02.
%%     ~~Describe \code{thresh} here~~
}
}
\details{
We use the decomposition \eqn{\bold{\beta}_i(s)=\sum_{k=1}^K
  \alpha_{ki}(s) \bold{w}_k}, \eqn{1\le i \le p}, based on the KL
expansion of \eqn{\sum_{i=1}^p \int X_i(s)\bold{\beta}_i(s)ds}. Let
\eqn{\bold{Y}_{\ell}=(Y_{\ell,1},...,Y_{\ell,m})^T} and
\eqn{\bold{X}_{\ell}(t)=(X_{\ell,1}(s), ...,X_{\ell,p}(s))^T}, \eqn{1
  \le \ell \le n}, denote \eqn{n} independent samples.  We estimate \eqn{\bold{\alpha}_k=(\alpha_{k1},...,\alpha_{kp})^T} for each \eqn{k} by solving the panelized generalized functional eigenvalue problem
\deqn{max_{\alpha} \frac{\int\int
    \bold{\alpha}(s)^T\hat{\bold{B}}(s,s')\bold{\alpha}(s')dsds'}{\int\int \bold{\alpha}(s)^T\hat{\bold{\Sigma}}(s,s')\bold{\alpha}(s')dsds'+P(\bold{\alpha})}}
\deqn{  {\rm{s.t.}} \quad \int\int \bold{\alpha}(s)^T \hat{\Sigma}(s,s')\bold{\alpha}(s')dsds'=1}
\deqn{  {\rm{and}} \quad  \int\int \bold{\alpha}(s)^T
  \hat{\bold{\Sigma}}(s,s')\bold{\alpha}_{k'}(s')dsds'=0 \quad {\rm{for}} \quad k'<k}
where \eqn{\hat{\bold{B}}(s,s')=\sum_{\ell=1}^n\sum_{\ell'=1}^n
  \{\bold{X}_{\ell}(s)-\bar{\bold{X}}(s)\}\{\bold{Y}_{\ell}-\bar{\bold{Y}}\}^T\{\bold{Y}_{\ell'}-\bar{\bold{Y}}\}
  \{\bold{X}_{\ell'}(s')-\bar{\bold{X}}(s')\}^T/n^2}, \eqn{\hat{\bold{\Sigma}}(s,s')=\sum_{\ell=1}^n
  \{\bold{X}_{\ell}(s)-\bar{\bold{X}}(s)\}\{\bold{X}_{\ell}(s')-\bar{\bold{X}}(s')\}^T/n},
and penalty
  \deqn{P(\alpha)=\tau \{(1-\lambda)\sum_{i=1}^p\|\alpha_i\|_{\eta}^2 +
    \lambda (\sum_{i=1}^p \|\alpha_i\|_{\eta})^2\}}
    where  \eqn{\|\alpha_i\|_{\eta}^2=\|\alpha_i\|^2+\eta\|\alpha_i''\|^2}.
Then we estimate \eqn{\{w_{k}(t), k>0\}}{{w_{k}(t), k>0}} by regressing
\eqn{\{\bold{Y}_{\ell}\}} on
\eqn{\{\hat{z}_{\ell,1},... \hat{z}_{\ell,K}\}} using least
square method. Here \eqn{\hat{z}_{\ell,k}= \int (\bold{X}_{\ell}(s)-\bar{\bold{X}}(s))^T\hat{\alpha}_{k}(s)ds}.
}
\value{An object of the ``cv.msof.hd'' class, which is used in the function    \code{\link{pred.msof.hd}} for prediction.
  \item{opt.K}{optimal number of components to be selected.}
  \item{opt.tau}{optimal value for \eqn{tau}.}
  \item{opt.lambda}{optimal value for \eqn{lambda}.}
  \item{opt.eta}{optimal value for \eqn{lambda}.}
  \item{min.error}{minimum CV error.}
  \item{errors}{list for CV errors. }
  \item{...}{for internal use}
}
\references{Ruiyan Luo and Xin Qi (Submitted)
%% ~put references to the literature/web site here ~
}
\author{ Ruiyan Luo and Xin Qi
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
# Example:  
#########################################################################

#toy example

ptm <- proc.time()
library(FRegSigCom)
data(air)
Y=matrix(0, nrow(air[[1]]), 3)
Y[,1]=apply(air[[1]][,1:8],1,mean)
Y[,2]=apply(air[[1]][,(1:8)+8],1,mean)
Y[,3]=apply(air[[1]][,(1:8)+16],1,mean)

X=list()
for(i in 1:(length(air)-1))
{
  X[[i]]=air[[i+1]]
}
 
ntrain=100 # in paper, we use 80 observations as training data
t.x=seq(0,1,length.out=ncol(X[[1]]))
train.index=sample(1:nrow(X[[1]]), ntrain)
X.train <- X.test <- list()
t.x.list=list()
for(i in 1:length(X))
{
  t.x.list[[i]]=t.x
  X.train[[i]]=X[[i]][train.index,]
  X.test[[i]]=X[[i]][-(train.index),]
}
Y.train <- Y[train.index,]
Y.test <- Y[-(train.index),]

fit.cv=cv.msof.hd(X.train, Y.train, t.x.list, upper.comp=5)
# in practice, use the default values (or larger) for
# "upper.comp" and "n.basis".

Y.pred=pred.msof.hd(fit.cv, X.test)
pred.error=mean((Y.pred-Y.test)^2)
print(c("pred.error=",pred.error))

print(proc.time()-ptm)


}