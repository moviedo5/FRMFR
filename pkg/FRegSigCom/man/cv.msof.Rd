\name{cv.msof}
\alias{cv.msof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for linear multivariate scalar-on-function regression
%%  ~~function to do ... ~~
}
\description{This function is used to perform cross-validation and build the final model using the signal compression approach for the following linear multivariate scalar-on-function regression model:
\deqn{\bold{Y}= \bold{\mu}+\sum_{i=1}^p\int_{a_i}^{b_i}X_i(s)\bold{\beta}_i(s)ds+\bold{\epsilon},}{Y= \mu+\int X_1(s)\beta_1(s)ds+...+\int X_p(s)\beta_p(s)ds+\epsilon,}
where \eqn{\bold{Y}}{Y} is an \eqn{m}-dimensional multivariate response variable, \eqn{\bold{\mu}}{\mu} is the  \eqn{m}-dimensional intercept vector. The \eqn{\{X_i(s),1\le i\le p\}}{{X_i(s),1\le i\le p}}
are \eqn{p} functional predictors and \eqn{\{\bold{\beta}_i(s),1\le i\le
  p\}}{{\beta_i(s),1\le i\le p}} are their corresponding \eqn{m}-dimensional vector of coefficient
functions, where \eqn{p} is a positive integer. The \eqn{\bold{\epsilon}}{\epsilon} is the random noise vector.

We require that all the sample curves of each functional predictor are observed in a common dense grid of time points, but the grid can be different for different predictors. All the sample curves of the functional response are observed in a common dense grid.
}
\usage{
cv.msof(X, Y, t.x.list, nbasis = 50, K.cv = 5, upper.comp = 10,
        thresh = 0.001)
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
  \item{nbasis}{the number of basis functions used for estimating the vector of functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s} (see the reference for details). Default is 50.
%%     ~~Describe \code{s.n.basis} here~~
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}
 \item{upper.comp}{the upper bound for the maximum number of components to be calculated. Default is 10.
%%     ~~Describe \code{upper.comp} here~~
}
  \item{thresh}{a number between 0 and 1 used to determine the maximum number of components we need to calculate. The maximum number is between one and the "upp.comp" above. The optimal number of components will be chosen between 1 and this maximum number, together with other tuning parameters by cross-validation. A smaller thresh value leads to a larger maximum number of components and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.001.
%%     ~~Describe \code{thresh} here~~
}
}
\details{
We use the decomposition \eqn{\bold{\beta}_i(s)=\sum_{k=1}^K
  \alpha_{ki}(s) \bold{w}_k}, \eqn{1\le i \le p}, based on the KL
expansion of \eqn{\sum_{i=1}^p \int X_i(s)\bold{\beta}_i(s)ds}. Let
\eqn{\bold{Y}_{\ell}=(Y_{\ell,1},...,Y_{\ell,m})^T} and
\eqn{\bold{X}_{\ell}(s)=(X_{\ell,1}(s), ...,X_{\ell,p}(s))^T}, \eqn{1
  \le \ell \le n}, denote \eqn{n} independent samples.  We estimate \eqn{\bold{\alpha}_k(s)=(\alpha_{k1}(s),...,\alpha_{kp}(s))^T} for each \eqn{k} by solving the panelized generalized functional eigenvalue problem
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
\deqn{P(\bold{\alpha})=\lambda\sum_{i=1}^p \{|| \alpha_i||^2+ \tau ||\alpha_i''||^2\}.}{\lambda{|| \alpha||^2} + \lambda\tau{|| \alpha''||^2}.}
Then we estimate \eqn{\{w_{k}(t), k>0\}}{w_{k}(t), k>0} by regressing
\eqn{\{\bold{Y}_{\ell}\}} on
\eqn{\{\hat{z}_{\ell,1},... \hat{z}_{\ell,K}\}} using least
square method. Here \eqn{\hat{z}_{\ell,k}= \int (\bold{X}_{\ell}(s)-\bar{\bold{X}}(s))^T\hat{\bold{\alpha}}_{k}(s)ds}.
}
\value{An object of the ``cv.msof'' class, which is used in the function   \code{\link{pred.msof}} for prediction.
\item{fitted_model}{a list containing information about fitted model.}
\item{is_Y_vector}{a logic value indicating whether \code{Y} is a
  vector.}
\item{Y}{input data \code{Y}.}
\item{x.smooth.params}{a list for internal use.}
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
# Example: multiple scalar-on-function regression
#########################################################################


ptm <- proc.time()
library(FRegSigCom)
data(corn)
X=corn$X
Y=corn$Y
ntrain=60 # in paper, we use 80 observations as training data
xtrange=c(0,1) # the range of t in x(t).
t.x.list=list(seq(0,1,length.out=ncol(X)))
train.index=sample(1:nrow(X), ntrain)
X.train <- X.test <- list()
X.train[[1]]=X[train.index,]
X.test[[1]]=X[-(train.index),]
Y.train <- Y[train.index,]
Y.test <- Y[-(train.index),]

fit.cv.1=cv.msof(X.train, Y.train, t.x.list)# the cv procedure for our method
Y.pred=pred.msof(fit.cv.1, X.test) # make prediction on the test data

pred.error=mean((Y.pred-Y.test)^2)
print(c("pred.error=",pred.error))

print(proc.time()-ptm)

}