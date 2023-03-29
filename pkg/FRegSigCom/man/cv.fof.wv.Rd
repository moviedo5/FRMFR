\name{cv.fof.wv}
\alias{cv.fof.wv}
%- Also NEED an '\alias' for EACH other topic documented here.
%ff.wave.sigcom.R
\title{Cross-validation for wavelet-based linear function-on-function regression method
}
\description{This function performs cross-validation and builds
  the final model for the following linear
  function-on-function regression model:
\deqn{Y(t)= \mu(t)+\sum_{i=1}^p\int_{a_i}^{b_i}Z_i(s)\beta_i(s,t)ds+\varepsilon(t)}{Y(t)= \mu(t)+\int Z_1(s)\beta_1(s,t)ds+...+\int Z_p(s)\beta_p(s,t)ds+\epsilon(t),}
where \eqn{\mu(t)}{\mu(t)} is the intercept function. The \eqn{\{Z_i(s),1\le i\le p\}}{{Z_i(s),1\le i\le p}}
are \eqn{p} functional predictors and \eqn{\{\beta_i(s,t),1\le i\le
  p\}}{{\beta_i(s,t),1\le i\le p}} are their corresponding coefficient
functions, where \eqn{p} is a positive integer. The
\eqn{\epsilon(t)}{\epsilon(t)} is the noise function. This method first
applies the fast wavelet transformation (FWT) to the predictor curves,  and transforms the original model
to a function-on-scalar model with wavelet coefficients as scalar
predictors. Then applying a dimension reduction based on signal compression to
approximate the regression function,  we have a
function-on-scalar model with a small number of uncorrelated scalar
predictors. %Assume sparsity (for example, when \eqn{Z_i(s)}{Z_i(s)}'s
%are smooth or smooth except a few jumps, the wavelet coefficient ), and
%we impose a sparsity penalty to fit the function-on-scalar regression model.
}
\usage{
cv.fof.wv(X, Y, t.y, K.cv = 5, upp.comp=10, thresh=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the \eqn{n\times p}{n*p} matrx of the wavelet coefficients of the predictor curves, where \eqn{n} is the sample size and \eqn{q} is the total number of wavelet coefficients for original functional predictors.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{the \eqn{n\times m}{n*m} data matrix for the functional response \eqn{Y(t)}{Y(t)}, where  \eqn{n}  is the sample size and  \eqn{m}  is the number of the observation time points for \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{Y} here~~
}
 \item{t.y}{the vector of obesrvation time points of the functional response \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{t.y} here~~
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}
  \item{upp.comp}{the upper bound for the maximum number of components to be calculated. Default is 10.
%%     ~~Describe \code{upper.comp} here~~
}
  \item{thresh}{a number between 0 and 1 used to determine the maximum number of components we need to calculate. The maximum number is between one and the "upp.comp" above. The optimal number of components will be chosen between 1 and this maximum number, together with other tuning parameters by cross-validation. A smaller thresh value leads to a larger maximum number of components and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.01.
%%     ~~Describe \code{thresh} here~~
}
}
\details{
This method first expresses the functional predictors and \eqn{\beta_i(.,t)} using wavelet expansion. Let \eqn{\bold{X}} and \eqn{\bold{\beta}(t)}  denote the vectors of all concatenated wavelet coefficients. Then the original model is transformed to
\deqn{Y(t)=\mu(t)+\bold{X}^T \bold{\beta}(t)+\varepsilon(t).}
We use the decomposition \eqn{\bold{\beta}(t)=\sum_{k=1}^K \bold{\alpha}_k w_k(t)} based on the KL expansion of \eqn{\bold{X}^T \bold{\beta}(t)}, where \eqn{\bold{\alpha}_k}'s are vectors of the same length as \eqn{\bold{\beta}(t)}. We estimate \eqn{\bold{\alpha}_k} for each \eqn{k} by solving the panelized generalized eigenvalue problem
\deqn{max_{\bold{\alpha}} \frac{\bold{\alpha}^T\hat{\bold{B}}\bold{\alpha}}{ \bold{\alpha}\hat{\bold{\Sigma}}\bold{\alpha}+P(\bold{\alpha})}}
\deqn{ {\rm{ s.t. }}\quad \bold{\alpha} \hat{\bold{\Sigma}}\bold{\alpha}=1}
\deqn{ {\rm{ and }}\quad  \bold{\alpha} \hat{\bold{\Sigma}}\bold{\alpha}_{k'}=0 \quad{\rm{for}}\quad k'<k}
where \eqn{\hat{\bold{B}}=\sum_{\ell=1}^n\sum_{\ell'=1}^n
  \{x_{\ell}-\bar{x}\}\int\{y_{\ell}(t)-\bar{y}(t)\}\{y_{\ell'}(t)-\bar{y}(t)\}dt
  \{x_{\ell'}-\bar{x}\}^T/n^2}, \eqn{\hat{\bold{\Sigma}}=\sum_{\ell=1}^n
  \{x_{\ell}-\bar{x}\}\{x_{\ell}-\bar{x}\}^T/n},
and penalty
\deqn{P(\bold{\alpha})=\tau\{(1-\lambda)||\bold{\alpha}||_2^2+ \lambda ||\bold{\alpha}||_1^2\}.}
Then we estimate \eqn{w_{k}(t), k>0}{w_{k}(t), k>0} by regressing
\eqn{Y(t)} on \eqn{\{\hat{z}_{1},... \hat{z}_{K}\}} with penalty \eqn{\kappa \sum_{k=0}^K
  \|w''_k\|^2} tuned by the smoothness parameter
\eqn{\kappa}. Here \eqn{\hat{z}_{k}=(\bold{X}-\bar{\bold{X}})^T\hat{\bold{\alpha}}_{k}}.
}
\value{An object of the ``cv.fof.wv'' class, which is used in the function \code{\link{pred.fof.wv}} for prediction. % and \code{\link{getcoef.fof.wv}} for extracting the estimated coefficient functions.
\item{min.error}{minimum CV error.}
\item{X}{input X.}
\item{Y}{input Y.}
\item{errors}{list for CV errors. }
\item{opt.K}{optimal number of components to be selected.}
\item{opt.smooth}{optimal smooth tuning parameter \eqn{\kappa}.}
\item{min.tau}{optimal tuning parameter \eqn{\tau}.}
\item{min.lambda}{optimal tuning parameter \eqn{\lambda}.}
\item{params.set}{set of tuning parameters using in CV.}
\item{...}{other output for internal use.}
}
\references{Ruiyan Luo, Xin Qi and Yanhong Wang. (2016) Functional wavelet regression for function-on-function linear models. Electronic Journal of Statistics. 10(2): 3179-3216.
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

########################################################################
## Example: wavelet function-on-function regresion
#########################################################################
ptm <- proc.time()
library(FRegSigCom)
library(wavethresh)
library(refund)
data(DTI)

I=which(is.na(apply(DTI$cca,1,mean)))
Y=DTI$cca[-I,] # functional response
X=DTI$rcst[-I,21:52] #functional predictor
n.wv=5

diagmat <- diag(2^n.wv)
W.x <- diagmat
for(i in 1:2^n.wv){
  tmp <- wd(diagmat[i,])
  tmp.cof <- accessC(tmp, level=0)
  for(j in 0:(n.wv-1))
    tmp.cof <- c(tmp.cof, accessD(tmp, level=j))
  W.x[,i] <- tmp.cof
}
X.wv=X%*% t(W.x)

t.y <- seq(0,1,length=dim(Y)[2])
# randomly split all the observations into a training set with 200 observations
# and a test set.
train.id=sample(1:nrow(Y), 50)
X.wv.train <- X.wv[train.id,]
Y.train <- Y[train.id, ]
X.wv.test <- X.wv[-(train.id),]
Y.test <- Y[-(train.id), ]

fit.cv=cv.fof.wv(X.wv.train, Y.train, t.y, upp.comp=5) # use default upp.comp or larger
Y.pred=pred.fof.wv(fit.cv, X.wv.test)
error<- mean((Y.pred-Y.test)^2)
print(c(" prediction error=", error))
print(proc.time()-ptm)

}
