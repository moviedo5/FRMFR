\name{cv.fof.spike}
\alias{cv.fof.spike}
%ff.spike.R
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for linear function-on-function regression on
  highly densely observed spiky data
%%  ~~function to do ... ~~
}
\description{This function is used to perform cross-validation and build the final model for highly densely observed spiky data using the signal compression approach for the following linear function-on-function regression model:
\deqn{Y(t)= \mu(t)+\sum_{i=1}^p\int_{a_i}^{b_i}X_i(s)\beta_i(s,t)ds+\varepsilon(t),}{Y(t)= \mu(t)+\int X_1(s)\beta_1(s,t)ds+...+\int X_p(s)\beta_p(s,t)ds+\epsilon(t),}
where \eqn{\mu(t)}{\mu(t)} is the intercept function. The \eqn{\{X_i(s),1\le i\le p\}}{{X_i(s),1\le i\le p}}
are \eqn{p} functional predictors and \eqn{\{\beta_i(s,t),1\le i\le
  p\}}{{\beta_i(s,t),1\le i\le p}} are their corresponding coefficient
functions, where \eqn{p} is a positive integer. The \eqn{\epsilon(t)}{\epsilon(t)} is the noise function.

We require that all the sample curves of each functional predictor are observed in a common dense grid of time points, but the grid can be different for different predictors. All the sample curves of the functional response are observed in a common dense grid.
}
\usage{
cv.fof.spike(X, Y, t.x, t.y, K.cv = 5, upper.comp = 10, thresh = 0.0001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of length \eqn{p}, the number of functional predictors. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of observation time points for  \eqn{X_i(s)}{X_i(s)}.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{the \eqn{n\times m}{n*m} data matrix for the functional response \eqn{Y(t)}{Y(t)}, where  \eqn{n}  is the sample size and  \eqn{m}  is the number of the observation time points for \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{Y} here~~
}
  \item{t.x}{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of observation time points of the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}.
%%     ~~Describe \code{t.x} here~~
}
  \item{t.y}{the vector of obesrvation time points of the functional response \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{t.y} here~~
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}
   \item{upper.comp}{the upper bound for the maximum number of components to be calculated. Default is 10.
%%     ~~Describe \code{upper.comp} here~~
}
  \item{thresh}{a number between 0 and 1 used to determine the maximum number of components we need to calculate. The maximum number is between one and the "upp.comp" above. The optimal number of components will be chosen between 1 and this maximum number, together with other tuning parameters by cross-validation. A smaller thresh value leads to a larger maximum number of components and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.0001.
%%     ~~Describe \code{thresh} here~~
}
}
\details{
  %Here we give the details for the model with \eqn{p=1}.
  We use the decomposition of the coefficient function:
\deqn{\beta_i(s,t)=\sum_{k=1}^\infty\psi_{ik}(s)w_k(t),  1\le i\le p}{\beta_i(s,t)=\psi_{i1}(s)w_1(t)+...+\psi_{iK}(s)w_K(t),   1\le i\le p,}
and estimate \eqn{\{\psi_{ik}(s), 1\le i\le p\}}{{\psi_{ik}(s), 1\le i\le p}} for each \eqn{k>0} by solving a generalized functional eigenvalue problem as given in \code{\link{cv.sigcom}} but with penalty replaced by
\deqn{\lambda\sum_{i=1}^p \sum_{j=0}^{J_1}\{e^{\tau(j-J_1)}2^{2\alpha
    j}||b_{ik,j}||^2+ e^{-\tau J_1} ||a_{ik,0}||^2\},}
where \eqn{a_{ik,0}} and \eqn{b_{ik,j}} \eqn{(0\le j \le J_1)} are vectors of
wavelet coefficients for \eqn{\psi_{ik}}.
Then we estimate \eqn{\{w_{k}(t), k>0\}}{{w_{k}(t), k>0}} by regressing
\eqn{Y(t)} on \eqn{\{\hat{z}_{1},... \hat{z}_{K}\}} using penalized
least square method, where \eqn{\hat{z}_{k}= \sum_{i=1}^p \int \{X_i(s)-\bar{X}_i(s)\}\hat{\psi}_{ik}(s)ds}.
%and then estimate \eqn{\{w_{k}(t), k>0\}}{{w_{k}(t), k>0}} by penalized least square with smoothness tuning parameter \eqn{\kappa}.
}
\value{An object of the ``cv.fof.spike'' class, which is used in the function \code{\link{pred.fof.spike}} for prediction. % and \code{\link{getcoef.fof.spike}} for extracting the estimated coefficient functions.
\item{mu}{a vector containing the estimated values of the intercept function at \code{t.y}.}
\item{Beta}{a list of \eqn{p} matrices, where the \eqn{i}-th matrix
  contains the estimated values of the slope coefficient function
  \eqn{beta_i(s,t)} at grid \code{t.x} * \code{t.y}.}
\item{c_fit_cv}{a list for internal use.}
}
\references{Xin Qi and Ruiyan Luo,  (manuscript) Functional regression for highly densely observed functional data with novel regularity.
%% ~put references to the literature/web site here ~
}
\author{Xin Qi and Ruiyan Luo
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

###########################################################################
# Example: spiky function-on-function regresion
###########################################################################
ptm <- proc.time()
library(FRegSigCom)
library(refund)
data(DTI)
I=which(is.na(apply(DTI$cca,1,mean)))
X=DTI$cca[-I,] # functional response
Y=DTI$rcst[-I,-(1:12)] #functional predictor
t.x <- list(seq(0,1,length=dim(X)[2]))
t.y <- seq(0,1,length=dim(Y)[2])
# randomly split all the observations into a training set with 200 observations
# and a test set.
train.id=sample(1:nrow(Y), 30)
X.train <- list(X[train.id,])
Y.train <- Y[train.id, ]
X.test <- list(X[-(train.id),])
Y.test <- Y[-(train.id), ]
fit.cv=cv.fof.spike(X.train, Y.train, t.x, t.y)
Y.pred=pred.fof.spike(fit.cv, X.test)
error<- mean((Y.pred-Y.test)^2)
print(c("prediction error=", error))

print(proc.time()-ptm)
}
