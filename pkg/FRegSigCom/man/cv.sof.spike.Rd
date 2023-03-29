\name{cv.sof.spike}
\alias{cv.sof.spike}
%spiky.scalar.on.function.R
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for linear scalar-on-function regression for highly densely observed spiky functional data
%%  ~~function to do ... ~~
}
\description{This function is used to perform cross-validation and build the final model for highly densely observed spiky data using the signal compression approach for the following linear scalar-on-function regression model:
\deqn{Y= \mu+\sum_{i=1}^p\int_{a_i}^{b_i}X_i(s)\beta_i(s)ds+\varepsilon}{Y= \mu+\int X_1(s)\beta_1(s)ds+...+\int X_p(s)\beta_p(s)ds+\epsilon,}
where \eqn{\mu}{\mu} is the intercept. The \eqn{\{X_i(s),1\le i\le p\}}{{X_i(s),1\le i\le p}}
are \eqn{p} functional predictors and \eqn{\{\beta_i(s),1\le i\le
  p\}}{{\beta_i(s),1\le i\le p}} are their corresponding coefficient
functions, where \eqn{p} is a positive integer. The \eqn{\epsilon}{\epsilon} is the random noise.

We require that all the sample curves of each functional predictor are observed in a common dense grid of time points, but the grid can be different for different predictors. All the sample curves of the functional response are observed in a common dense grid.
}
\usage{
cv.sof.spike(X, Y, t.x, K.cv = 10, upper.level = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of length \eqn{p}, the number of functional predictors. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of observation time points for  \eqn{X_i(s)}{X_i(s)}.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{an  \eqn{n} dimensional vector of the observed values for the response, where  \eqn{n}  is the sample size.
%%     ~~Describe \code{Y} here~~
}
  \item{t.x}{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of obesrvation time points of the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}.
%%     ~~Describe \code{t.x} here~~
}
  \item{K.cv}{the number of CV folds. Default is 10.
}
\item{upper.level}{the upper bound of the maximum resolution level. The optimal maximum resolution level is chosen between 1 and "upper.level", together with other tuning parameters, by cross-validation.
}
}
\details{
This method uses wavelet basis to expand \eqn{X_i(s)} and
\eqn{\beta_i(s)}, (\eqn{1 \le i \le p}), and estimates the expansion coefficients of
\eqn{\beta_i(s)}'s by penalized least squares method with penalty
\deqn{\lambda\sum_{i=1}^p \{\sum_{j=0}^{J_1}\{2^{-2\alpha e^{-(j-\tau)/\alpha}}2^{2\alpha
    j}||b_{ij}||^2+ \kappa ||b_i||^2\}\},}
where \eqn{b_{ij}} denotes the vector of wavelet coefficient for
\eqn{\beta_i(s)} at the \eqn{j}th level, and \eqn{b_{i}} is the vector
concatenating all \eqn{b_{ij}}, \eqn{(0\le j \le J_1)}.
}
\value{An object of the ``cv.sof.spike'' class, which is used in the function \code{\link{pred.sof.spike}} for prediction.
\item{mu}{the estimated intercept.}
\item{coef}{a list of \eqn{p} vectors, where the \eqn{i}-th vector contains the estimated values of the slope coefficient function \eqn{beta_i(s)} at \eqn{t.x}.}
\item{...}{optimal tuning parameters}
}
\references{Xin Qi and Ruiyan Luo,  (manuscript) Functional regression for highly densely observed functional data with novel regularity.
}
\author{Xin Qi and Ruiyan Luo,
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


##########################################################################
# Example: scalar-on-function for highly-densely observed curves
##########################################################################


ptm <- proc.time()
library(FRegSigCom)
data(Pork)
X=Pork$X
Y=Pork$Y
ntrain=40 # in paper, we use 80 observations as training data
xtrange=c(0,1) # the range of t in x(t).
t.x.list=list(seq(0,1,length.out=dim(X)[2]))
train.index=sample(1:dim(X)[1], ntrain)
X.train <- X.test <- list()

X.train[[1]]=X[train.index,]
X.test[[1]]=X[-(train.index),]
Y.train <- Y[train.index]
Y.test <- Y[-(train.index)]

fit.cv=cv.sof.spike(X.train, Y.train, t.x.list)
Y.pred=pred.sof.spike(fit.cv, X.test)
pred.error=mean((Y.pred-Y.test)^2)

print(c("pred.error=",pred.error))

print(proc.time()-ptm)


}