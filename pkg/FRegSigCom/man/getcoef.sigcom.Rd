\name{getcoef.sigcom}
\alias{getcoef.sigcom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the estimated intercept and coefficient functions for linear function-on-function models
%%  ~~function to do ... ~~
}
\description{This function is used to calculate the estimates for \eqn{\mu(t), \alpha(t), \beta_i(s,t)}{\mu(t), \alpha(t), \beta_i(s,t)} based on the object obtained from \code{\link{cv.sigcom}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getcoef.sigcom(fit.obj, t.x.coef=NULL, t.y.coef=NULL)
}
\arguments{
  \item{fit.obj}{the object obtained from \code{\link{cv.sigcom}}.}
  \item{t.x.coef}{a list of length \eqn{p} of vectors providing the observation time
    points of predictors on which coefficient functions will be
    evaluated.  If \code{t.x.coef}=NULL (default), \code{t.x} in
    \code{\link{cv.sigcom}} will be used.}
  \item{t.y.coef}{a vector of observation time points of response
    function on which the coefficient functions will be evaluated. If
    \code{t.y.coef}=NULL (default), \code{t.y} in
    \code{\link{cv.sigcom}} will be used.}
}
\value{a list containing
 \item{mu }{the vector of the estimated values of \eqn{\mu(t)} at the grid of observation points of the response function (\eqn{t.y.coef} for \code{\link{cv.sigcom}}).}
\item{beta }{ a list of length \eqn{p}, the number of functional
  predictors. Its \eqn{i}-th component is a matrix of the estimated
  values of coefficient functions \eqn{\beta_i(s,t)} of functional
  predictors at the full grid of observaion time points created by arguments \code{t.x.coef[[i]]} and \code{t.y} for \code{\link{cv.sigcom}}. The columns correspond to different observation points for response varible (\eqn{t.y.coef}). }
}
\references{Ruiyan Luo and Xin Qi, (2017) Function-on-Function Linear Regression by Signal Compression, Journal of the American Statistical Association. \url{http://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1164053}
%% ~put references to the literature/web site here ~
}
\author{Ruiyan Luo and Xin Qi
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{cv.sigcom}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{#See the examples in cv.sigcom().
}
