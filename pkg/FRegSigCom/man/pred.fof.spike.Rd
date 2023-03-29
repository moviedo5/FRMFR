\name{pred.fof.spike}
\alias{pred.fof.spike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for linear function-on-function regression on highly sensely observed spiky data
%%  ~~function to do ... ~~
}
\description{Make predition for functional response from the CV object obtained by \code{\link{cv.fof.spike}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pred.fof.spike(fit.cv, X.test)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.cv }{ the CV object obtained by \code{\link{cv.fof.spike}}.
%%     ~~Describe \code{fit.cv} here~~
}
  \item{X.test }{new observations for the functional predictors. It is a list of length \eqn{p},  the number of functional predcitors. Each element is the observed  matrix from a functional predictor, with rows repsenting observation vectors and columns corresponding to the observation time points.
%%     ~~Describe \code{X.test} here~~
}
%  \item{t.y.test}{a vector of observation time points where values of predicted response curves are to be calculated. If \code{t.y.test}=NULL (default), \code{t.y} in \code{\link{cv.sigcom}} will be used.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A matrix containing the predicted response for the new
  observations. The number of rows is equal to  the sample size of the
  new data set, and the number of columns is equal to the length of \code{t.y}.
%  \code{t.y.test} or \code{t.y} when \code{t.y.test}=NULL.
} 
\author{Xin Qi and Ruiyan Luo
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{cv.fof.spike}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{#See the examples in cv.fof.spike().
}
