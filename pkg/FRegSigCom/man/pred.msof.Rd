\name{pred.msof}
\alias{pred.msof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for linear multivariate scalar-on-function regression
%%  ~~function to do ... ~~
}
\description{Make predition for multivariate scalar response based on new observations of predictor curves from the CV object obtained by \code{\link{cv.msof}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pred.msof(fit.obj, X.test)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.obj}{ the CV object obtained by \code{\link{cv.msof}}.
%%     ~~Describe \code{fit.cv} here~~
}
 \item{X.test}{new observations of functional predictors. It is a list of length \eqn{p}, the number of functional predcitors. Each element is the observed  matrix from a functional predictor, with rows repsenting observation vectors and columns corresponding to the observation time points.
%%     ~~Describe \code{X.test} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A matrix which contains the predicted values of the multivariate
  response. Each row is the prediction for each new observation.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Ruiyan Luo and Xin Qi
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{cv.msof}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{#See the examples in cv.msof().}
