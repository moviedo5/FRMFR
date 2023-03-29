\name{pred.fof.wv}
\alias{pred.fof.wv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for linear function-on-function regression using signal compression
%%  ~~function to do ... ~~
}
\description{Make predition for functional response from the CV object obtained by \code{\link{cv.fof.wv}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pred.fof.wv(cv.obj, X.test)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cv.obj }{ the CV object obtained by \code{\link{cv.fof.wv}}.
%%     ~~Describe \code{fit.cv} here~~
}
  \item{X.test }{new observations for the functional predictors. It is a list of length \eqn{p},  the number of functional predcitors. Each element is the observed  matrix from a functional predictor, with rows repsenting observation vectors and columns corresponding to the observation time points.
%%     ~~Describe \code{X.test} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A matrix containing the predicted response for the new
  observations. The number of rows is equal to  the sample size of the
  new data set, and the number of columns is equal to the length of
   \code{t.y}.
}
\author{Ruiyan Luo and Xin Qi
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{cv.fof.wv}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{#See the examples in cv.fof.wv().
}
