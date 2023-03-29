\name{getcoef.ff.interaction}
\alias{getcoef.ff.interaction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the estimated coefficient functions for function-on-function interaction model
}
\description{This function is used to calculate the estimates of \eqn{\mu(t), \beta_i(s,t), \gamma_{ij}(u,v,t)}{\mu(t), \beta_i(s,t), \gamma_{ij}(u,v,t)} for function-on-function interaction model (see the description in \code{\link{cv.ff.interaction}}) based on the output object of \code{\link{cv.ff.interaction}}, or  \code{\link{step.ff.interaction}}.
}
\usage{
getcoef.ff.interaction(fit.obj, t.x.coef=NULL, t.y.coef=NULL)
%\method{coef}{hd}(fit.cv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.obj}{the output object of \code{\link{cv.ff.interaction}}, or \code{\link{step.ff.interaction}}}.
  \item{t.x.coef}{a list of length \eqn{p} of vectors providing the observation time
    points of predictors on which coefficient functions will be
    evaluated.  If \code{t.x.coef}=NULL (default), \code{t.x} in
    \code{\link{cv.ff.interaction}} or \code{\link{step.ff.interaction}}
    will be used.}
  \item{t.y.coef}{a vector of observation time points of response
    function on which the coefficient functions will be evaluated. If
    \code{t.y.coef}=NULL (default), \code{t.y} in
    \code{\link{cv.ff.interaction}} or \code{\link{step.ff.interaction}}
    will be used.}
}
\value{a list providing the given or selected main effects and interactions, together with the corresponding estimated coefficient functions.
 \item{intercept }{the vector of estimated \eqn{\mu(t)} evaluated at the vector \code{t.y.coef} of the observation points for the response function \eqn{y(t)}.}
\item{main_effects}{the index vector of the input main_effects for \code{\link{cv.ff.interaction}} or the selected main effects by \code{\link{step.ff.interaction}}.}
\item{coef_main }{a list of matrices of the estimated values of the coefficient functions of main effect specified by \code{main_effects}. Each matrix gives the estimated values of \eqn{\beta_i(s,t)} at the two-dimensional grid created by the observation point vectors \code{t.x.coef[[i]]} and \code{t.y.coef}, where \eqn{i} is an index in  \code{main_effects}.}
\item{inter_effects}{a matrix of two columns showing the input interactions for \code{\link{cv.ff.interaction}} or the selected interactions by \code{\link{step.ff.interaction}}. Each row shows the indices of the pair of functional variables in an interaction or quadratic effect.}
\item{coef_inter}{a list of three-dimensional arrays  of estimated values of the coefficient functions of interaction or quadratic effects specified by \code{inter_effects}. Each array gives the estimated values of \eqn{\gamma_{ij}(u,v,t)} at the three-dimensional grid created by the observation point vectors \code{t.x.coef[[i]]}, \code{t.x.coef[[j]]}  and \code{t.y.coef}, where the pair \eqn{i, j} is in  \code{inter_effects}.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Xin Qi and Ruiyan Luo
%%  ~~who you are~~
}
\references{Ruiyan Luo and Xin Qi (2018) Interaction model and model
  selection for function-on-function regression, Journal of
  Computational and Graphical Statistics.
  \url{https://doi.org/10.1080/10618600.2018.1514310}
  }
\seealso{\link{cv.ff.interaction}, \link{step.ff.interaction}.}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{ #See the examples in cv.ff.interaction() and step.ff.interaction().}
