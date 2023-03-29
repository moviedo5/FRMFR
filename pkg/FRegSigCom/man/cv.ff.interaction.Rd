\name{cv.ff.interaction}
\alias{cv.ff.interaction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for function-on-function regression models with specified main effects and two-way
  interaction terms
%%  ~~function to do ... ~~
}
\description{This function is used to perform cross-validation and build
  the final model using the signal compression approach for the
  following linear function-on-function regression model with main
  effects, two-way interaction effects and quadratic effects. Let
  \eqn{\{X_i(s), 1\le i\le p\}}{{X_i(s), 1\le i \le p}} be \eqn{p} potential functional predictors. The model is given by
\deqn{Y(t)= \mu(t)+\sum_{i \in M}\int_{a_i}^{b_i} X_i(s)\beta_i(s,t)ds+\sum_{(i,j) \in I}\int_{a_i}^{b_i}\int_{a_j}^{b_j} X_i(u)X_j(v)\gamma_{ij}(u,v,t)dudv+\epsilon(t),}{Y(t)= \mu(t)+\sum_{i \in M}\int_{a_i}^{b_i} X_i(s)\beta_i(s,t)ds+\sum_{(i,j) \in I}\int_{a_i}^{b_i}\int_{a_j}^{b_j} X_i(u)X_j(v)\gamma_{ij}(u,v,t)dudv+\epsilon(t),}
where \eqn{\mu(t)}{\mu(t)} is the intercept function. The index set
\eqn{M} of main effects is a subset of \eqn{\{1,...,p\}}{{1,...,p}}, and the index
set \eqn{I} of interactions and quadratic effects is a subset of the
collection of all possible pairs \eqn{\{(i,j), 1\le i\le j\le
    p\}}{{(i,j), 1\le i\le j\le p}}. The \eqn{\{\beta_i(s,t), i \in
  M\}}{{\beta_i(s,t), i\in M}} and  \eqn{\{\gamma_{ij}(u,v,t),(i,j)\in I\}}{{\gamma_{ij}(u,v,t),(i,j)\in I}} are the corresponding coefficient functions. The \eqn{\epsilon(t)}{\epsilon(t)} is the noise function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cv.ff.interaction(X, Y, t.x, t.y, main.effect, interaction.effect=NULL,
       adaptive=FALSE, s.n.basis=40, t.n.basis=40, inter.n.basis=20,
       basis.type.x="Bspline", basis.type.y="Bspline", K.cv=5, 
       upper.comp=8, thresh=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of \eqn{p} potential functional predictors. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th potential functional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of observation time points for  \eqn{X_i(s)}{X_i(s)}.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{the \eqn{n\times m}{n*m} data matrix for the functional response \eqn{Y(t)}{Y(t)}, where  \eqn{n}  is the sample size and  \eqn{m}  is the number of the observation time points for \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{Y} here~~
}
  \item{t.x}{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of observation time points of the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}.
%%     ~~Describe \code{t.x} here~~
}
  \item{t.y}{the vector of observation time points of the functional response \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{t.y} here~~
}
  \item{main.effect}{a vector of indices for main effects. It
    is a subset of \eqn{\{1,2,...,p\}}{{1,2,...,p}}.
}
  \item{interaction.effect}{a matrix of two columns.  Each row of this matrix specifies the index of a two-way interaction or a quadratic effect. Default is NULL.
}
  \item{adaptive}{a logic value indicating whether using adaptive
    penalty that has different smoothness tuning parameters for
    different target functions (see Details). Default is FALSE.
}
  \item{s.n.basis}{the number of basis functions used for estimating the functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s} (see Details). Default is 40.
%%     ~~Describe \code{s.n.basis} here~~
}
  \item{t.n.basis}{the number of basis functions used for estimating the functions \eqn{w_{k}(t)}{w_{k}(t)'s}. Default is 40.
%%     ~~Describe \code{t.n.basis} here~~
}
  \item{inter.n.basis}{the number of one-dimensional basis functions
    used to construct the tensor product basis functions for estimating the functions
    \eqn{\phi_{ijk}(u,v)}{\phi_{ijk}(u,v)'s}. Default is 20.
%%     ~~Describe \code{t.n.basis} here~~
}
   \item{basis.type.x}{the type of basis functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s}. Only "BSpline" (default) and "Fourier" are supported.
}
   \item{basis.type.y}{the type of basis functions \eqn{w_{k}(t)}{w_{k}(t)'s}. Only "BSpline" (default) and "Fourier" are supported.
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}
 \item{upper.comp}{the upper bound for the maximum number of components to be calculated. Default is 10.
%%     ~~Describe \code{upper.comp} here~~
}
  \item{thresh}{a number between 0 and 1 used to determine the maximum number of components we need to calculate. The maximum number is between one and the "upp.comp" above. The optimal number of components will be chosen between 1 and this maximum number, together with other tuning parameters by cross-validation. A smaller thresh value leads to a larger maximum number of components and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.01.
%%     ~~Describe \code{thresh} here~~
}
}
\details{
This method uses the decomposition of the coefficient functions
\deqn{\beta_i(s,t)=\sum_{k=1}^\infty\psi_{ik}(s)w_k(t),  i\in M}{\beta_i(s,t)=\psi_{i1}(s)w_1(t)+...+\psi_{iK}(s)w_K(t),   i\in M,}
and
\deqn{\gamma_{ij}(u,v,t)=\sum_{k=1}^\infty\phi_{ij,k}(s)w_k(t),
  (i,j)\in I}{\gamma_i(s,t)=\phi_{ij,1}(s)w_1(t)+...+\phi_{ij,K}(s)w_K(t),
  (i, j) \in I.}
where for each \eqn{k>0}, \eqn{\{\{\psi_{ik}(s),  i\in M\},
  \{\phi_{ij,k}(s),  (i,j)\in I\}\}}{{\psi_{ik}(s),  i\in M},
  {\phi_{ij,k}(s),  (i,j)\in I}}
 are estimated by solving a generalized functional
eigenvalue problem with the nonadaptive penalty
\deqn{\lambda\sum_{i\in M} \{|| \psi_{ik}||^2+ \tau ||\psi''_{ik}||^2\}
  + \lambda\sum_{(i,j)\in I} \{|| \phi_{ij,k}||^2+ \tau
  (||\partial_{uu}\psi_{ij,k}||^2+||\partial_{uv}\psi_{ij,k}||^2+||\partial_{vv}\psi_{ij,k}||^2)\}
} %{\lambda{|| \psi_{1k}||^2+...+  ||\psi_{pk}^2||^2} + \lambda\tau{|| \psi''_{1k}||^2+...+||\psi''_{pk}||^2},}
or the adaptive penalty
\deqn{\lambda\sum_{i\in M} \{\omega_{ik}^{(0)}||\psi_{ik}||^2+ \tau \omega_{ik}^{(2)}||\psi''_{ik}||^2\}+}
\deqn{+ \lambda\sum_{(i,j)\in I} \{\omega_{ij,k}^{(00)}|| \phi_{ij,k}||^2+ \tau
  (\omega_{ij,k}^{(20)}||\partial_{uu}\psi_{ij,k}||^2+\omega_{ij,k}^{(11)}||\partial_{uv}\psi_{ij,k}||^2+\omega_{ij,k}^{(02)}||\partial_{vv}\psi_{ij,k}||^2)\}
} %{\lambda{|| \psi_{1k}||^2+...+  ||\psi_{pk}^2||^2} + \lambda\tau{|| \psi''_{1k}||^2+...+||\psi''_{pk}||^2},}
and then  \eqn{\{w_{k}(t), k>0\}}{\{w_{k}(t), k>0\}} are estimate by regressing
\eqn{Y(t)} on \eqn{\{\hat{z}_{1},... \hat{z}_{K}\}} with nonadaptive penalty \eqn{\kappa \sum_{k=0}^K
  \|w''_k\|^2} or adaptive penalty \eqn{\kappa \sum_{k=0}^K
\omega_k^{(t)}\|w''_k\|^2} tuned by the smoothness parameter
\eqn{\kappa}. Here \eqn{\hat{z}_{k}= \sum_{i \in M}\int_{a_i}^{b_i}
  X_{i}(s)\hat{\psi}_{ik}(s)ds +\sum_{(i,j) \in
    I}\int_{a_i}^{b_i}\int_{a_j}^{b_j}
  X_{i}(u)X_{j}(v)\hat{\phi}_{ij,k}(u,v)dudv} and then centered around
its sample mean.
}
\value{An object of the ``cv.ff.interaction'' class, which is used in the function \code{\link{pred.ff.interaction}} for prediction and \code{\link{getcoef.ff.interaction}} for extracting the estimated coefficient functions.
\item{fitted_model}{a list for interval use.}
\item{y_penalty_inv}{a list for internal use.}
\item{X}{the input \code{X}.}
\item{Y}{the input \code{Y}.}
\item{x.smooth.params}{a list for internal use.}
\item{y.smooth.params}{a list for internal use.}
\item{basis.types}{a vector including \code{basis.type.x} and \code{basis.type.y}.}
}
\references{Ruiyan Luo and Xin Qi (2018) Interaction model and model
  selection for function-on-function regression, Journal of
  Computational and Graphical Statistics.
  \url{https://doi.org/10.1080/10618600.2018.1514310}
  }
\author{Xin Qi and Ruiyan Luo}
\seealso{ \code{\link{pred.ff.interaction}}, \code{\link{getcoef.ff.interaction}}}
\examples{
 
################################################################################# 
# Example: interaction function-on-function model with 
#          specified effects 
############################################################################### 


ptm <- proc.time()
library(FRegSigCom)
data(ocean)
Y=ocean[[1]]
Y.train=Y[1:50,]
Y.test=Y[-(1:50),]
t.y=seq(0,1, length.out = ncol(Y))
X.list=list()
X.train.list=list()
X.test.list=list()
t.x.list=list()
for(i in 1:4)
{
  X.list[[i]]=ocean[[i+1]]
  X.train.list[[i]]=X.list[[i]][1:50,]
  X.test.list[[i]]=X.list[[i]][-(1:50),]
  t.x.list[[i]]=seq(0,1, length.out = ncol(X.list[[i]]))
}
main.effect=1:2
inter.effect=rbind(c(1,1), c(1,2), c(2,2))
fit.fix.adaptive=cv.ff.interaction(X.train.list, Y.train, t.x.list, t.y,
           adaptive=TRUE, main.effect, inter.effect)
Y.pred=pred.ff.interaction(fit.fix.adaptive,  X.test.list)

error<- mean((Y.pred-Y.test)^2) 
print(c(" prediction error=", error))
print(proc.time()-ptm)

}
