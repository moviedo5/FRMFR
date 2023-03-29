\name{step.ff.interaction}
\alias{step.ff.interaction}
\title{Stepwise variable selection procedure for FOF regression models with two-way
  interactions
}
\description{This function conducts the stepwise procedure to select
  main effects, two-way interaction and quadratic effects for the
  following family of linear function-on-function interaction models. Let \eqn{\{X_i(s),1\le i\le p\}}{{X_i(s),1\le i\le p}} be \eqn{p} potential functional predictors. The family of models is given by
\deqn{Y(t)= \mu(t)+\sum_{i \in M}\int_{a_i}^{b_i} X_i(s)\beta_i(s,t)ds+\sum_{(i,j) \in I}\int_{a_i}^{b_i}\int_{a_j}^{b_j} X_i(u)X_j(v)\gamma_{ij}(u,v,t)dudv+\epsilon(t)}{Y(t)= \mu(t)+\sum_{i \in M}\int_{a_i}^{b_i} X_i(s)\beta_i(s,t)ds+\sum_{(i,j) \in I}\int_{a_i}^{b_i}\int_{a_j}^{b_j} X_i(u)X_j(v)\gamma_{ij}(u,v,t)dudv+\epsilon(t),}
where \eqn{\mu(t)}{\mu(t)} is the intercept function. The index set
\eqn{M} of the main effects is a subset of \eqn{\{1,...,p\}}{{1,...,p}}, and the index
set \eqn{I} of the interactions and quadratic effects is a subset of the
collection of all possible pairs \eqn{\{(i,j), 1\le i\le j\le p\}}{{(i,j),
  1\le i\le j\le p}}. We require that the models in each step satisfy the
hierarchy principle: if the interaction  \eqn{X_iX_j} is
included in the model, both the main effects \eqn{X_i} and \eqn{X_j} are included. The \eqn{\{\beta_i(s,t),i\in M\}}{{\beta_i(s,t),i \in M}} and  \eqn{\{\gamma_{ij}(u,v,t),(i,j)\in I\}}{{\gamma_{ij}(u,v,t), (i,j) \in
  I}} are the corresponding coefficient functions. The \eqn{\epsilon(t)}{\epsilon(t)} is the noise function. When the final model is selected, this function also fits the selected model.
}
\usage{
step.ff.interaction(X, Y, t.x, t.y, adaptive=FALSE, s.n.basis=40, t.n.basis=40,
          inter.n.basis=20, basis.type.x="Bspline", basis.type.y="Bspline",
          K.cv=5, upper.comp=8, thresh=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of \eqn{p} potential functional predictors. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th potential functional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of observation time points for  \eqn{X_i(s)}{X_i(s)}.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{the \eqn{n\times m}{n*m} data matrix for the functional response \eqn{Y(t)}{Y(t)}, where  \eqn{n}  is the sample size and  \eqn{m}  is the number of the observation time points for \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{Y} here~~
}
  \item{t.x}{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of obesrvation time points of the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}.
%%     ~~Describe \code{t.x} here~~
}
  \item{t.y}{the vector of observation time points of the functional response \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{t.y} here~~
}
  \item{adaptive}{a logic value indicating whether using adaptive
    penalty that uses different smoothness tuning parameters for
    different target functions. Default is FALSE.
}
  \item{s.n.basis}{the number of basis functions used for estimating the
    functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s} (see details in \code{\link{cv.ff.interaction}}). Default is 40.
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
   \item{basis.type.x}{the type of basis functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s} . Only "BSpline" (default) and "Fourier" are supported.
}
   \item{basis.type.y}{the type of basis functions \eqn{w_{k}(t)}{w_{k}(t)'s}. Only "BSpline" (default) and "Fourier" are supported.
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}
  \item{upper.comp}{the upper bound for the maximum number of components to be calculated. Default is 8.
%%     ~~Describe \code{upper.comp} here~~
}
  \item{thresh}{a number between 0 and 1 used to determine the maximum number of components we need to calculate. The maximum number is between one and the "upp.comp" above. The optimal number of components will be chosen between 1 and this maximum number, together with other tuning parameters by cross-validation. A smaller thresh value leads to a larger maximum number of components and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.01.
%%     ~~Describe \code{thresh} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{An object of the ``step.ff.interaction'' class, which is used in the function \code{\link{pred.ff.interaction}} for prediction and \code{\link{getcoef.ff.interaction}} for extracting the estimated coefficient functions.
\item{opt.main.effects}{a vector of indices of the selected main effects.}
\item{opt.interaction.effects}{a matrix of two columns.  Each row specifies the indices of a selected two-way interaction or quadratic effect.}
\item{fitted_model}{a list of the output of the fitted selected model and only for internal use.}
\item{y_penalty_inv}{a list for interval use.}
\item{X}{the input \code{X}.}
\item{Y}{the input \code{Y}.}
\item{x.smooth.params}{a list for internal use.}
\item{y.smooth.params}{a list for internal use.}
}
%\value{Y}
\references{Ruiyan Luo and Xin Qi (2018) Interaction model and model
  selection for function-on-function regression, Journal of
  Computational and Graphical Statistics.
  \url{https://doi.org/10.1080/10618600.2018.1514310}
}
\author{Xin Qi and Ruiyan Luo
}
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

library(FRegSigCom)
data(ocean)

Y=ocean$Salinity
X=list()
X[[1]]=ocean$Potential.density
X[[2]]=ocean$Temperature
X[[3]]=ocean$Oxygen
n.curves=length(X)
ntot=dim(Y)[1]
ntrain=50
ntest=ntot-ntrain
X.uncent=X
for(i in 1:n.curves){
  X[[i]]=scale(X.uncent[[i]],center=TRUE, scale=FALSE)
}
lengthX=dim(X[[1]])[2]
lengthY=dim(Y)[2]
t.x=seq(0,1,length=lengthX)
t.y=seq(0,1,length=lengthY)
I.train=sample(1:ntot, ntrain)
X.train=list()
X.test=list()
t.x.all=list()
for(j in 1:n.curves){
  X.train[[j]]=X[[j]][I.train,]
  X.test[[j]]=X[[j]][-I.train,]
  t.x.all[[j]]=t.x
}
Y.train=Y[I.train, ]
Y.test=Y[-I.train, ]


###############################
#model selection
###############################

fit.step=step.ff.interaction(X.train, Y.train, t.x.all, t.y)
Y.pred=pred.ff.interaction(fit.step,  X.test)
error.selected=mean((Y.pred-Y.test)^2)
print(c("error.selected=", error.selected))
#coef.obj=getcoef.ff.interaction(fit.step)
#str(coef.obj)
}
