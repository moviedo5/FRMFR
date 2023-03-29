\name{cv.nonlinear}
\alias{cv.nonlinear}
%ff.nonlinear.R
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for nonlinear function-on-function regression
%%  ~~function to do ... ~~
}
\description{This function is used to perform cross-validation and build the final model using signal compression approach for the following nonlinear function-on-function regression model:
\deqn{Y(t)= \mu(t)+\sum_{i=1}^p\int_{a_i}^{b_i}F_i(X_i(s), s,t)ds+\varepsilon(t)}{Y(t)= \mu(t)+\int_{a_1}^{b_1}F_1(X_1(s), s,t)ds+...+\int_{a_p}^{b_p}F_p(X_p(s), s,t)ds+\epsilon(t),}
where \eqn{\mu(t)}{\mu(t)} is the intercept function, \eqn{\{F_i(x,
  s,t),1\le i\le p\}}{{F_i(x,s,t), 1\le i\le p}} are all unspecified
nonlinear functions of \eqn{x, s, t},  \eqn{\{X_i(s),1\le i\le p\}}{{X_i(s), 1\le i\le p}} are functional predictors, and \eqn{\epsilon(t)}{\epsilon(t)} is the noise function.

In this method, we require that all the sample curves of each functional predictor be observed in a common dense set, but the observation points can be different for different functional predictors. All the sample curves of the functional response are observed in a common dense set.
}
\usage{
cv.nonlinear(X, Y, t.x.list, t.y, s.n.basis = 40, x.n.basis = 40,
        t.n.basis = 40, K.cv = 5, upper.comp = 10, thresh = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of length \eqn{p}, the number of functional predictors. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of observation time points for  \eqn{X_i(s)}{X_i(s)}.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{the \eqn{n\times m}{n*m} data matrix for the functional response \eqn{Y(t)}{Y(t)}, where  \eqn{n}  is the sample size and  \eqn{m}  is the number of the observation time points for \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{Y} here~~
}
  \item{t.x.list}{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of observation time points of the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}.
%%     ~~Describe \code{t.x} here~~
}
  \item{t.y}{the vector of observation time points of the functional response \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{t.y} here~~
}
  \item{s.n.basis}{the number of B-spline basis functions for the argument \eqn{s}. Default  is 40.
%%     ~~Describe \code{s.n.basis} here~~
}
  \item{x.n.basis}{the number of B-spline basis functions for \eqn{x}. Default is 40.
%%     ~~Describe \code{x.n.basis} here~~
}
  \item{t.n.basis}{the number of B-spline basis functions for \eqn{t}. Default is 4.
%%     ~~Describe \code{t.n.basis} here~~
}
  \item{K.cv}{the number of CV folds. Default is 5.
}
 \item{upper.comp}{the upper bound for the maximum number of components to be calculated. Default is 10.
%%     ~~Describe \code{upper.comp} here~~
}
  \item{thresh}{a number between 0 and 1 used to determine the maximum number of components we need to calculate. The maximum number is between one and the "upp.comp" above. The optimal number of components will be chosen between 1 and this maximum number, together with other tuning parameters by cross-validation. A smaller thresh value leads to a larger maximum number of components and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.01.
%%     ~~Describe \code{thresh} here~~
}
}
\details{
This method estimates a special decomposition:
\deqn{F_i(x,s,t)=\sum_{k=1}^\infty G_{ik}(x,s)w_k(t),  1\le i\le p.}{F_i(x,s,t)=G_{i1}(x,s)w_1(t)+...+G_{iK}(x,s)w_K(t),   1\le i\le p.}
We first estimate \eqn{\bold{G}_k=(G_{1k}(x_1,s),...,G_{pk}(x_p,s))^T} for each \eqn{k>0}  by solving a generalized penalized functional
eigenvalue problem
\deqn{max_{\bold{G}} \frac{\hat{\bold{\Lambda}}(\bold{G},\bold{G})}{\hat{\bold{\Sigma}}(\bold{G},\bold{G})+P(\bold{G})}}
\deqn{ {\rm{s.t.}} \quad \hat{\bold{\Sigma}}(\bold{G}_{k'},\bold{G})=0}
\deqn{  {\rm{and}} \quad  \int\int \bold{G}(s)^T
  \hat{\Sigma}(s,s')\bold{G}_{k'}(s')dsds'=0 \quad {\rm{for}}
  \quad k'<k}
where
\deqn{\hat{\bold{\Lambda}}(\bold{G},\bold{G})=\int \left[\sum_{\ell=1}^n\sum_{i=1}^p \int \{G_i(x_{\ell,i}(s),s)-\bar{G}_i(s)\}ds \{y_{\ell}(t)-\bar{y}(t)\}  \right]^2 dt /n^2,}
\deqn{\hat{\bold{\Sigma}}(\bold{G},\tilde{\bold{G}})=\sum_{\ell=1}^n \left[\sum_{i=1}^p \int \{G_i(x_{\ell,i}(s),s)-\bar{G}_i(s)\}ds\right]\left[\sum_{j=1}^p \int \{\tilde{G}_j(x_{\ell,j}(s),s)-\bar{\tilde{G}}_j(s)\}ds\right]/n,}
 penlaty
\deqn{P(G)=\lambda \sum_{j=1}^p \{\|G_j\|^2+\tau(\|\partial_{xx}G_j\|^2+\|\partial_{xs}G_j\|^2+\|\partial_{ss}G_j\|^2)\},}
and \eqn{
  <G,\tilde{G}>_{H^2}=<G,\tilde{G}>_{L^2}+<\partial_{xx}G,\partial_{xx}\tilde{G}>_{L^2}+<\partial_{xs}G,\partial_{xs}\tilde{G}>_{L^2}+<\partial_{ss}G,\partial_{ss}\tilde{G}>_{L^2}
} with \eqn{<G,\tilde{G}>_{L^2}=\int\int G(x,s)\tilde{G}(x,s) dxds}.
Then we estimate \eqn{\{w_{k}(t), k>0\}}{{w_{k}(t), k>0}} by regressing
\eqn{\{y_{\ell}(t)\}_{\ell=1}^n} on \eqn{\{\hat{z}_{\ell,1},... \hat{z}_{\ell,K}\}_{\ell=1}^n} with penalty \eqn{\kappa \sum_{k=0}^K
  \|w''_k\|^2} tuned by the smoothness parameter
\eqn{\kappa}. Here \eqn{\hat{z}_{\ell,k}= \sum_{i=1}^p \int (G_{ik}(x_{\ell,i}(s),s)-\bar{G}_{ik}(s))ds}, and \eqn{\bar{G}_{ik}(s)=\sum_{\ell=1}^n G_{ik}(X_{\ell,i}(s),s)/n}.
}

\value{A fitted CV-object, which is used in the function \code{\link{pred.nonlinear}} for prediction and \code{\link{getcoef.nonlinear}} for  extracting the estimated coefficient functions.
\item{opt.K}{optimal number of components to be selected.}
\item{opt.lambda}{optimal value for \eqn{\lambda}.} %, the tuning parameter  for \eqn{\{G_{ik}(x,s), 1 \le i \le p\}}{{G_{ik}(x,s), 1 \le i \le p}}.}
\item{opt.tau}{optimal value for \eqn{\tau}.}%, the smoothness tuning parameter for \eqn{w_k(t)}{w_k(t)}.}
\item{opt.kappa}{optimal value for \eqn{\kappa}, the smoothness tuning parameter for \eqn{w_k(t)}{w_k(t)}.}
\item{...}{other output for internal use.}
%\item{y.params}{a list for internal use.}
}
%\item{opt.index}{index of the optimal lambda.}
%\item{min.error}{minimum CV error.}
%\item{errors}{list for CV errors. }
%\item{lambda.set}{the set of \eqn{\lambda} values used in CV.}
%\item{kappa.set}{the set of \eqn{\kappa} values used in CV.}
%\item{opt.K}{optimal number of components to be selected.}
%\item{opt.lambda}{optimal value for \eqn{\lambda}, the tuning parameter
%  for \eqn{\{G_{ik}(x,s), 1 \le i \le p\}}{{G_{ik}(x,s), 1 \le i \le p}}.}
%\item{opt.kappa}{optimal value for \eqn{\kappa}, the smoothness tuning parameter for \eqn{w_k(t)}{w_k(t)}.}
%\item{fit.1}{a list for internal use.}
%\item{y.params}{a list for internal use.}
%}
\references{
   Xin Qi and Ruiyan Luo. (Accepted) Nonlinear function on function additive model
   with multiple predictor curves. Stistica Sinica.
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

#################################################################################
# Example: Nonlinear function-on-function model
###############################################################################

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
X.train.list <- list(X[train.id,])
Y.train <- Y[train.id, ]
X.test.list <- list(X[-(train.id),])
Y.test <- Y[-(train.id), ]
fit.cv=cv.nonlinear(X.train.list, Y.train, t.x, t.y, upper.comp=3,
s.n.basis=20, x.n.basis=20,t.n.basis=20) 
           # in practice, use the default values (or larger) for
           # "upper.comp", "s.n.basis","x.n.basis", and "t.n.basis".
Y.pred=pred.nonlinear(fit.cv, X.test.list)
nonlinear.error= mean((Y.pred-Y.test)^2)
print(c("prediction error=",nonlinear.error))
print(proc.time()-ptm)
}