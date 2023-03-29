---
editor_options: 
  markdown: 
    wrap: 72
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{=html}
<!--
# fda.usc.devel
README.md is generated from README.Rmd. Please edit that file
# [![Travis-CI Build Status](https://api.travis-ci.org/maierhofert/classiFunc.svg?branch=master)](https://travis-ci.org/maierhofert/classiFunc)
# [![packageversion](https://img.shields.io/badge/Package%20version-0.1.1-orange.svg?style=flat-square)](https://CRAN.R-project.org/package=classiFunc)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/fda.usc.devel)](https://cran.r-project.org/package=fda.usc.devel)
[![](https://cranlogs.r-pkg.org/badges/fda.usc.devel)](https://cran.r-project.org/package=fda.usc.devel)
# [![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
Incluir reference card
[![Travis-CI Build Status](https://api.travis-ci.org/moviedo5/fda.usc.devel.svg?branch=master)](https://travis-ci.org/moviedo5/fda.usc.devel)
## fda.usc.devel: Functional Data Analysis and Utilities for Statistical Computing
-->
```
# fda.usc.devel: Functional Data Analysis and Utilities for Statistical Computing <img src="inst/figures/fda.usc.png" alt="Logo fda.usc" align="right" width="140"/>

[![Licence](https://img.shields.io/badge/licence-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)



## Package overview (devel version)

```{=html}
<!--

C:\Program Files\R\Rlibrary
The **fda.usc.devel** package implements methods  for exploratory and descriptive analysis of functional data such as depth measurements, atypical curves detection, regression models, supervised classification, unsupervised classification and functional analysis of variance.
-->
```
**fda.usc.devel** package carries out exploratory and descriptive
analysis of functional data exploring its most important features such
as depth measurements or functional outliers detection, among others. It
also helps to explain and model the relationship between a dependent
variable and independent (regression models) and make predictions.
Methods for supervised or unsupervised classification of a set of
functional data regarding a feature of the data are also included. It
can perform functional ANOVA, hypothesis testing, functional response
models and many others.

## Installation

You can install the devel fda.usc.devel version from Github with:

``` r
library(devtools)
devtools::install_github("moviedo5/fda.usc.devel")
```

## Issues & Feature Requests

For issues, bugs, feature requests etc. please use the [Github
Issues](https://github.com/moviedo5/fda.usc.devel/issues). Input is
always welcome.

## Documentation

A hands on introduction to can be found in the reference
[vignette](https://www.jstatsoft.org/index.php/jss/article/view/v051i04/v51i04.pdf).

Cheatsheet [fda.usc.devel reference
card](https://zenodo.org/record/3386752/files/RefCard_fda.usc_v1.pdf?download=1).

## References

Febrero-Bande, M. and Oviedo de la Fuente, M. (2012). Statistical
Computing in Functional Data Analysis: The R Package fda.usc. *Journal
of Statistical Software*, 51(4):1-28.
<http://www.jstatsoft.org/v51/i04/>


Darbalaei, M., Amini, M., Febrero-Bande, M., & la Fuente, M. O. D. (2022). Functional Regression Models with Functional Response: New Approaches and a Comparative Study. arXiv preprint arXiv:2207.04773. <https://doi.org/10.48550/arXiv.2207.04773>


```{=html}
<!--
# install.packages("fda")
# library(fda)
# version minima de R 3.5.0
# pruebas wavelets en create.fdata.basis
# se incluye ORCID

library(roxygen2)
# setwd("C:/Users/moviedo/github/fda.usc.devel/")
#setwd("C:/Users/manue/Nextcloud/github/fda.usc.devel")
getwd()

pkgbuild::compile_dll()

roxygenize()
#
devtools::document()

library(devtools)
# devtools::install("C:/Users/moviedo/github/fda.usc.devel")
# devtools::install("/home/manuel.oviedo/Nextcloud/github/fda.usc.devel")
# devtools::install("C:/Users/moviedo/github/fda.usc.devel_2.0.4")
# devtools::install("C:/Users/moviedo/github/fda.usc.devel")
# devtools::install("D:/Users/moviedo/github/fda.usc.devel") # AreaCientifica
#' 3
# devtools::install("C:/Users/manue/Nextcloud/github/fda.usc.devel"
 #  , auth_token = "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"
, upgrade = FALSE  )
# token del  20220607

#devtools::check("C:/Users/moviedo/github/fda.usc.devel",args="--as-cran",manual=TRUE,remote=TRUE)

devtools::check()

devtools::install()

#' 4

# devtools::check("C:/Users/moviedo/github/fda.usc.devel_2.0.3",args="--as-cran",manual=TRUE)

# https://github.com/moviedo5/fda.usc.devel.git

# args(devtools::install_svn)
#devtools::install_github("https://github.com/moviedo5/fda.usc.devel.git")

devtools::install_github("https://github.com/moviedo5/fda.usc.devel"
,auth_token ="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")


devtools::install_github("moviedo5/fda.usc.devel"
,auth_token ="d515a42a964561386415028c2040d1b4281c14e2")


devtools::install_github(
  repo = "moviedo5/fda.usc.devel",
 # ref = "devel",
  #subdir = "Library/PackageName",
  auth_token = "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"
  , upgrade = FALSE)

# old token ,auth_token ="8a58c906d7eb2b319c51b1fee4e4eb6f7152e796")

# git@github.com:moviedo5/fda.usc.devel.git

# devtools::install("C:/Users/moviedo/github/fda.usc.devel")
devtools::install()

R CMD build fda.usc.devel
R CMD check fda.usc.devel_2.0.4.tar.gz --as-cran
R CMD INSTALL fda.usc.devel_2.0.4.tar.gz --build
R CMD build fda.usc.devel
R CMD check fda.usc.devel_2.0.2.tar.gz --as-cran
R CMD INSTALL fda.usc.devel_2.0.2.tar.gz --build

R CMD INSTALL fda.usc.devel_2.0.4.tar.gz --build

** checking loading without being on the library search path ... OK
** checking use of S3 registration ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... NOTE
fregre.lm.fr.pruebas: possible error in fdata2model.fr(vfunc = vfunc,
  vnf = vnf, response = response, XX = XX, data = data, basis.x =
  basis.x, basis.y = basis.y, pf = pf, tf = tf, lambda =
lambda, P =
  P): unused argument (basis.y = basis.y)
fregre.mlm.fr: possible error in fdata2model.fr(vfunc = vfunc, vnf =
  vnf, response = response, XX = XX, data = data, basis.x = basis.x,
  basis.y = basis.y, pf = pf, tf = tf, lambda = lambda, P = P): unused
  argument (basis.y = basis.y)
predict.fregre.lm.fr.recu: no visible binding for global
variable 'bsp'
Undefined global functions or variables:
  bsp
  
Manuel Oviedo PhD thesis [Advances in functional regression and classification models](http://hdl.handle.net/10347/18236)
-->
```
<!-- ```{r example} -->

```{=html}
<!-- ## basic example code 
Falta copiar careta de scripts
   
-->
```
```{=html}
<!-- 
<img src="inst/figures/baseplot.png" height="300"/>
-->
```
