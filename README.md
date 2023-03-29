---
output: github_document
bibliography: bibliografia.bib
---

# FRMFR
Supplementary codes and data used in the paper:
Darbalaei, M., Amini, M., Febrero-Bande, M., & Oviedo-de-la-Fuente, M. O. D. (2022). *Functional Regression Models with Functional Response: New Approaches and a Comparative Study*. arXiv preprint [arXiv:2207.04773](https://doi.org/10.48550/arXiv.2207.04773).


Please cite this paper as:

```
@article{darbalaei2022functional,
  title={Functional Regression Models with Functional Response: New Approaches and a Comparative Study},
  author={Darbalaei, Mohammad and Amini, Morteza and Febrero-Bande, Manuel and Oviedo de-la Fuente, Manuel},
  journal={arXiv preprint arXiv:2207.04773},
  year={2022}
}
```
https://github.com/moviedo5/FRMFR/blob/main/pkg/fda.usc.devel_2.1.1.tar.gz

# Installation  
In order to use paper implementation and run all files (numerical and real examples), the following prerequisites are needed:

<!--Necessary packages for running the paper codes.-->


## 1. `fda.usc` package
To install `fda.usc.devel` package (devel version of `fda.usc`) from Github with (2023/03/29):

```{r , eval = FALSE}
# install.packages("devtools")
require(devtools)
devtools::install_github("moviedo5/FRMFR/pkg/fda.usc.devel")
```


<!--
```{r , eval = FALSE}
packageurl <- "https://github.com/moviedo5/FRMFR/blob/main/pkg/fda.usc.devel_2.1.1.tar.gz"
install.packages(packageurl, repos=NULL)
url <- "https://github.com/rubenfcasal/simres/releases/download/v0.1/simres_0.1.3.zip"
install.packages(url, repos = NULL)
```
-->

To compares our proposed methods (namely, FLMFR, FSAMFR, and FKAMFR, which are available in the `fda.usc` package (devel version) through the commands `fregre.lm.fr`, `fregre.sam.fr` and `fregre.kam.fr`) with the four mentioned competitor methods (namely, PFR, FAMM, LSC, and DISC).


## 2. `refund` package

PFR and FAMM methods are available in the `refund` package through the command `pffr`, where the argument formula allows us to include linear `ffpc`, `ff` or nonlinear term `sff`. 

To install refund package from CRAN or Github.

```{r , eval = FALSE}
# install.packages("refund")
# install_version("refund", version = "0.1-30", repos = "http://cran.us.r-project.org")
devtools::install_github("refunders/refund") # latest patched version directly from Github
```

##  3. `FRegSigCom` package

The authors considered the latter as an experimental feature. LSC and DISC methods are available in the `FRegSigCom` package through the commands `cv.sigcom` and `cv.nonlinear`.  

This package is not currently maintained and its latest version was published in November 2018 but, anyway, it can be downloaded and installed from the Packages/Archive section of CRAN.


```{r , eval = FALSE}
devtools::install_github("moviedo5/FRMFR/pkg/fda.usc.devel")

devtools::install_github("moviedo5/FRMFR/pkg/FRegSigCom")
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/FRegSigCom/FRegSigCom_0.3.0.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")

```

<!--The package FRegSigCom, also needed for the previous codes, must be downloaded from the Archive section of CRAN (https://cran.r-project.org/rc/contrib/Archive/FRegSigCom/) and manually installed.-->

