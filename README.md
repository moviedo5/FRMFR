---
output: github_document
bibliography: bibliografia.bib
---

# FRMFR
 Supplementary Codes and Data of "Functional Regression Models with Functional Response: New Approaches and a Comparative Study" paper

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


# Installation  
In order to use paper implementation and run all files (numerical and real examples), the following prerequisites are needed:

<!--Necessary packages for running the paper codes.-->


## 1. `fda.usc` package
To install `fda.usc.devel` package (devel version of `fda.usc`) from Github with (2023/03/29):

```{r , eval = FALSE}
# install.packages("devtools")
# require(devtools)
# devtools::install_github("moviedo5/fda.usc.devel")
remotes::install_github("moviedo5/FRMFR/pkg/fda.usc.devel_2.1.1.tar.gz, INSTALL_opts = "--with-keep.source")
```