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
# install_version("refund", version = "0.1-30", 
#                 repos = "http://cran.us.r-project.org")
# latest patched version directly from Github
devtools::install_github("refunders/refund") 
```

##  3. `FRegSigCom` package

The authors considered the latter as an experimental feature. LSC and DISC methods are available in the `FRegSigCom` package through the commands `cv.sigcom` and `cv.nonlinear`.  

This package is not currently maintained and its latest version was published in November 2018 but, anyway, it can be downloaded and installed from the Packages/Archive section of CRAN.


```{r , eval = FALSE}
devtools::install_github("moviedo5/FRMFR/pkg/FRegSigCom")
# url <- "https://cran.r-project.org/src/contrib/Archive/FRegSigCom/FRegSigCom_0.3.0.tar.gz"
# install.packages(url, repos=NULL, type="source")
```

<!--The package FRegSigCom, also needed for the previous codes, must be downloaded from the Archive section of CRAN (https://cran.r-project.org/rc/contrib/Archive/FRegSigCom/) and manually installed.-->



# Simulation (Numerical Studies)

<!--Consult a detailed documentation of the code and examples of use in-->

+ `./inst/script/Simulation.R`: Code for main simulation. Scenarios 1--4.

  - Linear smooth (LS)
  - Linear non-smooth (LNS)
  - Nonlinear smooth (NLS)
  - Nonlinear non{smooth (NLNS)
  
  
```{r , eval =T,message=F,warning=F}
library(fda.usc.devel)
# source("./inst/script/Simulation.R")
```


 
# Real Data Applications

Consult a detailed documentation of the data examples and R code of used.

## 1. Air Quality Data
Our last example is the Air Quality dataset (AQI) available from the UCI machine learning repository @Qi2019. AQI is a popular dataset consisting of five metal oxide chemical sensors embedded into an air quality multisensor device. The column names in the dataset begin with *PT*. The sensors are labeled with:

+ Carbon monoxide (`CO`), 
+ Non-methane hydrocarbons (`NMHC`), 
+ total Nitrogen Oxides (`NOx`), 
+ Ozone (`O3`) 
because it is supposed that its measures are related with the respective pollutants. 

The goal of this study is to predict the content of Benzene (`C6H6`) obtained through an independent analyzer considered the Ground Truth. 
These sensors were collected as 24 hourly averaged concentration values each day jointly with the relative humidity (`rH`) as an external factor. 

```{r , eval =T,message=F,warning=F}
data("AirQuality")
names(AirQuality)
par(mfrow=c(3,2))
plot(AirQuality$NO2)
plot(AirQuality$CO)
plot(AirQuality$NMHC)
plot(AirQuality$NOx)
plot(AirQuality$C6H6)
plot(AirQuality$Temp)

# source("BikeSharing.R")
```

+ `./inst/script/data-real-aqi.R`: Code for example Air Quality.


<!--+ AirQualityUCI.xlsx: Air Quality Data.-->
 
## 2. Bike-sharing Data

To illustrate how our proposed function-on-function methods work, we use the Bike-sharing data [@Fanaee-T2014] as our first example. This dataset is collected by [Capital Bikeshare System (CBS), Washington D.C., USA](https://ride.capitalbikeshare.com/system-data).
The number of casual bike rentals (`NCBR`) is considered as our functional response variable


Functional predictors: 
+ Temperature (`T`), 
+ Humidity (`H`), 
+ Wind Speed (`WS`) and
+ Feeling Temperature (`FT`)


The corresponding plots are displayed in Figure XXX. 

```{r , eval = T}
data("AirQuality")
names(AirQuality)

listfinal <- AirQuality
dim(listfinal$df)
al=lapply(listfinal[-1],function(v) which(fda.usc.devel:::is.na.fdata(v)))
al=unique(rapply(al,drop))
listfinal=subset(listfinal,-al)
dim(listfinal$df)

par(mfrow=c(2,3))
plot(listfinal$NMHC,col="grey")
lines(func.mean(listfinal$NMHC))
plot(listfinal$CO,col="grey")
lines(func.mean(listfinal$CO))
plot(listfinal$NOx,col="grey")
lines(func.mean(listfinal$NMHC))
plot(listfinal$O3,col="grey")
lines(func.mean(listfinal$O3))
plot(AirQuality$RH,col="grey")
lines(func.mean(listfinal$RH))
plot(listfinal$C6H6,col="grey")
lines(func.mean(listfinal$C6H6))

# source("AirQuality.R")
```

These variables are recorded each hour from January 1, 2011, to December 31, 2012. Similar to @Kim2018, we only consider the data for Saturday trajectories, and `NBCR` is log--transformed to avoid its natural heteroskedasticity. Ignoring three curves with missing values, the dataset contains 102 trajectories, each with 24 data points (hourly) for all variables. 

+ `./inst/script/bike-sharing2.R`: Code for Bike--sharing data example.


<!--+ bike-sharing2.R: Code for Bike--sharing data example.
+ hour.csv: Bike--sharing data.-
-->

## 3. Electricity Demand and Price Data

Daily profiles of Electricity Price  and Demand, both measured hourly, are obtained from two biannual periods separated by ten years: 2008-2009 and 2018-2019 (source:omie.es). 

```{r , eval = T}
data(omel2008)
par(mfrow=c(2,2))
plot(omel2008$precio,col="grey")
lines(func.mean(omel2008$precio),lwd=2)
plot(omel2008$energia,col="grey")
lines(func.mean(omel2008$energia),lwd=2)
data(omel2018)
plot(omel2018$precio,col="grey")
lines(func.mean(omel2018$precio),lwd=2)
plot(omel2018$energia,col="grey")
lines(func.mean(omel2018$energia),lwd=2)
# library(lubridate)
# plot(omel2008$precio,col=year(omel2008$df$ifecha)-2007)
# plot(omel2008$energia,col=year(omel2008$df$ifecha)-2007)
# data(omel2018)
# omel2018 <- ldat
# plot(omel2018$precio,col=year(omel2018$df$ifecha)-2015)
# plot(omel2018$energia,col=year(omel2018$df$ifecha)-2015)
```

Profiles for Electricity Demand (first row) and Electricity Price (second row) for the periods 2008-09 (first column) and 2018-19 (second column). The black line corresponds to the functional mean of each dataset.

+ `./inst/script/Exampleomel.R`: Code for Electricity Demand and Price example.

<!--
+ Exampleomel.R: Code for Electricity Demand and Price example.

+ omel2008-09.rda: Electricity data for 2008-09 period.

+ omel2018-19.rda: Electricity data for 2018-19 period.

@darbalaei2022functional

-->



## References


<!--Darbalaei, M., Amini, M., Febrero-Bande, M., & Oviedo-de-la-Fuente, M. O. D. (2022). Functional Regression Models with Functional Response: New Approaches and a Comparative Study. arXiv preprint [arXiv:2207.04773](https://doi.org/10.48550/arXiv.2207.04773).

Fanaee-T, H., & Gama, J. (2014). Event labeling combining ensemble detectors and background knowledge. Progress in Artificial Intelligence, 2, 113-127.

Febrero-Bande, M., González-Manteiga, W. & Oviedo de la Fuente, M. Variable selection in functional additive regression models. Comput Stat 34, 469–487 (2019). https://doi.org/10.1007/s00180-018-0844-5

Kim, J. S., Staicu, A. M., Maity, A., Carroll, R. J., & Ruppert, D. (2018). Additive function-on-function regression. Journal of Computational and Graphical Statistics, 27(1), 234-244.

Qi, X., & Luo, R. (2019). Nonlinear function-on-function additive model with multiple predictor curves. Statistica Sinica, 29(2), 719-739.

-->
