#' @title omel dataset
#'
#' @description Electricity Demand and Price Data from the Iberian Energy Market 
#' 
#' @details
#' 730  daily curves for price and volume from 2018-01-01 to 2019-12-31.
#'
#' @format 
#' \describe{
#' 
#' The variables are as follows:
#' \itemize{
#' \item{\code{df}:}{ \code{data.frame} with the following variables:}
#' \itemize{
#' \item{\code{index}:}{ Index of day }
#' \item{\code{ifecha}:}{ Date, from 2018-01-01 to 2019-12-31. }
#' \item{\code{dias}:}{ Day of week.}
#' }   
#' \item{\code{precio}:}{ \code{fdata} class object. The daily profiles of Electricity Price measured hourly}
#' \item{\code{energia}:}{ \code{fdata} class object. The daily profiles of Electricity Demand measured hourly}
#' }
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name omel2018
#' @usage data(omel2018)
#' @source \url{https://www.omie.es/es/market-results/daily/daily-market/daily-hourly-price}
#' @references 
#' Febrero-Bande, M., Gonzalez-Manteiga, W., Oviedo De La Fuente, M. (2019). 
#' Variable selection in functional additive regression models.
#' Computational Statistics, 34(2), 469-487.
NULL
