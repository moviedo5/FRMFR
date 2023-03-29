#' @title Air Quality Data
#'
#' @description 
#' This dataset contain a \code{ldata} with a data frame and 8 functional
#' \code{fdata} objects with 389 trajectories each with 24 data points 
#' (measured hourly) for all variables.
#' 
#' Air Quality is a popular dataset that consists of a list of four metal oxide
#'  chemical sensors embedded in an air quality chemical multisensor device 
#'  placed in an Italian city with a significant pollution indicator.
#'  
#'  These pollutant factors are, Nitrogen Dioxide (NO2), Carbon monoxide (CO), 
#'  Non-methane hydrocarbons(NMHC), total Nitrogen Oxides (NOx), Ozone (O3) and Benzene (C6H6). 
#'  These factors were collected as a 24 hourly averaged concentration values in each day.
#' 
#' @format 
#' \describe{
#' \code{ldata} objects: 
#' \itemize{
#' \item{\code{df}:}{ \code{data.frame} with:}
#' \itemize{
#' \item{\code{iday:}}{ Index.}
#' }
#' }
#' \itemize{
#'  8 functional \code{fdata} objects with with 389 curves
#'  (per row) observed hourly (per column).
#' \itemize{
#' \item{\code{C6H6}:}{  (log transformation of) daily curves of Benzene concentrations measured every hour (in microg/m^3) }
#' \item{\code{CO}:}{ (log transformation of) daily curves of CO concentrations measured every hour (in mg/m^3) }
#' \item{\code{NMHC}:}{ (log transformation of) daily curves of Non Metanic Hydrocarbons concentrations measured every hour (in ppm) }
#' \item{\code{NOx}:}{  (log transformation of) daily curves of Total Nitrogen Oxides concentrations measured every hour (in mg/m^3) }
#' \item{\code{NO2}:}{ (log transformation of) daily curves of Nitrogen Dioxide concentrations measured every hour  (in microg/m^3).}
#' \item{\code{O3}:}{ (log transformation of) daily curves of Ozone concentrations measured every hour  (in microg/m^3).}
#' \item{\code{Temp}:}{   (log transformation of) daily curves of temperatura measured every hour (in Celsius) }
#' \item{\code{RH}:}{  daily curves of relative humidity measured every hour (in percentage) }
#' } 
#' }
#' }
#' @details
#' \describe{
#' The dataset contains 389 daily curves of NO2, CO, NMHC, NOx, C6H6, 
#' temperature and humidity measured in a polluted area within
#'  an Italian city.  Data were recorded from March 2004 to
# February 2005. 34 days have missing values.
#' }
#' @docType data
#' @keywords datasets
#' @name AirQuality
#' @usage data(AirQuality)
#' @source The UCI Machine Learning Repository: 
#' \url{https://archive.ics.uci.edu/ml/datasets/air+quality}
#' @references 
#' De Vito, S. et al (2008) Sensors and Actuators B: Chemical, 129: 50-757.
NULL

