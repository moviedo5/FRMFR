#' 
#' @title Bike Sharing Data
#'
#' @description 
#' This dataset contain a \code{ldata} with a data frame and 5 functional
#' \code{fdata} objects with 102 trajectories each with 24 data points 
#' (measured hourly) for all variables.
#' 
#' @format 
#' \describe{
#' \code{ldata} objects: 
#' \itemize{
#' \item{\code{df}:}{ \code{data.frame} with 2 variables:}
#' \itemize{
#' \item{\code{iday:}}{(numeric) Index of the day (from 1 to 729).}
#' \item{\code{date}:}{ (Date) variable from 2011-01-01 to 2012-12-29.}
#' }
#' }
#' \itemize{
#'  5 functional \code{fdata} objects with with n=102 curves
#'  (per row) observed hourly (per column).
#' \itemize{
#' \item{\code{logNBCR}:}{ Logarithm of number of casual bike rentals.}
#' \item{\code{temp}:}{ Normalized temperature in Celsius.}
#' \item{\code{feeltemp}:}{ Normalized feeling temperature in Celsius.}
#' \item{\code{humidity}:}{ Normalized humidity.}
#' \item{\code{windspeed}:}{ Normalized wind speed.}
#' } 
#' }
#' These functional variables are recorded each hour from Jan 1st, 2011 to Dec 
#' 1st, 2012. Similar to Kim et al. (2018), we only consider the data for Saturday
#' trajectories and three curves with missing values have ignored. 
#' }
#' @details
#' \describe{
#' Data manipulation 
#' \itemize{
#' \item{\code{logNBCR}:}{ NBCR is log-transformed to avoid its natural 
#' heteroskedasticity.}
#' \item{\code{temp}:}{ The normalized values are derived via 
#' \eqn{(t-t_{min})/(t_max-t_min)}, \eqn{t_min=-8}, \eqn{t_max=+39}.}
#' \item{\code{feeltemp}:}{ The normalized values are derived via
#'  \eqn{(t-t_min)/(t_max-t_min)}, \eqn{t_min=-16}, \eqn{t_max=+50}.}
#' \item{\code{humidity}:}{ The normalized values are divided to 100 (max).}
#' \item{\code{windspeed}:}{ The normalized values are divided by 67 (max).}
#' }
#' Amini et al. 2023 consider the number of casual bike rentals (\code{logNBCR}) 
#' as our functional  response, and Temperature (\code{temp}), Feeling 
#' Temperature (\code{feeltemp}), Humidity (\code{humidity}) and Wind Speed 
#' (\code{windspeed}) as the functional covariates. 
#' 
#' This dataset is originally collected by Capital Bikeshare System (CBS), 
#' Washington D.C., USA. 

#' }
#' @docType data
#' @keywords datasets
#' @name BikeSharing
#' @usage data(BikeSharing)
#' @source The UCI Machine Learning Repository: 
#' \url{https://archive.ics.uci.edu/ml/datasets/bike+sharing+dataset}
#' @references 
#' Darbalaei, M., Amini, M., Febrero-Bande, M., & la Fuente, M. O. D. (2022).
#'  Functional Regression Models with Functional Response: New Approaches and 
#'  a Comparative Study. arXiv preprint \url{https://arxiv.org/abs/2207.04773}.
#' 
NULL

