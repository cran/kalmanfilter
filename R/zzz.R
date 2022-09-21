#' Kalman Filter
#' 
#' \emph{kalmanfilter} Rcpp implementation of the multivariate Kalman filter for 
#' state space models that can handle missing values and exogenous data in the 
#' observation and state equations. See the package vignette using 
#' \code{browseVignettes("kalmanfilter")} to view it in your browser.
#'  
#' @docType package
#' @author Alex Hubbard
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib kalmanfilter, .registration=TRUE
#' @name kalmanfilter
NULL