% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kalman_filter.R
\name{kalman_filter}
\alias{kalman_filter}
\title{Kalman Filter}
\usage{
kalman_filter(ssm, yt, Xo = NULL, Xs = NULL, weight = NULL, smooth = FALSE)
}
\arguments{
\item{ssm}{list describing the state space model, must include names
B0 - N_b x 1 matrix (or array of length yt), initial guess for the unobserved components 
P0 - N_b x N_b matrix (or array of length yt), initial guess for the covariance matrix of the unobserved components
Dm - N_b x 1 matrix (or array of length yt), constant matrix for the state equation
Am - N_y x 1 matrix (or array of length yt), constant matrix for the observation equation
Fm - N_b X p matrix (or array of length yt), state transition matrix
Hm - N_y x N_b matrix (or array of length yt), observation matrix
Qm - N_b x N_b matrix (or array of length yt), state error covariance matrix
Rm - N_y x N_y matrix (or array of length yt), state error covariance matrix
betaO - N_y x N_o matrix (or array of length yt), coefficient matrix for the observation exogenous data
betaS - N_b x N_s matrix (or array of length yt), coefficient matrix for the state exogenous data}

\item{yt}{N x T matrix of data}

\item{Xo}{N_o x T matrix of exogenous observation data}

\item{Xs}{N_s x T matrix of exogenous state}

\item{weight}{column matrix of weights, T x 1}

\item{smooth}{boolean indication whether to run the backwards smoother}
}
\value{
list of cubes and matrices output by the Kalman filter
}
\description{
Kalman Filter
}
\examples{
\dontrun{
#Stock and Watson Markov switching dynamic common factor
library(kalmanfilter)
library(data.table)
data(sw_dcf)
data = sw_dcf[, colnames(sw_dcf) != "dcoinc", with = FALSE]
vars = colnames(data)[colnames(data) != "date"]

#Set up the state space model
ssm = list()
ssm[["Fm"]] = rbind(c(0.8760, -0.2171, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                  c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                  c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                  c(0, 0, 0, 0, 0.0364, -0.0008, 0, 0, 0, 0, 0, 0), 
                  c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                  c(0, 0, 0, 0, 0, 0, -0.2965, -0.0657, 0, 0, 0, 0), 
                  c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                  c(0, 0, 0, 0, 0, 0, 0, 0, -0.3959, -0.1903, 0, 0),
                  c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0), 
                  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2436, 0.1281), 
                  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
ssm[["Fm"]] = array(ssm[["Fm"]], dim = c(dim(ssm[["Fm"]]), 2))
ssm[["Dm"]] = matrix(c(-1.5700, rep(0, 11)), nrow = nrow(ssm[["Fm"]]), ncol = 1)
ssm[["Dm"]] = array(ssm[["Dm"]], dim = c(dim(ssm[["Dm"]]), 2))
ssm[["Dm"]][1,, 2] = 0.2802
ssm[["Qm"]] = diag(c(1, 0, 0, 0, 0.0001, 0, 0.0001, 0, 0.0001, 0, 0.0001, 0))
ssm[["Qm"]] = array(ssm[["Qm"]], dim = c(dim(ssm[["Qm"]]), 2))
ssm[["Hm"]] = rbind(c(0.0058, -0.0033, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), 
                  c(0.0011, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
                  c(0.0051, -0.0033, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0), 
                  c(0.0012, -0.0005, 0.0001, 0.0002, 0, 0, 0, 0, 0, 0, 1, 0))
ssm[["Hm"]] = array(ssm[["Hm"]], dim = c(dim(ssm[["Hm"]]), 2))
ssm[["Am"]] = matrix(0, nrow = nrow(ssm[["Hm"]]), ncol = 1)
ssm[["Am"]] = array(ssm[["Am"]], dim = c(dim(ssm[["Am"]]), 2))
ssm[["Rm"]] = matrix(0, nrow = nrow(ssm[["Am"]]), ncol = nrow(ssm[["Am"]]))
ssm[["Rm"]] = array(ssm[["Rm"]], dim = c(dim(ssm[["Rm"]]), 2))
ssm[["B0"]] = matrix(c(rep(-4.60278, 4), 0, 0, 0, 0, 0, 0, 0, 0)) 
ssm[["B0"]] = array(ssm[["B0"]], dim = c(dim(ssm[["B0"]]), 2))
ssm[["B0"]][1:4,, 2] = rep(0.82146, 4)
ssm[["P0"]] = rbind(c(2.1775, 1.5672, 0.9002, 0.4483, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(1.5672, 2.1775, 1.5672, 0.9002, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(0.9002, 1.5672, 2.1775, 1.5672, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(0.4483, 0.9002, 1.5672, 2.1775, 0, 0, 0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0.0001, 0, 0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0.0001,  0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 0.0001, -0.0001, 0, 0, 0, 0),
                    c(0, 0, 0, 0, 0, 0, -0.0001, 0.0001, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, 0.0001, -0.0001, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, -0.0001, 0.0001, 0, 0), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0001, -0.0001), 
                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.0001, 0.0001))
ssm[["P0"]] = array(ssm[["P0"]], dim = c(dim(ssm[["P0"]]), 2))
  
#Log, difference and standardize the data
data[, c(vars) := lapply(.SD, log), .SDcols = c(vars)]
data[, c(vars) := lapply(.SD, function(x){
  x - shift(x, type = "lag", n = 1)
}), .SDcols = c(vars)]
data[, c(vars) := lapply(.SD, scale), .SDcols = c(vars)]
  
#Convert the data to an NxT matrix
yt = t(data[, c(vars), with = FALSE])
kf = kalman_filter(ssm, yt, smooth = TRUE)
}
}
