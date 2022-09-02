#' Kalman Filter
#' @param ssm list describing the state space model
#' @param yt matrix of data
#' @param Xo matrix of exogenous observation data
#' @param Xs matrix of exogenous state data
#' @param smooth boolean indication whether to run the backwards smoother
#' @useDynLib kalmanfilter, .registration=TRUE
#' @return list
#' @examples
#' #Nelson-Siegel Dynamic Factor Yield Curve Model
#' library(kalmanfilter)
#' library(data.table)
#' data(treasuries)
#' tau = unique(treasuries$maturity)
#' 
#' #Set up the state space model
#' ssm = list()
#' ssm[["Fm"]] = rbind(c(0.97, -0.03, -0.01), 
#'                     c(0.08, 0.79, -0.15), 
#'                     c(-0.15, 0.04, 0.88))
#' ssm[["Dm"]] = matrix(c(0.14, -0.1, 0.25), nrow = nrow(ssm[["Fm"]]), ncol = 1)
#' ssm[["Qm"]] = rbind(c(0.13, 0.12, -0.05), 
#'                     c(0.12, 0.25, -0.07), 
#'                     c(-0.05, -0.07, 1.02))
#' ssm[["Hm"]] = cbind(rep(1, 11),
#'                     -(1 - exp(-tau*0.04))/(tau*0.04), 
#'                     (1 - exp(-tau*0.04))/(tau*0.04) - exp(-tau*0.04))
#' ssm[["Am"]] = matrix(0, nrow = length(tau), ncol = 1)
#' ssm[["Rm"]] = diag(c(0.01, 0.00, 0.012, 0.01, 0.01, 0.00, 
#'                      0.01, 0.01, 0.01, 0.03, 0.04))
#' ssm[["betaO"]] = matrix(0, nrow = length(tau), ncol = 1)
#' ssm[["betaS"]] = matrix(0, nrow = nrow(ssm[["Fm"]]), ncol = 1)
#' ssm[["B0"]] = matrix(c(5.72, -0.80, 1.51), nrow = nrow(ssm[["Fm"]]), ncol = 1)
#' ssm[["P0"]] = diag(rep(0.01, nrow(ssm[["Fm"]])))
#' 
#' #Convert to an NxT matrix
#' yt = dcast(treasuries, "date ~ maturity", value.var = "value")
#' yt = t(yt[, 2:ncol(yt)])
#' kalman_filter(ssm, yt)
#' @export
kalman_filter = function(ssm, yt, Xo = NULL, Xs = NULL, smooth = FALSE){
  if(!is.list(ssm)){
    stop("ssm must be a list")
  }else{
    if(!all(c("Fm", "Dm", "Qm", "Hm", "Am", "Rm", "betaO", "betaS", "B0", "P0") %in% names(ssm))){
      stop("Fm, Dm, Qm, Hm, Am, Rm, betaO, betaS, B0, P0 must be the names of the matrices in ssm")
    }
    if(any(sapply(ssm, is.matrix) == FALSE)){
      stop("All elements of ssm must be matrices.")
    }
    if(length(unique(sapply(c("Fm", "Dm", "Qm", "betaS", "B0", "P0"), function(x){nrow(ssm[[x]])}))) != 1){
      stop("Fm, Dm, Qm, betaS, B0, P0 must have the same nrow")
    }
    if(length(unique(sapply(c("Hm", "Am", "Rm", "betaO"), function(x){nrow(ssm[[x]])}))) != 1){
      stop("Hm, Am, Rm, and betaO must have the same nrow")
    }
    if(any(sapply(c("P0", "Rm", "Qm"), function(x){nrow(ssm[[x]]) == ncol(ssm[[x]])}) == FALSE)){
      stop("nrow and ncol of P0, Rm, and Qm must be the same")
    }
    if(any(sapply(c("B0", "Dm", "Am"), function(x){ncol(ssm[[x]])}) != 1)){
      stop("ncol of B0, Am, and Dm must be 1")
    }
    if(ncol(ssm[["Hm"]]) != nrow(ssm[["Fm"]])){
      stop("ncol of Hm must equal nrow of Fm")
    }
  }
  if(!is.matrix(yt)){
    stop("yt must be a matrix")
  }
  if(is.null(Xo)){
    Xo = matrix(0, nrow = 1, ncol = ncol(yt))
  }else{
    if(!is.matrix(Xo)){
      stop("Xo must be a matrix.")
    }
    if(ncol(Xo) != ncol(yt)){
      stop("ncol of Xo must equal ncol yt")
    }
  }
  if(is.null(Xs)){
    Xs = matrix(0, nrow = 1, ncol = ncol(yt))
  }else{
    if(!is.matrix(Xs)){
      stop("Xs must be a matrix.")
    }
    if(ncol(Xs) != ncol(yt)){
      stop("ncol of Xs must equal ncol yt")
    }
  }
  if(!is.logical(smooth)){
    stop("smooth must be logical")
  }
  
  #Solve for the best initial starting values of the factors
  objective = function(par, yt){
    B1 = ssm[["Dm"]] + ssm[["Fm"]] %*% par + ssm[["betaS"]] %*% Xs[, 1]
    return(sum((yt[, 1] - (ssm[["Am"]] + ssm[["Hm"]]%*% B1 + ssm[["betaO"]] %*% Xo[, 1]))^2, na.rm = TRUE))
  }
  out = stats::optim(ssm[["B0"]], fn = objective, method = "BFGS", 
              control = list(trace = FALSE), yt = yt)
  ssm[["B0"]] = out$par
  return(filter(ssm, yt, Xo, Xs, smooth))
}

#' Kalman Likelihood
#' @param ssm list describing the state space model
#' @param yt matrix of data
#' @param Xo matrix of exogenous observation data
#' @param Xs matrix of exogenous state data
#' @param w matrix of weights
#' @useDynLib kalmanfilter, .registration=TRUE
#' @return numeric
#' @examples
#' #Nelson-Siegel Dynamic Factor Yield Curve Model
#' library(kalmanfilter)
#' library(data.table)
#' data(treasuries)
#' tau = unique(treasuries$maturity)
#' 
#' #Set up the state space model
#' ssm = list()
#' ssm[["Fm"]] = rbind(c(0.97, -0.03, -0.01), 
#'                     c(0.08, 0.79, -0.15), 
#'                     c(-0.15, 0.04, 0.88))
#' ssm[["Dm"]] = matrix(c(0.14, -0.1, 0.25), nrow = nrow(ssm[["Fm"]]), ncol = 1)
#' ssm[["Qm"]] = rbind(c(0.13, 0.12, -0.05), 
#'                     c(0.12, 0.25, -0.07), 
#'                     c(-0.05, -0.07, 1.02))
#' ssm[["Hm"]] = cbind(rep(1, 11),
#'                     -(1 - exp(-tau*0.04))/(tau*0.04), 
#'                     (1 - exp(-tau*0.04))/(tau*0.04) - exp(-tau*0.04))
#' ssm[["Am"]] = matrix(0, nrow = length(tau), ncol = 1)
#' ssm[["Rm"]] = diag(c(0.01, 0.00, 0.012, 0.01, 0.01, 0.00, 
#'                      0.01, 0.01, 0.01, 0.03, 0.04))
#' ssm[["betaO"]] = matrix(0, nrow = length(tau), ncol = 1)
#' ssm[["betaS"]] = matrix(0, nrow = nrow(ssm[["Fm"]]), ncol = 1)
#' ssm[["B0"]] = matrix(c(5.72, -0.80, 1.51), nrow = nrow(ssm[["Fm"]]), ncol = 1)
#' ssm[["P0"]] = diag(rep(0.01, nrow(ssm[["Fm"]])))
#' 
#' #Convert to an NxT matrix
#' yt = dcast(treasuries, "date ~ maturity", value.var = "value")
#' yt = t(yt[, 2:ncol(yt)])
#' kalman_lik(ssm, yt)
#' @export
kalman_lik = function(ssm, yt, Xo = NULL, Xs = NULL, w = NULL){
  if(!is.list(ssm)){
    stop("ssm must be a list")
  }else{
    if(!all(c("Fm", "Dm", "Qm", "Hm", "Am", "Rm", "betaO", "betaS", "B0", "P0") %in% names(ssm))){
      stop("Fm, Dm, Qm, Hm, Am, Rm, betaO, betaS, B0, P0 must be the names of the matrices in ssm")
    }
    if(any(sapply(ssm, is.matrix) == FALSE)){
      stop("All elements of ssm must be matrices.")
    }
    if(length(unique(sapply(c("Fm", "Dm", "Qm", "betaS", "B0", "P0"), function(x){nrow(ssm[[x]])}))) != 1){
      stop("Fm, Dm, Qm, betaS, B0, P0 must have the same nrow")
    }
    if(length(unique(sapply(c("Hm", "Am", "Rm", "betaO"), function(x){nrow(ssm[[x]])}))) != 1){
      stop("Hm, Am, Rm, and betaO must have the same nrow")
    }
    if(any(sapply(c("P0", "Rm", "Qm"), function(x){nrow(ssm[[x]]) == ncol(ssm[[x]])}) == FALSE)){
      stop("nrow and ncol of P0, Rm, and Qm must be the same")
    }
    if(any(sapply(c("B0", "Dm", "Am"), function(x){ncol(ssm[[x]])}) != 1)){
      stop("ncol of B0, Am, and Dm must be 1")
    }
    if(ncol(ssm[["Hm"]]) != nrow(ssm[["Fm"]])){
      stop("ncol of Hm must equal nrow of Fm")
    }
  }
  if(!is.matrix(yt)){
    stop("yt must be an NxT matrix")
  }
  if(is.null(Xo)){
    Xo = matrix(0, nrow = 1, ncol = ncol(yt))
  }else{
    if(!is.matrix(Xo)){
      stop("Xo must be a matrix.")
    }
    if(ncol(Xo) != ncol(yt)){
      stop("ncol of Xo must equal ncol yt")
    }
  }
  if(is.null(Xs)){
    Xs = matrix(0, nrow = 1, ncol = ncol(yt))
  }else{
    if(!is.matrix(Xs)){
      stop("Xs must be a matrix.")
    }
    if(ncol(Xs) != ncol(yt)){
      stop("ncol of Xs must equal ncol yt")
    }
  }
  if(is.null(w)){
    w = matrix(1, nrow = ncol(yt), ncol = 1)
  }else{
    if(!is.matrix(w)){
      stop("w must be a Tx1 matrix")
    }
    if(nrow(w) != ncol(yt)){
      stop("nrow of w must equal ncol of yt")
    }
    if(ncol(w) != 1){
      stop("ncol of w must equal 1")
    }
  }
  return(likelihood(ssm, yt, Xo, Xs, w))
}

