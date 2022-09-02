## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE, include = TRUE, message = FALSE, warning = FALSE, eval = FALSE
)

## -----------------------------------------------------------------------------
#  lnl = kalman_lik(ssm, yt, Xo, Xs, w)

## -----------------------------------------------------------------------------
#  f = kalman_filter(ssm, yt, Xo, Xs, smooth)

## -----------------------------------------------------------------------------
#  library(kalmanfilter)
#  library(maxLik)
#  library(ggplot2)
#  library(data.table)
#  library(gridExtra)
#  data(treasuries)
#  tau = unique(treasuries$maturity)
#  
#  #State space model for the Nelson-Siegel dynamic factor yield curve
#  yc_ssm = function(par, tau){
#    #Set factor names
#    tau = unique(tau)
#    factors = c("l", "s", "c")
#  
#    #Loading on the slope factor
#    Ls = function(par, tau){
#      return(-(1 - exp(-tau*par["lambda"]))/(tau*par["lambda"]))
#    }
#  
#    #Loading on the curvature factor
#    Lc = function(par, tau){
#      return(-Ls(par["lambda"], tau) - exp(-tau*par["lambda"]))
#    }
#  
#    ########## Define the state space model
#    #Transition matrix
#    Fm = matrix(par[grepl("phi_", names(par))], nrow = 3, ncol = 3,
#                dimnames = list(paste0(factors, "_t0"),
#                                paste0(factors, "_t1")))
#    Fm[is.na(Fm)] = 0
#  
#    #Transition exogenous matrix
#    betaS = matrix(par[grepl("betaS_", names(par))], nrow = 3,
#                   dimnames = list(rownames(Fm),
#                                   unique(sapply(names(par)[grepl("betaS_", names(par))], function(x){
#                                     s = strsplit(x, "\\.")[[1]]
#                                     return(paste(s[2:length(s)], collapse = "."))
#                                   }))))
#  
#    #Transition intercept matrix
#    Dm = matrix(par[grepl("D_", names(par))], ncol = 1,
#                dimnames = list(rownames(Fm), NULL))
#  
#    #Transition error covariance matrix
#    Qm = matrix(par[unlist(lapply(factors, function(x){paste0("sig_", x, factors)}))], nrow = 3, ncol = 3)
#    Qm[lower.tri(Qm)] = Qm[upper.tri(Qm)]
#    dimnames(Qm) = list(rownames(Fm), rownames(Fm))
#  
#    #Observation equation matrix
#    n = length(tau)
#    Hm = matrix(NA, nrow = n, ncol = 3,
#                dimnames = list(as.character(tau), rownames(Fm)))
#    tau = as.numeric(tau)
#    if("l" %in% factors){
#      Hm[, "l_t0"] = 1
#    }
#    if("s" %in% factors){
#      Hm[, "s_t0"] = Ls(par, tau)
#    }
#    if("c" %in% factors){
#      Hm[, "c_t0"] = Lc(par, tau)
#    }
#  
#    betaO = matrix(par[grepl("betaO_", names(par))], nrow = n,
#                   dimnames = list(rownames(Hm),
#                                   unique(sapply(names(par)[grepl("betaO_", names(par))], function(x){
#                                     s = strsplit(x, "\\.")[[1]]
#                                     return(paste(s[2:length(s)], collapse = "."))
#                                   }))))
#  
#    #Observation error variance covariance matrix
#    Rm = diag(sum(grepl("sig_\\d+", names(par))))*par[grepl("sig_\\d+", names(par))]^2
#    dimnames(Rm) = list(rownames(Hm), rownames(Hm))
#  
#    #Observation intercept matrix
#    Am = matrix(par[grepl("A_", names(par))], nrow = nrow(Hm), ncol = 1,
#                dimnames = list(rownames(Hm), NULL))
#  
#    #Initial guess for the unobserved vector
#    B0 = matrix(par[names(par) %in% c("level", "slope", "curvature")], ncol = 1, nrow = nrow(Fm),
#                dimnames = list(rownames(Fm), NULL))
#  
#    #Initial guess for the variance of the unobserved vector
#    P0 = diag(par["P0"]^2, nrow = nrow(Qm), ncol = ncol(Qm))
#  
#    return(list(Fm = Fm, Dm = Dm, Qm = Qm,
#                Hm = Hm, Am = Am, Rm = Rm,
#                betaO = betaO, betaS = betaS,
#                B0 = B0, P0 = P0))
#  }
#  
#  #Set the initial values
#  init = c(level = 5.73, slope = -0.46, curvature = 0.67, lambda = 0.04,
#           D_l = 0.14, D_s = -0.08, D_c = 0.24,
#           phi_ll = 0.97, phi_ls = -0.02, phi_lc = 0.01,
#           phi_sl = 0.082, phi_ss = 0.79, phi_sc = -0.14,
#           phi_cl = -0.15, phi_cs = 0.05, phi_cc = 0.88,
#           betaS_l.X = 0, betaS_s.X = 0, betaS_c.X = 0,
#           sig_ll = 0.13,
#           sig_sl = 0.12, sig_ss = 0.24,
#           sig_cl = -0.05, sig_cs = -0.08, sig_cc = 1.03,
#           A_1 = 0, A_3 = 0, A_6 = 0, A_12 = 0, A_24 = 0,
#           A_36 = 0, A_60 = 0, A_84 = 0, A_120 = 0, A_240 = 0, A_360 = 0,
#           betaO_1.X = 0, betaO_3.X = 0, betaO_6.X = 0, betaO_12.X = 0,
#           betaO_24.X = 0, betaO_36.X = 0, betaO_60.X = 0, betaO_84.X = 0,
#           betaO_120.X = 0, betaO_240.X = 0, betaO_360.X = 0,
#           sig_1 = 0.08, sig_3 = 0.01, sig_6 = 0.11, sig_12 = 0.12,
#           sig_24 = 0.07, sig_36 = 0.01, sig_60 = 0.06, sig_84 = 0.07,
#           sig_120 = 0.04, sig_240 = 0.18, sig_360 = 0.2,
#           P0 = 0.01)
#  
#  #Define the fixed values: not using exogenous data or observation intercept
#  fixed = c("betaS_l.x", "betaS_s.X", "betaS_c.X",
#            "A_1", "A_3", "A_6", "A_12", "A_24", "A_36", "A_60",
#            "A_84", "A_120", "A_240", "A_360",
#            "betaO_1.X", "betaO_3.X", "betaO_6.X", "betaO_12.X",
#            "betaO_24.X", "betaO_36.X", "betaO_60.X", "betaO_84.X",
#            "betaO_120.X", "betaO_240.X", "betaO_360.X")
#  
#  #Set up the constraints: lambda and all standard deviation parameters must be positive
#  ineqA = matrix(0, nrow = 15, ncol = length(init), dimnames = list(NULL, names(init)))
#  diag(ineqA[, c("lambda", "sig_ll", "sig_ss", "sig_cc", paste0("sig_", tau))]) = 1
#  ineqB = matrix(0, nrow = nrow(ineqA), ncol = 1)
#  
#  #Convert to an NxT matrix
#  yt = dcast(treasuries, "date ~ maturity", value.var = "value")
#  yt = t(yt[, 2:ncol(yt)])
#  
#  #Set the objective function
#  objective = function(par, data){
#    ssm = yc_ssm(par, unique(data$maturity))
#    return(kalman_lik(ssm, yt))
#  }
#  
#  #Solve the model
#  solve = maxLik(logLik = objective,
#                 start = init, method = "BFGS",
#                 finalHessian = FALSE, hess = NULL,
#                 control = list(printLevel = 2, iterlim = 10000),
#                 constraints = list(ineqA = ineqA, ineqB = ineqB),
#                 data = treasuries, fixed = fixed)
#  
#  #Get the estimated state space model
#  ssm = yc_ssm(solve$estimate, tau)
#  
#  #Filter the data with the model
#  filter = kalman_filter(ssm, yt, smooth = TRUE)
#  
#  #Get the estimated unobserved factors
#  B_tt = data.table(t(filter[["B_tt"]]))[, "date" := unique(treasuries$date)]
#  colnames(B_tt)[1:3] = c("level", "slope", "curvature")
#  
#  #Get the estimated yields
#  y_tt = data.table(t(filter[["y_tt"]]))[, "date" := unique(treasuries$date)]
#  colnames(y_tt)[1:(ncol(y_tt) - 1)] = tau
#  
#  #Plot the data
#  g1 = ggplot(treasuries, id.vars = "date") +
#    geom_line(aes(x = date, y = value, group = factor(maturity), color = factor(maturity))) +
#    theme_minimal() + theme(legend.position = "bottom") +
#    labs(title = "Actual Treasury Yields", x = "", y = "value") +
#    guides(color = guide_legend(title = NULL))
#  g2 = ggplot(melt(y_tt, id.vars = "date")) +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    theme_minimal() + theme(legend.position = "bottom") +
#    labs(title = "Estimated Treasury Yields", x = "", y = "value") +
#    guides(color = guide_legend(title = NULL))
#  g3 = ggplot(melt(B_tt, id.vars = "date")) +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    theme_minimal() + theme(legend.position = "bottom") +
#    labs(title = "Estimated Factors", x = "", y = "value") +
#    guides(color = guide_legend(title = NULL))
#  grid.arrange(g1, g2, g3, layout_matrix = rbind(c(1, 2), c(3, 3)))

