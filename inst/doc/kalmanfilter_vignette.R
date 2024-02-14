## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE, include = TRUE, message = FALSE, warning = FALSE, eval = FALSE
)

## -----------------------------------------------------------------------------
#  kf = kalman_filter(ssm, yt, Xo, Xs, w, smooth)

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
#    #Observation error variance covariance matrix
#    Rm = diag(sum(grepl("sig_\\d+", names(par))))*par[grepl("sig_\\d+", names(par))]^2
#    dimnames(Rm) = list(rownames(Hm), rownames(Hm))
#  
#    #Observation intercept matrix
#    Am = matrix(0, nrow = nrow(Hm), ncol = 1,
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
#                B0 = B0, P0 = P0))
#  }
#  
#  #Set the initial values
#  init = c(level = 5.9030, slope = -0.7090, curvature = 0.8690, lambda = 0.0423,
#           D_l = 0.1234, D_s = -0.2285, D_c = 0.2020,
#           phi_ll = 0.9720, phi_ls = 0.1009, phi_lc = -0.1226,
#           phi_sl = -0.0209, phi_ss = 0.8189, phi_sc = 0.0192,
#           phi_cl = -0.0061, phi_cs = -0.1446, phi_cc = 0.8808,
#           sig_ll = 0.1017,
#           sig_sl = 0.0937, sig_ss = 0.2267,
#           sig_cl = 0.0303, sig_cs = 0.0351, sig_cc = 0.7964,
#           sig_1 = 0.0934, sig_3 = 0.0001, sig_6 = 0.1206, sig_12 = 0.1525,
#           sig_24 = 0.1328, sig_36 = 0.0855, sig_60 = 0.0001, sig_84 = 0.0397,
#           sig_120 = 0.0595, sig_240 = 0.1438, sig_360 = 0.1450,
#           P0 = 0.0001)
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
#  objective = function(par){
#    ssm = yc_ssm(par, tau)
#    return(kalman_filter(ssm, yt,)$lnl)
#  }
#  
#  #Solve the model
#  solve = maxLik(logLik = objective, start = init, method = "BFGS",
#                 finalHessian = FALSE, hess = NULL,
#                 control = list(printLevel = 2, iterlim = 10000),
#                 constraints = list(ineqA = ineqA, ineqB = ineqB))
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
#    labs(title = "Actual Treasury Yields", x = "", y = "Value") +
#    guides(color = guide_legend(title = NULL))
#  g2 = ggplot(melt(y_tt, id.vars = "date")) +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    theme_minimal() + theme(legend.position = "bottom") +
#    labs(title = "Estimated Treasury Yields", x = "", y = "Value") +
#    guides(color = guide_legend(title = NULL))
#  g3 = ggplot(melt(B_tt, id.vars = "date")) +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    theme_minimal() + theme(legend.position = "bottom") +
#    labs(title = "Estimated Factors", x = "", y = "Value") +
#    guides(color = guide_legend(title = NULL))
#  grid.arrange(g1, g2, g3, layout_matrix = rbind(c(1, 2), c(3, 3)))

## -----------------------------------------------------------------------------
#  #Set the objective function
#  objective = function(par){
#    ssm = yc_ssm(par, tau)
#    for(x in names(ssm)){
#      if(!x %in% c("B0", "P0")){
#        ssm[[x]] = array(ssm[[x]], dim = c(dim(ssm[[x]]), ncol(yt)))
#      }
#    }
#    return(kalman_filter(ssm, yt)$lnl)
#  }
#  
#  #Solve the model
#  solve = maxLik(logLik = objective, start = init, method = "BFGS",
#                 finalHessian = FALSE, hess = NULL,
#                 control = list(printLevel = 2, iterlim = 10000),
#                 constraints = list(ineqA = ineqA, ineqB = ineqB))

## -----------------------------------------------------------------------------
#  library(kalmanfilter)
#  library(data.table)
#  library(maxLik)
#  library(ggplot2)
#  library(gridExtra)
#  data(sw_dcf)
#  data = sw_dcf[, colnames(sw_dcf) != "dcoinc", with = F]
#  vars = colnames(data)[colnames(data) != "date"]
#  
#  #State space model for the Stock and Watson Dynamic Common Factor model
#  dcf_ssm = function(par, yt){
#    #Get the parameters
#    vars = dimnames(yt)[which(unlist(lapply(dimnames(yt), function(x){!is.null(x)})))][[1]]
#    phi = par[grepl("phi", names(par))]
#    names(phi) = gsub("phi", "", names(phi))
#    gamma = par[grepl("gamma_", names(par))]
#    names(gamma) = gsub("gamma_", "", names(gamma))
#    psi = par[grepl("psi_", names(par))]
#    names(psi) = gsub("psi_", "", names(psi))
#    sig = par[grepl("sigma_", names(par))]
#    names(sig) = gsub("sigma_", "", names(sig))
#  
#    #Build the transition matrix
#    phi_dim = max(c(length(phi)), length(unique(sapply(strsplit(names(gamma), "\\."), function(x){x[2]}))))
#    psi_dim = sapply(unique(sapply(strsplit(names(psi), "\\."), function(x){x[1]})), function(x){
#      max(as.numeric(sapply(strsplit(names(psi)[grepl(paste0("^", x), names(psi))], "\\."), function(x){x[2]})))
#    })
#    Fm = matrix(0, nrow = phi_dim + length(psi), ncol = phi_dim + length(psi),
#                dimnames = list(
#                  c(paste0("ct.", 0:(phi_dim - 1)),
#                    unlist(lapply(names(psi_dim), function(x){paste0("e_", x, ".", 0:(psi_dim[x] - 1))}))),
#                  c(paste0("ct.", 1:phi_dim),
#                    unlist(lapply(names(psi_dim), function(x){paste0("e_", x, ".", 1:psi_dim[x])})))
#                ))
#    Fm["ct.0", paste0("ct.", names(phi))] = phi
#    for(i in 1:length(vars)){
#      Fm[paste0("e_", i, ".0"),
#         paste0("e_", names(psi)[grepl(paste0("^", i), names(psi))])] = psi[grepl(paste0("^", i), names(psi))]
#    }
#    diag(Fm[intersect(rownames(Fm), colnames(Fm)), intersect(rownames(Fm), colnames(Fm))]) = 1
#  
#    #Build the observation matrix
#    Hm = matrix(0, nrow = nrow(yt), ncol = nrow(Fm), dimnames = list(rownames(yt), rownames(Fm)))
#    for(i in 1:length(vars)){
#      Hm[i, paste0("ct.", sapply(strsplit(names(gamma)[grepl(paste0("^", i), names(gamma))], "\\."), function(x){x[2]}))] =
#        gamma[grepl(paste0("^", i), names(gamma))]
#    }
#    diag(Hm[, paste0("e_", 1:length(vars), ".0")]) = 1
#  
#    #Build the state covariance matrix
#    #Set the dynamic common factor standard deviation to 1
#    Qm = matrix(0, ncol = ncol(Fm), nrow = nrow(Fm), dimnames = list(rownames(Fm), rownames(Fm)))
#    Qm["ct.0", "ct.0"] = 1
#    for(i in 1:length(vars)){
#      Qm[paste0("e_", i, ".0"), paste0("e_", i, ".0")] = sig[names(sig) == i]^2
#    }
#  
#    #Build the observation equation covariance matrix
#    Rm = matrix(0, ncol = nrow(Hm), nrow = nrow(Hm), dimnames = list(vars, vars))
#  
#    #Transition equation intercept matrix
#    Dm = matrix(0, nrow = nrow(Fm), ncol = 1, dimnames = list(rownames(Fm), NULL))
#  
#    #Observation equation intercept matrix
#    Am = matrix(0, nrow = nrow(Hm), ncol = 1)
#  
#    #Initialize the filter for each state
#    B0 = matrix(0, nrow(Fm), 1)
#    P0 = diag(nrow(Fm))
#    dimnames(B0) = list(rownames(Fm), NULL)
#    dimnames(P0) = list(rownames(Fm), rownames(Fm))
#  
#    B0 = solve(diag(ncol(Fm)) - Fm) %*% Dm
#    VecP_ll = solve(diag(prod(dim(Fm))) - kronecker(Fm, Fm)) %*% matrix(as.vector(Qm), ncol = 1)
#    P0 = matrix(VecP_ll[, 1], ncol = ncol(Fm))
#  
#    return(list(B0 = B0, P0 = P0, Am = Am, Dm = Dm, Hm = Hm, Fm = Fm, Qm = Qm, Rm = Rm))
#  }
#  
#  #Log the data
#  data.log = copy(data)
#  data.log[, c(vars) := lapply(.SD, log), .SDcols = c(vars)]
#  
#  #Difference the data
#  data.logd = copy(data.log)
#  data.logd[, c(vars) := lapply(.SD, function(x){
#    x - shift(x, type = "lag", n = 1)
#  }), .SDcols = c(vars)]
#  
#  #Standardize the data
#  data.logds = copy(data.logd)
#  data.logds[, c(vars) := lapply(.SD, scale, scale = FALSE), .SDcols = c(vars)]
#  
#  #Transpose the data
#  yt = t(data.logds[, c(vars), with = FALSE])
#  
#  #Set the initial values
#  init = c(phi1 = 0.8588, phi2 = -0.1526,
#           psi_1.1 = -0.1079, psi_1.2 = -0.1415,
#           psi_2.1 = -0.3166, psi_2.2 = -0.0756,
#           psi_3.1 = -0.3994, psi_3.2 = -0.2028,
#           psi_4.1 = -0.0370, psi_4.2 = 0.3646,
#           gamma_1.0 = 0.0064, gamma_1.1 = -0.0020,
#           gamma_2.0 = 0.0018,
#           gamma_3.0 = 0.0059, gamma_3.1 = -0.0027,
#           gamma_4.0 = 0.0014, gamma_4.1 = -0.0004, gamma_4.2 = 0.0001, gamma_4.3 = 0.0004,
#           sigma_1 = 0.0049, sigma_2 = 0.0056, sigma_3 = 0.0077, sigma_4 = 0.0013)
#  
#  #Set the constraints
#  ineqA = matrix(0, nrow = 14,
#                 ncol = length(init), dimnames = list(NULL, names(init)))
#  #Stationarity constraints
#  ineqA[c(1, 2), c("phi1", "phi2")] = rbind(c(1, 1), c(-1, -1))
#  ineqA[c(3, 4), grepl("psi_1", colnames(ineqA))] = rbind(c(1, 1), c(-1, -1))
#  ineqA[c(5, 6), grepl("psi_2", colnames(ineqA))] = rbind(c(1, 1), c(-1, -1))
#  ineqA[c(7, 8), grepl("psi_3", colnames(ineqA))] = rbind(c(1, 1), c(-1, -1))
#  ineqA[c(9, 10), grepl("psi_4", colnames(ineqA))] = rbind(c(1, 1), c(-1, -1))
#  #Non-negativity constraints
#  diag(ineqA[c(11, 12, 13, 14), grepl("sigma_", colnames(ineqA))]) = 1
#  ineqB = matrix(c(rep(1, 10),
#                   rep(0, 4)), nrow = nrow(ineqA), ncol = 1)
#  
#  #Define the objective function
#  objective = function(par, yt){
#    ssm = dcf_ssm(par, yt)
#    return(kalman_filter(ssm, yt, smooth = FALSE)$lnl)
#  }
#  
#  #Solve the model
#  solve = maxLik(logLik = objective, start = init, method = "BFGS",
#                 finalHessian = FALSE, hess = NULL,
#                 control = list(printLevel = 2, iterlim = 10000),
#                 constraints = list(ineqA = ineqA, ineqB = ineqB),
#                 yt = yt)
#  
#  #Get the estimated state space model
#  ssm = dcf_ssm(solve$estimate, yt)
#  
#  #Get the column means and standard deviations
#  M = matrix(unlist(data.logd[, lapply(.SD, mean, na.rm = TRUE), .SDcols = c(vars)]),
#                 ncol = 1, dimnames = list(vars, NULL))
#  
#  #Get the coefficient matrices
#  Hm = ssm[["Hm"]]
#  Fm = ssm[["Fm"]]
#  
#  #Final K_t is approximation to steady state K
#  filter = kalman_filter(ssm, yt, smooth = T)
#  K = filter$K_t[,, dim(filter$K_t)[3]]
#  W = solve(diag(nrow(K)) - (diag(nrow(K)) - K %*% Hm) %*% Fm) %*% K
#  d = (W %*% M)[1, 1]
#  
#  #Get the intercept terms
#  gamma = Hm[, grepl("ct", colnames(Hm))]
#  D = M - gamma %*% matrix(rep(d, ncol(gamma)))
#  
#  #Initialize first element of the dynamic common factor
#  Y1 = t(data.log[, c(vars), with = F][1, ])
#  initC = function(par){
#    return(sum((Y1 - D - gamma %*% par)^2))
#  }
#  C10 = optim(par = Y1, fn = initC, method = "BFGS", control = list(trace = FALSE))$par[1]
#  Ctt = rep(C10, ncol(yt))
#  
#  #Build the rest of the level of the dynamic common factor
#  ctt = filter$B_tt[which(rownames(Fm) == "ct.0"), ]
#  for(j in 2:length(Ctt)){
#    Ctt[j] = ctt[j] + Ctt[j - 1] + c(d)
#  }
#  Ctt = data.table(date = data$date, dcf = Ctt, d.dcf = ctt)
#  
#  #Plot the outputs
#  g1 = ggplot(melt(data.log, id.vars = "date")[, "value" := scale(value), by = "variable"]) +
#    ggtitle("Data", subtitle = "Log Levels (Rescaled)") +
#    scale_y_continuous(name = "Value") +
#    scale_x_date(name = "") +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    theme_minimal() + theme(legend.position = "bottom") + guides(color = guide_legend(title = NULL))
#  
#  g2 = ggplot( melt(data.logds, id.vars = "date")) +
#    ggtitle("Data", subtitle = "Log Differenced & Standardized") +
#    scale_y_continuous(name = "Value") +
#    scale_x_date(name = "") +
#    geom_hline(yintercept = 0, color = "black") +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    theme_minimal() + theme(legend.position = "bottom") + guides(color = guide_legend(title = NULL))
#  
#  g3 = ggplot(melt(Ctt, id.vars = "date")[variable == "dcf", ]) +
#    ggtitle("Dynamic Common Factor", subtitle = "Levels") +
#    scale_x_date(name = "") +
#    geom_hline(yintercept = 0, color = "grey") +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    theme_minimal() + theme(legend.position = "bottom") + scale_color_manual(values = "black") +
#    guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
#    scale_y_continuous(name = "Value", limits = range(Ctt$dcf, na.rm = TRUE))
#  
#  g4 = ggplot(melt(Ctt, id.vars = "date")[variable == "d.dcf", ]) +
#    ggtitle("Dynamic Common Factor", subtitle = "Differenced") +
#    scale_x_date(name = "") +
#    geom_hline(yintercept = 0, color = "grey") +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    theme_minimal() + theme(legend.position = "bottom") + scale_color_manual(values = "black") +
#    guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
#    scale_y_continuous(name = "Value", limits = range(Ctt$d.dcf, na.rm = TRUE))
#  
#  grid.arrange(g1, g2, g3, g4, layout_matrix = matrix(c(1, 3, 2, 4), nrow = 2))

