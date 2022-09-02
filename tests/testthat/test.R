#Nelson-Siegel Dynamic Factor Yield Curve Model
library(kalmanfilter)
library(data.table)
data(treasuries)
tau = unique(treasuries$maturity)

#Set up the state space model
ssm = list()
ssm[["Fm"]] = rbind(c(0.97, -0.03, -0.01), 
                   c(0.08, 0.79, -0.15), 
                   c(-0.15, 0.04, 0.88))
ssm[["Dm"]] = matrix(c(0.14, -0.1, 0.25), nrow = nrow(ssm[["Fm"]]), ncol = 1)
ssm[["Qm"]] = rbind(c(0.13, 0.12, -0.05), 
                   c(0.12, 0.25, -0.07), 
                   c(-0.05, -0.07, 1.02))
ssm[["Hm"]] = cbind(rep(1, 11),
                   -(1 - exp(-tau*0.04))/(tau*0.04), 
                   (1 - exp(-tau*0.04))/(tau*0.04) - exp(-tau*0.04))
ssm[["Am"]] = matrix(0, nrow = length(tau), ncol = 1)
ssm[["Rm"]] = diag(c(0.01, 0.00, 0.012, 0.01, 0.01, 0.00, 
                    0.01, 0.01, 0.01, 0.03, 0.04))
ssm[["betaO"]] = matrix(0, nrow = length(tau), ncol = 1)
ssm[["betaS"]] = matrix(0, nrow = nrow(ssm[["Fm"]]), ncol = 1)
ssm[["B0"]] = matrix(c(5.72, -0.80, 1.51), nrow = nrow(ssm[["Fm"]]), ncol = 1)
ssm[["P0"]] = diag(rep(0.01, nrow(ssm[["Fm"]])))

#Convert to an NxT matrix
yt = dcast(treasuries, "date ~ maturity", value.var = "value")
yt = t(yt[, 2:ncol(yt)])

test_that("kalman likelihood", {
  expect_equal(class(kalman_lik(ssm, yt)), "numeric")
})

test_that("kalman filter", {
  expect_equal(class(kalman_filter(ssm, yt)), "list")
})

