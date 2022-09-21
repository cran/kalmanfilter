#Nelson-Siegel Dynamic Factor Yield Curve Model
library(kalmanfilter)
library(data.table)
data(treasuries)
tau = unique(treasuries$maturity)

#Set up the state space model
ssm = list()
ssm[["Fm"]] = rbind(c(0.9720, -0.0209, -0.0061), 
                    c(0.1009 , 0.8189, -0.1446), 
                    c(-0.1226, 0.0192, 0.8808))
ssm[["Dm"]] = matrix(c(0.1234, -0.2285, 0.2020), nrow = nrow(ssm[["Fm"]]), ncol = 1)
ssm[["Qm"]] = rbind(c(0.1017, 0.0937, 0.0303), 
                    c(0.0937, 0.2267, 0.0351), 
                    c(0.0303, 0.0351, 0.7964))
ssm[["Hm"]] = cbind(rep(1, 11),
                     -(1 - exp(-tau*0.0423))/(tau*0.0423), 
                     (1 - exp(-tau*0.0423))/(tau*0.0423) - exp(-tau*0.0423))
ssm[["Am"]] = matrix(0, nrow = length(tau), ncol = 1)
ssm[["Rm"]] = diag(c(0.0087, 0, 0.0145, 0.0233, 0.0176, 0.0073, 
                     0, 0.0016, 0.0035, 0.0207, 0.0210))
ssm[["B0"]] = matrix(c(5.9030, -0.7090, 0.8690), nrow = nrow(ssm[["Fm"]]), ncol = 1)
ssm[["P0"]] = diag(rep(0.0001, nrow(ssm[["Fm"]])))

#Convert to an NxT matrix
yt = dcast(treasuries, "date ~ maturity", value.var = "value")
yt = t(yt[, 2:ncol(yt)])

test_that("kalman filter", {
  expect_equal(class(kalman_filter(ssm, yt)), "list")
})

ssm[["betaO"]] = matrix(0, nrow = nrow(yt), ncol = 1)
ssm[["betaS"]] = matrix(0, nrow = nrow(ssm[["Fm"]]), ncol = 1)
Xo = matrix(rnorm(ncol(yt)), ncol = ncol(yt), nrow = 1)
Xs = matrix(rnorm(ncol(yt)), ncol = ncol(yt), nrow = 1)

test_that("kalman filter exo", {
  expect_equal(class(kalman_filter(ssm, yt, Xo, Xs)), "list")
})

weights = matrix(rnorm(ncol(yt)), ncol = 1, nrow = ncol(yt))

test_that("kalman filter weighted", {
  expect_equal(class(kalman_filter(ssm, yt, weight = weights)), "list")
})
