################################################################################

context("GIBBS_SAMPLER")

skip_if_not_installed("tmvtnorm")

################################################################################

test_that("Gibbs sampler works", {
  
  cov <- matrix(c(1, 0.2, 0.2, 0.5), 2)
  lower <- c(-Inf, 0)
  upper <- c(0, Inf)
  
  true <- tmvtnorm::rtmvnorm(100e3, mean = rep(0, 2), H = solve(cov),
                             lower = lower, upper = upper)
  
  test <- rtmvnorm.gibbs(100e3, sigma = cov, lower = lower, upper = upper, ind = 1:2)
  
  expect_true(all(test[, 1] <= 0))
  expect_true(all(test[, 2] >= 0))
  expect_equal(colMeans(test), colMeans(true), tolerance = 0.01)
  expect_equal(cov(test), cov(true), tolerance = 0.01)
})

################################################################################

test_that("Gibbs sampler works with any cov", {
  
  cov <- matrix(c(1, 0.2, 0.2, 0.5), 2) * 100
  lower <- c(-Inf, 10)
  upper <- c(0, 20)
  
  true <- tmvtnorm::rtmvnorm(100e3, mean = rep(0, 2), H = solve(cov),
                             lower = lower, upper = upper, algorithm = "gibbs")
  
  test <- rtmvnorm.gibbs(100e3, sigma = cov, lower = lower, upper = upper, ind = 1:2)
  
  expect_true(all(test[, 1] <= 0))
  expect_true(all(test[, 2] >= 10))
  expect_true(all(test[, 2] <= 20))
  expect_equal(colMeans(test), colMeans(true), tolerance = 0.01)
  expect_equal(cov(test), cov(true), tolerance = 0.1)
})

################################################################################

test_that("Gibbs sampler with fixed values works", {
  
  cov <- matrix(c(1, 0.2, 0.2, 0.5), 2)
  lower <- c(-Inf, 0)
  upper <- c(0, 1e-5)
  
  true <- tmvtnorm::rtmvnorm(100e3, mean = rep(0, 2), H = solve(cov),
                             lower = lower, upper = upper, algorithm = "gibbs")
  
  test <- rtmvnorm.gibbs(100e3, sigma = cov, ind = 1:2, 
                         lower = lower, upper = upper)
  
  expect_true(all(test[, 1] <= 0))
  val <- test[1, 2]
  expect_gt(val, 0)
  expect_lt(val, 1e-5)
  expect_true(all(test[, 2] == val))
  expect_equal(colMeans(test), colMeans(true), tolerance = 0.01)
  expect_equal(cov(test), cov(true), tolerance = 0.01)
})

################################################################################
