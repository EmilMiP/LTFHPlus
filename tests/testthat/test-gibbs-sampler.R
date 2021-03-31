context("GIBBS_SAMPLER")

test_that("Gibbs sampler works", {
  
  skip_if_not_installed("tmvtnorm")
  
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
