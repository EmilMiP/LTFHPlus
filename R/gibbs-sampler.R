

#'
#' Gibbs Sampler for the truncated multivariate normal distribution
#'
#' @param n_sim number of draws after burn in.
#' @param sigma covariance matrix for the normal distribution
#' @param lower lower truncation for the normal distribution. It can be -Inf
#' @param upper upper truncation for the normal distribution. It can be Inf
#' @param fixed Are any entries fixed ? Logical vector with an entry for each dimension
#' @param ind Indices to return from the gibbs sampler. 1 corresponds to the genetic liability of the first phenotype. c(1,5) corresponds to the genetic liability of the first two phenotypes provided, but with no siblings in the model.
#' @param burn_in Number of iterations that count as burn in for the gibbs sampler.
#'
#' @return Returns the estimated genetic liabilities.
#'
#' @export
#' 

rtmvnorm.gibbs <- function(n_sim, sigma, lower, upper, fixed, ind = 1,  burn_in = 1000) {

  # Start value from support region,
  # may be lower or upper bound if they are finite,
  # if both are infinite, we take 0.
  x0 <- ifelse(is.finite(lower), lower, ifelse(is.finite(upper), upper, 0))

  # Pre-computations
  d <- nrow(sigma)
  sd <- rep(NA, d)
  P <- matrix(NA, d, d)
  for (j in 1:d) {
    P_j <- rep(0, d)
    P_j[-j] <- solve(sigma[-j, -j] , sigma[-j, j])
    P[, j] <- P_j
    sd[j] <- drop(sqrt(sigma[j, j] - crossprod(P_j, sigma[, j])))
  }
  # vector entries to return. Baseline is genetic liability for first phenotype. e.g. ind = c(1,5) for two phenotypes with no siblings.
  to_return <- rep(-1, d)
  to_return[ind] <- seq_along(ind) - 1L
  
  
#  print(P)
  # Actual Gibbs sampler
  rtmvnorm_gibbs_cpp(P, sd, lower, upper, fixed, to_return, x0, n_sim, burn_in)
}

if (FALSE) {

  sd(replicate(100, {
    colMeans(tmvtnorm::rtmvnorm(10e3, mean = rep(0, 4), H = inv_cov,
                                lower = lims[, "lower"], upper = lims[, "upper"],
                                algorithm = "gibbs", burn.in.samples = 1000))[1]
  }))
  # 0.01322415 / 0.01394558

  test <- rtmvnorm.gibbs(10e3, sigma = cov, lower = lims[, "lower"], upper = lims[, "upper"])
  plot(test)

  sd(replicate(100, {
    mean(rtmvnorm.gibbs(10e3, sigma = cov, lower = lims[, "lower"], upper = lims[, "upper"]))
  }))
  # 0.008360882 / 0.008005332

  replicate(10, {
    rtmvnorm.gibbs(10e3, sigma = cov, lower = lims[, "lower"], upper = lims[, "upper"])
  })

  hist(replicate(100000, 2 + 3 * qnorm(runif(1, 0.2, 0.7))))
  hist(replicate(100000, qnorm(runif(1, 0.2, 0.7), mean = 2, sd = 3)))

  microbenchmark::microbenchmark(
    tmvtnorm::rtmvnorm(10e3, mean = rep(0, 4), H = inv_cov,
                       lower = lims[, "lower"], upper = lims[, "upper"],
                       algorithm = "gibbs", burn.in.samples = 1000),
    rtmvnorm.gibbs(10e3, sigma = cov, lower = lims[, "lower"], upper = lims[, "upper"]),
    times = 10
  )

}

