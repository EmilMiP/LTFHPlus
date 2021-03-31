################################################################################

#' Gibbs Sampler for the truncated multivariate normal distribution
#'
#' @param n_sim Number of draws from the Gibbs sampler (after burn-in).
#' @param sigma Covariance matrix for the multivariate normal distribution.
#' @param lower Lower truncation for the normal distribution. It can be `-Inf`.
#' @param upper Upper truncation for the normal distribution. It can be `Inf`.
#' @param fixed Should variables be fixed? Logical vector with an entry for each dimension.
#' @param ind Indices of variables to return from the Gibbs sampler. 
#'   Default is `1` and corresponds to the first variable (usually the genetic 
#'   liability of the first phenotype). Use `c(1, 5)` to return both the first and 
#'   fifth variables.
#' @param burn_in Number of iterations that count as burn in for the Gibbs sampler.
#'
#' @return The sampling values from the Gibbs sampler for the variables `ind`.
#' 
#' @references
#' Kotecha, J. H., & Djuric, P. M. (1999, March). Gibbs sampling approach for 
#' generation of truncated multivariate gaussian random variables. In 1999 IEEE 
#' International Conference on Acoustics, Speech, and Signal Processing. 
#' Proceedings. ICASSP99 (Cat. No. 99CH36258) (Vol. 3, pp. 1757-1760). IEEE.
#' \doi{10.1109/ICASSP.1999.756335}
#' 
#' Wilhelm, S., & Manjunath, B. G. (2010). tmvtnorm: A package for the truncated 
#' multivariate normal distribution. The R Journal. \doi{10.32614/RJ-2010-005}
#'
#' @export
#' 
rtmvnorm.gibbs <- function(n_sim, sigma, lower, upper, fixed, ind = 1, burn_in = 1000) {
  
  # Vector entries to return
  to_return <- rep(-1, d)
  to_return[ind] <- seq_along(ind) - 1L
  
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
  
  # Actual Gibbs sampler
  rtmvnorm_gibbs_cpp(P, sd, lower, upper, fixed, to_return, x0, n_sim, burn_in)
}

################################################################################

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

################################################################################
