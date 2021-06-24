################################################################################

#' Gibbs Sampler for the truncated multivariate normal distribution
#'
#' @param n_sim Number of draws from the Gibbs sampler (after burn-in).
#' @param sigma Covariance matrix for the multivariate normal distribution.
#' @param lower Lower truncation for the normal distribution. It can be `-Inf`.
#' @param upper Upper truncation for the normal distribution. It can be `Inf`.
#' @param fixed Should variables be fixed? Logical vector with an entry for each 
#'   dimension. Default is `TRUE` when `upper` is equal to `lower`.
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
#' @examples
#' samp <- rtmvnorm.gibbs(10e3, sigma = matrix(c(1, 0.2, 0.2, 0.5), 2),
#'                        lower = c(-Inf, 0), upper = c(0, Inf), ind = 1:2)
#' 
rtmvnorm.gibbs <- function(n_sim, sigma, lower, upper, 
                           fixed = (lower == upper), 
                           ind = 1, burn_in = 1000) {
  
  #stopifnot(all(upper >= lower))
  
  # Start with medians of univariate distributions
  sd0 <- sqrt(diag(sigma))
  sum_p0 <- stats::pnorm(lower, sd = sd0) + stats::pnorm(upper, sd = sd0)
  x0 <- stats::qnorm(sum_p0 / 2, sd = sd0) 

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
  
  # Vector entries to return
  to_return <- rep(-1, d)
  to_return[ind] <- seq_along(ind) - 1L
  
  # Actual Gibbs sampler
  rtmvnorm_gibbs_cpp(P, sd, lower, upper, fixed, to_return, x0, n_sim, burn_in)
}

################################################################################
