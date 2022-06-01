#' Gibbs Sampler for the truncated multivariate normal distribution
#' 
#' \code{rtmvnorm.gibbs} implements Gibbs sampler for the truncated 
#' multivariate normal distribution with covariance matrix \code{covmat}.
#' 
#' Given a covariance matrix \code{covmat} and lower and upper cutoff points,
#' the function \code{rtmvnorm.gibbs} can be used to perform Gibbs sampler on a truncated 
#' multivariable normal distribution. It is possible to specify which variables 
#' to return from the Gibbs sampler, making it convenient to use when estimating
#' only the full liability or the genetic component of the full liability.
#'
#' @param n_sim A positive number representing the number of draws from the 
#' Gibbs sampler after burn-in.. Defaults to 1e+05.
#' @param covmat A symmetric and numeric matrix representing the covariance
#' matrix for the multivariate normal distribution.
#' @param lower A number or numeric vector representing the lower cutoff point(s) for the 
#' truncated normal distribution. The length of lower must be 1 or equal 
#' to the dimension of the multivariable normal distribution.
#' Defaults to \code{-Inf}.
#' @param upper A number or numeric vector representing the upper cutoff point(s) for the 
#' truncated normal distribution. Must be greater or equal to lower.
#' In addition the length of upper must be 1 or equal to the dimension 
#' of the multivariable normal distribution.
#' Defaults to \code{Inf}.
#' @param fixed A logical scalar or a logical vector indicating which
#' variables to fix. If \code{fixed} is a vector, it must have the same length as 
#' lower and upper. Defaults to \code{TRUE} when \code{lowert} is equal to 
#' \code{upper} and \code{FALSE} otherwise.
#' @param out An integer or numeric vector indicating which variables should be returned
#' from the Gibbs sampler. If out = c(1), the first variable (usually the genetic 
#' component of the full liability of the first phenotype) is estimated and returned. 
#' If out = c(2), the second variable (usually full liability) is estimated and returned. 
#' If out = c(1,2), both the first and the second variable are estimated and returned. 
#' Defaults to c(1).
#' @param burn_in A number of iterations that count as burn in for the Gibbs sampler.
#' Must be non-negative. Defaults to 1000.
#'
#' @return If \code{covmat} is a symmetric and numeric matrix, if \code{n_sim} and 
#' \code{burn_in} are positive/non-negative numbers, if \code{out} is a numeric vector and 
#' \code{lower}, \code{upper} and \code{fixed} are numbers or vectors of the same length 
#' and the required format, \code{rtmvnorm.gibbs} returns the sampling values 
#' from the Gibbs sampler for all variables specified in \code{out}.
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
#' @examples
#' samp <- rtmvnorm.gibbs(10e3, covmat = matrix(c(1, 0.2, 0.2, 0.5), 2),
#'                        lower = c(-Inf, 0), upper = c(0, Inf), out = 1:2)
#' @export
rtmvnorm.gibbs <- function(n_sim = 1e+05, covmat, lower = -Inf, upper, fixed = (lower == upper), out = c(1), burn_in = 1000) {
  
  # Checking that n_sim is a number
  if(!is.numeric(n_sim)) stop("The number of simulations n_sim must be numeric!")
  # Checking that n_sim is strictly positive
  if(n_sim <=0)stop("n_sim must be a positive number!")
  # Checking that covmat is symmetric
  if(!isSymmetric.matrix(covmat)) stop("The covariance matrix covmat must be symmetric!")
  # and numeric
  if(!is.numeric(covmat)) stop("The covariance matrix covmat must be numeric!")
  # Checking that the lower and upper cutoff points are valid
  if(!is.numeric(lower)) stop("The lower cutoff point(s) must be numeric!")
  if(!is.numeric(upper)) stop("The upper cutoff point(s) must be numeric!")
  if(length(lower)!= length(upper)) stop("The lower and upper cutoff point(s) differ in length!")
  if(length(lower)!= 1 && length(lower)!= ncol(covmat)) stop("The length of the lower and upper cutoff point(s) must be 1 or equal to the dimension of the mutlivariable normal distribution!")
  if(any(upper < lower)){
    cat("Warning message: \n Some lower cutoff points are larger than the corresponding upper cutoff points! \n
        The lower and upper cutoff points will be swapped...")
    
    swapping_indx <- which(upper < lower)
    
    lower[swapping_indx] <- lower[swapping_indx] + upper[swapping_indx]
    upper[swapping_indx] <- lower[swapping_indx] - upper[swapping_indx]
    lower[swapping_indx] <- lower[swapping_indx] - upper[swapping_indx]
  }
  
  # Checking that fixed is valid
  fixed <- as.logical(fixed)
  if(length(fixed) != 1 && length(fixed) != length(lower)) stop("fixed must be of length 1 or length(lower)!")
  # Checking whether out is valid
  if(is.numeric(out)){
    
    out <- intersect(out, c(1:ncol(covmat)))
    
  }else{
    stop("out must be a numeric vector!")
  }
  
  # Checking whether out is empty
  if(length(out) == 0){
    
    cat("Warning message: \n out is not of the required format! \n The function will return the first estimated variable!")
    out <- c(1)
  }
  # Sorting out
  out <- sort(out)
  # Checking that burn_in is valid
  if(!is.numeric(burn_in)) stop("burn_in must be numeric!")
  if(burn_in < 0) stop("burn_in must be non-negative!")
  if(burn_in == Inf) stop("burn_in must be finite!")
  
  
  # Start with medians of univariate distributions
  sd0 <- sqrt(diag(covmat))
  sum_p0 <- stats::pnorm(lower, sd = sd0) + stats::pnorm(upper, sd = sd0)
  x0 <- stats::qnorm(sum_p0 / 2, sd = sd0) 

  # Pre-computations
  d <- nrow(covmat)
  sd <- rep(NA, d)
  P <- matrix(NA, d, d)
  for (j in 1:d) {
    P_j <- rep(0, d)
    P_j[-j] <- solve(covmat[-j, -j] , covmat[-j, j])
    P[, j] <- P_j
    sd[j] <- drop(sqrt(covmat[j, j] - crossprod(P_j, covmat[, j])))
  }
  
  # Vector entries to return
  to_return <- rep(-1, d)
  to_return[out] <- seq_along(out) - 1L
  
  # Actual Gibbs sampler
  rtmvnorm_gibbs_cpp(P, sd, lower, upper, fixed, to_return, x0, n_sim, burn_in)
}
