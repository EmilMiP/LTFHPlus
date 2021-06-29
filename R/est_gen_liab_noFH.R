

#'
#' Estimate genetic liability similar to LT-FH
#'
#' @param h2 Liability scale heritability of the trait being analysed.
#' @param phen tibble or data.frame with status of the individual with columns ids, lower, and upper (thresholds).
#' @param tol Convergence criteria of the gibbs sampler. Default is 0.01, meaning a standard error of the mean below 0.01
#'
#' @return Returns the estimated genetic liabilities based on prevalence information.
#'
#' @examples 
#' # See R/Example/example_nosib.R for an example of use and input.
#' @importFrom dplyr %>%
#' @export
#' 


estimate_gen_liability_noFH = function(phen,
                                       h2,
                                       tol = 0.01) {
  

  phen$post_gen_no_fam <- NA
  phen$post_gen_no_fam_se <- NA
  
  cov_mat = matrix(h2, ncol = 2, nrow = 2)
  cov_mat[2,2] = 1
  
  ph = future.apply::future_lapply(X = 1:nrow(phen), FUN = function(i){ 

    cur_lower = c(-Inf, phen$lower[i])
    cur_upper = c( Inf, phen$upper[i])
    
    # estimate genetic liability & covergence check
    se = NULL 
    vals = list() #store simulated values
    vals.ctr = 1
    while (is.null(se) || se > tol) {
      gen_liabs = rtmvnorm.gibbs(50e3, 
                                 burn_in = 1000,
                                 sigma = cov_mat,
                                 lower = cur_lower, 
                                 upper = cur_upper,
                                 fixed = cur_upper - cur_lower < 1e-4)
      vals[[vals.ctr]] = gen_liabs
      se = batchmeans::bm(unlist(vals))$se
      vals.ctr =  vals.ctr + 1
    }
    #calculate the final values
    batchmeans::bm(unlist(vals))
    
  }, future.seed = TRUE)
  
  phen$post_gen_no_fam    = sapply(ph, function(x) x$est)
  phen$post_gen_no_fam_se = sapply(ph, function(x) x$se)
  
  return(phen)
}
