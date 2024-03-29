utils::globalVariables("i")
utils::globalVariables("string")

#'
#' Estimate genetic liability based off of family history and prevalence information.
#' 
#' This function has been replaced by \code{\link{estimate_liability}}. This function is kept for legacy reasons.
#'
#' @param h2 Liability scale heritability of the trait being analysed.
#' @param phen tibble or data.frame with the IDs and status of the genotyped individual and the parents and siblings.
#' @param thr tibble or data.frame with a row of each individual in the provided phen. Each row should contain the ID and threshold value needed for the model.
#' @param ind Indices to return from the gibbs sampler. 1 corresponds to the genetic liability of the first phenotype. c(1,5) corresponds to the genetic liability of the first two phenotypes provided, but with no siblings in the model.
#' @param id_col Column name for IDs of family members.
#' @param tol Convergence criteria of the gibbs sampler. Default is 0.01, meaning a standard error of the mean below 0.01
#'
#' @return Returns the estimated genetic liabilities.
#' @noRd

estimate_gen_liability = function(h2,
                                  phen, 
                                  thr, 
                                  ind = c(1),
                                  id_col = "ids",
                                  tol = 0.01) {
  
  warning("'estimate_gen_liability()' was deprecated in LTFHPlus v1.0.0. It is only kept for legacy reasons.\n
          Please consider using 'estimate_liability()' instead.\
          The signature and semantics have changed, see '?estimate_liability'.")
  
  phen$post_gen_liab <- NA
  phen$post_gen_liab_se <- NA
  
  
  ph = future.apply::future_lapply(X = 1:nrow(phen), FUN = function(i){
    
    full_fam = phen[[id_col]][[i]]
    n_sib = length(full_fam) - 3
    
    cov = suppressWarnings(get_cov(h2 = h2, n_sib = n_sib))
    
    cur_indiv = thr[match(full_fam, thr[[1]]), ]
    lower = c(-Inf, cur_indiv$lower)
    upper = c(Inf , cur_indiv$upper)
    fixed <- (upper - lower) < 1e-4
    
    #covergence check
    se = NULL 
    vals = list() #store simulated values
    vals.ctr = 1
    while (is.null(se) || se > tol) {
      gen_liabs = rtmvnorm.gibbs(1e5,
                                 burn_in = 1000,
                                 covmat  = cov,
                                 lower   = lower, 
                                 upper   = upper,
                                 out     = ind,
                                 fixed   = fixed)
      vals[[vals.ctr]] = gen_liabs
      se = batchmeans::bm(unlist(vals))$se
      vals.ctr =  vals.ctr + 1
    }

    batchmeans::bm(unlist(vals))
  }, future.seed = TRUE)
  
  phen$post_gen_liab      = sapply(ph, FUN = function(x) x$est)
  phen$post_gen_liab_se   = sapply(ph, FUN = function(x) x$se)
  
  return(phen)
}