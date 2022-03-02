
#' Estimate genetic liability when multiple traits and prevalence information for each trait are available, but family history is not
#'
#' @param h2_vec Liability scale heritability of the trait being analysed.
#' @param gen_cor_vec vector of genetic correlations, in order of appearance in phenotype file
#' @param phen tibble or data.frame with status of the individual in order of h2_vec.
#' @param prev_vec Vector of prevalences in order of appearance in h2_vec.
#' @param tol Convergence criteria of the Gibbs sampler. Default is 0.01, meaning a standard error of the mean below 0.01
#'
#' @return Returns the estimated genetic liabilities.
#'
#' @examples 
#' # See R/Example/example_nosib.R for an example of use and input.
#' @importFrom dplyr %>%
#' @export
#' 


estimate_gen_liability_noFH_corr_traits = function(phen, 
                                                   h2_vec, 
                                                   gen_cor_vec,
                                                   prev_vec,
                                                   tol = .01) {
  ntraits = length(h2_vec)
  
  reduced = phen %>% # we assume first column is ID, rest are phenotype columns in order of h2_vec
    dplyr::select(-1) %>% 
    dplyr::group_by_all() %>% 
    dplyr::slice_head() %>%
    dplyr::ungroup()
  
  reduced$post_gen_no_fam <- NA
  reduced$post_gen_no_fam_se <- NA
  
  cov_mat = construct_covmat()
  
  for (i in 1:nrow(reduced)) {
    cur_config = unlist(reduced[i,1:ntraits])
    
    lower = rep(-Inf, 1 + ntraits)
    upper = rep(Inf,  1 + ntraits)
    
    lower[-1][!is.na(cur_config) & (cur_config == 1)] <- stats::qnorm(prev_vec[!is.na(cur_config) & (cur_config == 1)], lower.tail = FALSE)
    upper[-1][!is.na(cur_config) & (cur_config == 0)] <- stats::qnorm(prev_vec[!is.na(cur_config) & (cur_config == 0)], lower.tail = FALSE)
    
    fixed = rep(FALSE, 1 + ntraits)
    
    # estimate genetic liability & covergence check
    se = NULL 
    vals = list() #store simulated values
    vals.ctr = 1
    while (is.null(se) || se > tol) {
      gen_liabs = rtmvnorm.gibbs(n_sim = 50e3, 
                                 burn_in = 1000,
                                 covmat  = cov_mat,
                                 lower = lower, 
                                 upper = upper,
                                 fixed = fixed)
      vals[[vals.ctr]] = gen_liabs
      se = batchmeans::bm(unlist(vals))$se
      vals.ctr =  vals.ctr + 1
    }
    #calculate the final values
    est = batchmeans::bm(unlist(vals))
    reduced$post_gen_no_fam[i]    = est$est
    reduced$post_gen_no_fam_se[i] = est$se
    
  }
  return(dplyr::left_join(phen, reduced))
}
