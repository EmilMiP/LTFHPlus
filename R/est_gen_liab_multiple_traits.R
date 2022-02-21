utils::globalVariables("i")

#'
#' Estimates genetic liability based off of multiple phenotypes, family history, and prevalence information for each trait.
#'
#' @param corr.mat Correlation matrix for each of the traits being analysed. 
#' @param phen.list List of ids for family members. The order is assumed to be OS, P1, P2, and any siblings.
#' @param thr.list list of tibble or data.frame with a row of each individual in the provided phen. Each row should contain the ID and threshold value needed for the model.
#' @param id_col Column names of IDs for family members. 
#' @param ind Indices to return from the Gibbs sampler. 1 corresponds to the genetic liability of the first phenotype. c(1,5) corresponds to the genetic liability of the first two phenotypes provided, but with no siblings in the model.
#' @param tol Convergence criteria of the Gibbs sampler. Default is 0.01, meaning a standard error of the mean below 0.01
#'
#' @return Returns the estimated genetic liabilities.
#'
#' @examples 
#' # See See vignettes for example of use and input.
#'
#' @export
#' 

estimate_gen_liability_multi_trait = function(corr.mat,
                                              phen.list,
                                              thr.list, 
                                              ind = c(1),
                                              id_col = "ids",
                                              tol = 0.01) {
  
  ntrait = nrow(corr.mat)
  #get ids from the first phenotype tibble
  res = dplyr::tibble(IID = sapply(phen.list[[1]][[id_col]], function(x) x[1]),
                      post_gen_liab = NA,
                      post_gen_liab_se = NA)
  
  
  # p <- progressr::progressor(along = 1:nrow(phen))
  
  
  ph = future.apply::future_lapply(X = 1:nrow(res), FUN = function(i){
    thres_info = list()
    for(ii in 1:ntrait) {
    full_fam = phen.list[[ii]][[id_col]][[i]]
    thres_info[[ii]] = thr.list[[ii]][match(full_fam, thr.list[[ii]][[1]]), ] %>%
      dplyr::select(-1) %>%
      dplyr::bind_rows(dplyr::as_tibble(list(lower = -Inf, upper = Inf)), .)
   # fixed <- (upper - lower) < 1e-4
    }
    thres_info_comb = do.call("rbind", thres_info)
    
    if(any(is.na(thres_info_comb))) {
      thres_info_comb[which(is.na(thres_info_comb$lower)),1] = -Inf
      thres_info_comb[which(is.na(thres_info_comb$upper)),2] = Inf
    }
    
    n_sib = length(full_fam) - 3
    cov = construct_covmat(corr_mat = corr.mat, n_sib = n_sib)
    
    #covergence check
    se = NULL 
    vals = list() #store simulated values
    vals.ctr = 1
    while (is.null(se) || se > tol) {
      gen_liabs = rtmvnorm.gibbs(1e5,
                                 burn_in = 1000,
                                 sigma   = cov,
                                 lower   = thres_info_comb$lower, 
                                 upper   = thres_info_comb$upper,
                                 ind     = ind,
                                 fixed   = thres_info_comb$upper - thres_info_comb$lower < 1e-4)
      vals[[vals.ctr]] = gen_liabs
      se = batchmeans::bm(unlist(vals))$se
      vals.ctr =  vals.ctr + 1
    }
    #calculate the final values
    # p(sprintf("%g", i))
    
    batchmeans::bm(unlist(vals))
  }, future.seed = TRUE)
  
  if (length(ind) > 1) {
    tmp <- t(sapply(ph, FUN = function(x) x[1,]))
    for (ii in 1:ntrait) {
      res[[paste("post_gen_liab_", ii, sep = "")]] = unlist(tmp[,ii])
    }
    tmp <- t(sapply(ph, FUN = function(x) x[2,]))
    for (ii in 1:ntrait) {
      res[[paste("post_gen_liab_", ii,"_se", sep = "")]] = unlist(tmp[,ii])
    }
    
  } else {
    res$post_gen_liab      = unlist(sapply(ph, FUN = function(x) x[1]))
    res$post_gen_liab_se   = unlist(sapply(ph, FUN = function(x) x[2]))
  }
  return(res)
}

  

