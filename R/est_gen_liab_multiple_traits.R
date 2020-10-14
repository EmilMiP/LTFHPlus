utils::globalVariables("i")

#'
#' Estimates genetic liability based off of multiple phenotypes, family history, and temportal information.
#'
#' @param corr_mat Correlation matrix for each of the traits being analysed. 
#' @param phen.list list of tibble or data.frame with the IDs and status of the genotyped individual and the parents and siblings.
#' @param thr.list list of tibble or data.frame with a row of each individual in the provided phen. Each row should contain the ID and threshold value needed for the model.
#' @param status_cols Vector with the names of the columns that has the status of each family. default is c("child_stat", "father_stat", "mother_stat").
#' @param ids Column names of IDs for family members. 
#' @param ind Indices to return from the gibbs sampler. 1 corresponds to the genetic liability of the first phenotype. c(1,5) corresponds to the genetic liability of the first two phenotypes provided, but with no siblings in the model.
#' @param nthreads number of threads to use in estimating the genetic liabilities. Do not exceed the number of threads your CPU has available. 
#' @param tol Convergence criteria of the gibbs sampler. Default is 0.01, meaning a standard error of the mean below 0.01
#'
#' @return Returns the estimated genetic liabilities.
#'
#' @examples 
#' # See R/Example/example_multi.R for an example of use and input.
#'
#' @export
#' 
estimate_gen_liability_multi_trait = function(phen.list,
                                              thr.list, 
                                              corr_mat,
                                              status_cols = c("child_stat", "father_stat", "mother_stat"),
                                              ids = c("FID", "pid_f", "pid_m"),
                                              ind = c(1,5),
                                              nthreads = 10,
                                              tol = 0.01) {
  
  n_trait = length(phen.list)
  #get ids from the first phenotype tibble
  phen = phen.list[[1]]
  
  iterations = nrow(phen)
  
  cat("starting parallelization backend with", nthreads, "threads for generation of children:\n")
  cl =   parallel::makeCluster(nthreads, type = "SOCK")
  doParallel::registerDoParallel(cl)
  
  ph = foreach(i = 1:nrow(phen),
               .export = c("get_full_cov", "check_positive_definite","rtmvnorm.gibbs", "rtmvnorm_gibbs_cpp"),
               .inorder = T) %dopar% { 
                 fam = unlist(phen[i,ids])
                 n_sib = length(fam) - 3
                 full_cov = get_full_cov(corr_mat = corr_mat, n_sib = n_sib)
                 full_cov = check_positive_definite(full_cov = full_cov, corr_mat = corr_mat,
                                                    correction_val = 0.99, n_sib = n_sib)
                 cov_size = nrow(full_cov)
                 lower = rep(-Inf, cov_size)
                 upper = rep(Inf, cov_size) 
                 
                 for (k in 1:n_trait) {
                   cur_phen = phen.list[[k]]
                   cur_thr  = thr.list[[k]]
                   full_fam = cur_phen[i,]
                   cur_status = unlist(cur_phen[i, status_cols])
                   for (ii in 1:length(fam) + 1) {
                     indiv_thr = cur_thr[cur_thr[[1]] == fam[ii - 1], ]
                     if (is.na(cur_status[ii - 1])) {
                       #here to deal with NAs  
                     } else if (cur_status[ii - 1] == 1) {
                       lower[(4 + n_sib) * (k - 1) + ii] <- upper[(4 + n_sib) * (k - 1) + ii] <- indiv_thr$thr
                     } else {
                       upper[(4 + n_sib) * (k - 1) + ii] <- indiv_thr$thr
                     }
                     
                   }
                 }
                 fixed <- (upper - lower) < 1e-4
                 
                 #covergence check
                 se = NULL 
                 vals = list() #store simulated values
                 vals.ctr = 1
                 while (is.null(se) || se > tol) {
                   gen_liabs = rtmvnorm.gibbs(1e5, burn_in = 1000,
                                              sigma = full_cov,
                                              lower = lower, 
                                              upper = upper,
                                              ind = ind,
                                              fixed = fixed)
                   vals[[vals.ctr]] = gen_liabs
                   se = batchmeans::bm(unlist(vals))$se
                   vals.ctr =  vals.ctr + 1
                 }
                 #calculate the final values
                 vals = do.call("rbind", vals)
                 sapply(1:length(ind), FUN = function(n) {
                   batchmeans::bm(vals[,n])
                 })
                 #batchmeans::bm(unlist(vals))
               }
  parallel::stopCluster(cl)
  if (length(ind) > 1) {
    tmp <- t(sapply(ph, FUN = function(x) x[1,]))
    for (ii in 1:n_trait) {
      phen[[paste("post_gen_liab_", ii, sep = "")]] = unlist(tmp[,ii])
    }
    tmp <- t(sapply(ph, FUN = function(x) x[2,]))
    for (ii in 1:n_trait) {
      phen[[paste("post_gen_liab_", ii,"_se", sep = "")]] = unlist(tmp[,ii])
    }
    
  } else {
    phen$post_gen_liab      = unlist(sapply(ph, FUN = function(x) x[1,]))
    phen$post_gen_liab_se   = unlist(sapply(ph, FUN = function(x) x[2,]))
  }

  
  return(phen)
}
