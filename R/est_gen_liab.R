
#'
#' Estimate genetic liability based off of family history and temporal information.
#'
#' @param h2 Liability scale heritability of the trait being analysed.
#' @param phen tibble or data.frame with the IDs and status of the genotyped individual and the parents and siblings.
#' @param thr tibble or data.frame with a row of each individual in the provided phen. Each row should contain the ID and threshold value needed for the model.
#' @param status_cols Vector with the names of the columns that has the status of each family. default is c("child_stat", "father_stat", "mother_stat").
#' @param ids Column names of IDs for family members. 
#' @param ind Indices to return from the gibbs sampler. 1 corresponds to the genetic liability of the first phenotype. c(1,5) corresponds to the genetic liability of the first two phenotypes provided, but with no siblings in the model.
#' @param nthreads number of threads to use in estimating the genetic liabilities. Do not exceed the number of threads your CPU has available. 
#' @param tol Convergence criteria of the gibbs sampler. Default is 0.01, meaning a standard error of the mean below 0.1
#'
#' @return Returns the estimated genetic liabilities.
#'
#' @examples See the example.R for an example of use and input.
#'
#' @export

estimate_gen_liability = function(h2,
                                  phen, 
                                  thr, 
                                  ids = c("FID", "pid_f", "pid_m"),
                                  ind = c(1),
                                  status_cols = c("child_stat", "father_stat", "mother_stat"),
                                  nthreads = 10,
                                  tol = 0.01) {
  phen$post_gen_liab <- NA
  phen$post_gen_liab_se <- NA
  
  
  iterations = nrow(phen)
  
  cat("starting parallelization backend with", nthreads, "threads for generation of children:\n")
  cl = parallel::makeCluster(nthreads, type = "SOCK")
  doParallel::registerDoParallel(cl)
  
  ph = foreach::foreach(i = 1:nrow(phen)
                 .export = c("get_cov", "rtmvnorm.gibbs", "rtmvnorm_gibbs_cpp"),
                 .inorder = T) %dopar% {


                 full_fam = phen[i,]
                 fam = unlist(full_fam[,ids])
                 n_sib = length(fam) - 3
                 
                 cov = get_cov(h2 = h2, n_sib = n_sib)
                 
                 cov_size = nrow(cov)
                 
                 lower = rep(-Inf, cov_size)
                 upper = rep(Inf, cov_size) 
                 cur_status = unlist(phen[i, status_cols])
                for (ii in 1:length(fam) + 1) {
                    cur_indiv = thr[thr$ids == fam[ii - 1], ]
                    
                    if (is.na(cur_status[ii - 1])) {
                      #here to deal with NAs for now 
                    } else if (cur_status[ii - 1] == 1) {
                      lower[ii] <- upper[ii] <- cur_indiv$thr
                    } else {
                      upper[ii] <- cur_indiv$thr
                    }
                    
                  }
                 fixed <- (upper - lower) < 1e-4
                 
                 #covergence check
                se = NULL 
                vals = list() #store simulated values
                vals.ctr = 1
                while (is.null(se) || se > tol) {
                  gen_liabs = rtmvnorm.gibbs(10e4, burn_in = 1000,
                                             sigma = cov,
                                             lower = lower, 
                                             upper = upper,
                                             fixed = fixed)
                   vals[[vals.ctr]] = gen_liabs
                   se = batchmeans::bm(unlist(vals))$se
                   vals.ctr =  vals.ctr + 1
                 }
                #calculate the final values
                batchmeans::bm(gen_liabs)
               }
  parallel::stopCluster(cl)
  phen$post_gen_liab      = sapply(ph, FUN = function(x) x$est)
  phen$post_gen_liab_se   = sapply(ph, FUN = function(x) x$se)

  return(phen)
}
