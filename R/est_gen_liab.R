library(tidyverse)
library(doSNOW)
library(progress)

estimate_gen_liability = function(h2,
                                  phen, 
                                  thr, 
                                  ids = c("FID", "pid_f", "pid_m"),
                                  status_cols = c("child_stat", "father_stat", "mother_stat"),
                                  nthreads = 10,
                                  tol = 0.01) {
  phen$post_gen_liab <- NA
  phen$post_gen_liab_se <- NA
  
  
  iterations = nrow(phen)
  
  cat("starting parallelization backend with", nthreads, "threads for generation of children:\n")
  cl = makeCluster(nthreads, type = "SOCK")
  registerDoSNOW(cl)
  
  pb = progress_bar$new(
    format = "[:bar] :percent",
    total = iterations,
    width = 100)
  
  progress_num = 1:iterations
  progress = function(n){
    pb$tick(tokens = list(letter = progress_num[n]))
  }
  
  opts = list(progress = progress)
  
  ph = foreach(i = 1:nrow(phen),
                 .options.snow = opts,
                 .inorder = T) %dopar% {


                 full_fam = phen[i,]
                 fam = unlist(full_fam[,ids])
                 n_sib = length(fam) - 3
                 
                 cov = get_cov(h2 = h2, n_sib = n_sib)
                 
                 cov_size = nrow(cov)
                 
                 lower = rep(-Inf, cov_size)
                 upper = rep(Inf, cov_size) 
                 
                for (ii in 1:length(fam) + 1) {
                    cur_indiv = thr[thr$ids == fam[ii - 1], ]
                    
                    if (is.na(full_fam[[status_cols[ii - 1]]])) {
                      #here to deal with NAs for now 
                    } else if (full_fam[[status_cols[ii - 1]]] == 1) {
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
                   se = bm(unlist(vals))$se
                   vals.ctr =  vals.ctr + 1
                 }
                #calculate the final values
                batchmeans::bm(gen_liabs)
               }
  stopCluster(cl)
  phen$post_gen_liab      = sapply(ph, FUN = function(x) x$est)
  phen$post_gen_liab_se   = sapply(ph, FUN = function(x) x$se)

  return(phen)
}
