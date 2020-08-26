
estimate_gen_liability_multi_trait = function(phen.list,
                                              thr.list, 
                                              nthreads = 10,
                                              corr_mat,
                                              ids = c("FID", "pid_f", "pid_m"),
                                              status_cols = c("child_stat", "father_stat", "mother_stat"),
                                              tol = 0.01) {
  
  n_trait = length(phen.list)
  #get ids from the first phenotype tibble
  phen = phen.list[[1]]
  
  
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
               .packages = "LTFHPlus",
               .inorder = T) %dopar% { 
                 #fam = c(phen$FID[i], phen$pid_f[i], phen$pid_m[i], phen$sib_ids[[i]])
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
                   for (ii in 1:length(fam) + 1) {
                     indiv_thr = cur_thr[cur_thr$ids == fam[ii - 1], ]
                     if (is.na(full_fam[[status_cols[ii - 1]]])) {
                       #here to deal with NAs  
                     } else if (full_fam[[status_cols[ii - 1]]] == 1) {
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
                   gen_liabs = rtmvnorm.gibbs(10e4, burn_in = 1000,
                                              sigma = full_cov,
                                              lower = lower, 
                                              upper = upper,
                                              fixed = fixed)
                   vals[[vals.ctr]] = gen_liabs
                   se = batchmeans::bm(unlist(vals))$se
                   vals.ctr =  vals.ctr + 1
                 }
                 #calculate the final values
                 batchmeans::bm(unlist(vals))
               }
  
  phen$post_gen_liab      = sapply(ph, FUN = function(x) x$est)
  phen$post_gen_liab_se   = sapply(ph, FUN = function(x) x$se)
  
  return(phen)
}
