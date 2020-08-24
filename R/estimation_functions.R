
estimate_gen_liability_multi_trait = function(phen.list , thr.list, nthreads = 20) {
  
    n_trait = length(phen.list)
    phen = all_phenotypes[[1]]
    
    phen$post_gen_liab <- NA
    phen$post_gen_liab_se <- NA

    
    ph = foreach(i = 1:nrow(phen), .packages = c("batchmeans", "LTFHPlus"),
                 .options.snow = opts,
                 .inorder = T) %dopar% { 
                   fam = c(phen$FID[i], phen$pid_f[i], phen$pid_m[i], phen$sib_ids[[i]])
                   n_sib = length(fam) - 3
                   full_cov = get_full_cov(corr_mat = corr_mat, n_sib = n_sib)
                   full_cov = check_positive_definite(full_cov = full_cov, corr_mat = corr_mat,
                                                      correction_val = 0.99, n_sib = n_sib)
                   cov_size = nrow(full_cov)
                   lower = rep(-Inf, cov_size)
                   upper = rep(Inf, cov_size) 
                   
                   for (k in 1:n_trait) {
                     cur_phen = all_phenotypes[[k]]
                     cur_thr  = all_thr[[k]]
                     
                     for (ii in 1:length(fam) + 1) {
                       cur_fam = cur_thr[cur_thr$pid == fam[ii - 1], ]
                       if (is.na(cur_fam$status)) {
                         #here to deal with NAs  
                       } else if (cur_fam$status == 1) {
                         lower[(4 + n_sib) * (k - 1) + ii] <- upper[(4 + n_sib) * (k - 1) + ii] <- curfam$thr
                       } else {
                         upper[(4 + n_sib) * (k - 1) + ii] <- cur_fam$thr
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
                                                              sigma = cov,
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
    
    means = sapply(ph, FUN = function(x) x$est)
    ses   = sapply(ph, FUN = function(x) x$se)
    
    return(list("post_gen_liab"    = means,
                "post_gen_liab_se" = ses  ))
  
}
