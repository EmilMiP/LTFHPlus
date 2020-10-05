library(LTFHPlus)
library(dplyr)
library(stringr)
library(ggplot2) 
library(gridExtra)
library(doSNOW)
library(progress)


N = 5000
h2_1 = .5
h2_2 = .5
gen_cor = .5
corr_mat = diag(c(h2_1, h2_2))
corr_mat[1,2] <- corr_mat[2,1] <- gen_cor
nthreads = 5  # number of threads to use for ltfh++
nsib = 3
tol = 0.01

#calculates the thresholds used to determine status:
multiplier = 1
prev = c(0.08, .02) * multiplier
#(thr = qnorm(1 - prev))


#### THE NEXT SECTION REQUIRES YOU TO HAVE THE SOURCE CODE FOR LT-FH LOADED OR SOURCING IT ####
source("C:/Code/LTFH/assign_ltfh.R")
## download from here: https://alkesgroup.broadinstitute.org/UKBB/LTFH/


est_cir = function(data, 
                   indivs = c("child", "father", "mother"), 
                   ids = c("FID", "pid_f", "pid_m")) {
  
  res = tibble()
  for(i in seq_along(indivs)) {
    indiv = indivs[i]
    stat_col = paste(indiv, "_stat", sep = "")
    age_col  = paste(indiv, "_age", sep = "")
    
    if (indivs[i] %in% c("father", "mother")) {
      ph = data %>% 
        arrange(!!as.name(age_col)) %>%
        mutate(cir = (cumsum(!!as.name(stat_col)) + 1)/n()) %>%
        mutate(thr = qnorm(cir, lower.tail = FALSE)) %>% 
        select(!!as.name(ids[i]), thr)
    } else {
      sex_col = paste(indiv, "_sex", sep = "")
      ph = data %>% 
        group_by(!!as.name(sex_col)) %>%
        arrange(!!as.name(age_col)) %>%
        mutate(cir = (cumsum(!!as.name(stat_col)) + 1)/n()) %>%
        mutate(thr = qnorm(cir, lower.tail = FALSE)) %>% 
        ungroup() %>%
        select(!!as.name(ids[i]), thr)
    }
    colnames(ph) = c("FID", "thr")
    res = bind_rows(res, ph)
  }
  return(res)
}

#simulate liabilities
liabs = MASS::mvrnorm(n = N, mu = rep(0, 2*(4 + nsib)), Sigma = get_full_cov(corr_mat = corr_mat, n_sib = nsib))

child_age = runif(N, 10, 60)
father_age  = child_age + runif(N, 20, 35)
mother_age  = child_age + runif(N, 20, 35)

all_phen = list()
for (i in 1:2) {
  simu_liab = tibble(
    FID         = as.character(1:N),
    IID         = 1:N,
    child_gen   = liabs[,(i - 1) * 4 + 1],
    child_full  = liabs[,(i - 1) * 4 + 2],
    father_full = liabs[,(i - 1) * 4 + 3],
    mother_full = liabs[,(i - 1) * 4 + 4],
    child_sex   = sample(1:2, size = N, replace = TRUE),
    child_stat  = (child_full  > qnorm(prev[child_sex], lower.tail = F)) + 0L,
    father_stat = (father_full > qnorm(prev[1], lower.tail = F)) + 0L,
    mother_stat = (mother_full > qnorm(prev[2], lower.tail = F)) + 0L,
    child_age   = child_age,
    father_age  = father_age,
    mother_age  = mother_age,
    pid_f       = paste(FID, "_f", sep = ""),
    pid_m       = paste(FID, "_m", sep = "")
  )
  
  if (nsib > 0) {
    for (ii in 1:nsib) {
      sib_stat = paste("sib", ii, "_stat", sep = "")
      sib_sex  = paste("sib", ii, "_sex", sep = "")
      sib_age  = paste("sib", ii, "_age", sep = "")
      sib_full = paste("sib", ii, "_full", sep = "")
      sib_id   = paste("pid_s", ii, sep = "")
      simu_liab[[sib_full]] = liabs[, 4 + ii]
      simu_liab[[sib_sex]]  = sample(1:2, size = N, replace = TRUE)
      simu_liab[[sib_age]]  = sample(-5:5, size = N, replace = TRUE) + simu_liab$child_age
      simu_liab[[sib_stat]] = (simu_liab[[sib_full]] > qnorm(prev[simu_liab[[sib_sex]]], lower.tail = F)) + 0L
      simu_liab[[sib_id]]   = paste(simu_liab$FID, "_s", ii, sep = "")
    }
  }
  all_phen[[i]] = simu_liab
}

all_thr = list() 
for (ii in 1:2) { #The two tables are identical here, but in real data, we would see differences depending on age of onset, cohort effects etc.
  simu_liab = all_phen[[ii]]
  indivs = c("child", "father", "mother", if (nsib > 0) paste("sib", 1:nsib, sep = ""))
  ids = c("FID", "pid_f", "pid_m", if(nsib > 0) paste("pid_s", 1:nsib, sep = ""))
  
  all_thr[[ii]] = est_cir(data = simu_liab,
                          indivs = indivs,
                          ids = ids)
}


# data = estimate_gen_liability_multi_trait(corr_mat = corr_mat,
#                                           phen.list = all_phen,
#                                           thr.list = all_thr,
#                                           ids = ids,
#                                           ind = c(1,4 + nsib), 
#                                           status_cols = paste(indivs, "_stat", sep = ""),
#                                           nthreads = nthreads)
data = all_phen[[1]]
ind = c(1,4 + nsib)
status_cols = paste(indivs, "_stat", sep = "")
n_trait = length(all_phen)
cat("starting parallelization backend with", nthreads, "threads for generation of children:\n")
cl = makeCluster(nthreads, type = "SOCK")
registerDoSNOW(cl)
iterations = nrow(simu_liab)

pb = progress_bar$new(
  format = "[:bar] :percent",
  total = iterations,
  width = 100)

progress_num = 1:iterations
progress = function(n){
  pb$tick(tokens = list(letter = progress_num[n]))
}

opts = list(progress = progress)

ph = foreach(i = 1:nrow(data),
             .options.snow = opts,
             .export = c("get_full_cov", "check_positive_definite","rtmvnorm.gibbs", "rtmvnorm_gibbs_cpp"),
             .inorder = T) %dopar% { 
               fam = unlist(data[i,ids])
               full_cov = get_full_cov(corr_mat = corr_mat, n_sib = nsib)
               full_cov = check_positive_definite(full_cov = full_cov, corr_mat = corr_mat,
                                                  correction_val = 0.99, n_sib = nsib)
               cov_size = nrow(full_cov)
               lower = rep(-Inf, cov_size)
               upper = rep(Inf, cov_size) 
               
               for (k in 1:n_trait) {
                 cur_phen = all_phen[[k]]
                 cur_thr  = all_thr[[k]]
                 full_fam = cur_phen[i,]
                 cur_status = unlist(cur_phen[i, status_cols])
                 cur_liab   = unlist(cur_phen[i, paste(indivs, "_full", sep = "")])
                 for (ii in 1:length(fam) + 1) {
                   indiv_thr = cur_thr[cur_thr[[1]] == fam[ii - 1], ]
                   if (is.na(cur_status[ii - 1])) {
                     #here to deal with NAs  
                   } else if (cur_status[ii - 1] == 1) {
                     lower[(4 + nsib) * (k - 1) + ii] <- cur_liab[[ii - 1]]
                   } else {
                     upper[(4 + nsib) * (k - 1) + ii] <- indiv_thr$thr
                   }
                   
                 }
               }
               fixed <- (upper - lower) < 1e-4
               
               #covergence check
               se = NULL 
               vals = list() #store simulated values
               vals.ctr = 1
               while (is.null(se) || se > tol) {
                 gen_liabs = rtmvnorm.gibbs(5e4, burn_in = 1000,
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
             }
stopCluster(cl)
if (length(ind) > 1) {
  tmp <- t(sapply(ph, FUN = function(x) x[1,]))
  for (ii in 1:n_trait) {
    data[[paste("post_gen_liab_", ii, sep = "")]] = unlist(tmp[,ii])
  }
  tmp <- t(sapply(ph, FUN = function(x) x[2,]))
  for (ii in 1:n_trait) {
    data[[paste("post_gen_liab_", ii,"_se", sep = "")]] = unlist(tmp[,ii])
  }
  
} else {
  data$post_gen_liab      = unlist(sapply(ph, FUN = function(x) x[1,]))
  data$post_gen_liab_se   = unlist(sapply(ph, FUN = function(x) x[2,]))
}









res = as_tibble(as.data.frame(matrix(NA, nrow = nrow(data), ncol = 7)))
res[,1] = as.double(data$FID)
res[,2] = as.double(data$IID)
colnames(res) = c("FID", "IID", "CHILD_STATUS", "P1_STATUS", "P2_STATUS", "NUM_SIBS", "SIB_STATUS")
res$CHILD_STATUS = data$child_stat
res$P1_STATUS = data$father_stat
res$P2_STATUS = data$mother_stat
res$NUM_SIBS = nsib
res$SIB_STATUS = 0
if (nsib > 0) res$SIB_STATUS = (rowSums(simu_liab[,paste("sib", 1:nsib, "_stat", sep = "")]) > 0) + 0L

ltfh = create_pheno(data = as.data.frame(res),
                    trait_h2 = h2_1,
                    T_val_child = qnorm(mean(prev), lower.tail = F),
                    T_val_parent = qnorm(mean(prev), lower.tail = F),
                    relevant_trait_child = "CHILD_STATUS",
                    relevant_trait_dad = "P1_STATUS",
                    relevant_trait_mom = "P2_STATUS",
                    number_siblings_col = "NUM_SIBS",
                    relevant_trait_sib = "SIB_STATUS",
                    maximum_siblings_to_compute = nsib)

data = left_join(data, as.data.frame(ltfh))


with(data, c( "primary"      = cov(child_gen, post_gen_liab_1), 
              "secondary"    = cov(child_gen, post_gen_liab_2), 
              "LT-FH"        = cov(child_gen, ltfh), 
              "LT-FH_nosib"  = cor(child_gen, ltfh_nosib), 
              "Case-Control" = cov(child_gen, child_stat)))
with(data, c( "primary"      = cor(child_gen, post_gen_liab_1), 
              "secondary"    = cor(child_gen, post_gen_liab_2), 
              "LT-FH"        = cor(child_gen, ltfh), 
              "LT-FH_nosib"  = cor(child_gen, ltfh_nosib), 
              "Case-Control" = cor(child_gen, child_stat)))




p1 = ggplot(data, aes(x = post_gen_liab_1, y = child_gen, color = rowSums(all_phen[[2]][,c("child_stat", "father_stat", "mother_stat")]) > 0)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Any case of \n2nd Phenotype \nin Family") +
  xlab("Estimated Genetic Liability (LTFH++) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


p2 = ggplot(data, aes(x = ltfh, y = child_gen, color = rowSums(all_phen[[2]][,c("child_stat", "father_stat", "mother_stat")]) > 0)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Any case of \n2nd Phenotype \nin Family") +
  xlab("Estimated Genetic Liability (LTFH) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2)

