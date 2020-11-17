library(LTFHPlus)
library(dplyr)
library(stringr)
library(ggplot2) 
library(gridExtra)
library(tidyverse)
library(doSNOW)
library(progressr)

N = 5000 
h2 = .5 
nsib = 3
nthreads = 5  # number of threads to use for ltfh++
tol = 0.01

#calculates the thresholds used to determine status:
multiplier = 1
prev = c(0.05, 0.05) * multiplier

#### THE NEXT SECTION REQUIRES YOU TO HAVE THE SOURCE CODE FOR LT-FH LOADED OR SOURCING IT ####
#source("C:/Code/LTFH/assign_ltfh.R")
## download from here: https://alkesgroup.broadinstitute.org/UKBB/LTFH/

#plan(multisession)
handlers("progress")


#simulate liabilities
liabs = MASS::mvrnorm(n = N, mu = rep(0, 4 + nsib), Sigma = LTFHPlus:::get_cov(h2, n_sib = nsib))

simu_liab = tibble(
  FID         = as.character(1:N),
  IID         = 1:N,
  child_gen   = liabs[,1],
  child_full  = liabs[,2],
  father_full = liabs[,3],
  mother_full = liabs[,4],
  child_sex   = sample(1:2, size = N, replace = TRUE),
  child_stat  = (child_full  > qnorm(prev[child_sex], lower.tail = F)) + 0L,
  father_stat = (father_full > qnorm(prev[1], lower.tail = F)) + 0L,
  father_sex  = rep(1, N),
  mother_stat = (mother_full > qnorm(prev[2], lower.tail = F)) + 0L,
  mother_sex  = rep(2, N),
  child_age   = runif(N, 10, 60),
  father_age  = child_age + runif(N, 20, 35),
  mother_age  = child_age + runif(N, 20, 35),
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



# LTFH++_fast -------------------------------------------------------------
## more elegant implementaion of estimate_gen_liability_ltfh is pending ##
simu_liab$NUM_SIBS = nsib
simu_liab$SIB_STATUS = 0
if (nsib > 0) simu_liab$SIB_STATUS = (rowSums(simu_liab[,paste("sib", 1:nsib, "_stat", sep = "")]) > 0) + 0L  

with_progress({data = estimate_gen_liability_ltfh(h2 = h2,
                                   phen = simu_liab,
                                   child_threshold = qnorm(mean(prev), lower.tail = F),
                                   parent_threshold = qnorm(mean(prev), lower.tail = F),
                                   status_col_offspring = "child_stat",
                                   status_col_father    = "father_stat",
                                   status_col_mother    = "mother_stat",
                                   status_col_siblings  = "SIB_STATUS",
                                   number_of_siblings_col = "NUM_SIBS")})

colnames(data)[ncol(data) - 1:0 ] = paste(colnames(data)[ncol(data) - 1:0 ], "_fast", sep = "")



#  LT-FH++  ---------------------------------------------------------------

cov = get_cov(h2 = h2, n_sib = nsib)

indivs = c("child", "father", "mother", if(nsib > 0) paste("sib", 1:nsib, sep = ""))
ids = c("FID", "pid_f", "pid_m", if(nsib > 0) paste("pid_s", 1:nsib, sep = ""))
status_cols = paste(indivs, "_stat", sep = "")
age_cols = paste(indivs, "_age", sep = "")
sex_cols = paste(indivs, "_sex", sep = "")
liab_cols = paste(indivs, "_full", sep = "")

thr = tibble::tibble(ids = character(0), lower = numeric(0), upper = numeric(0))
for (jj in seq_along(status_cols)) {
  lower = rep(-Inf, nrow(simu_liab))
  upper = rep(Inf, nrow(simu_liab))
  cases = simu_liab[[status_cols[jj]]] == 1 
  lower[cases] <- upper[cases] <- age_to_thres(age = liab_to_aoo(simu_liab[[liab_cols[jj]]][cases], pop_prev = prev[simu_liab[[sex_cols[jj]]]][cases]), 
                                                                           pop_prev = prev[simu_liab[[sex_cols[jj]]]][cases])
  upper[!cases] <- age_to_thres(age = simu_liab[[age_cols[jj]]][!cases], pop_prev = prev[simu_liab[[sex_cols[jj]]]][!cases])
  thr = bind_rows(thr, as_tibble(list("ids" = simu_liab[[ids[jj]]],"lower" = lower, "upper" = upper)))
}  

estimate_gen_liability(h2 = h2, phen = simu_liab,thr = thr, status_cols = status_cols, ids = ids,nthreads = nthreads)

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
ph = foreach(i = 1:nrow(simu_liab),
             .options.snow = opts) %dopar% {
               full_fam = simu_liab[i,]
               fam = unlist(full_fam[,ids])
               n_sib = length(fam) - 3
               lia_cols = names(full_fam)[grep("full",names(full_fam))]
               cov = LTFHPlus::get_cov(h2 = h2, n_sib = n_sib)
               
               cov_size = nrow(cov)
               
               lower = rep(-Inf, cov_size)
               upper = rep(Inf, cov_size) 
               cur_status = unlist(simu_liab[i, status_cols])
               for (ii in 1:(3 + n_sib)) {
                 #cur_indiv = thr[thr[[1]] == fam[ii], ]
                 
                 new_a = qnorm(prev[full_fam[sex_cols[[ii]]][[1]]], lower.tail = F)
                 if (is.na(cur_status[ii])) {
                   #here to deal with NAs for now 
                 } else if (cur_status[ii] == 1) {
                   lower[ii + 1] <- upper[[ii + 1]] <- LTFHPlus::age_to_thres(age = LTFHPlus::liab_to_aoo(full_fam[liab_cols[[ii]]][[1]], 
                                                                                                          pop_prev = prev[full_fam[sex_cols[[ii]]][[1]]]), 
                                                                              pop_prev = prev[full_fam[sex_cols[[ii]]][[1]]])
                 } else {
                   upper[ii + 1] <- LTFHPlus::age_to_thres(age = full_fam[age_cols[ii]][[1]], 
                                                           pop_prev = prev[full_fam[sex_cols[ii]][[1]]])
                 }
                 
               }
               fixed <- (upper - lower) < 1e-4
               
               #covergence check
               se = NULL 
               vals = list() #store simulated values
               vals.ctr = 1
               while (is.null(se) || se > tol) {
                 gen_liabs = LTFHPlus::rtmvnorm.gibbs(1e4,
                                                      burn_in = 1000,
                                                      sigma   = cov,
                                                      lower   = lower, 
                                                      upper   = upper,
                                                      fixed   = fixed)
                 vals[[vals.ctr]] = gen_liabs
                 se = batchmeans::bm(unlist(vals))$se
                 vals.ctr =  vals.ctr + 1
               }
               #calculate the final values
               batchmeans::bm(unlist(vals))
             }
stopCluster(cl)
simu_liab[["post_gen_liab"]]    <- sapply(ph, FUN = function(x) x$est)
simu_liab[["post_gen_liab_se"]] <- sapply(ph, FUN = function(x) x$se)


data = left_join(data, simu_liab)



#  LT-FH ------------------------------------------------------------------


res = as_tibble(as.data.frame(matrix(NA, nrow = nrow(simu_liab), ncol = 7)))
res[,1] = as.double(simu_liab$FID)
res[,2] = as.double(simu_liab$IID)
colnames(res) = c("FID", "IID", "CHILD_STATUS", "P1_STATUS", "P2_STATUS", "NUM_SIBS", "SIB_STATUS")
res$CHILD_STATUS = simu_liab$child_stat
res$P1_STATUS = simu_liab$father_stat
res$P2_STATUS = simu_liab$mother_stat
res$NUM_SIBS = nsib
res$SIB_STATUS = 0
if (nsib > 0) res$SIB_STATUS = (rowSums(simu_liab[,paste("sib", 1:nsib, "_stat", sep = "")]) > 0) + 0L

ltfh = create_pheno(data = as.data.frame(res),
                    trait_h2 = h2,
                    T_val_child = qnorm(mean(prev), lower.tail = F),
                    T_val_parent = qnorm(mean(prev), lower.tail = F),
                    relevant_trait_child = "CHILD_STATUS",
                    relevant_trait_dad = "P1_STATUS",
                    relevant_trait_mom = "P2_STATUS",
                    number_siblings_col = "NUM_SIBS",
                    relevant_trait_sib = "SIB_STATUS",
                    maximum_siblings_to_compute = nsib)

simu_liab = left_join(data, as.data.frame(ltfh))



# Plots & Results ---------------------------------------------------------



#Morale of this example: The estimatation methods all perform similarly, but with a threshold for each group, estimate_gen_liability_ltfh, is the fastest.

with(simu_liab, c("LTFH++" = cov(child_gen, post_gen_liab),"LTFH_gibbs" = cov(child_gen, post_gen_liab_fast), "LTFH" = cov(child_gen, ltfh)))
with(simu_liab, c("LTFH++" = cor(child_gen, post_gen_liab),"LTFH_gibbs" = cor(child_gen, post_gen_liab_fast), "LTFH" = cor(child_gen, ltfh)))

#comparing LTFH to LTFH++_fast
simu_liab = simu_liab %>% mutate(difference = ltfh - post_gen_liab_fast)

comp1 = ggplot(simu_liab, aes(x = post_gen_liab_fast, y = ltfh, color = as.factor(child_stat))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

comp2 = ggplot(simu_liab %>% arrange(ltfh) %>% select(difference) %>% unique() %>% mutate(grp = 1:n()), aes(y = difference, x = as_factor(grp))) +
  geom_bar(stat = "identity")

grid.arrange(comp1, comp2)



#plotting estimating genetic liability for the different methods against the true genetic liabilities:
p1 = ggplot(simu_liab, aes(x = post_gen_liab, y = child_gen, color = rowSums(data[,c("child_stat", "father_stat", "mother_stat", if(nsib > 0) paste("sib", 1:nsib, "_stat", sep = ""))]) > 0)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status in Family") +
  xlab("Estimated Genetic Liability (LTFH++) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")


p2 = ggplot(simu_liab, aes(x = post_gen_liab_fast, y = child_gen, color = rowSums(data[,c("child_stat", "father_stat", "mother_stat", if(nsib > 0) paste("sib", 1:nsib, "_stat", sep = ""))]) > 0)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status in Family") +
  xlab("Estimated Genetic Liability (LTFH_gibbs) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

p3 = ggplot(simu_liab, aes(x = ltfh, y = child_gen, color = rowSums(data[,c("child_stat", "father_stat", "mother_stat", if(nsib > 0) paste("sib", 1:nsib, "_stat", sep = ""))]) > 0)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status in Family") +
  xlab("Estimated Genetic Liability (LTFH) ") +
  ylab("True Genetic Liability") + 
#  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top")

grid.arrange(p1, p2, p3)

grid.arrange(p1, p3)


