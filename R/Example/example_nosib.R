library(LTFHPlus)
library(dplyr)
library(stringr)
library(ggplot2) 
library(gridExtra)
library(tidyverse)
library(doSNOW)
library(progress)

N = 5000 
h2 = .5 
nsib = 1
nthreads = 20  # number of threads to use for ltfh++

#calculates the thresholds used to determine status:
K = .05
multiplier = 1
prev = c(0.08, .02) * multiplier

#### THE NEXT SECTION REQUIRES YOU TO HAVE THE SOURCE CODE FOR LT-FH LOADED OR SOURCING IT ####
source("D:/Work/Project1/LTFH/software v2/assign_ltfh.R")
## download from here: https://alkesgroup.broadinstitute.org/UKBB/LTFH/

#age of onset to liability. simulated age is age of onset if indiv is a case.

est_cir = function(data, indivs = c("child", "father", "mother"), ids = c("FID", "pid_f", "pid_m")) {
  
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
  mother_stat = (mother_full > qnorm(prev[2], lower.tail = F)) + 0L,
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


simu_liab$NUM_SIBS = nsib
simu_liab$SIB_STATUS = 0
if (nsib > 0) simu_liab$SIB_STATUS = (rowSums(simu_liab[,paste("sib", 1:nsib, "_stat", sep = "")]) > 0) + 0L  


## more elegant implementaion of estimate_gen_liability_ltfh is pending ##

# LTFH++_fast -------------------------------------------------------------


data = estimate_gen_liability_ltfh(h2 = h2,
                                   phen = simu_liab,
                                   child_threshold = qnorm(mean(prev), lower.tail = F),
                                   parent_threshold = qnorm(mean(prev), lower.tail = F),
                                   status_col_offspring = "child_stat",
                                   status_col_father    = "father_stat",
                                   status_col_mother    = "mother_stat",
                                   status_col_siblings  = "SIB_STATUS",
                                   number_of_siblings_col = "NUM_SIBS")

colnames(data)[ncol(data) - 1:0 ] = paste(colnames(data)[ncol(data) - 1:0 ], "_fast", sep = "")



#  LT-FH++  ---------------------------------------------------------------

cov = get_cov(h2 = h2, n_sib = nsib)

indivs = c("child", "father", "mother", if(nsib > 0) paste("sib", 1:nsib, sep = ""))
ids = c("FID", "pid_f", "pid_m", if(nsib > 0) paste("pid_s", 1:nsib, sep = ""))
thr = est_cir(data = simu_liab,
              indivs = indivs,
              ids = ids)
status_cols = paste(indivs, "_stat", sep = "")


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
h2 = .5
tol = 0.01
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
                 cur_indiv = thr[thr[[1]] == fam[ii], ]
                 
                 if (is.na(cur_status[ii])) {
                   #here to deal with NAs for now 
                 } else if (cur_status[ii] == 1) {
                   lower[ii + 1] <- full_fam[[lia_cols[ii]]] #upper[[ii + 1]] <-  
                 } else {
                   upper[ii + 1] <- cur_indiv$thr
                 }
                 
               }
               fixed <- (upper - lower) < 1e-4
               
               #covergence check
               se = NULL 
               vals = list() #store simulated values
               vals.ctr = 1
               while (is.null(se) || se > tol) {
                 gen_liabs = GibbsSampler::rtmvnorm.gibbs(1e4,
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

with(simu_liab, c("LTFH++" = cov(child_gen, post_gen_liab),"LTFH++_fast" = cov(child_gen, post_gen_liab_fast), "LTFH" = cov(child_gen, ltfh)))^2
with(simu_liab, c("LTFH++" = cor(child_gen, post_gen_liab),"LTFH++_fast" = cor(child_gen, post_gen_liab_fast), "LTFH" = cor(child_gen, ltfh)))

q1 = ggplot(simu_liab, aes(x = post_gen_liab_fast, y = ltfh, color = as.factor(child_stat))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

q2 = ggplot(simu_liab, aes(x = post_gen_liab, y = ltfh, color = as.factor(child_stat))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
grid.arrange(q1, q2)

p1 = ggplot(simu_liab, aes(x = post_gen_liab, y = child_gen, color = as.factor(child_stat))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) 

p2 = ggplot(simu_liab, aes(x = post_gen_liab_fast, y = child_gen, color = as.factor(child_stat))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) 


p3 = ggplot(simu_liab, aes(x = ltfh, y = child_gen, color = as.factor(child_stat))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

grid.arrange(p1, p2, p3)

simu_liab = simu_liab %>% mutate(difference = ltfh - post_gen_liab_fast)

comp1 = ggplot(simu_liab, aes(x = post_gen_liab_fast, y = ltfh, color = as.factor(child_stat))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

comp2 = ggplot(simu_liab %>% arrange(ltfh) %>% select(difference) %>% unique() %>% mutate(grp = 1:n()), aes(y = difference, x = as_factor(grp))) +
  geom_bar(stat = "identity")

grid.arrange(comp1, comp2)

