library(LTFHPlus)
library(dplyr)
library(stringr)
library(ggplot2) 
library(gridExtra)
library(tidyverse)

N = 5000 
h2 = .5 
nsib = 2
nthreads = 6  # number of threads to use for ltfh++

#calculates the thresholds used to determine status:
K = .05
multiplier = 1
prev = c(0.08, .02) * multiplier

#### THE NEXT SECTION REQUIRES YOU TO HAVE THE SOURCE CODE FOR LT-FH LOADED OR SOURCING IT ####
source("C:/Code/LTFH/assign_ltfh.R")
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

indivs = c("child", "father", "mother", if(nsib > 0) paste("sib", 1:nsib, sep = ""))
ids = c("FID", "pid_f", "pid_m", if(nsib > 0) paste("pid_s", 1:nsib, sep = ""))

if (nsib > 0) {
  sib_vec = paste("sib", 1:nsib, "_stat", sep = "")
} else {
  sib_vec = ""
}
ids = c("FID", "pid_f", "pid_m", if(nsib > 0) paste("pid_s", 1:nsib, sep = ""))
thr = simu_liab[,ids] %>% 
  pivot_longer(cols = ids,names_to = "IDs") %>% 
  select(IDs) %>% 
  mutate(thr = qnorm(mean(prev), lower.tail = F))
  

## more elegant implementaion of estimate_gen_liability_ltfh is pending ##
data = estimate_gen_liability_ltfh(h2 = h2,
                               phen = simu_liab,
                               thr = thr,
                               ids = ids,
                               status_col_offspring = "child_stat",
                               status_col_parents = c("father_stat", "mother_stat"),
                               status_col_siblings = sib_vec)

colnames(data)[ncol(data) - 1:0 ] = paste(colnames(data)[ncol(data) - 1:0 ], "_fast", sep = "")
# 
# thr2 = est_cir(data = simu_liab, 
#                indivs = indivs,
#                ids = ids)
# 
# 
# data2 = LTFHPlus::estimate_gen_liability(h2 = h2,
#                                          phen = simu_liab,
#                                          thr = thr2,
#                                          ids = ids,
#                                          status_cols = paste(indivs, "_stat", sep = ""),
#                                          nthreads = nthreads)
# data = left_join(data, data2)

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

0.000964305 / sqrt(105979)
qnorm(0.7805/2)

se = sqrt(1 / (36306.88 ) * 0.8470 * 2 * (1 - 0.8470))
se * -0.279
P = 0.0002707
CHISQ = 13.26
SE = 0.02414
1 / sqrt(2 * 0.152364 * (1 - 0.152364) * (105979 + (qnorm(6.8e-1/2))^2) )
