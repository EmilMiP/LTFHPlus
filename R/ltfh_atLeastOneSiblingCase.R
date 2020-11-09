library(dplyr)
library(stringr)


N = 10000 
h2 = .5 
nsib = 3
nthreads = 5  # number of threads to use for ltfh++
tol = 0.01

#calculates the thresholds used to determine status:
multiplier = 1
prev = c(0.05, 0.05) * multiplier


trans = function(x) qbinom(x, size = 1, prob = .5)


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

ages = simu_liab %>% select(contains("age"))
full = simu_liab %>% select(contains("full"))
stat = simu_liab %>% select(contains("stat"))
sex  = simu_liab %>% select(contains("sex"))

simu_liab = mutate(simu_liab, "NUM_SIB" = nsib,
                   "SIB_STAT" = ((select(simu_liab, paste0("sib", 1:nsib, "_stat")) %>% rowSums) > 0) + 0) 
#, paste0("sib", 1:nsib, "_stat")
configs = simu_liab %>% 
  select(c("child_stat", "father_stat", "mother_stat", "NUM_SIB", "SIB_STAT")) %>%
  distinct() %>% arrange(select(., contains("stat", ignore.case = F)))

# dealing with no sibling cases -------------------------------------------

nosibs_configs = configs %>% filter(SIB_STAT == 0)

nosibs_configs[["est"]] = c()

for (i in 1:nrow(nosibs_configs)) {

  cur_stat = c(unlist(nosibs_configs[i,1:3]),rep(0, nosibs_configs$NUM_SIB[i]))
  lower = rep(-Inf, 4 + nosibs_configs$NUM_SIB[i])
  upper = rep(Inf,  4 + nosibs_configs$NUM_SIB[i])
  upper[-1][cur_stat == 0] = qnorm(mean(prev), lower.tail = F)
  lower[-1][cur_stat == 1] = qnorm(mean(prev), lower.tail = F)
  
  fixed = upper - lower < 1e-3
  
  nosibs_configs$est[i] = mean(LTFHPlus::rtmvnorm.gibbs(n_sim = 100e3,
                                             sigma = LTFHPlus:::get_cov(h2, n_sib = nosibs_configs$NUM_SIB[i]),
                                             lower = lower,
                                             upper = upper,
                                             fixed = fixed,
                                             ind = 1,
                                             burn_in = 1000))
}
# Dealing with "at least one sibling is a case" ---------------------------

sib_configs = configs %>% filter(NUM_SIB > 1, SIB_STAT == 1)

sib_configs[["est_per_sib"]] = list(c())
for (i in 1:nrow(sib_configs)) {
  combs = expand.grid(lapply(1:sib_configs$NUM_SIB[i], FUN = function(x) c(0,1)))
  res = c()
  for (j in 2:nrow(combs)) {
    cur_stat = c(unlist(sib_configs[i,1:3]),unlist(combs[j,]))
    lower = rep(-Inf, 4 + ncol(combs))
    upper = rep(Inf,  4 + ncol(combs))
    upper[-1][cur_stat == 0] = qnorm(mean(prev), lower.tail = F)
    lower[-1][cur_stat == 1] = qnorm(mean(prev), lower.tail = F)

    fixed = upper - lower < 1e-3
    
    res[j - 1] = mean(LTFHPlus::rtmvnorm.gibbs(n_sim = 100e3,
                                  sigma = LTFHPlus:::get_cov(h2, n_sib = ncol(combs)),
                                  lower = lower,
                                  upper = upper,
                                  fixed = fixed,
                                  ind = 1,
                                  burn_in = 1000))
  }
  sib_configs$est_per_sib[[i]] =  as_tibble(cbind(combs[-1,], res)) %>%
    mutate(string = rowSums(select(., 1:ncol(combs)))) %>%
    group_by(string) %>%
    summarise(grp_est = mean(res), .groups = "drop") %>% arrange(string) %>% pull(grp_est)
}

sib_configs$probs = lapply(1:nrow(combs), FUN = function(n) {
  lower = rep(-Inf, 4 + configs$NUM_SIB[n]) 
  upper = rep(Inf, 4 + configs$NUM_SIB[n])
  
  lower[2:4 * (unlist(configs[n, 1:3]) == 1)] = qnorm(mean(prev), lower.tail = F)
  upper[2:4 * (unlist(configs[n, 1:3]) != 1)] = qnorm(mean(prev), lower.tail = F)
    
  tmp <- tmvtnorm::rtmvnorm(n = 10e3, mean = rep(0, 4 + nsib), sigma = LTFHPlus:::get_cov(h2, n_sib = nsib),
                            lower = lower,
                            upper = upper, 
                            algorithm = "gibbs") %>% as_tibble()
  colnames(tmp) = c("child_gen", paste0(c("child", "father", "mother", paste0("sib", 1:nsib)), "_full"))
  
  cur_string = paste(configs[n,1:3], collapse = " ")
  
  mutate_all(tmp, .funs =  pnorm ) %>%
    mutate_all(.funs = trans) %>% 
    select(-child_gen) %>%
    group_by_all() %>%
    summarise(n = n(), .groups = 'drop') %>% 
    mutate("cases" = (select(., contains("sib")) %>% rowSums)) %>% 
    filter(cases != 0) %>%
    select(-contains("sib")) %>%
    group_by(select(., -n)) %>%
    summarise(n = sum(n), .groups = 'drop') %>% 
    mutate(string = do.call(paste, select(.,1:3))) %>% 
    select(-contains("full")) %>% 
    filter(string == cur_string) %>%
    mutate(n_tot = sum(n)) %>%
    mutate(prob = n / n_tot) %>% 
    arrange(cases) %>% 
    pull(prob) 
}) 

sib_configs = sib_configs %>% mutate("est" = sapply((Map("*", est_per_sib, probs)), sum))

ltfh_data_2 = mutate(ltfh_data, "NUM_SIB" = nsib,
       "SIB_STAT" = ((select(ltfh_data, paste0("sib", 1:nsib, "_stat")) %>% rowSums) > 0) + 0) %>% distinct(select(., contains("stat")))

bind_rows(nosibs_configs, sib_configs) %>% 
  select(-est_per_sib, -probs) %>% 
  left_join(ltfh_data_2) %>% 
  group_by(child_stat, father_stat, mother_stat, NUM_SIB, SIB_STAT) %>%
  summarise(mean(est), mean(ltfh)) %>% distinct()

ltfh_data[["SIB_STAT"]] = select(ltfh_data, contains(ltfh_data, "sib"))
ltfh_data %>% group_by(ltfh) %>% summarise(mean(ltfh))
ltfh_data[["NUM_SIB"]] = 3
#simulate liabilities

probs = MASS::mvrnorm(n = 1e6, mu = rep(0, 4 + nsib), Sigma = LTFHPlus:::get_cov(h2, n_sib = nsib)) %>% 
  as_tibble()
colnames(probs)[-1] = colnames(stat)
#  rename_at(vars(starts_with("v")), funs(colnames(stat)))

 # rename(select(. ,-1),  colnames(stat)) colnames(probs)[-1] = colnames(stat)
probs[["SIB_STAT"]] = (rowSums(probs[, ncol(probs) - 0:1] > qnorm(mean(prev), lower.tail = F)) > 0) + 0
probs[["NUM_SIB"]]  = nsib

probs2 = lapply(select(probs, str_subset(colnames(probs), "stat")), FUN = function(n) (n > qnorm(mean(prev), lower.tail = F)) + 0) %>% 
  as_tibble() %>% 
  group_by_all() %>% 
  summarise(prob = n() / 10e3) %>%
  ungroup() %>% 
  mutate("SIB_STAT" = (select(.,contains("sib")) %>% rowSums) ) %>%
  select(-contains("sib",ignore.case = F)) %>%
  group_by(select(.,contains("stat"))) %>% 
  summarise(probs = sum(prob)) %>% 
  ungroup() %>%
  mutate(string = do.call(paste, select(.,1:3))) %>% 
  filter(string == paste(cur_stat[1:3], collapse = " ") & SIB_STAT != 0) %>% 
  left_join(all_est, by = c("SIB_STAT" = "string")) %>%
  mutate(probs = probs / sum(probs)) %>%
  mutate(comb_est = sum(grp_est * probs)) %>% 
  pull(comb_est) %>% 
  unique() 
probs2



# prev 0.05, nsib = 3, h2=.5
# simu_liab %>% 
#   mutate("SIB_STAT" = (select(.,contains("sib")) %>% rowSums) ) %>% 
#   group_by(SIB_STAT) %>%
#   summarise(ltfh = mean(ltfh), .groups = "drop") %>% pull(ltfh) %>% unique()
#gives estimates below:

#-0.1493510  0.3103226  0.3008392  0.6329396  0.7251343  0.9683431  0.8948305  1.1549660  1.1384194  1.3746923  1.6224080
#-0.1493510  0.3103226  0.3008392  0.6329396             0.9683431  0.8948305  1.1549660  1.1384194  1.3746923  1.6224080