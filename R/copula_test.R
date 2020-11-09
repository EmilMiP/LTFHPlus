

h2_1 = .5
h2_2 = .5
corr = .5
nsib = 3

trans = function(x) qbinom(x, size = 1, prob = .5)

#simulate liabilities
liabs = MASS::mvrnorm(n = 100e3, mu = rep(0, 4 + nsib), Sigma = LTFHPlus:::get_cov(h2, n_sib = nsib)) %>% as_tibble()
colnames(liabs) = c("child_gen", paste0(c("child", "father", "mother", paste0("sib", 1:nsib)), "_full"))
#select(liabs, contains("sib")) %>% 

#simulate from a truncated normal, so we do not need to "bulk sample"

tmp <- tmvtnorm::rtmvnorm(n = 10e3, mean = rep(0, 4 + nsib), sigma = LTFHPlus:::get_cov(h2, n_sib = nsib),
                   lower = c(-Inf, rep(qnorm(mean(prev), lower.tail = F), 3), rep(-Inf, nsib)),
                   upper = rep(Inf, 4 + nsib), 
                   algorithm = "gibbs") %>% as_tibble()
colnames(tmp) = c("child_gen", paste0(c("child", "father", "mother", paste0("sib", 1:nsib)), "_full"))

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
  filter(string == paste(cur_stat[1:3], collapse = " ")) %>%
  mutate(n_tot = sum(n)) %>%
  mutate(prob = n / n_tot) %>% 
  select(cases, prob) %>% 
  left_join(all_est, by = c("cases" = "string")) %>%
  mutate(comb_est = sum(grp_est * prob)) %>% 
  pull(comb_est) %>% 
  unique()
  

  group_by(cases) %>% 
  summarise(n_tot = sum(n))
paste(cur_stat[1:3], collapse = " ")
 cor(stat_tmp)
sigma = matrix(corr / 2, ncol = nsib, nrow = nsib)
diag(sigma) = 1
cop = MASS::mvrnorm(n = 1e6, mu = rep(0, nsib), Sigma = sigma)
cor(cop)

unifs = pnorm(cop)
qbinom(unifs, size = 1, prob = .005) %>% 
  as_tibble() %>% 
  group_by_all() %>% 
  summarise(prob = n() / 1e6) %>%
  ungroup() %>% 
  mutate("SIB_STAT" = (select(.,contains("v")) %>% rowSums) ) %>%
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
cor(tmp)



stats = cop > qnorm(.9, sd = sqrt(1 - corr^2))
table(stats[,1], stats[,2]) / 1e6
cor(stats)
corr
chisq = qchisq(pnorm(cop), df = 1)
cor(chisq)

diff_fct = function(x, target_corr = corr,  prob_dat = cop) {
  prob1 = mean(prob_dat[,1] >= x)#pnorm(x, lower.tail = F)
  prob2 = mean(prob_dat[,2] >= x) #pnorm(x, lower.tail = F)
  prob12 = mean(rowSums(prob_dat >= x) == 2)
  return(target_corr - ( prob12 - prob1 * prob2) / (sqrt(prob1 - prob1^2) * sqrt(prob2 - prob2^2)))
}

vals = numeric(1000)
xs = seq(-4,4, length.out = 1000)
for(i in 1:1000) {
  if(i %% 100 == 0) print(i)
  vals[i] = diff_fct(x = xs[i])
}
plot(x = xs, y = vals)
curve(diff_fct, from = -4, to = 4)
plot(x = seq(-1,1,.01),y = diff_fct(seq(-1,1,.01)))

diff_fct2 = function(x, target_corr = corr, prob_dat = cop) {
  prob1 = mean(prob_dat[,1] >= x)#pnorm(x, lower.tail = F)
  prob2 = mean(prob_dat[,2] >= x) #pnorm(x, lower.tail = F)
  prob12 = mean(rowSums(prob_dat >= x) == 2)
  return(( prob12 - prob1 * prob2) / (sqrt(prob1 - prob1^2) * sqrt(prob2 - prob2^2)))
}
diff_fct2(x = 4)
plot(x = seq(-1,1,.01),y = diff_fct2(seq(-1,1,.01)))

diff_fct(x = 1)
curve(diff_fct2, from = -1, to = 1, xname = "x")

abline(h = 0)
uniroot(diff_fct, c(-5, 5))
diff_fct(-1)
diff_fct(1)
diff_fct(0)
diff_fct(4)
diff_fct(.5)

cor(liabs[,6 - 1:0])
cor(simu_liab[, c("sib1_stat", "sib2_stat")])




mean(rowSums(stats) == 1)
mean(rowSums(stats) == 2)


10e3 * .05

expand.grid(lapply(1:nsib, FUN = function(x) c(0,1)))
