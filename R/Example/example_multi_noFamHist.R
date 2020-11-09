library(LTFHPlus)
library(dplyr)
library(stringr)
library(ggplot2) 
library(gridExtra)
library(doSNOW)
library(progress)

N = 5000 
h2 = .5 
gen_cor = .5
ntraits = 4
nthreads = 5


#calculates the thresholds used to determine status:
multiplier = 2
prev = c(0.08, 0.02) * multiplier

#including gen liab for all disorders:
cor_mat = matrix(gen_cor * h2, ncol = 2*ntraits, nrow = 2*ntraits)
for (i in 1:ntraits) {
  cor_mat[(i-1)*2 + 1:2,(i-1)*2 + 1:2] = matrix(c(h2,h2,h2,1), ncol = 2, nrow =2)
}


cor_mat = matrix(gen_cor * h2, ncol = 1 + ntraits, nrow = 1 + ntraits)
diag(cor_mat) = 1
cor_mat[1:2,1:2] = matrix(c(h2,h2,h2,1), ncol = 2, nrow =2)


#simulate liabilities
liabs = MASS::mvrnorm(n = N, mu = rep(0, 2*ntraits), Sigma = cor_mat)

simu_liab = list()

for (i in 1:ntraits) {
  simu_liab[[paste0("child_gen_", i)]]  = liabs[,2*(i - 1) + 1]
  simu_liab[[paste0("child_full_", i)]] = liabs[,2*(i - 1) + 2]
  simu_liab[[paste0("child_sex_", i)]]  = sample(1:2, size = N, replace = TRUE)
  simu_liab[[paste0("child_stat_", i)]] = (simu_liab[[paste0("child_full_", i)]]  > qnorm(prev[simu_liab[[paste0("child_sex_", i)]]], lower.tail = F)) + 0L
  simu_liab[[paste0("child_age_", i)]]  = runif(N, 10, 100)
}
simu_liab = as_tibble(simu_liab)

select(simu_liab, str_subset(colnames(simu_liab), "stat")) %>% colMeans()

ages = select(simu_liab, str_subset(colnames(simu_liab), "age"))
full = select(simu_liab, str_subset(colnames(simu_liab), "full"))
stat = select(simu_liab, str_subset(colnames(simu_liab), "stat"))
sex  = select(simu_liab, str_subset(colnames(simu_liab), "sex"))

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
               
               
               lower = rep(-Inf, 1 + ntraits)
               upper = rep(Inf,  1 + ntraits)
               
               upper[1 + 1:ntraits] = LTFHPlus::age_to_thres(age = unlist(ages[i,]), pop_prev = prev[unlist(sex[i,])])
               for (ii in 1:ncol(stat)) {
                 if(stat[i,ii][[1]] == 1) {
                   lower[1 + ii] = LTFHPlus::age_to_thres(age = unlist(ages[i,ii]), pop_prev = prev[unlist(sex[i,ii])])
                 }
               }
               
               fixed = upper - lower < 1e-3
               
               
               mean(LTFHPlus::rtmvnorm.gibbs(n_sim = 10e3,
                                             sigma = cor_mat,
                                             lower = lower,
                                             upper = upper,
                                             fixed = fixed,
                                             ind = 1,
                                             burn_in = 1000))    



}
stopCluster(cl)

simu_liab[["post_gen_no_fam"]] = unlist(ph)

with(simu_liab, c(cov(child_stat_1, child_gen_1), cov(post_gen_no_fam, child_gen_1))) # 0.07809969 0.12120399 #with all gen liabs
with(simu_liab, c(cor(child_stat_1, child_gen_1), cor(post_gen_no_fam, child_gen_1))) # 0.3762194 0.4285388  #with all gen liabs
#it does not appear to matter if the gen liab is included for all disorders. cov and cor were identical to 3rd or 4th decimal point.

#plotting estimating genetic liability for the different methods against the true genetic liabilities:

simu_liab = mutate(simu_liab, nb_stats = select(simu_liab, str_subset(colnames(simu_liab), "stat")) %>% rowSums())

ggplot(simu_liab, aes(x = post_gen_no_fam, y = child_gen_1, color = as.factor(nb_stats))) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status") +
  xlab("Estimated Genetic Liability (LTFH++) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) #+ 
  #coord_equal()

ggplot(simu_liab, aes(x = child_stat_1, y = child_gen_1, color = as.factor(nb_stats))) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status") +
  xlab("Estimated Genetic Liability (LTFH++) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) #+ 
  #coord_equal()
