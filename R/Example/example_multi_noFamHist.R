library(LTFHPlus)
library(dplyr)
library(stringr)
library(ggplot2) 
library(gridExtra)

N = 10000 
h2 = .5 
gen_cor = .5
ntraits = 4
nthreads = 5
prev = .5

#uncomment to and fill with other values - if not positive semidefinite use Matrix::nearPD() for approximation 
# h2_vec  = c(.6, .6, .4, .3)
# gen_cor_vec = c(.3, .4, .3,
#                     .3, .5,
#                         .2)

h2_vec  = rep(h2, ntraits) 
gen_cor_vec = rep(gen_cor, ntraits*(ntraits - 1)/2 )

prev_vec = rep(mean(prev), ntraits)

cov_mat = generate_cov_matrix_noFH(h2_vec = h2_vec,
                                   gen_cor_vec = gen_cor_vec)

#calculates the thresholds used to determine status:
multiplier = 1
prev = c(0.08, 0.02) * multiplier


#simulate liabilities & make tibble
liabs = MASS::mvrnorm(n = N, mu = rep(0, 1 + ntraits), Sigma = cov_mat)

simu_liab = list()
simu_liab[[paste0("child_gen_", 1)]]  = liabs[,1]
for (i in 1:ntraits) {
  simu_liab[[paste0("child_full_", i)]] = liabs[,1 + i ]
  simu_liab[[paste0("child_sex_", i)]]  = sample(1:2, size = N, replace = TRUE)
  simu_liab[[paste0("child_stat_", i)]] = (simu_liab[[paste0("child_full_", i)]]  > qnorm(prev[simu_liab[[paste0("child_sex_", i)]]], lower.tail = F)) + 0L
  simu_liab[[paste0("child_age_", i)]]  = runif(N, 10, 100)
}
simu_liab = as_tibble(simu_liab)

#extract relevant information
phen = simu_liab %>% 
  select(., 1, contains("stat"))



phen = LTFHPlus::estimate_gen_liability_noFH_corr_traits(phen = phen,
                                                  h2_vec = h2_vec,
                                                  gen_cor_vec = gen_cor_vec,
                                                  prev_vec = prev_vec)
simu_liab = left_join(simu_liab, phen)

with(simu_liab, c(cov(child_stat_1, child_gen_1), cov(post_gen_no_fam, child_gen_1))) # 0.07809969 0.12120399 #with all gen liabs
with(simu_liab, c(cor(child_stat_1, child_gen_1), cor(post_gen_no_fam, child_gen_1))) # 0.3762194 0.4285388  #with all gen liabs
#it does not appear to matter if the gen liab is included for all disorders. cov and cor were identical to 3rd or 4th decimal point.

#plotting estimating genetic liability for the different methods against the true genetic liabilities:

ggplot(simu_liab, aes(x = post_gen_no_fam, y = child_gen_1, color = as.factor(child_stat_1))) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status") +
  xlab("Estimated Genetic Liability (LTMT) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_smooth(color = "black", method = "lm")#+ 
  #coord_equal()


