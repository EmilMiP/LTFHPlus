library(LTFHPlus)
library(tidyverse)
library(doSNOW)
library(progress)


N = 5000
h2_1 = .5
h2_2 = .5
gen_cor = .5
corr_mat = diag(c(h2_1, h2_2))
corr_mat[1,2] <- corr_mat[2,1] <- gen_cor
nthreads = 6  # number of threads to use for ltfh++

#calculates the thresholds used to determine status:
K = .05
#multiplier = 1
#prev = c(0.08, .02) * multiplier
#(thr = qnorm(1 - prev))

#covariate matrix
cov = get_full_cov(corr_mat = corr_mat)

#age of onset to liability. simulated age is age of onset if indiv is a case.
aoo_to_liab = function(age) qnorm( age / 500, lower.tail = FALSE)



#simulate liabilities
liabs = MASS::mvrnorm(n = N, mu = rep(0, 2*(4 + nsib)), Sigma = cov)

child_age = runif(N, 10, 60)
father_age  = child_age + runif(N, 20, 35)
mother_age  = child_age + runif(N, 20, 35)

all_phen = list()
for (i in 1:2) {
  all_phen[[i]] = tibble(
    FID         = 1:N,
    IID         = 1:N,
    child_gen   = liabs[,(i - 1) * 4 + 1],
    child_full  = liabs[,(i - 1) * 4 + 2],
    father_full = liabs[,(i - 1) * 4 + 3],
    mother_full = liabs[,(i - 1) * 4 + 4],
    child_stat  = (child_full  > qnorm(K, lower.tail = F)) + 0L,
    father_stat = (father_full > qnorm(K, lower.tail = F)) + 0L,
    mother_stat = (mother_full > qnorm(K, lower.tail = F)) + 0L,
    child_age   = child_age,
    father_age  = father_age,
    mother_age  = mother_age,
    pid_f       = paste(FID, "_f", sep = ""),
    pid_m       = paste(FID, "_m", sep = "")
  )
}

all_thr = list() 
for (i in 1:2) { #The two tables are identical here, but in real data, we would see differences depending on age of onset, cohort effects etc.
  all_thr[[i]] = tibble(
    FID = c(all_phen[[i]]$FID, all_phen[[i]]$pid_f, all_phen[[i]]$pid_m),
    thr = c(aoo_to_liab(all_phen[[i]]$child_age), aoo_to_liab(all_phen[[i]]$father_age), aoo_to_liab(all_phen[[i]]$mother_age))
  )
}


data = estimate_gen_liability_multi_trait(corr_mat = corr_mat,
                                          phen.list = all_phen,
                                          thr.list = all_thr,
                                          ids = c("FID", "pid_f", "pid_m"),
                                          ind = c(1,5), 
                                          status_cols = c("child_stat", "father_stat", "mother_stat"),
                                          nthreads = nthreads)

cov(data[,c("child_gen", "post_gen_liab_1", "post_gen_liab_2", "child_stat")])
library(ggplot2)
ggplot(data, aes(x = post_gen_liab_1, y = child_gen, color = rowSums(all_phen[[2]][,c("child_stat", "father_stat", "mother_stat")]) > 0)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Any case of \n2nd Phenotype \nin Family") +
  xlab("Estimated Genetic Liability (LTFH++) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
