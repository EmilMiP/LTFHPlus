library(LTFHPlus)
library(dplyr)
library(stringr)
library(ggplot2) 
library(gridExtra)


N = 10000
h2 = .5 
nsib = 0
nthreads = 20  # number of threads to use for ltfh++

#calculates the thresholds used to determine status:
#K = .05
multiplier = 1
prev = c(0.08, .02) * multiplier


#### THE NEXT SECTION REQUIRES YOU TO HAVE THE SOURCE CODE FOR LT-FH LOADED OR SOURCING IT ####
source("D:/Work/Project1/LTFH/software v2/assign_ltfh.R")
## download from here: https://alkesgroup.broadinstitute.org/UKBB/LTFH/



#covariate matrix
cov = get_cov(h2)

##age of onset to liability. simulated age is age of onset if indiv is a case.
#aoo_to_liab = function(age) qnorm( age / 500, lower.tail = FALSE)



#simulate liabilities
liabs = MASS::mvrnorm(n = N, mu = rep(0, 4 + nsib), Sigma = cov)

child_sex   = sample(1:2, size = N, replace = TRUE)
child_age   = runif(N, 10, 60)
father_age  = child_age + runif(N, 20, 35)
mother_age  = child_age + runif(N, 20, 35)

simu_liab = tibble(
  FID         = 1:N,
  IID         = 1:N,
  child_gen   = liabs[,1],
  child_full  = liabs[,2],
  father_full = liabs[,3],
  mother_full = liabs[,4],
  child_sex   = child_sex,
  child_stat  = (child_full  > qnorm(prev[child_sex], lower.tail = F)) + 0L,
  father_stat = (father_full > qnorm(prev[1], lower.tail = F)) + 0L,
  mother_stat = (mother_full > qnorm(prev[2], lower.tail = F)) + 0L,
  child_age   = child_age,
  father_age  = father_age,
  mother_age  = mother_age,
  pid_f       = paste(FID, "_f", sep = ""),
  pid_m       = paste(FID, "_m", sep = "")
)

indivs = c("child", "father", "mother")
all_liabs = c()
all_stat = c()
for (i in seq_along(indivs)) {
  cur_stat = simu_liab[[paste(indivs[i], "_stat", sep ="")]] == 0
  cur_thr = rep(NA, length(cur_stat))
  if (!(indivs[i] %in% c("father", "mother"))) {
    cur_thr[cur_stat]  = qnorm(prev, lower.tail = F)[simu_liab[[paste(indivs[i], "_sex", sep = "")]][cur_stat]]
    
  }
  if (indivs[i] == "father"){
    cur_thr[cur_stat] = qnorm(prev[1], lower.tail = F)
  }
  if (indivs[i] == "mother"){
    cur_thr[cur_stat] = qnorm(prev[2], lower.tail = F)
  }
  cur_thr[!cur_stat] = simu_liab[[paste(indivs[i], "_full", sep = "")]][!cur_stat]
  all_liabs = c(all_liabs, cur_thr)
  all_stat  = c(all_stat, cur_stat)
}


thr = tibble(
  ids = c(simu_liab$FID, simu_liab$pid_f, simu_liab$pid_m),
  thr = all_liabs
)


data = LTFHPlus::estimate_gen_liability(h2 = h2,
                                         phen = simu_liab,
                                         thr = thr,
                                         ids = c("FID", "pid_f", "pid_m"),
                                         status_cols = c("child_stat", "father_stat", "mother_stat"),
                                         nthreads = nthreads)

res = as_tibble(as.data.frame(matrix(NA, nrow = nrow(simu_liab), ncol = 7)))
res[,1] = as.double(simu_liab$FID)
res[,2] = as.double(simu_liab$IID)
colnames(res) = c("FID", "IID", "CHILD_STATUS", "P1_STATUS", "P2_STATUS", "NUM_SIBS", "SIB_STATUS")
res$CHILD_STATUS = simu_liab$child_stat
res$P1_STATUS = simu_liab$father_stat
res$P2_STATUS = simu_liab$mother_stat
res$NUM_SIBS = 0
res$SIB_STATUS = 0

ltfh = create_pheno(data = as.data.frame(res),
                    trait_h2 = h2,
                    T_val_child = qnorm(mean(prev), lower.tail = F),
                    T_val_parent = qnorm(mean(prev), lower.tail = F),
                    relevant_trait_child = "CHILD_STATUS",
                    relevant_trait_dad = "P1_STATUS",
                    relevant_trait_mom = "P2_STATUS",
                    number_siblings_col = "NUM_SIBS",
                    relevant_trait_sib = "SIB_STATUS",
                    maximum_siblings_to_compute = 0)

data = left_join(data, as.data.frame(ltfh))


with(data, c("LTFH++" = cov(child_gen, post_gen_liab), "LTFH" = cov(child_gen, ltfh), "Case-Control" = cov(child_gen, child_stat)))^2


summary(data$post_gen_liab)
summary(data$ltfh)




p1 = ggplot(data, aes(x = post_gen_liab, y = child_gen, color = rowSums(data[,c("child_stat", "father_stat", "mother_stat")]) > 0)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status in Family") +
  xlab("Estimated Genetic Liability (LTFH++) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
#  xlim(-0.4, 2) +
  theme(plot.title = element_text(hjust = 0.5))


p2 = ggplot(data, aes(x = ltfh, y = child_gen, color = rowSums(data[,c("child_stat", "father_stat", "mother_stat")]) > 0)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status in Family") +
  xlab("Estimated Genetic Liability (LTFH++) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
 # xlim(-0.4, 2) +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2)


