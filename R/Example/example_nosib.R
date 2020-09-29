library(LTFHPlus)
library(dplyr)
library(stringr)
library(ggplot2) 
library(gridExtra)


N = 5000 
h2 = .5 
nsib = 0
nthreads = 6  # number of threads to use for ltfh++

#calculates the thresholds used to determine status:
K = .05
multiplier = 1
prev = c(0.08, .02) * multiplier
(thr = qnorm(1 - prev))

#### THE NEXT SECTION REQUIRES YOU TO HAVE THE SOURCE CODE FOR LT-FH LOADED OR SOURCING IT ####
source("C:/Code/LTFH/assign_ltfh.R")
## download from here: https://alkesgroup.broadinstitute.org/UKBB/LTFH/



#covariate matrix
cov = LTFHPlus:::get_cov(h2)

#age of onset to liability. simulated age is age of onset if indiv is a case.
aoo_to_liab = function(age) qnorm( age / 500, lower.tail = FALSE)



#simulate liabilities
liabs = MASS::mvrnorm(n = N, mu = rep(0, 4 + nsib), Sigma = cov)

simu_liab = tibble(
  FID         = 1:N,
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

thr2 = tibble(
  ids = c(simu_liab$FID, simu_liab$pid_f, simu_liab$pid_m),
  thr = mean(qnorm(K, lower.tail = F)))


## more elegant implementaion of estimate_gen_liability_ltfh is pending ##
data = estimate_gen_liability_ltfh(h2 = h2,
                               phen = simu_liab,
                               thr = thr2,
                               ids = c("FID", "pid_f", "pid_m"),
                               status_cols = c("child_stat", "father_stat", "mother_stat"))

colnames(data)[ncol(data) - 1:0 ] = paste(colnames(data)[ncol(data) - 1:0 ], "_fast", sep = "")

data2 = LTFHPlus::estimate_gen_liability(h2 = h2,
                                         phen = simu_liab,
                                         thr = thr2,
                                         ids = c("FID", "pid_f", "pid_m"),
                                         status_cols = c("child_stat", "father_stat", "mother_stat"),
                                         nthreads = nthreads)
data = left_join(data, data2)

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
                    T_val_child = qnorm(K, lower.tail = F),
                    T_val_parent = qnorm(K, lower.tail = F),
                    relevant_trait_child = "CHILD_STATUS",
                    relevant_trait_dad = "P1_STATUS",
                    relevant_trait_mom = "P2_STATUS",
                    number_siblings_col = "NUM_SIBS",
                    relevant_trait_sib = "SIB_STATUS",
                    maximum_siblings_to_compute = 0)

simu_liab = left_join(data, as.data.frame(ltfh))




#Morale of this example: The estimatation methods all perform similarly, but with a threshold for each group, estimate_gen_liability_ltfh, is the fastest.

with(simu_liab, c(cov(child_gen, post_gen_liab), cov(child_gen, post_gen_liab_fast), cov(child_gen, ltfh)))^2

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
