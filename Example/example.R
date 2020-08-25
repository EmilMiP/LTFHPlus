#library(LTFHPlus)
library(tidyverse)
library(doSNOW)
library(progress)

N = 5000
nsib = 0 # doesnt actaully do anthing, shows we have 0 siblings in the below example.
h2 = .5
nthreads = 8  # number of threads to use for ltfh++

#calculates the thresholds used to determine status:
K = .05
multiplier = 1
prev = c(0.08, .02) * multiplier
(thr = qnorm(1 - prev))

#covariate matrix
cov = get_cov(h2)

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
  child_stat  = (child_full  > qnorm(K, lower.tail = F)) + 0L,
  father_stat = (father_full > qnorm(K, lower.tail = F)) + 0L,
  mother_stat = (mother_full > qnorm(K, lower.tail = F)) + 0L,
  child_age   = runif(N, 10, 60),
  father_age  = child_age + runif(N, 20, 35),
  mother_age  = child_age + runif(N, 20, 35),
  pid_f       = paste(FID, "_f", sep = ""),
  pid_m       = paste(FID, "_m", sep = "")
)

thr = tibble(
  ids = c(simu_liab$FID, simu_liab$pid_f, simu_liab$pid_m),
  thr = c(aoo_to_liab(simu_liab$child_age), aoo_to_liab(simu_liab$father_age), aoo_to_liab(simu_liab$mother_age))
)


data = estimate_gen_liability(h2 = h2,
                              phen = simu_liab,
                              thr = thr,
                              ids = c("FID", "pid_f", "pid_m"),
                              status_cols = c("child_stat", "father_stat", "mother_stat"))
