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

#constructs covariance matrix with a baseling of 2 parents and n_sib siblings (with-in disorder):
get_cov = function(h2, n_sib = 0) {
  cov <- matrix(h2/2, 4 + n_sib, 4 + n_sib)
  diag(cov) <- 1
  cov[3,4] <- cov[4,3] <- 0
  cov[1:2, 1] <- cov[1, 1:2] <- h2
  cov
} 
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
  child_stat  = (child_full  > qnorm(K, lower.tail = F)) + 0L,
  father_stat = (father_full > qnorm(K, lower.tail = F)) + 0L,
  mother_stat = (mother_full > qnorm(K, lower.tail = F)) + 0L,
  child_age   = runif(N, 10, 60),
  father_age  = child_age + runif(N, 20, 35),
  mother_age  = child_age + runif(N, 20, 35),
)

simu_liab[["post_gen_liab"]]    <- NA
simu_liab[["post_gen_liab_se"]] <- NA


stat = paste(c("child", "father", "mother"), "_stat", sep = "")
age  = paste(c("child", "father", "mother"), "_age" , sep = "")
liab = paste(c("child", "father", "mother"), "_full", sep = "")



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
             .options.snow = opts,
             .inorder = TRUE) %dopar% {
               #  if (i %% 500 == 0) print(i)
               lower = rep(-Inf, ncol(cov))
               upper = rep(Inf, ncol(cov))
               x <- simu_liab[i, ]
               for (ii in 1:3) {
                 if (x[[stat[ii]]] == 1) {
                   lower[1 + ii] <- upper[1 + ii] <- x[[liab[ii]]]
                 } else {
                   upper[1 + ii] <- aoo_to_liab(x[[age[ii]]])
                 }
                 
               }
               
               fixed <- (upper - lower) < 1e-4
               
               gen_liabs <- LTFHPlus::rtmvnorm.gibbs(10e3, burn_in = 1000, sigma = cov, ### Using GibbsSampler Package ####
                                                         lower = lower, upper = upper, fixed = fixed)
               batchmeans::bm(gen_liabs)
             }
stopCluster(cl)
simu_liab[["post_gen_liab"]]    = sapply(ph, FUN = function(x) x$est)
simu_liab[["post_gen_liab_se"]] = sapply(ph, FUN = function(x) x$se)

#### THE NEXT SECTION REQUIRES YOU TO HAVE THE SOURCE CODE FOR LT-FH LOADED OR SOURCING IT ####
source("../Project1/LTFH/software v2/assign_ltfh.R")

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

simu_liab = left_join(simu_liab, as.data.frame(ltfh))


library(ggplot2) 
library(gridExtra)


with(simu_liab, c(cov(child_gen, post_gen_liab), cov(child_gen, ltfh)))^2

p1 = ggplot(simu_liab, aes(x = post_gen_liab, y = child_gen, color = as.factor(child_stat))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) 


p2 = ggplot(simu_liab, aes(x = ltfh, y = child_gen, color = as.factor(child_stat))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

grid.arrange(p1, p2)
