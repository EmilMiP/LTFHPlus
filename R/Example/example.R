library(LTFHPlus)
library(dplyr)
library(stringr)
library(ggplot2) 
library(gridExtra)
library(progressr)
library(future)
library(doFuture)


N = 5000 
h2 = .5 
nsib = 3
nthreads = 5  # number of threads to use for ltfh++
tol = 0.01

#calculates the thresholds used to determine status:
multiplier = 1
prev = c(0.08, 0.02) * multiplier

#### THE NEXT SECTION REQUIRES YOU TO HAVE THE SOURCE CODE FOR LT-FH LOADED OR SOURCING IT ####
source("C:/Code/LTFH/assign_ltfh.R")
## download from here: https://alkesgroup.broadinstitute.org/UKBB/LTFH/

#plan(multisession)
handlers("progress")


#simulate liabilities
liabs = MASS::mvrnorm(n = N, mu = rep(0, 4 + nsib), Sigma = LTFHPlus:::get_cov(h2, n_sib = nsib))

#initialize phenotype tibble with to fill out
phen = tibble(
  ids = character(N),
  sex = numeric(N),
  status = numeric(N),
  age = numeric(N),
  aoo = numeric(N)
)
#filling out tibble with list entries, each row is a family. assuming order: os, father, mother, siblings
for(i in 1:N) {
  id_vec     = c(i , paste0(i, c("_f", "_m")), if(nsib > 0) paste0(i, "_s", 1:nsib) )
  sex_vec    = c(sample(1:2, size = 1, replace = TRUE), 1:2,
                 if(nsib > 0) sample(1:2, size = nsib, replace = TRUE))
  status_vec = sapply(seq_along(sex_vec), function(ii) {
    (liabs[,-1][i,ii] > qnorm(prev[sex_vec[ii]], lower.tail = F)) + 0L
  })
  cur_child_age = runif(1, 10, 60)
  age_vec       = cur_child_age + c(0 , runif(2, 20, 35), if (nsib > 0) runif(nsib, -15, 15))
  
  aoo_vec = rep(NA, length(id_vec))
  cases = status_vec == 1
  aoo_vec[cases] = liab_to_aoo(liab = liabs[,-1][i, cases],
                                         pop_prev = prev[sex_vec[cases]])
  
  phen$ids[i]    = list(id_vec)
  phen$sex[i]    = list(sex_vec)
  phen$status[i] = list(status_vec)
  phen$age[i]    = list(age_vec)
  phen$aoo[i]    = list(aoo_vec)
}

#constructing input for LT-FH
res = tibble(
  ids          = sapply(phen$ids,    function(x) x[1]),
  CHILD_STATUS = sapply(phen$status, function(x) x[1]),
  P1_STATUS    = sapply(phen$status, function(x) x[2]),
  P2_STATUS    = sapply(phen$status, function(x) x[3]),
  NUM_SIBS     = nsib,
  SIB_STATUS   = sapply(phen$status, function(x) (sum(x[-(1:3)]) > 0) + 0L) 
)


# LTFH++_fast -------------------------------------------------------------

with_progress({data = estimate_gen_liability_ltfh(h2 = h2,
                                   phen = res,
                                   child_threshold = qnorm(mean(prev), lower.tail = F),
                                   parent_threshold = qnorm(mean(prev), lower.tail = F),
                                   status_col_offspring = "CHILD_STATUS",
                                   status_col_father    = "P1_STATUS",
                                   status_col_mother    = "P2_STATUS",
                                   status_col_siblings  = "SIB_STATUS",
                                   number_of_siblings_col = "NUM_SIBS")
})

colnames(data)[ncol(data) - 1:0 ] = paste(colnames(data)[ncol(data) - 1:0 ], "_fast", sep = "")



#  LT-FH++  ---------------------------------------------------------------

#assigning thresholds for each family
tmp = lapply(1:nrow(phen), function(n) {
  ids   = phen$ids[[n]]
  cases = phen$status[[n]] == 1
  lower = rep(-Inf, length(ids))
  upper = rep( Inf, length(ids))
  
  lower[cases] <- upper[cases] <- age_to_thres(age = phen$aoo[[n]][cases], pop_prev = prev[phen$sex[[n]][cases]])
  
  upper[!cases] <- age_to_thres(age = phen$age[[n]][!cases], pop_prev = prev[phen$sex[[n]][!cases]])
  return(tibble(ids = ids, lower = lower, upper = upper))
})
#combining all values into one tibble
thr = do.call("rbind", tmp)


#performs LT-FH++ analysis with progress bar
with_progress({simu_liab = estimate_gen_liability(h2 = h2, 
                       phen = phen, 
                       thr = thr,  
                       id_col = "ids",
                       nthreads = nthreads)
})

#combining data
simu_liab$ids = as.character(1:N)
data = left_join(data, simu_liab %>% select(post_gen_liab, post_gen_liab_se, ids))



#  LT-FH ------------------------------------------------------------------
#running LT-FH
colnames(res)[1] = "IID"
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
#combining data
ltfh = apply(ltfh, MARGIN = 2, FUN = as.numeric)
data$IID = 1:N
simu_liab = left_join(data, as.data.frame(ltfh))
simu_liab$child_gen = liabs[,1]



# Plots & Results ---------------------------------------------------------

simu_liab[["case-control"]] = sapply(phen$status, function(x) x[1])

#Morale of this example: The estimatation methods all perform similarly, but with a threshold for each group, estimate_gen_liability_ltfh, is the fastest.

with(simu_liab, c("LTFH++" = cov(child_gen, post_gen_liab), "LTFH_gibbs" = cov(child_gen, post_gen_liab_fast), "LTFH" = cov(child_gen, ltfh)))
with(simu_liab, c("LTFH++" = cor(child_gen, post_gen_liab), "LTFH_gibbs" = cor(child_gen, post_gen_liab_fast), "LTFH" = cor(child_gen, ltfh)))

#comparing LTFH to LTFH++_fast
simu_liab = simu_liab %>% mutate(difference = ltfh - post_gen_liab_fast)

comp1 = ggplot(simu_liab, aes(x = post_gen_liab_fast, y = ltfh, color = as.factor(CHILD_STATUS))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

comp2 = ggplot(simu_liab %>% arrange(ltfh) %>% select(difference) %>% unique() %>% mutate(grp = 1:n()), aes(y = difference, x = as.factor(grp))) +
  geom_bar(stat = "identity")

grid.arrange(comp1, comp2)



#plotting estimating genetic liability for the different methods against the true genetic liabilities:
p1 = ggplot(simu_liab, aes(x = post_gen_liab, y = child_gen, color = sapply(phen$status, function(x) sum(x) > 0 ))) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status in Family") +
  xlab("Estimated Genetic Liability (LTFH++) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")


p2 = ggplot(simu_liab, aes(x = post_gen_liab_fast, y = child_gen, color = sapply(phen$status, function(x) sum(x) > 0 ))) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status in Family") +
  xlab("Estimated Genetic Liability (LTFH_gibbs) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

p3 = ggplot(simu_liab, aes(x = ltfh, y = child_gen, color = sapply(phen$status, function(x) sum(x) > 0 ))) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Status in Family") +
  xlab("Estimated Genetic Liability (LTFH) ") +
  ylab("True Genetic Liability") + 
#  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top")


#grid.arrange(p1, p2, p3)

grid.arrange(p1, p3)


