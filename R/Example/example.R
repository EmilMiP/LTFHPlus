library(LTFHPlus)
library(ggplot2) 
library(gridExtra)
library(dplyr)



N = 10000 #number of individuals
h2 = .5 #heritability
nsib = 0 #number of siblings
nthreads = 5  # number of threads to use for ltfh++


#calculates the sex-specific prevalences used to determine status:
multiplier = 1
prev = c(0.08, 0.02) * multiplier

#### THE NEXT SECTION REQUIRES YOU TO HAVE THE SOURCE CODE FOR LT-FH LOADED OR SOURCING IT ####
source("D:/Work/LTFH/assign_ltfh.R")
## download from here: https://alkesgroup.broadinstitute.org/UKBB/LTFH/

#start sessions to parallelize over later
plan(multisession)



#simulate liabilities
liabs = MASS::mvrnorm(n = N, mu = rep(0, 4 + nsib), Sigma = LTFHPlus:::get_cov(h2, n_sib = nsib))

#Age of children, parents, and siblings
child_age = runif(N, 10, 60)
father_age  = child_age + runif(N, 20, 35)
mother_age  = child_age + runif(N, 20, 35)
children_sex = sample(1:2, size = N, replace = TRUE)
if(nsib > 0) {
  sibling_sex = replicate(nsib, sample(1:2, size = N, replace = TRUE))
  sibling_age = replicate(nsib, runif(N, -15, 15) + child_age)
}

#Constructing tibble with list entries, each row is a family. assuming order: os, father, mother, siblings (if any)
phen = lapply(1:N, function(i) {
  cur_liabs  = liabs[i, 1:(4 + nsib)][-1]
  id_vec     = c(i , paste0(i, c("_f", "_m")), if(nsib > 0) paste0(i, "_s", 1:nsib) )
  sex_vec    = c(children_sex[i], 1:2, if(nsib > 0) sibling_sex[i,])
  age_vec    = c(child_age[i], father_age[i], mother_age[i], if (nsib > 0) sibling_age[i,])
  status_vec = sapply(seq_along(sex_vec), function(ii) {
    (cur_liabs[ii] > qnorm(prev[sex_vec[ii]], lower.tail = F)) + 0L
  })
  aoo_vec = rep(NA, length(id_vec))
  
  cases = status_vec == 1
  aoo_vec[cases] = LTFHPlus::liab_to_aoo(liab = cur_liabs[cases],
                                         pop_prev = prev[sex_vec[cases]])
  c(list(id_vec),
    list(sex_vec),
    list(status_vec),
    list(age_vec),
    list(aoo_vec))
}) %>% do.call("rbind", . ) %>% #outputs list, combining with rbind 
  as_tibble(.name_repair = 'unique') %>% # converting to tibble
  rename("ids" = "V1", "sex" = "V2", "status" = "V3", "age" = "V4", "aoo" = "V5") # renaming columns


#constructing input for LT-FH
res = tibble(
  ids          = sapply(phen$ids,    function(x) x[1]),
  CHILD_STATUS = sapply(phen$status, function(x) x[1]),
  P1_STATUS    = sapply(phen$status, function(x) x[2]),
  P2_STATUS    = sapply(phen$status, function(x) x[3]),
  NUM_SIBS     = nsib,
  SIB_STATUS   = sapply(phen$status, function(x) (sum(x[-(1:3)]) > 0) + 0L) 
)


# LTFH with gibbs sampler --------------------------------------------------

data = LTFHPlus::estimate_gen_liability_ltfh(h2 = h2,
                                   phen = res,
                                   child_threshold = qnorm(mean(prev), lower.tail = F),
                                   parent_threshold = qnorm(mean(prev), lower.tail = F),
                                   status_col_offspring = "CHILD_STATUS",
                                   status_col_father    = "P1_STATUS",
                                   status_col_mother    = "P2_STATUS",
                                   status_col_siblings  = "SIB_STATUS",
                                   number_of_siblings_col = "NUM_SIBS")


colnames(data)[ncol(data) - 1:0 ] = paste(colnames(data)[ncol(data) - 1:0 ], "_fast", sep = "")



#  LT-FH++  ---------------------------------------------------------------

#assigning thresholds for each family
thr = lapply(1:nrow(phen), function(n) {
  res = list()
  ids   = phen$ids[[n]]
  cases = phen$status[[n]] == 1
  lower = rep(-Inf, length(ids))
  upper = rep( Inf, length(ids))
  
  lower[cases] <- upper[cases] <- LTFHPlus::age_to_thres(age = phen$aoo[[n]][cases], pop_prev = prev[phen$sex[[n]][cases]])
  
  upper[!cases] <- LTFHPlus::age_to_thres(age = phen$age[[n]][!cases], pop_prev = prev[phen$sex[[n]][!cases]])
  res$ids = ids
  res$lower = lower
  res$upper = upper
  res
}) %>% 
  lapply(., as_tibble) %>% 
  bind_rows()

#performs LT-FH++ analysis 
simu_liab = estimate_gen_liability(h2 = h2, 
                       phen = phen, 
                       thr = thr,  
                       id_col = "ids")


#adding ids from the offspring
simu_liab$ids = as.character(1:N)
#combining data based on offspring's id
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
                    maximum_siblings_to_compute = nsib) %>% as_tibble()
ltfh[-1] = apply(ltfh[-1], MARGIN = 2, FUN = as.numeric) %>% as_tibble()

#combining data
simu_liab = left_join(data, as.data.frame(ltfh), by = c("ids" = "IID"))
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

