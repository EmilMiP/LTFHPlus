library(LTFHPlus)
library(ggplot2)
library(gridExtra)


#number of individuals
N = 5000
#constructing correlation matrix for LT-FH++
h2_1 = .5
h2_2 = .5
gen_cor = .5
corr_mat = diag(c(h2_1, h2_2))
corr_mat[1,2] <- corr_mat[2,1] <- gen_cor * sqrt(h2_1 * h2_2)

nthreads = 5  # number of threads to use for ltfh++
nsib = 0 #number of siblings

#start sessions to parallelize over later
plan(multisession, workers = nthreads)

#This example will rely on the LT-FH implementation in LTFHPlus. 
#For Hujoel et al's implementation, see example.R



#calculates the sex-specific prevalences used to determine status:
multiplier = 1
prev = c(0.08, .02) * multiplier


#simulate liabilities
liabs = MASS::mvrnorm(n = N, mu = rep(0, 2*(4 + nsib)), Sigma = get_full_cov(corr_mat = corr_mat, n_sib = nsib))

#Age of children, parents, and siblings
child_age = runif(N, 10, 60)
father_age  = child_age + runif(N, 20, 35)
mother_age  = child_age + runif(N, 20, 35)
children_sex = sample(1:2, size = N, replace = TRUE)
if(nsib > 0) {
  sibling_sex = replicate(nsib, sample(1:2, size = N, replace = TRUE))
  sibling_age = replicate(nsib, runif(N, -15, 15) + child_age)
}

#constructing a list of phenotypes
all_phen = list()
for (j in 1:nrow(corr_mat)) {

  tmp = lapply(1:N, function(i) {
    cur_liabs  = liabs[i, (4 + nsib) * (j-1) + 1:(4 + nsib)][-1]
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
  })

  
  all_phen[[j]] = do.call("rbind", tmp) %>% 
    as_tibble() %>%
    rename("ids" = "V1", "sex" = "V2", "status" = "V3", "age" = "V4", "aoo" = "V5")
}

#calculating threshold for every individual
all_thr = list()
for (i in seq_along(all_phen)) {
  phen = all_phen[[i]]
  #assigning thresholds for each family
  tmp = lapply(1:nrow(phen), function(n) {
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
#    return(tibble(ids = ids, lower = lower, upper = upper))
    res
  })
  #combining all values into one tibble
  all_thr[[i]] = lapply(tmp, as_tibble) %>% bind_rows()
}
#### from this point onwards, only ids in all_phen per phenotype will be used to match the thresholds


#performing the LT-FH++ analysis
multi = estimate_gen_liability_multi_trait(corr_mat = corr_mat,
                                           phen.list = all_phen,
                                           thr.list = all_thr)
plan(sequential) # removing the sessions created by multisession.

#constructing input for LT-FH
res = tibble(
  IID          = sapply(all_phen[[1]]$ids,    function(x) x[1]),
  CHILD_STATUS = sapply(all_phen[[1]]$status, function(x) x[1]),
  P1_STATUS    = sapply(all_phen[[1]]$status, function(x) x[2]),
  P2_STATUS    = sapply(all_phen[[1]]$status, function(x) x[3]),
  NUM_SIBS     = nsib,
  SIB_STATUS   = sapply(all_phen[[1]]$status, function(x) (sum(x[-(1:3)]) > 0) + 0L) 
)

#performing LT-FH article with gibbs sampling implementation
ltfh = estimate_gen_liability_ltfh(h2 = h2_1,
                                   phen = res,
                                   child_threshold = qnorm(mean(prev), lower.tail = F),
                                   parent_threshold = qnorm(mean(prev), lower.tail = F),
                                   status_col_offspring = "CHILD_STATUS",
                                   status_col_father    = "P1_STATUS",
                                   status_col_mother    = "P2_STATUS",
                                   status_col_siblings  = "SIB_STATUS",
                                   number_of_siblings_col = "NUM_SIBS")
colnames(ltfh)[7:8] = c("ltfh", "ltfh_se")


data = left_join(multi, ltfh %>% select(IID, ltfh, ltfh_se)) %>%
  left_join(res[,1:2])
data$child_gen = liabs[,1]

with(data, c( "primary"      = cov(child_gen, post_gen_liab), 
             # "secondary"    = cov(child_gen, post_gen_liab_2), 
              "LT-FH"        = cov(child_gen, ltfh), 
#              "LT-FH_nosib"  = cor(child_gen, ltfh_nosib), 
              "Case-Control" = cov(child_gen, CHILD_STATUS)))
with(data, c( "primary"      = cor(child_gen, post_gen_liab), 
            #  "secondary"    = cor(child_gen, post_gen_liab_2), 
              "LT-FH"        = cor(child_gen, ltfh), 
#             "LT-FH_nosib"  = cor(child_gen, ltfh_nosib), 
              "Case-Control" = cor(child_gen, CHILD_STATUS)))




p1 = ggplot(data, aes(x = post_gen_liab, y = child_gen, color = sapply(all_phen[[2]]$status, sum) > 0)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Any case of \n2nd Phenotype \nin Family") +
  xlab("Estimated Genetic Liability (LT-FH++) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


p2 = ggplot(data, aes(x = ltfh, y = child_gen, color = sapply(all_phen[[2]]$status, sum) > 0)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Any case of \n2nd Phenotype \nin Family") +
  xlab("Estimated Genetic Liability (LT-FH) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2)

