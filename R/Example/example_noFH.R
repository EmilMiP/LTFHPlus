library(future.apply)
library(dplyr)
library(LTFHPlus)
library(ggplot2)



N = 10000 
h2 = .5 

multiplier = 1
prev = c(0.02, 0.08) * multiplier

cov_mat = matrix(h2, ncol = 2, nrow = 2)
cov_mat[2,2] = 1

#simulate liabilities & make tibble
liabs = MASS::mvrnorm(n = N, mu = rep(0, 2), Sigma = cov_mat) %>%
  as_tibble() %>%
  rename("child_gen"  = "V1",
         "child_full" = "V2")

phen = tibble(
  ids = 1:N,
  child_gen  = liabs$child_gen,
  child_full = liabs$child_full,
  age  = runif(N, 10, 60),
  sex        = sample(1:2, size = N, replace = T), # 1 and 2 to avoid "+1" through out the example
  status     = (child_full > qnorm(prev[sex], lower.tail = F)) + 0L,
  lower      = -Inf,
  upper      = Inf
) 

#who are the cases=
cases = which(phen$status == 1) 

phen$aoo = NA
phen$aoo[cases] = LTFHPlus::liab_to_aoo(liab = phen$child_full[cases], pop_prev = prev[phen$sex[cases]])

phen$upper[-cases] = LTFHPlus::age_to_thres(age = phen$age[-cases], pop_prev = prev[phen$sex[-cases]])
phen$lower[cases] = LTFHPlus::age_to_thres(age = phen$aoo[cases], pop_prev = prev[phen$sex[cases]])


phen = select(phen, ids, status, lower, upper, child_gen)


plan(multisession)

ltfhpp_nofh = estimate_gen_liability_noFH(phen = phen, h2 = h2)

#hist(ltfhpp_nofh$post_gen_no_fam)

with(ltfhpp_nofh, c("cc" = cor(child_gen, status), "ltfh_noFH" = cor(child_gen, post_gen_no_fam)))


ggplot(ltfhpp_nofh, aes(x = post_gen_no_fam, y = child_gen, color = status)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0) + 
  labs(color = "Case-Control Status") +
  xlab("Estimated Genetic Liability (LTFH++) ") +
  ylab("True Genetic Liability") + 
  ggtitle("True vs Estimated Genetic Liability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")


ltfhpp_nofh_ds = group_by(ltfhpp_nofh, status) %>%
  sample_n(length(cases)) %>%
  ungroup()
ltfhpp_nofh_ds %>% count(status)  

with(ltfhpp_nofh_ds, c("cc" = cor(child_gen, status), "ltfh_noFH" = cor(child_gen, post_gen_no_fam)))
