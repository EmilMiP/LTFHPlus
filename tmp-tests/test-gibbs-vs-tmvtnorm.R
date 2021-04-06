sd(replicate(100, {
  colMeans(tmvtnorm::rtmvnorm(10e3, mean = rep(0, 4), H = inv_cov,
                              lower = lims[, "lower"], upper = lims[, "upper"],
                              algorithm = "gibbs", burn.in.samples = 1000))[1]
}))
# 0.01322415 / 0.01394558

test <- rtmvnorm.gibbs(10e3, sigma = cov, lower = lims[, "lower"], upper = lims[, "upper"])
plot(test)

sd(replicate(100, {
  mean(rtmvnorm.gibbs(10e3, sigma = cov, lower = lims[, "lower"], upper = lims[, "upper"]))
}))
# 0.008360882 / 0.008005332

replicate(10, {
  rtmvnorm.gibbs(10e3, sigma = cov, lower = lims[, "lower"], upper = lims[, "upper"])
})

hist(replicate(100000, 2 + 3 * qnorm(runif(1, 0.2, 0.7))))
hist(replicate(100000, qnorm(runif(1, 0.2, 0.7), mean = 2, sd = 3)))

microbenchmark::microbenchmark(
  tmvtnorm::rtmvnorm(10e3, mean = rep(0, 4), H = inv_cov,
                     lower = lims[, "lower"], upper = lims[, "upper"],
                     algorithm = "gibbs", burn.in.samples = 1000),
  rtmvnorm.gibbs(10e3, sigma = cov, lower = lims[, "lower"], upper = lims[, "upper"]),
  times = 10
)