################################################################################
# Testing truncated_normal_cdf
################################################################################

# Positive Tests
test_that("truncated_normal_cdf works with default parameters", {
  liability <- 0.5
  test <- truncated_normal_cdf(liability)
  expected <- -0.2585375387
  expect_equal(test, expected)
})

test_that("truncated_normal_cdf works with non-default parameters", {
  liability <- 0.5
  lower <- 1.0
  upper <- 2.0
  test <- truncated_normal_cdf(liability, lower, upper)
  expected <- -0.1498822848
  expect_equal(test, expected)
})

test_that("truncated_normal_cdf works with swapped lower and upper values", {
  liability <- 0.5
  test <- truncated_normal_cdf(0.5, lower = Inf, upper = 1.644853627)
  expected <- -0.2585375387
  expect_equal(test, expected)
})

# Negative Tests

test_that("truncated_normal_cdf rejects non-numeric arguments", {
  expect_error(truncated_normal_cdf("liability"), "liability must be numeric")
  expect_error(truncated_normal_cdf(0.5, lower = "lower"), "lower cutoff point must be numeric")
  expect_error(truncated_normal_cdf(0.5, upper = "upper"), "upper cutoff point must be numeric")
})


################################################################################
# Testing convert_age_to_cir
################################################################################

# Positive Tests
test_that("convert_age_to_cir works with default parameters", {
  age <- 40
  expected_cir <- 0.007585818
  test <- convert_age_to_cir(age)
expect_equal(test, expected_cir)
})

test_that("convert_age_to_cir works with non-default parameters", {
  age <- 50
  pop_prev <- 0.2
  mid_point <- 70
  slope <- 1/10
  expected_cir <- 0.02384058
  test <- convert_age_to_cir(age, pop_prev, mid_point, slope)
  expect_equal(test, expected_cir)
})

test_that("convert_age_to_cir works with boundary values", {
  age_boundary <- 0
  pop_prev_boundary <- 1
  mid_point_boundary <- 1
  slope_boundary <- 0
  expected_cir_boundary <- 0.5
  test <- convert_age_to_cir(age_boundary, pop_prev_boundary, mid_point_boundary, slope_boundary)
  expect_equal(test, expected_cir_boundary)
})

test_that("convert_age_to_cir works with negative slope", {
  test <- convert_age_to_cir(age = 40, slope = -1/12)
  expect_equal(test, 0.0841131)
})

# Negative Tests
test_that("convert_age_to_cir rejects negative age", {
  expect_error(convert_age_to_cir(-1), "age must be non-negative")
})

test_that("convert_age_to_cir warns of unrealistic age", {
  expect_warning(convert_age_to_cir(150), "it is unrealistic to be of age 150 or older")
})

test_that("convert_age_to_cir rejects invalid pop_prev", {
  expect_error(convert_age_to_cir(30, pop_prev = -0.1), "pop_prev must be positive")
  expect_error(convert_age_to_cir(30, pop_prev = 1.1), "pop_prev must be smaller or equal to 1")
})

test_that("convert_age_to_cir rejects non-positive mid_point", {
  expect_error(convert_age_to_cir(30, mid_point = 0), "mid_point must be positive")
})

test_that("convert_age_to_cir rejects non-numeric arguments", {
  expect_error(convert_age_to_cir("age"), "age must be numeric")
  expect_error(convert_age_to_cir(30, pop_prev = "pop_prev"), "pop_prev must be numeric")
  expect_error(convert_age_to_cir(30, mid_point = "mid_point"), "mid_point must be numeric")
  expect_error(convert_age_to_cir(30, slope = "slope"), "slope must be numeric")
})

################################################################################
# Testing convert_age_to_thresh
################################################################################

# Positive Tests
test_that("convert_age_to_thresh works with default parameters in logistic mode", {
  age <- 40
  expected_value <- 2.428255857
  test <- convert_age_to_thresh(age)
  expect_equal(test, expected_value)
})

test_that("convert_age_to_thresh defaults to logistic with invalid dist", {
  age <- 40
  expected_value <- 2.428255857
  test <- convert_age_to_thresh(age, dist = "invalid")
  expect_equal(test, expected_value)
})

test_that("convert_age_to_thresh uses first valid dist", {
  age <- 40
  expected_value <- 2.1280452341849840359
  test <- convert_age_to_thresh(age, dist = c("invalid", "normal", "logistic"))
  expect_equal(test, expected_value)
})

test_that("convert_age_to_thresh works with default parameters in normal mode", {
  age <- 40
  expected_value <- 2.1280452341849840359
  test <- convert_age_to_thresh(age, dist = "normal")
  expect_equal(test, expected_value)
})

test_that("convert_age_to_thresh works in normal mode with min and max age swapped", {
  age <- 40
  min_age <- 90
  max_age <- 10
  expected_value <- 2.1280452341849840359
  test <- convert_age_to_thresh(age, dist = "normal", min_age = min_age, max_age = max_age)
  expect_equal(test, expected_value)
})

test_that("convert_age_to_thresh works in normal mode with lower and upper swapped", {
  age <- 40
  lower <- Inf
  upper <- 1.644853626951472636
  expected_value <- 2.1280452341849840359
  test <- convert_age_to_thresh(age, dist = "normal", lower = lower, upper = upper)
  expect_equal(test, expected_value)
})

test_that("convert_age_to_thresh works with non-default parameters", {
  age <- 50
  pop_prev <- 0.2
  mid_point <- 70
  slope <- 1/10
  min_age <- 20
  max_age <- 80
  lower <- 1.96
  upper <- 3
  expected_value_logistic <- 1.98019903
  expected_value_normal <- 2.318248972
  test_logistic <- convert_age_to_thresh(age, dist = "logistic", pop_prev = pop_prev, mid_point = mid_point, slope = slope)
  test_normal <- convert_age_to_thresh(age, dist = "normal", min_age = min_age, max_age = max_age, lower = lower, upper = upper)
  expect_equal(test_logistic, expected_value_logistic)
  expect_equal(test_normal, expected_value_normal)
})

# Negative Tests
test_that("convert_age_to_thresh rejects non-positive age", {
  expect_error(convert_age_to_thresh(-1), "age must be non-negative")
})

test_that("convert_age_to_thresh rejects invalid logistic parameters", {
  expect_error(convert_age_to_thresh(30, dist = "logistic", pop_prev = -0.1), "pop_prev must be positive")
  expect_error(convert_age_to_thresh(30, dist = "logistic", pop_prev = 1.1), "pop_prev must be smaller or equal to 1")
  expect_error(convert_age_to_thresh(30, dist = "logistic", mid_point = -60), "mid_point must be positive")
  
})

test_that("convert_age_to_thresh rejects invalid normal parameters", {
  expect_error(convert_age_to_thresh(30, dist = "normal", min_age = -10), "min_age must be positive")
  expect_error(convert_age_to_thresh(30, dist = "normal", max_age = -10), "max_age must be positive")
})

test_that("convert_age_to_thresh rejects invalid dist", {
  expect_error(convert_age_to_thresh(30, dist = 123), "dist must be a string")
})

test_that("convert_age_to_thresh rejects non-numeric arguments", {
  expect_error(convert_age_to_thresh("age"), "age must be numeric")
  expect_error(convert_age_to_thresh(30, pop_prev = "pop_prev"), "pop_prev must be numeric")
  expect_error(convert_age_to_thresh(30, dist = "normal", min_age = "min_age"), "min_age must be numeric")
  expect_error(convert_age_to_thresh(30, dist = "normal", max_age = "max_age"), "max_age must be numeric")
  expect_error(convert_age_to_thresh(30, dist = "normal", lower = "lower"), "lower cutoff point must be numeric")
  expect_error(convert_age_to_thresh(30, dist = "normal", upper = "upper"), "upper cutoff point must be numeric")
  expect_error(convert_age_to_thresh(30, dist = "logistic", slope = "invalid_slope"), "slope must be numeric")
  expect_error(convert_age_to_thresh(30, dist = "logistic", mid_point = "invalid_mid_point"), "mid_point must be numeric")
  expect_error(convert_age_to_thresh(30, dist = "logistic", pop_prev = "invalid_pop_prev"), "pop_prev must be numeric")
  })

################################################################################
# Testing convert_cir_to_age
################################################################################

# Positive Tests
test_that("convert_cir_to_age works with default parameters", {
  cir <- 0.08
  expected_age <- 71.09035489
  test <- convert_cir_to_age(cir)
  expect_equal(test, expected_age)
})

test_that("convert_cir_to_age works with non-default parameters", {
  cir <- 0.08
  pop_prev <- 0.2
  mid_point <- 70
  slope <- 1/10
  expected_age <- 65.94534892
  test <- convert_cir_to_age(cir, pop_prev, mid_point, slope)
  expect_equal(test, expected_age)
})

test_that("convert_cir_to_age works with boundary values", {
  cir <- 0.0001
  pop_prev <- 1
  mid_point <- 0.0001
  slope <- -1/8
  expected_age <- 73.68202294
  test <- convert_cir_to_age(cir, pop_prev, mid_point, slope)
  expect_equal(test, expected_age)
})

test_that("convert_cir_to_age returns 0 instead of negative value", {
  cir <- 0.1
  mid_point <- 1
  pop_prev <- 1
  test <- convert_cir_to_age(cir, pop_prev, mid_point)
  expect_equal(test, 0)
})

# Negative Tests
test_that("convert_cir_to_age rejects non-positive cir", {
  expect_error(convert_cir_to_age(0), "cir must be positive")
  expect_error(convert_cir_to_age(-1), "cir must be positive")
})

test_that("convert_cir_to_age rejects invalid pop_prev", {
  expect_error(convert_cir_to_age(0.05, pop_prev = 0), "pop_prev must be positive")
  expect_error(convert_cir_to_age(0.05, pop_prev = -0.1), "pop_prev must be positive")
  expect_error(convert_cir_to_age(0.05, pop_prev = 1.1), "pop_prev must be smaller or equal to 1")
})

test_that("convert_cir_to_age returns NA when cir is greater than pop_prev", {
  expect_equal(convert_cir_to_age(0.5, pop_prev = 0.1), NA)
})

test_that("convert_cir_to_age rejects non-positive mid_point", {
  expect_error(convert_cir_to_age(0.05, mid_point = 0), "mid_point must be positive")
  expect_error(convert_cir_to_age(0.05, mid_point = -60), "mid_point must be positive")
})

test_that("convert_cir_to_age rejects non-numeric arguments", {
  expect_error(convert_cir_to_age("a"), "cir must be numeric")
  expect_error(convert_cir_to_age(0.05, pop_prev = "a"), "pop_prev must be numeric")
  expect_error(convert_cir_to_age(0.05, mid_point = "a"), "mid_point must be numeric")
  expect_error(convert_cir_to_age(0.05, slope = "a"), "slope must be numeric")
})

################################################################################
# Testing convert_liability_to_aoo
################################################################################

# Positive Tests
test_that("convert_liability_to_aoo works with default parameters in logistic mode", {
  liability <- 2
  expected_age <- 50.2202059
  test <- convert_liability_to_aoo(liability)
  expect_equal(test, expected_age)
})

test_that("convert_liability_to_aoo works with default parameters in normal mode", {
  liability <- 2
  expected_age <- 97.547511875336127218
  test <- convert_liability_to_aoo(liability, dist = "normal")
  expect_equal(test, expected_age)
})

test_that("convert_liability_to_aoo works with non-default parameters", {
  liability <- 2
  pop_prev <- 0.2
  mid_point <- 70
  slope <- 1/10
  min_aoo <- 20
  max_aoo <- 80
  lower <- 1.96
  upper <- 3
  
  expected_age_logistic <- 49.47010524
  expected_age_normal <- 99.82017894
  
  test_logistic <- convert_liability_to_aoo(liability, dist = "logistic", pop_prev = pop_prev, mid_point = mid_point, slope = slope)
  test_normal <- convert_liability_to_aoo(liability, dist = "normal", min_aoo = min_aoo, max_aoo = max_aoo, lower = lower, upper = upper)
  expect_equal(test_logistic, expected_age_logistic)
  expect_equal(test_normal, expected_age_normal)
})

test_that("convert_liability_to_aoo works with boundary values in logistic mode", {
  liability <- 2
  pop_prev <- 1
  expected_age <- 29.91862861
  test <- convert_liability_to_aoo(liability, pop_prev = pop_prev)
  expect_equal(test, expected_age)
})

test_that("convert_liability_to_aoo defaults to logistic with invalid dist", {
  liability <- 2
  expected_age <- 50.2202059
  test <- convert_liability_to_aoo(liability, dist = "invalid")
  expect_equal(test, expected_age)
})

test_that("convert_liability_to_aoo uses first valid dist", {
  liability <- 2
  expected_age <- 97.54751188
  test <- convert_liability_to_aoo(liability, dist = c("invalid", "normal", "logistic"))
  expect_equal(test, expected_age)
})

test_that("convert_liability_to_aoo can swap lower and upper in normal mode", {
  liability <- 2
  expected_age <- 97.547511875336127218
  test <- convert_liability_to_aoo(liability, dist = "normal", lower = Inf, upper = 1.644853626951472636)
  expect_equal(test, expected_age)
})

test_that("convert_liability_to_aoo can swap min_aoo and max_aoo in normal mode", {
  liability <- 2
  expected_age <- 97.547511875336127218
  test <- convert_liability_to_aoo(liability, dist = "normal", min_aoo = 90, max_aoo = 10)
  expect_equal(test, expected_age)
})

test_that("convert_liability_to_aoo returns 0 instead of negative values in logistic mode", {
  liability <- 10
  expected_age <- 0
  test <- convert_liability_to_aoo(liability)
  expect_equal(test, expected_age)
})

test_that("convert_liability_to_aoo returns 0 instead of negative values in normal mode", {
  liability <- 10
  expected_age <- 0
  test <- convert_liability_to_aoo(liability)
  expect_equal(test, expected_age)
})

# Negative Tests
test_that("convert_liability_to_aoo defaults to logistic when given invalid dist value", {
  expect_equal(convert_liability_to_aoo(2, dist = "invalid_dist"), 50.2202059)
})

test_that("convert_liability_to_aoo returns NA when liability is too low in logistic mode", {
  liability <- 1
  test <- convert_liability_to_aoo(liability)
  expect_equal(test, NA)
})


test_that("convert_liability_to_aoo rejects invalid logistic parameters", {
  expect_error(convert_liability_to_aoo(0.05, dist = "logistic", pop_prev = -0.1), "pop_prev must be positive")
  expect_error(convert_liability_to_aoo(0.05, dist = "logistic", pop_prev = 0), "pop_prev must be positive")
  expect_error(convert_liability_to_aoo(0.05, dist = "logistic", pop_prev = 1.1), "pop_prev must be smaller or equal to 1")
  expect_error(convert_liability_to_aoo(0.05, dist = "logistic", mid_point = -60), "mid_point must be positive")
  expect_error(convert_liability_to_aoo(0.05, dist = "logistic", mid_point = 0), "mid_point must be positive")
})

test_that("convert_liability_to_aoo rejects invalid normal parameters", {
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", min_aoo = -10), "min_aoo must be positive")
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", min_aoo = 0), "min_aoo must be positive")
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", max_aoo = -10), "max_aoo must be positive")
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", max_aoo = 0), "max_aoo must be positive")
})

test_that("convert_liability_to_aoo rejects non-numeric arguments", {
  expect_error(convert_liability_to_aoo("liability"), "liability must be numeric")
  expect_error(convert_liability_to_aoo(0.05, dist = 123), "dist must be a string")
  expect_error(convert_liability_to_aoo(0.05, mid_point = "mid_point"), "mid_point must be numeric")
  expect_error(convert_liability_to_aoo(0.05, slope = "slope"), "slope must be numeric")
  expect_error(convert_liability_to_aoo(0.05, pop_prev = "pop_prev"), "pop_prev must be numeric")
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", min_aoo = "min_aoo"), "min_aoo must be numeric")
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", max_aoo = "max_aoo"), "max_aoo must be numeric")
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", lower = "lower"), "lower cutoff point must be numeric")
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", upper = "upper"), "upper cutoff point must be numeric")
})

################################################################################
# Testing convert_observed_to_liability_scale
################################################################################

# Positive Tests
test_that("convert_observed_to_liability_scale works with default parameters", {
  test <- convert_observed_to_liability_scale()
  expected <- 0.4242283384
  expect_equal(test, expected)
})

test_that("convert_observed_to_liability_scale works with single non-default parameters", {
  obs_h2 <- 0.6
  pop_prev <- 0.1
  prop_cases <- 0.6
  test <- convert_observed_to_liability_scale(obs_h2, pop_prev, prop_cases)
  expected <- 0.657474694
  expect_equal(test, expected)
})

test_that("convert_observed_to_liability_scale works with vector parameters", {
  obs_h2 <- c(0.6, 0.7)
  pop_prev <- c(0.1, 0.2)
  prop_cases <- c(0.6, 0.7)
  test <- convert_observed_to_liability_scale(obs_h2, pop_prev, prop_cases)
  expected <- c(0.657474694, 1.088731486)
  expect_equal(test, expected)
})

test_that("convert_observed_to_liability_scale works with NULL prop_cases", {
  obs_h2 <- 0.6
  pop_prev <- 0.1
  test <- convert_observed_to_liability_scale(obs_h2, pop_prev, NULL)
  expected <- 1.753265851
  expect_equal(test, expected)
})

# Negative Tests
test_that("convert_observed_to_liability_scale rejects negative inputs", {
  expect_error(convert_observed_to_liability_scale(-0.1, 0.05, 0.5), "observed heritability\\(ies\\) must be non-negative")
  expect_error(convert_observed_to_liability_scale(0.5, -0.1, 0.5), "population prevalence\\(s\\) must be non-negative")
  expect_error(convert_observed_to_liability_scale(0.5, 0.05, -0.1), "proportion\\(s\\) of cases must be non-negative")
})

test_that("convert_observed_to_liability_scale rejects inputs larger than 1", {
  expect_error(convert_observed_to_liability_scale(1.1, 0.05, 0.8), "observed heritability\\(ies\\) must be smaller than or equal to one")
  expect_error(convert_observed_to_liability_scale(0.5, 1.1, 0.8), "population prevalence\\(s\\) must be smaller than or equal to one")
  expect_error(convert_observed_to_liability_scale(0.5, 0.05, 1.1), "proportion\\(s\\) of cases must be smaller than or equal to one")
})

test_that("convert_observed_to_liability_scale rejects non-numeric inputs", {
  expect_error(convert_observed_to_liability_scale("obs_h2", 0.05, 0.5), "observed heritability\\(ies\\) must be numeric")
  expect_error(convert_observed_to_liability_scale(0.5, "pop_prev", 0.5), "population prevalence\\(s\\) must be numeric")
  expect_error(convert_observed_to_liability_scale(0.5, 0.05, "prop_cases"), "proportion\\(s\\) of cases must be numeric")
})