testing_tolerance <- 1e-4

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

# Negative Tests
test_that("convert_age_to_cir rejects negative age", {
  expect_error(convert_age_to_cir(-1), "age must be non-negative")
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
  expected_value <- 2.428256
  test <- convert_age_to_thresh(age)
  expect_equal(test, expected_value, tolerance = testing_tolerance)
})

test_that("convert_age_to_thresh defaults to logistic with invalid dist", {
  age <- 40
  expected_value <- 2.428256
  test <- convert_age_to_thresh(age, dist = "invalid")
  expect_equal(test, expected_value, tolerance = testing_tolerance)
})

test_that("convert_age_to_thresh uses first dist", {
  age <- 40
  expected_value <- 2.128167
  test <- convert_age_to_thresh(age, dist = c("invalid", "normal", "dist"))
  expect_equal(test, expected_value, tolerance = testing_tolerance)
})

test_that("convert_age_to_thresh works with default parameters in normal mode", {
  age <- 40
  expected_value <- 2.128167
  test <- convert_age_to_thresh(age, dist = "normal")
  expect_equal(test, expected_value, tolerance = testing_tolerance)
})

test_that("convert_age_to_thresh works in normal mode with min and max age swapped", {
  age <- 40
  min_age <- 90
  max_age <- 10
  expected_value <- 2.128167
  test <- convert_age_to_thresh(age, dist = "normal", min_age = min_age, max_age = max_age)
  expect_equal(test, expected_value, tolerance = testing_tolerance)
})

test_that("convert_age_to_thresh works in normal mode with lower and upper swapped", {
  age <- 40
  lower <- Inf
  upper <- 1.645
  expected_value <- 2.128167
  test <- convert_age_to_thresh(age, dist = "normal", lower = lower, upper = upper)
  expect_equal(test, expected_value, tolerance = testing_tolerance)
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
  expected_value_logistic <- 1.980199
  expected_value_normal <- 2.318249
  test_logistic <- convert_age_to_thresh(age, dist = "logistic", pop_prev = pop_prev, mid_point = mid_point, slope = slope)
  test_normal <- convert_age_to_thresh(age, dist = "normal", min_age = min_age, max_age = max_age, lower = lower, upper = upper)
  expect_equal(test_logistic, expected_value_logistic, tolerance = testing_tolerance)
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

test_that("convert_age_to_thresh rejects non-numeric arguments", {
  expect_error(convert_age_to_thresh("age"), "age must be numeric")
  expect_error(convert_age_to_thresh(30, dist = 123), "dist must be a string")
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
  cir <- 0.05
  expected_age <- 60 - log(0.1/cir - 1) * 1/(1/8)
  test <- convert_cir_to_age(cir)
  expect_equal(test, expected_age)
})

test_that("convert_cir_to_age works with non-default parameters", {
  cir <- 0.05
  pop_prev <- 0.2
  mid_point <- 70
  slope <- 1/10
  expected_age <- mid_point - log(pop_prev/cir - 1) * 1/slope
  test <- convert_cir_to_age(cir, pop_prev, mid_point, slope)
  expect_equal(test, expected_age)
})

test_that("convert_cir_to_age works with boundary values", {
  cir <- 0.0001
  pop_prev <- 1
  mid_point <- 0.0001
  slope <- -1/8
  expected_age <- mid_point - log(pop_prev/cir - 1) * 1/slope
  test <- convert_cir_to_age(cir, pop_prev, mid_point, slope)
  expect_equal(test, expected_age)
})

test_that("convert_cir_to_age returns 0 instead of negative value", {
  cir <- 0.05
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
  liability <- 0.05
  cir <- liability  # Assuming liability represents cir in this context
  expected_age <- 60 - log(0.1/cir - 1) * 1/(1/8)
  test <- convert_liability_to_aoo(liability)
  expect_equal(test, expected_age)
})

test_that("convert_liability_to_aoo works with default parameters in normal mode", {
  liability <- 0.05
  expected_age <- (1 - truncated_normal_cdf(liability, lower=1.645, upper=Inf)) * 90 + 10
  test <- convert_liability_to_aoo(liability, dist = "normal")
  expect_equal(test, expected_age)
})

test_that("convert_liability_to_aoo works with non-default parameters", {
  liability <- 0.05
  pop_prev <- 0.2
  mid_point <- 70
  slope <- 1/10
  min_aoo <- 20
  max_aoo <- 80
  lower <- 1.96
  upper <- 3
  
  expected_age_logistic <- mid_point - log(pop_prev/liability - 1) * 1/slope
  expected_age_normal <- (1 - truncated_normal_cdf(liability, lower, upper)) * (max_aoo - min_aoo) + min_aoo
  
  test_logistic <- convert_liability_to_aoo(liability, dist = "logistic", pop_prev = pop_prev, mid_point = mid_point, slope = slope)
  test_normal <- convert_liability_to_aoo(liability, dist = "normal", min_aoo = min_aoo, max_aoo = max_aoo, lower = lower, upper = upper)
  expect_equal(test_logistic, expected_age_logistic)
  expect_equal(test_normal, expected_age_normal)
})

# Negative Tests
test_that("convert_liability_to_aoo rejects invalid dist value", {
  expect_error(convert_liability_to_aoo(0.05, dist = "invalid_dist"), "dist must be either 'logistic' or 'normal'")
})

test_that("convert_liability_to_aoo rejects invalid logistic parameters", {
  expect_error(convert_liability_to_aoo(0.05, dist = "logistic", pop_prev = -0.1), "pop_prev must be positive and at most 1")
  expect_error(convert_liability_to_aoo(0.05, dist = "logistic", mid_point = -60), "mid_point must be positive")
  expect_error(convert_liability_to_aoo(0.05, dist = "logistic", slope = -1/8), "slope must be a valid number")
})

test_that("convert_liability_to_aoo rejects invalid normal parameters", {
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", min_aoo = -10), "min_aoo must be positive")
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", max_aoo = 5), "max_aoo must be greater than min_aoo")
  expect_error(convert_liability_to_aoo(0.05, dist = "normal", lower = 2, upper = 1), "upper must be greater or equal to lower")
})

test_that("convert_liability_to_aoo rejects non-numeric arguments", {
  expect_error(convert_liability_to_aoo("liability"), "liability must be numeric")
  expect_error(convert_liability_to_aoo(0.05, dist = 123), "dist must be a string")
  expect_error(convert_liability_to_aoo(0.05, pop_prev = "pop_prev"), "pop_prev must be numeric")
  expect_error(convert_liability_to_aoo(0.05, min_aoo = "min_aoo"), "min_aoo must be numeric")
})