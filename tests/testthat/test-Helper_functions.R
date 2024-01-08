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

