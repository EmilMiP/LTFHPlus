################################################################################

test_that("convert_age_to_cir works with defaults", {
  test <- convert_age_to_cir(age = 40)
  expect_equal(test, 0.00758582)
})

################################################################################

test_that("convert_age_to_cir works with boundary values", {
  test <- convert_age_to_cir(age = 0, pop_prev = 1)
  expect_equal(test, 0.00055277)
})

################################################################################

test_that("convert_age_to_cir works with negative slope", {
  test <- convert_age_to_cir(age = 40, slope = -1/12)
  expect_equal(test, 0.0841131)
})

################################################################################

test_that("convert_age_to_cir rejects negative arguments", {
  expect_error(convert_age_to_cir(age = -1), "age must be non-negative")
  expect_error(convert_age_to_cir(age = 0, pop_prev = -1), "pop_prev must be positive")
  expect_error(convert_age_to_cir(age = 0, mid_point = 0), "mid_point must be positive")
})

################################################################################

test_that("convert_age_to_cir warns when age is big", {
  expect_warning(convert_age_to_cir(age = 150), "it is unrealistic to be of age 150 or older")
})

################################################################################