# Convert age to cumulative incidence rate

`convert_age_to_cir` computes the cumulative incidence rate from a
person's age.

## Usage

``` r
convert_age_to_cir(age, pop_prev = 0.1, mid_point = 60, slope = 1/8)
```

## Arguments

- age:

  A non-negative number representing the individual's age.

- pop_prev:

  A positive number representing the overall population prevalence. Must
  be at most 1. Defaults to 0.1.

- mid_point:

  A positive number representing the mid point logistic function.
  Defaults to 60.

- slope:

  A number holding the rate of increase. Defaults to 1/8.

## Value

If age and mid_point are positive numbers, if pop_prev is a positive
number between 0 and 1 and if slope is a valid number, then
`convert_age_to_cir` returns a number, which is equal to the cumulative
incidence rate.

## Details

Given a person's age, `convert_age_to_cir` can be used to compute the
cumulative incidence rate (cir), which is given by the formula \$\$pop\\
prev / (1 + exp((mid\\ point - age) \* slope))\$\$

## Examples

``` r
curve(sapply(age, convert_age_to_cir), from = 10, to = 110, xname = "age")
```
