# Convert cumulative incidence rate to age

`convert_cir_to_age` computes the age from a person's cumulative
incidence rate.

## Usage

``` r
convert_cir_to_age(cir, pop_prev = 0.1, mid_point = 60, slope = 1/8)
```

## Arguments

- cir:

  A positive number representing the individual's cumulative incidence
  rate.

- pop_prev:

  A positive number representing the overall population prevalence. Must
  be at most 1 and must be larger than cir. Defaults to 0.1.

- mid_point:

  A positive number representing the mid point logistic function.
  Defaults to 60.

- slope:

  A number holding the rate of increase. Defaults to 1/8.

## Value

If cir and mid_point are positive numbers, if pop_prev is a positive
number between 0 and 1 and if slope is a valid number, then
`convert_cir_to_age` returns a number, which is equal to the current
age.

## Details

Given a person's cumulative incidence rate (cir), `convert_cir_to_age`
can be used to compute the corresponding age, which is given by
\$\$mid\\ point - \log(pop\\ prev/cir - 1) \* 1/slope\$\$

## Examples

``` r
curve(sapply(cir, convert_cir_to_age), from = 0.001, to = 0.099, xname = "cir") 
```
