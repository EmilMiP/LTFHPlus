# CDF for truncated normal distribution.

`truncated_normal_cdf` computes the cumulative density function for a
truncated normal distribution.

## Usage

``` r
truncated_normal_cdf(
  liability,
  lower = stats::qnorm(0.05, lower.tail = FALSE),
  upper = Inf
)
```

## Arguments

- liability:

  A number representing the individual's true underlying liability.

- lower:

  A number representing the lower cutoff point for the truncated normal
  distribution. Defaults to 1.645 (stats::qnorm(0.05, lower.tail =
  FALSE)).

- upper:

  A number representing the upper cutoff point of the truncated normal
  distribution. Must be greater or equal to lower. Defaults to Inf.

## Value

If liability is a number and the lower and upper cutoff points are
numbers satisfying lower \<= upper, then `truncated_normal_cdf` returns
the probability that the liability will take on a value less than or
equal to `liability`.

## Details

This function can be used to compute the value of the cumulative density
function for a truncated normal distribution given an individual's true
underlying liability.

## Examples

``` r
curve(sapply(liability, truncated_normal_cdf), from = qnorm(0.05, lower.tail = FALSE), to = 3.5,
 xname = "liability")

 
```
