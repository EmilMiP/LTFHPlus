# Convert liability to age of onset

`convert_liability_to_aoo` computes the age of onset from an
individual's true underlying liability using either the logistic
function or the truncated normal distribution.

## Usage

``` r
convert_liability_to_aoo(
  liability,
  dist = "logistic",
  pop_prev = 0.1,
  mid_point = 60,
  slope = 1/8,
  min_aoo = 10,
  max_aoo = 90,
  lower = stats::qnorm(0.05, lower.tail = FALSE),
  upper = Inf
)
```

## Arguments

- liability:

  A number representing the individual's true underlying liability.

- dist:

  A string indicating which distribution to use. If dist = "logistic",
  the logistic function will be used to compute the age of onset. If
  dist = "normal", the truncated normal distribution will be used
  instead. Defaults to "logistic".

- pop_prev:

  Only necessary if dist = "logistic". A positive number representing
  the overall population prevalence. Must be at most 1. Defaults to 0.1.

- mid_point:

  Only necessary if dist = "logistic". A positive number representing
  the mid point logistic function. Defaults to 60.

- slope:

  Only necessary if dist = "logistic". A number holding the rate of
  increase. Defaults to 1/8.

- min_aoo:

  Only necessary if dist = "normal". A positive number representing the
  individual's earliest age of onset. Defaults to 10.

- max_aoo:

  Only necessary if dist = "normal". A positive number representing the
  individual's latest age of onset. Must be greater than min_aoo.
  Defaults to 90.

- lower:

  Only necessary if dist = "normal". A number representing the lower
  cutoff point for the truncated normal distribution. Defaults to 1.645
  (stats::qnorm(0.05, lower.tail = FALSE)).

- upper:

  Only necessary if dist = "normal". A number representing the upper
  cutoff point of the truncated normal distribution. Must be greater or
  equal to lower. Defaults to Inf.

## Value

If liability is a number and all other necessary arguments are valid,
then `convert_liability_to_aoo` returns a positive number, which is
equal to the age of onset.

## Details

Given a person's cumulative incidence rate (cir),
`convert_liability_to_aoo` can be used to compute the corresponding age.
Under the logistic function, the age is given by \$\$mid\\ point -
log(pop\\ prev/cir - 1) \* 1/slope\$\$, while it is given by \$\$(1 -
truncated\\ normal\\ cdf(liability = liability, lower = lower , upper =
upper)) \* max\\ aoo + min\\ aoo\$\$ under the truncated normal
distribution.

## Examples

``` r
curve(sapply(liability, convert_liability_to_aoo), from = 1.3, to = 3.5, xname = "liability") 

curve(sapply(liability, convert_liability_to_aoo, dist = "normal"),
 from = qnorm(0.05, lower.tail = FALSE), to = 3.5, xname = "liability") 

```
