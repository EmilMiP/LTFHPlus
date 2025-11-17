# Gibbs Sampler for the truncated multivariate normal distribution

`rtmvnorm.gibbs` implements Gibbs sampler for the truncated multivariate
normal distribution with covariance matrix `covmat`.

## Usage

``` r
rtmvnorm.gibbs(
  n_sim = 1e+05,
  covmat,
  lower = -Inf,
  upper,
  fixed = (lower == upper),
  out = c(1),
  burn_in = 1000
)
```

## Arguments

- n_sim:

  A positive number representing the number of draws from the Gibbs
  sampler after burn-in.. Defaults to `1e+05`.

- covmat:

  A symmetric and numeric matrix representing the covariance matrix for
  the multivariate normal distribution.

- lower:

  A number or numeric vector representing the lower cutoff point(s) for
  the truncated normal distribution. The length of lower must be 1 or
  equal to the dimension of the multivariable normal distribution.
  Defaults to `-Inf`.

- upper:

  A number or numeric vector representing the upper cutoff point(s) for
  the truncated normal distribution. Must be greater or equal to lower.
  In addition the length of upper must be 1 or equal to the dimension of
  the multivariable normal distribution. Defaults to `Inf`.

- fixed:

  A logical scalar or a logical vector indicating which variables to
  fix. If `fixed` is a vector, it must have the same length as lower and
  upper. Defaults to `TRUE` when `lower` is equal to `upper` and `FALSE`
  otherwise.

- out:

  An integer or numeric vector indicating which variables should be
  returned from the Gibbs sampler. If `out = c(1)`, the first variable
  (usually the genetic component of the full liability of the first
  phenotype) is estimated and returned. If `out = c(2)`, the second
  variable (usually full liability) is estimated and returned. If
  `out = c(1,2)`, both the first and the second variable are estimated
  and returned. Defaults to `c(1)`.

- burn_in:

  A number of iterations that count as burn in for the Gibbs sampler.
  Must be non-negative. Defaults to `1000`.

## Value

If `covmat` is a symmetric and numeric matrix, if `n_sim` and `burn_in`
are positive/non-negative numbers, if `out` is a numeric vector and
`lower`, `upper` and `fixed` are numbers or vectors of the same length
and the required format, `rtmvnorm.gibbs` returns the sampling values
from the Gibbs sampler for all variables specified in `out`.

## Details

Given a covariance matrix `covmat` and lower and upper cutoff points,
the function `rtmvnorm.gibbs()` can be used to perform Gibbs sampler on
a truncated multivariable normal distribution. It is possible to specify
which variables to return from the Gibbs sampler, making it convenient
to use when estimating only the full liability or the genetic component
of the full liability.

## References

Kotecha, J. H., & Djuric, P. M. (1999, March). Gibbs sampling approach
for generation of truncated multivariate gaussian random variables. In
1999 IEEE International Conference on Acoustics, Speech, and Signal
Processing. Proceedings. ICASSP99 (Cat. No. 99CH36258) (Vol. 3, pp.
1757-1760). IEEE.
[doi:10.1109/ICASSP.1999.756335](https://doi.org/10.1109/ICASSP.1999.756335)

Wilhelm, S., & Manjunath, B. G. (2010). tmvtnorm: A package for the
truncated multivariate normal distribution. The R Journal.
[doi:10.32614/RJ-2010-005](https://doi.org/10.32614/RJ-2010-005)

## Examples

``` r
samp <- rtmvnorm.gibbs(10e3, covmat = matrix(c(1, 0.2, 0.2, 0.5), 2),
                       lower = c(-Inf, 0), upper = c(0, Inf), out = 1:2)
```
