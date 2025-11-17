# Convert the heritability on the observed scale to that on the liability scale

`convert_observed_to_liability_scale` transforms the heritability on the
observed scale to the heritability on the liability scale.

## Usage

``` r
convert_observed_to_liability_scale(
  obs_h2 = 0.5,
  pop_prev = 0.05,
  prop_cases = 0.5
)
```

## Arguments

- obs_h2:

  A number or numeric vector representing the liability-scale
  heritability(ies)on the observed scale. Must be non-negative and at
  most 1. Defaults to 0.5

- pop_prev:

  A number or numeric vector representing the population prevalence(s).
  All entries must be non-negative and at most one. If it is a vector,
  it must have the same length as obs_h2. Defaults to 0.05.

- prop_cases:

  Either NULL or a number or a numeric vector representing the
  proportion of cases in the sample. All entries must be non-negative
  and at most one. If it is a vector, it must have the same length as
  obs_h2. Defaults to 0.5.

## Value

If `obs_h2`, `pop_prev` and `prop_cases` are non-negative numbers that
are at most one, the function returns the heritability on the liability
scale using Equation 23 from Sang Hong Lee, Naomi R. Wray, Michael E.
Goddard and Peter M. Visscher, "Estimating Missing Heritability for
Diseases from Genome-wide Association Studies", The American Journal of
Human Genetics, Volume 88, Issue 3, 2011, pp. 294-305,
[doi:10.1016/j.ajhg.2011.02.002](https://doi.org/10.1016/j.ajhg.2011.02.002)
. If `obs_h2`, `pop_prev` and `prop_cases` are non-negative numeric
vectors where all entries are at most one, the function returns a vector
of the same length as obs_h2. Each entry holds to the heritability on
the liability scale which was obtained from the corresponding entry in
obs_h2 using Equation 23. If `obs_h2` and `pop_prev` are non-negative
numbers that are at most one and `prop_cases` is `NULL`, the function
returns the heritability on the liability scale using Equation 17 from
Sang Hong Lee, Naomi R. Wray, Michael E. Goddard and Peter M. Visscher,
"Estimating Missing Heritability for Diseases from Genome-wide
Association Studies", The American Journal of Human Genetics, Volume 88,
Issue 3, 2011, pp. 294-305,
[doi:10.1016/j.ajhg.2011.02.002](https://doi.org/10.1016/j.ajhg.2011.02.002)
. If `obs_h2` and `pop_prev` are non-negative numeric vectors such that
all entries are at most one, while `prop_cases` is `NULL`,
`convert_observed_to_liability_scale` returns a vector of the same
length as obq_h2. Each entry holds to the liability-scale heritability
that was obtained from the corresponding entry in obs_h2 using Equation
17.

## Details

This function can be used to transform the heritability on the observed
scale to that on the liability scale.
`convert_observed_to_liability_scale` uses either Equation 17 (if
prop_cases = NULL) or Equation 23 from Sang Hong Lee, Naomi R. Wray,
Michael E. Goddard and Peter M. Visscher, "Estimating Missing
Heritability for Diseases from Genome-wide Association Studies", The
American Journal of Human Genetics, Volume 88, Issue 3, 2011, pp.
294-305,
[doi:10.1016/j.ajhg.2011.02.002](https://doi.org/10.1016/j.ajhg.2011.02.002)
to transform the heritability on the observed scale to the heritability
on the liability scale.

## References

Sang Hong Lee, Naomi R. Wray, Michael E. Goddard, Peter M. Visscher
(2011, March). Estimating Missing Heritability for Diseases from
Genome-wide Association Studies. In The American Journal of Human
Genetics (Vol. 88, Issue 3, pp. 294-305).
[doi:10.1016/j.ajhg.2011.02.002](https://doi.org/10.1016/j.ajhg.2011.02.002)

## Examples

``` r
convert_observed_to_liability_scale()
#> [1] 0.4242283
convert_observed_to_liability_scale(prop_cases=NULL)
#> [1] 2.232781
convert_observed_to_liability_scale(obs_h2 = 0.8, pop_prev = 1/44, 
                                    prop_cases = NULL)
#> [1] 6.105859
convert_observed_to_liability_scale(obs_h2 = c(0.5,0.8), 
                                    pop_prev = c(0.05, 1/44), 
                                    prop_cases = NULL)
#> [1] 2.232781 6.105859
```
