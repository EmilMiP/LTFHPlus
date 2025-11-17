# Relatedness between a pair of family members

`get_relatedness` returns the relatedness times the liability-scale
heritability for a pair of family members

## Usage

``` r
get_relatedness(s1, s2, h2 = 0.5, from_covmat = FALSE)
```

## Arguments

- s1, s2:

  Strings representing the two family members. The strings must be
  chosen from the following list of strings:

  - `g` (Genetic component of full liability)

  - `o` (Full liability)

  - `m` (Mother)

  - `f` (Father)

  - `c[0-9]*.[0-9]*` (Children)

  - `mgm` (Maternal grandmother)

  - `mgf` (Maternal grandfather)

  - `pgm` (Paternal grandmother)

  - `pgf` (Paternal grandfather)

  - `s[0-9]*` (Full siblings)

  - `mhs[0-9]*` (Half-siblings - maternal side)

  - `phs[0-9]*` (Half-siblings - paternal side)

  - `mau[0-9]*` (Aunts/Uncles - maternal side)

  - `pau[0-9]*` (Aunts/Uncles - paternal side).

- h2:

  A number representing the squared heritability on liability scale.
  Must be non-negative and at most 1. Defaults to 0.5

- from_covmat:

  logical variable. Only used internally. allows for skip of negative
  check.

## Value

If both `s1` and `s2` are strings chosen from the mentioned list of
strings and `h2` is a number satisfying \\0 \leq h2 \leq 1\\, then the
output will be a number that equals the percentage of shared DNA between
`s1` and `s2` times the squared heritability `h2`.

## Details

This function can be used to get the percentage of shared DNA times the
liability-scale heritability \\h^2\\ for two family members.

## Note

If you are only interested in the percentage of shared DNA, set
`h2 = 1`.

## Examples

``` r
get_relatedness("g","o")
#> [1] 0.5
get_relatedness("g","f", h2 = 1)
#> [1] 0.5
get_relatedness("o","s", h2 = 0.3)
#> [1] 0.15


# This will result in errors:
try(get_relatedness("a","b"))
#> Error in validate_relatives(s1) : 
#>   s1 contains invalid abbreviations! Use strings from the following list: 
#> 
#>   - g (Genetic component of full liability)
#> 
#>   - o (Full liability)
#> 
#>   - m (Mother)
#> 
#>   - f (Father)
#> 
#>   - c[0-9]*.[0-9]* (Children)
#> 
#>   - mgm (Maternal grandmother)
#> 
#>   - mgf (Maternal grandfather)
#> 
#>   - pgm (Paternal grandmother)
#> 
#>   - pgf (Paternal grandfather)
#> 
#>   - s[0-9]* (Full siblings)
#> 
#>   - mhs[0-9]* (Half-siblings - maternal side)
#> 
#>   - phs[0-9]* (Half-siblings - paternal side)
#> 
#>   - mau[0-9]* (Aunts/Uncles - maternal side)
#> 
#>   - pau[0-9]* (Aunts/Uncles - paternal side).
try(get_relatedness(m, mhs))
#> Error in eval(expr, envir) : object 'm' not found
```
