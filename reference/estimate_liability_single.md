# Estimating the genetic or full liability

`estimate_liability_single` estimates the genetic component of the full
liability and/or the full liability for a number of individuals based on
their family history.

## Usage

``` r
estimate_liability_single(
  .tbl = NULL,
  family_graphs = NULL,
  h2 = 0.5,
  pid = "PID",
  fam_id = "fam_ID",
  family_graphs_col = "fam_graph",
  role = NULL,
  out = c(1),
  tol = 0.01
)
```

## Arguments

- .tbl:

  A matrix, list or data frame that can be converted into a tibble. Must
  have at least five columns that hold the family identifier, the
  personal identifier, the role and the lower and upper thresholds. Note
  that the role must be one of the following abbreviations

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

  - `pau[0-9]*` (Aunts/Uncles - paternal side). Defaults to `NULL`.

- family_graphs:

  A tibble with columns pid and family_graph_col. See prepare_graph for
  construction of the graphs. The family graphs Defaults to NULL.

- h2:

  A number representing the heritability on liability scale for a single
  phenotype. Must be non-negative. Note that under the liability
  threshold model, the heritability must also be at most 1. Defaults to
  0.5.

- pid:

  A string holding the name of the column in `.tbl` (or `family` and
  `threshs`) that hold the personal identifier(s). Defaults to "PID".

- fam_id:

  A string holding the name of the column in `.tbl` or `family` that
  holds the family identifier. Defaults to "fam_ID".

- family_graphs_col:

  Name of column with family graphs in family_graphs. Defaults to
  "fam_graph".

- role:

  A string holding the name of the column in `.tbl` that holds the role.
  Each role must be chosen from the following list of abbreviations

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

  - `pau[0-9]*` (Aunts/Uncles - paternal side). Defaults to "role".

- out:

  A character or numeric vector indicating whether the genetic component
  of the full liability, the full liability or both should be returned.
  If `out = c(1)` or `out = c("genetic")`, the genetic liability is
  estimated and returned. If `out = c(2)` or `out = c("full")`, the full
  liability is estimated and returned. If `out = c(1,2)` or
  `out = c("genetic", "full")`, both components are estimated and
  returned. Defaults to `c(1)`.

- tol:

  A number that is used as the convergence criterion for the Gibbs
  sampler. Equals the standard error of the mean. That is, a tolerance
  of 0.2 means that the standard error of the mean is below 0.2.
  Defaults to 0.01.

## Value

If `family` and `threshs` are two matrices, lists or data frames that
can be converted into tibbles, if `family` has two columns named like
the strings represented in `pid` and `fam_id`, if `threshs` has a column
named like the string given in `pid` as well as a column named "lower"
and a column named "upper" and if the liability-scale heritability `h2`,
`out`, `tol` and `always_add` are of the required form, then the
function returns a tibble with either four or six columns (depending on
the length of out). The first two columns correspond to the columns
`fam_id` and `pid` ' present in `family`. If `out` is equal to `c(1)` or
`c("genetic")`, the third and fourth column hold the estimated genetic
liability as well as the corresponding standard error, respectively. If
`out` equals `c(2)` or `c("full")`, the third and fourth column hold the
estimated full liability as well as the corresponding standard error,
respectively. If `out` is equal to `c(1,2)` or `c("genetic","full")`,
the third and fourth column hold the estimated genetic liability as well
as the corresponding standard error, respectively, while the fifth and
sixth column hold the estimated full liability as well as the
corresponding standard error, respectively.

## Details

This function can be used to estimate either the genetic component of
the full liability, the full liability or both. It is possible to input
either

## See also

[`future_apply`](https://future.apply.futureverse.org/reference/future_apply.html),
[`estimate_liability_multi`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability_multi.md),
[`estimate_liability`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md)

## Examples

``` r
sims <- simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, 
add_ind = TRUE, h2 = 0.5, n_sim=10, pop_prev = .05)
#
estimate_liability_single(.tbl = sims$thresholds, 
h2 = 0.5, pid = "indiv_ID", fam_id = "fam_ID", role = "role", out = c(1), 
tol = 0.01)
#> The number of workers is 1
#> # A tibble: 10 × 4
#>    fam_ID    indiv_ID  genetic_est genetic_se
#>    <chr>     <chr>           <dbl>      <dbl>
#>  1 fam_ID_1  fam_ID_1     -0.0116     0.00446
#>  2 fam_ID_2  fam_ID_2     -0.0240     0.00446
#>  3 fam_ID_3  fam_ID_3     -0.0385     0.00390
#>  4 fam_ID_4  fam_ID_4     -0.0275     0.00419
#>  5 fam_ID_5  fam_ID_5      1.33       0.00166
#>  6 fam_ID_6  fam_ID_6     -0.0121     0.00438
#>  7 fam_ID_7  fam_ID_7     -0.00609    0.00433
#>  8 fam_ID_8  fam_ID_8     -0.0240     0.00435
#>  9 fam_ID_9  fam_ID_9     -0.0377     0.00409
#> 10 fam_ID_10 fam_ID_10    -0.0304     0.00426
# 
sims <- simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, 
h2 = 0.5, n_sim=10, pop_prev = .05)
#> Warning: Neither fam_vec nor n_fam is specified...
#
estimate_liability_single(.tbl = sims$thresholds, 
h2 = 0.5, pid = "indiv_ID", fam_id = "fam_ID", role = "role",
out = c("genetic"), tol = 0.01)
#> The number of workers is 1
#> # A tibble: 10 × 4
#>    fam_ID    indiv_ID  genetic_est genetic_se
#>    <chr>     <chr>           <dbl>      <dbl>
#>  1 fam_ID_1  fam_ID_1     0.000939    0.00365
#>  2 fam_ID_2  fam_ID_2    -0.00306     0.00374
#>  3 fam_ID_3  fam_ID_3    -0.000359    0.00393
#>  4 fam_ID_4  fam_ID_4    -0.00829     0.00363
#>  5 fam_ID_5  fam_ID_5    -0.00799     0.00386
#>  6 fam_ID_6  fam_ID_6     0.00318     0.00370
#>  7 fam_ID_7  fam_ID_7    -0.00404     0.00384
#>  8 fam_ID_8  fam_ID_8    -0.00414     0.00383
#>  9 fam_ID_9  fam_ID_9     0.00231     0.00403
#> 10 fam_ID_10 fam_ID_10   -0.00171     0.00403
```
