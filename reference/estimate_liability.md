# Estimating the genetic or full liability for a variable number of phenotypes

`estimate_liability` estimates the genetic component of the full
liability and/or the full liability for a number of individuals based on
their family history for one or more phenotypes. It is a wrapper around
[`estimate_liability_single`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability_single.md)
and
[`estimate_liability_multi`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability_multi.md).

## Usage

``` r
estimate_liability(
  .tbl = NULL,
  family_graphs = NULL,
  h2 = 0.5,
  pid = "PID",
  fam_id = "fam_ID",
  role = "role",
  family_graphs_col = "fam_graph",
  out = c(1),
  tol = 0.01,
  genetic_corrmat = NULL,
  full_corrmat = NULL,
  phen_names = NULL
)
```

## Arguments

- .tbl:

  A matrix, list or data frame that can be converted into a tibble. Must
  have at least five columns that hold the family identifier, the
  personal identifier, the role and the lower and upper thresholds for
  all phenotypes of interest. Note that the role must be one of the
  following abbreviations

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

  Either a number representing the heritability on liability scale for a
  single phenotype, or a numeric vector representing the liability-scale
  heritabilities for all phenotypes. All entries in h2 must be
  non-negative and at most 1.

- pid:

  A string holding the name of the column in `family` and `threshs` that
  hold the personal identifier(s). Defaults to `"PID"`.

- fam_id:

  A string holding the name of the column in `family` that holds the
  family identifier. Defaults to `"fam_ID"`.

- role:

  A string holding the name of the column in `.tbl` that holds the
  role.Each role must be chosen from the following list of abbreviations

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

- family_graphs_col:

  Name of column with family graphs in family_graphs. Defaults to
  "fam_graph".

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

- genetic_corrmat:

  Either `NULL` (if `h2` is a number) or a numeric matrix (if `h2` is a
  vector of length \> 1) holding the genetic correlations between the
  desired phenotypes. All diagonal entries must be equal to one, while
  all off-diagonal entries must be between -1 and 1. In addition, the
  matrix must be symmetric. Defaults to `NULL`.

- full_corrmat:

  Either `NULL` (if `h2` is a number) or a numeric matrix (if `h2` is a
  vector of length \> 1) holding the full correlations between the
  desired phenotypes. All diagonal entries must be equal to one, while
  all off-diagonal entries must be between -1 and 1. In addition, the
  matrix must be symmetric. Defaults to `NULL`.

- phen_names:

  Either `NULL` or a character vector holding the phenotype names. These
  names will be used to create the row and column names for the
  covariance matrix. If it is not specified, the names will default to
  phenotype1, phenotype2, etc. Defaults to NULL.

## Value

If `family` and `threshs` are two matrices, lists or data frames that
can be converted into tibbles, if `family` has two columns named like
the strings represented in `pid` and `fam_id`, if `threshs` has a column
named like the string given in `pid` as well as a column named "lower"
and a column named "upper" and if the liability-scale heritability `h2`
is a number (`length(h2)=1`), and `out`, `tol` and `always_add` are of
the required form, then the function returns a tibble with either four
or six columns (depending on the length of out). The first two columns
correspond to the columns `fam_id` and `pid` ' present in `family`. If
`out` is equal to `c(1)` or `c("genetic")`, the third and fourth column
hold the estimated genetic liability as well as the corresponding
standard error, respectively. If `out` equals `c(2)` or `c("full")`, the
third and fourth column hold the estimated full liability as well as the
corresponding standard error, respectively. If `out` is equal to
`c(1,2)` or `c("genetic","full")`, the third and fourth column hold the
estimated genetic liability as well as the corresponding standard error,
respectively, while the fifth and sixth column hold the estimated full
liability as well as the corresponding standard error, respectively. If
`h2` is a numeric vector of length greater than 1 and if
`genetic_corrmat`, `full_corrmat`, `out` and `tol` are of the required
form, then the function returns a tibble with at least six columns
(depending on the length of out). The first two columns correspond to
the columns `fam_id` and `pid` present in the tibble `family`. If `out`
is equal to `c(1)` or `c("genetic")`, the third and fourth columns hold
the estimated genetic liability as well as the corresponding standard
error for the first phenotype, respectively. If `out` equals `c(2)` or
`c("full")`, the third and fourth columns hold the estimated full
liability as well as the corresponding standard error for the first
phenotype, respectively. If `out` is equal to `c(1,2)` or
`c("genetic","full")`, the third and fourth columns hold the estimated
genetic liability as well as the corresponding standard error for the
first phenotype, respectively, while the fifth and sixth columns hold
the estimated full liability as well as the corresponding standard error
for the first phenotype, respectively. The remaining columns hold the
estimated genetic liabilities and/or the estimated full liabilities as
well as the corresponding standard errors for the remaining phenotypes.

## Details

This function can be used to estimate either the genetic component of
the full liability, the full liability or both for a variable number of
traits.

## See also

[`future_apply`](https://future.apply.futureverse.org/reference/future_apply.html),
[`estimate_liability_single`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability_single.md),
[`estimate_liability_multi`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability_multi.md)

## Examples

``` r
genetic_corrmat <- matrix(0.4, 3, 3)
diag(genetic_corrmat) <- 1
full_corrmat <- matrix(0.6, 3, 3)
diag(full_corrmat) <- 1
#
sims <- simulate_under_LTM(fam_vec = c("m","f"), n_fam = NULL, add_ind = TRUE, 
genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, h2 = rep(.5,3), 
n_sim = 1, pop_prev = rep(.1,3))
#> Warning: There was 1 warning in `mutate()`.
#> ℹ In argument: `max_age = purrr::invoke(pmax, select(., ends_with("_age")))`.
#> Caused by warning:
#> ! `invoke()` was deprecated in purrr 1.0.0.
#> ℹ Please use `exec()` instead.
#> ℹ The deprecated feature was likely used in the LTFHPlus package.
#>   Please report the issue at <https://github.com/EmilMiP/LTFHPlus/issues>.
estimate_liability(.tbl = sims$thresholds, h2 = rep(.5,3), 
genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat,
pid = "indiv_ID", fam_id = "fam_ID", role = "role", out = c(1), 
phen_names = paste0("phenotype", 1:3), tol = 0.01)
#> The number of workers is 1
#> # A tibble: 1 × 8
#>   fam_ID   indiv_ID genetic_phenotype1_est genetic_phenotype1_se
#>   <chr>    <chr>                     <dbl>                 <dbl>
#> 1 fam_ID_1 fam_ID_1                -0.0920               0.00720
#> # ℹ 4 more variables: genetic_phenotype2_est <dbl>,
#> #   genetic_phenotype2_se <dbl>, genetic_phenotype3_est <dbl>,
#> #   genetic_phenotype3_se <dbl>
```
