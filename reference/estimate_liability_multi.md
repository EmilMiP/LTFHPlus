# Estimating the genetic or full liability for multiple phenotypes

`estimate_liability_multi` estimates the genetic component of the full
liability and/or the full liability for a number of individuals based on
their family history for a variable number of phenotypes.

## Usage

``` r
estimate_liability_multi(
  .tbl = NULL,
  family_graphs = NULL,
  h2_vec,
  genetic_corrmat,
  full_corrmat,
  phen_names = NULL,
  pid = "PID",
  fam_id = "fam_ID",
  role = "role",
  family_graphs_col = "fam_graph",
  out = c(1),
  tol = 0.01
)
```

## Arguments

- .tbl:

  A matrix, list or data frame that can be converted into a tibble. Must
  have at least seven columns that hold the family identifier, the
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

- h2_vec:

  A numeric vector representing the liability-scale heritabilities for
  all phenotypes. All entries in h2_vec must be non-negative and at most
  1.

- genetic_corrmat:

  A numeric matrix holding the genetic correlations between the desired
  phenotypes. All diagonal entries must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric.

- full_corrmat:

  A numeric matrix holding the full correlations between the desired
  phenotypes. All diagonal entries must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric.

- phen_names:

  A character vector holding the phenotype names. These names will be
  used to create the row and column names for the covariance matrix. If
  it is not specified, the names will default to phenotype1, phenotype2,
  etc. Defaults to NULL.

- pid:

  A string holding the name of the column in `family` and `threshs` that
  hold the personal identifier(s). Defaults to "PID".

- fam_id:

  A string holding the name of the column in `family` that holds the
  family identifier. Defaults to "fam_ID".

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

## Value

If `family` and `threshs` are two matrices, lists or data frames that
can be converted into tibbles, if `family` has two columns named like
the strings represented in `pid` and `fam_id`, if `threshs` has a column
named like the string given in `pid` as well as a column named `"lower"`
and a column named `"upper"` and if the liability-scale heritabilities
in `h2_vec`, `genetic_corrmat`, `full_corrmat`, `out` and `tol` are of
the required form, then the function returns a tibble with at least six
columns (depending on the length of out). The first two columns
correspond to the columns `fam_id` and `pid` present in the tibble
`family`. If `out` is equal to `c(1)` or `c("genetic")`, the third and
fourth columns hold the estimated genetic liability as well as the
corresponding standard error for the first phenotype, respectively. If
`out` equals `c(2)` or `c("full")`, the third and fourth columns hold
the estimated full liability as well as the corresponding standard error
for the first phenotype, respectively. If `out` is equal to `c(1,2)` or
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
[`estimate_liability`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md)

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
estimate_liability_multi(.tbl = sims$thresholds, h2_vec = rep(.5,3), 
genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat,
pid = "indiv_ID", fam_id = "fam_ID", role = "role", out = c(1), 
phen_names = paste0("phenotype", 1:3), tol = 0.01)
#> The number of workers is 1
#> # A tibble: 1 × 8
#>   fam_ID   indiv_ID genetic_phenotype1_est genetic_phenotype1_se
#>   <chr>    <chr>                     <dbl>                 <dbl>
#> 1 fam_ID_1 fam_ID_1                  0.410               0.00686
#> # ℹ 4 more variables: genetic_phenotype2_est <dbl>,
#> #   genetic_phenotype2_se <dbl>, genetic_phenotype3_est <dbl>,
#> #   genetic_phenotype3_se <dbl>

```
