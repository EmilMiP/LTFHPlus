# Estimate genetic liability similar to LT-FH

Estimate genetic liability similar to LT-FH

## Usage

``` r
estimate_gen_liability_ltfh(
  h2,
  phen,
  child_threshold,
  parent_threshold,
  status_col_offspring = "CHILD_STATUS",
  status_col_father = "P1_STATUS",
  status_col_mother = "P2_STATUS",
  status_col_siblings = "SIB_STATUS",
  number_of_siblings_col = "NUM_SIBS",
  tol = 0.01
)
```

## Arguments

- h2:

  Liability scale heritability of the trait being analysed.

- phen:

  tibble or data.frame with status of the genotyped individual, parents
  and siblings.

- child_threshold:

  single numeric value that is used as threshold for the offspring and
  siblings.

- parent_threshold:

  single numeric value that is used as threshold for both parents

- status_col_offspring:

  Column name of status for the offspring

- status_col_father:

  Column name of status for the father

- status_col_mother:

  Column name of status for the mother

- status_col_siblings:

  Column name of status for the siblings

- number_of_siblings_col:

  Column name for the number of siblings for a given individual

- tol:

  Convergence criteria of the Gibbs sampler. Default is 0.01, meaning a
  standard error of the mean below 0.01

## Value

Returns the estimated genetic liabilities.

## Examples

``` r
phen <- data.frame(
CHILD_STATUS = c(0,0),
P1_STATUS = c(1,1),
P2_STATUS = c(0,1),
SIB_STATUS = c(1,0),
NUM_SIBS = c(2,0))

h2 <- 0.5
child_threshold <- 0.7
parent_threshold <- 0.8

estimate_gen_liability_ltfh(h2, phen, child_threshold, parent_threshold)
#>   CHILD_STATUS P1_STATUS P2_STATUS SIB_STATUS NUM_SIBS post_gen_liab
#> 1            0         1         0          1        2     0.1591523
#> 2            0         1         1          0        0     0.3500480
#>   post_gen_liab_se
#> 1      0.003150908
#> 2      0.003289368
```
