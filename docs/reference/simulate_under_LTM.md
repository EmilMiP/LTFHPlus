# Simulate under the liability threshold model.

`simulate_under_LTM` simulates families and thresholds under the
liability threshold model for a given family structure and a variable
number of phenotypes.Please note that it is not possible to simulate
different family structures.

## Usage

``` r
simulate_under_LTM(
  fam_vec = c("m", "f", "s1", "mgm", "mgf", "pgm", "pgf"),
  n_fam = NULL,
  add_ind = TRUE,
  h2 = 0.5,
  genetic_corrmat = NULL,
  full_corrmat = NULL,
  phen_names = NULL,
  n_sim = 1000,
  pop_prev = 0.1
)
```

## Arguments

- fam_vec:

  A vector of strings holding the different family members. All family
  members must be represented by strings from the following list:

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

  - `pau[0-9]*` (Aunts/Uncles - paternal side). Defaults to
    `c("m","f","s1","mgm","mgf","pgm","pgf")`.

- n_fam:

  A named vector holding the desired number of family members. See
  [`setNames`](https://rdrr.io/r/stats/setNames.html). All names must be
  picked from the list mentioned above. Defaults to `NULL`.

- add_ind:

  A logical scalar indicating whether the genetic component of the full
  liability as well as the full liability for the underlying target
  individual should be included in the covariance matrix. Defaults to
  `TRUE`.

- h2:

  Either a number or a numeric vector holding the liability-scale
  heritability(ies) for one or more phenotypes. All entries in `h2` must
  be non-negative. Note that under the liability threshold model, the
  heritabilities must also be at most 1. Defaults to 0.5.

- genetic_corrmat:

  Either `NULL` or a numeric matrix holding the genetic correlations
  between the desired phenotypes. Must be specified, if
  `length(h2)`\\\>0\\, and will be ignored if `h2` is a number. All
  diagonal entries in `genetic_corrmat` must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric. Defaults to `NULL`.

- full_corrmat:

  Either `NULL` or a numeric matrix holding the full correlations
  between the desired phenotypes. Must be specified, if
  `length(h2)`\\\>0\\, and will be ignored if `h2` is a number. All
  diagonal entries in `full_corrmat` must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric. Defaults to `NULL`.

- phen_names:

  Either `NULL` or character vector holding the phenotype names. These
  names will be used to create the row and column names for the
  covariance matrix. Must be specified, if `length(h2)` \\\> 0\\, and
  will be ignored if `h2` is a number. If it is not specified, the names
  will default to phenotype1, phenotype2, etc. Defaults to `NULL`.

- n_sim:

  A positive number representing the number of simulations. Defaults to
  1000.

- pop_prev:

  Either a number or a numeric vector holding the population
  prevalence(s), i.e. the overall prevalence(s) in the population. All
  entries in `pop_prev` must be positive and smaller than 1. Defaults to
  0.1.

## Value

If either `fam_vec` or `n_fam` is used as the argument, if it is of the
required format, if the liability-scale heritability `h2` is a number
satisfying \\0 \leq h^2\\, `n_sim` is a strictly positive number, and
`pop_prev` is a positive number that is at most one, then the output
will be a list containing two tibbles. The first tibble, `sim_obs`,
holds the simulated liabilities, the disease status and the current
age/age-of-onset for all family members in each of the `n_sim` families.
The second tibble, `thresholds`, holds the family identifier, the
personal identifier, the role (specified in fam_vec or n_fam) as well as
the lower and upper thresholds for all individuals in all families. Note
that this tibble has the format required in
[`estimate_liability`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md).
If either `fam_vec` or `n_fam` is used as the argument and if it is of
the required format, if `genetic_corrmat` and `full_corrmat` are two
numeric and symmetric matrices satisfying that all diagonal entries are
one and that all off-diagonal entries are between -1 and 1, if the
liability-scale heritabilities in `h2_vec` are numbers satisfying \\0
\leq h^2_i\\ for all \\i \in \\1,...,n_pheno\\\\, `n_sim` is a strictly
positive number, and `pop_prev` is a positive numeric vector such that
all entries are at most one, then the output will be a list containing
the following lists. The first outer list, which is named after the
first phenotype in `phen_names`, holds the tibble `sim_obs`, which holds
the simulated liabilities, the disease status and the current
age/age-of-onset for all family members in each of the `n_sim` families
for the first phenotype. As the first outer list, the second outer list,
which is named after the second phenotype in `phen_names`, holds the
tibble `sim_obs`, which holds the simulated liabilities, the disease
status and the current age/age-of-onset for all family members in each
of the `n_sim` families for the second phenotype. There is a list
containing `sim_obs` for each phenotype in `phen_names`. The last list
entry, `thresholds`, holds the family identifier, the personal
identifier, the role (specified in fam_vec or n_fam) as well as the
lower and upper thresholds for all individuals in all families and all
phenotypes. Note that this tibble has the format required in
[`estimate_liability`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md).
Finally, note that if neither `fam_vec` nor `n_fam` are specified, the
function returns the disease status, the current age/age-of-onset, the
lower and upper thresholds, as well as the personal identifier for a
single individual, namely the individual under consideration (called
`o`). If both `fam_vec` and `n_fam` are defined, the user is asked to '
decide on which of the two vectors to use.

## Details

This function can be used to simulate the case-control status, the
current age and age-of-onset as well as the lower and upper thresholds
for a variable number of phenotypes for all family members in each of
the `n_sim` families. If `h2` is a number, `simulate_under_LTM`
simulates the case- control status, the current age and age-of-onset as
well as thresholds for a single phenotype. However, if `h2` is a numeric
vector, if `genetic_corrmat` and `full_corrmat` are two symmetric
correlation matrices, and if `phen_names` and `pop_prev` are to numeric
vectors holding the phenotype names and the population prevalences,
respectively, then `simulate_under_LTM` simulates the case-control
status, the current age and age-of-onset as well as thresholds for two
or more (correlated) phenotypes. The family members can be specified
using one of two possible formats.

## See also

[`construct_covmat`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat.md)
[`simulate_under_LTM_single`](https://emilmip.github.io/LTFHPlus/reference/simulate_under_LTM_single.md)
[`simulate_under_LTM_multi`](https://emilmip.github.io/LTFHPlus/reference/simulate_under_LTM_multi.md)

## Examples

``` r
simulate_under_LTM()
#> $sim_obs
#> # A tibble: 1,000 × 26
#>    fam_ID        g      o      m       f      s1     mgm      mgf    pgm     pgf
#>    <chr>     <dbl>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl>    <dbl>  <dbl>   <dbl>
#>  1 fam_ID… -0.0519  0.852  0.335  0.276   1.05    0.0926  6.94e-1  1.34  -0.0264
#>  2 fam_ID…  0.511   1.44   0.539 -0.459   1.54    0.986   1.74e+0  1.57  -0.135 
#>  3 fam_ID…  0.895   0.356 -0.172  1.10   -0.697  -0.111  -8.26e-2  1.48  -1.32  
#>  4 fam_ID…  1.06    0.524  1.89   2.30    0.576  -0.277  -6.85e-2  1.27   1.14  
#>  5 fam_ID… -1.45   -1.13  -1.18  -0.126  -0.595  -0.958  -7.26e-2 -1.73  -0.680 
#>  6 fam_ID…  0.394  -0.191 -0.703  0.931   0.415   0.138  -3.81e-1  1.81   1.41  
#>  7 fam_ID…  0.0366 -0.220  1.02   0.152  -0.0920  1.90    1.36e-1 -1.49  -2.48  
#>  8 fam_ID… -0.562   0.425 -0.823 -0.119  -0.197  -1.26    9.60e-4 -0.195 -0.200 
#>  9 fam_ID… -0.0462  0.788 -0.155 -0.285   0.956   0.449   1.13e+0 -0.251 -0.787 
#> 10 fam_ID…  0.670   1.73   0.841 -0.0677 -0.521   0.824   5.34e-1  0.744 -0.191 
#> # ℹ 990 more rows
#> # ℹ 16 more variables: o_status <lgl>, m_status <lgl>, f_status <lgl>,
#> #   s1_status <lgl>, mgm_status <lgl>, mgf_status <lgl>, pgm_status <lgl>,
#> #   pgf_status <lgl>, o_aoo <dbl>, m_aoo <dbl>, f_aoo <dbl>, s1_aoo <dbl>,
#> #   mgm_aoo <dbl>, mgf_aoo <dbl>, pgm_aoo <dbl>, pgf_aoo <dbl>
#> 
#> $thresholds
#> # A tibble: 8,000 × 5
#>    fam_ID    indiv_ID    role    lower upper
#>    <chr>     <chr>       <chr>   <dbl> <dbl>
#>  1 fam_ID_1  fam_ID_1_1  o     -Inf     3.24
#>  2 fam_ID_2  fam_ID_2_1  o        1.44  1.44
#>  3 fam_ID_3  fam_ID_3_1  o     -Inf     2.63
#>  4 fam_ID_4  fam_ID_4_1  o     -Inf     2.83
#>  5 fam_ID_5  fam_ID_5_1  o     -Inf     3.03
#>  6 fam_ID_6  fam_ID_6_1  o     -Inf     3.45
#>  7 fam_ID_7  fam_ID_7_1  o     -Inf     2.95
#>  8 fam_ID_8  fam_ID_8_1  o     -Inf     3.06
#>  9 fam_ID_9  fam_ID_9_1  o     -Inf     3.48
#> 10 fam_ID_10 fam_ID_10_1 o        1.74  1.74
#> # ℹ 7,990 more rows
#> 

genetic_corrmat <- matrix(0.4, 3, 3)
diag(genetic_corrmat) <- 1
full_corrmat <- matrix(0.6, 3, 3)
diag(full_corrmat) <- 1

simulate_under_LTM(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2,2), 
c("m","mgm","mgf","s","mhs")))
#> $sim_obs
#> # A tibble: 1,000 × 26
#>    fam_ID        g      o       m    mgm     mgf       s1       s2   mhs1   mhs2
#>    <chr>     <dbl>  <dbl>   <dbl>  <dbl>   <dbl>    <dbl>    <dbl>  <dbl>  <dbl>
#>  1 fam_ID…  0.513   0.187 -0.0519  0.292 -0.772   1.40     0.868   -1.20   0.507
#>  2 fam_ID… -1.29   -1.91  -2.06   -1.14   0.0605  0.00188 -0.841   -0.475 -2.38 
#>  3 fam_ID… -0.0984 -0.642 -0.802   1.44  -2.50   -1.07     0.933    0.399 -2.62 
#>  4 fam_ID… -2.20   -1.99  -1.54   -0.535 -0.358  -1.94    -2.01    -2.52  -1.31 
#>  5 fam_ID… -0.323  -0.311 -1.40   -0.208 -0.146  -0.827    1.02    -0.537 -1.07 
#>  6 fam_ID…  0.483   3.10   1.20   -0.364 -0.805   0.350    0.550   -0.420  0.235
#>  7 fam_ID… -0.503  -1.27   0.730   1.66   1.91   -0.120   -0.352    0.926  0.465
#>  8 fam_ID… -0.274  -0.261 -1.32   -0.483 -2.05   -1.24    -1.18     0.939  1.12 
#>  9 fam_ID… -0.285  -0.339 -1.36   -2.32  -0.310  -1.22    -2.00    -0.250 -0.576
#> 10 fam_ID… -0.218  -0.238 -0.0574 -0.154 -0.578   1.91    -0.00965  0.408  0.921
#> # ℹ 990 more rows
#> # ℹ 16 more variables: o_status <lgl>, m_status <lgl>, mgm_status <lgl>,
#> #   mgf_status <lgl>, s1_status <lgl>, s2_status <lgl>, mhs1_status <lgl>,
#> #   mhs2_status <lgl>, o_aoo <dbl>, m_aoo <dbl>, mgm_aoo <dbl>, mgf_aoo <dbl>,
#> #   s1_aoo <dbl>, s2_aoo <dbl>, mhs1_aoo <dbl>, mhs2_aoo <dbl>
#> 
#> $thresholds
#> # A tibble: 8,000 × 5
#>    fam_ID    indiv_ID    role    lower upper
#>    <chr>     <chr>       <chr>   <dbl> <dbl>
#>  1 fam_ID_1  fam_ID_1_1  o     -Inf     3.35
#>  2 fam_ID_2  fam_ID_2_1  o     -Inf     2.79
#>  3 fam_ID_3  fam_ID_3_1  o     -Inf     3.28
#>  4 fam_ID_4  fam_ID_4_1  o     -Inf     2.91
#>  5 fam_ID_5  fam_ID_5_1  o     -Inf     3.03
#>  6 fam_ID_6  fam_ID_6_1  o        3.10  3.10
#>  7 fam_ID_7  fam_ID_7_1  o     -Inf     3.52
#>  8 fam_ID_8  fam_ID_8_1  o     -Inf     3.10
#>  9 fam_ID_9  fam_ID_9_1  o     -Inf     3.31
#> 10 fam_ID_10 fam_ID_10_1 o     -Inf     2.91
#> # ℹ 7,990 more rows
#> 

simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE, 
genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, n_sim = 200)
#> $sim_obs
#> # A tibble: 200 × 10
#>    fam_ID        m      f      s1 m_status f_status s1_status m_aoo f_aoo s1_aoo
#>    <chr>     <dbl>  <dbl>   <dbl> <lgl>    <lgl>    <lgl>     <dbl> <dbl>  <dbl>
#>  1 fam_ID…  0.559  -0.406 -0.346  FALSE    FALSE    FALSE        50    47     26
#>  2 fam_ID…  0.494  -1.46   1.57   FALSE    FALSE    TRUE         43    42     63
#>  3 fam_ID… -1.53    0.266 -2.39   FALSE    FALSE    FALSE        37    33     11
#>  4 fam_ID… -1.20    0.125  0.939  FALSE    FALSE    FALSE        47    55     29
#>  5 fam_ID… -0.966   0.736 -1.12   FALSE    FALSE    FALSE        52    49     23
#>  6 fam_ID…  1.73   -0.675 -0.0872 TRUE     FALSE    FALSE        57    52     32
#>  7 fam_ID…  0.488  -0.286 -0.668  FALSE    FALSE    FALSE        36    46     17
#>  8 fam_ID…  0.944  -0.266  0.133  FALSE    FALSE    FALSE        37    39     12
#>  9 fam_ID…  2.57   -0.126  0.0400 TRUE     FALSE    FALSE        37    46     26
#> 10 fam_ID…  0.0848 -1.34  -0.0793 FALSE    FALSE    FALSE        59    62     38
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 600 × 5
#>    fam_ID    indiv_ID    role    lower upper
#>    <chr>     <chr>       <chr>   <dbl> <dbl>
#>  1 fam_ID_1  fam_ID_1_1  m     -Inf     2.01
#>  2 fam_ID_2  fam_ID_2_1  m     -Inf     2.30
#>  3 fam_ID_3  fam_ID_3_1  m     -Inf     2.55
#>  4 fam_ID_4  fam_ID_4_1  m     -Inf     2.13
#>  5 fam_ID_5  fam_ID_5_1  m     -Inf     1.93
#>  6 fam_ID_6  fam_ID_6_1  m        1.74  1.74
#>  7 fam_ID_7  fam_ID_7_1  m     -Inf     2.59
#>  8 fam_ID_8  fam_ID_8_1  m     -Inf     2.55
#>  9 fam_ID_9  fam_ID_9_1  m        2.55  2.55
#> 10 fam_ID_10 fam_ID_10_1 m     -Inf     1.68
#> # ℹ 590 more rows
#> 

simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5, 
n_sim = 200, pop_prev = 0.05)
#> Warning: Neither fam_vec nor n_fam is specified...
#> $sim_obs
#> # A tibble: 200 × 5
#>    fam_ID          g       o o_status o_aoo
#>    <chr>       <dbl>   <dbl> <lgl>    <dbl>
#>  1 fam_ID_1   0.803   0.267  FALSE       30
#>  2 fam_ID_2  -0.709  -1.60   FALSE       15
#>  3 fam_ID_3   1.07    1.25   FALSE       24
#>  4 fam_ID_4  -0.787  -1.96   FALSE       23
#>  5 fam_ID_5   0.161   0.558  FALSE       21
#>  6 fam_ID_6   0.0272  0.0553 FALSE       20
#>  7 fam_ID_7  -0.475  -0.287  FALSE       29
#>  8 fam_ID_8   0.974   1.60   FALSE       26
#>  9 fam_ID_9  -0.591  -0.841  FALSE       34
#> 10 fam_ID_10 -0.772  -1.38   FALSE       16
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 200 × 5
#>    fam_ID    indiv_ID    role  lower upper
#>    <chr>     <chr>       <chr> <dbl> <dbl>
#>  1 fam_ID_1  fam_ID_1_1  o      -Inf  3.05
#>  2 fam_ID_2  fam_ID_2_1  o      -Inf  3.57
#>  3 fam_ID_3  fam_ID_3_1  o      -Inf  3.26
#>  4 fam_ID_4  fam_ID_4_1  o      -Inf  3.30
#>  5 fam_ID_5  fam_ID_5_1  o      -Inf  3.37
#>  6 fam_ID_6  fam_ID_6_1  o      -Inf  3.40
#>  7 fam_ID_7  fam_ID_7_1  o      -Inf  3.09
#>  8 fam_ID_8  fam_ID_8_1  o      -Inf  3.19
#>  9 fam_ID_9  fam_ID_9_1  o      -Inf  2.90
#> 10 fam_ID_10 fam_ID_10_1 o      -Inf  3.54
#> # ℹ 190 more rows
#> 
```
