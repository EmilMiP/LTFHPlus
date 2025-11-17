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
#>    fam_ID        g       o       m       f      s1     mgm    mgf    pgm     pgf
#>    <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
#>  1 fam_ID…  0.882   0.658   0.420   1.04    0.375   0.691   1.36   0.128 -0.105 
#>  2 fam_ID…  1.25    0.743  -0.252   1.76    1.32    1.61    1.55   0.145  1.87  
#>  3 fam_ID… -0.0854 -0.203   1.15   -1.09   -0.0995  0.0126  1.41  -1.63   0.793 
#>  4 fam_ID…  0.641   1.97    2.21   -0.0599 -0.149  -0.0414  1.31   0.929 -1.77  
#>  5 fam_ID… -0.646  -1.05    0.0164 -0.362  -0.891   0.119  -1.72  -0.754  0.524 
#>  6 fam_ID… -0.262  -0.616   0.527   0.297   0.313  -0.324   1.74   1.31   0.0850
#>  7 fam_ID…  0.0423  0.674   0.809  -0.328   1.61    0.301  -1.36  -2.44  -0.985 
#>  8 fam_ID…  0.177  -0.510  -0.135  -0.179  -1.15    0.0461 -0.223 -0.215 -0.683 
#>  9 fam_ID…  0.605  -0.0532 -0.0695  1.04    0.605   1.15   -0.285 -0.640  0.361 
#> 10 fam_ID…  1.28    1.15    0.298  -0.399   0.930   0.557   0.870 -0.143  0.919 
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
#>  1 fam_ID_1  fam_ID_1_1  o     -Inf     2.63
#>  2 fam_ID_2  fam_ID_2_1  o     -Inf     2.83
#>  3 fam_ID_3  fam_ID_3_1  o     -Inf     3.03
#>  4 fam_ID_4  fam_ID_4_1  o        1.97  1.97
#>  5 fam_ID_5  fam_ID_5_1  o     -Inf     2.95
#>  6 fam_ID_6  fam_ID_6_1  o     -Inf     3.06
#>  7 fam_ID_7  fam_ID_7_1  o     -Inf     3.48
#>  8 fam_ID_8  fam_ID_8_1  o     -Inf     2.76
#>  9 fam_ID_9  fam_ID_9_1  o     -Inf     2.47
#> 10 fam_ID_10 fam_ID_10_1 o     -Inf     2.99
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
#>    fam_ID          g       o      m    mgm    mgf      s1     s2   mhs1    mhs2
#>    <chr>       <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
#>  1 fam_ID_1  -0.261   0.0258 -0.465  0.130 -0.641 -0.572  -1.64  -2.22  -1.36  
#>  2 fam_ID_2   0.994   0.687  -0.630  0.627  0.690  0.888   1.54  -1.07   1.27  
#>  3 fam_ID_3  -0.978  -1.26   -1.85  -0.925 -1.66   0.0732  0.265  0.850  0.288 
#>  4 fam_ID_4  -0.575  -0.0307 -0.518  0.896  1.84  -1.83   -2.19  -2.20  -1.32  
#>  5 fam_ID_5   0.922   0.398  -0.550 -0.208  1.80  -0.201   0.518  2.19  -0.340 
#>  6 fam_ID_6   0.922   0.259   1.38   0.640 -0.549  0.395   1.39   0.890  0.137 
#>  7 fam_ID_7  -0.0667 -0.824  -0.373 -0.910  1.05   0.309  -0.708 -0.486 -2.02  
#>  8 fam_ID_8  -0.293   1.72    0.246 -1.24  -0.732  0.901   0.410 -0.135  0.562 
#>  9 fam_ID_9  -0.356  -1.85    0.575  0.674 -0.673  0.366   1.12   0.130  0.0916
#> 10 fam_ID_10 -1.55   -0.235  -0.866 -0.376  1.12  -1.08   -0.982 -0.970 -1.42  
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
#>  1 fam_ID_1  fam_ID_1_1  o     -Inf     2.91
#>  2 fam_ID_2  fam_ID_2_1  o     -Inf     3.03
#>  3 fam_ID_3  fam_ID_3_1  o     -Inf     2.55
#>  4 fam_ID_4  fam_ID_4_1  o     -Inf     3.52
#>  5 fam_ID_5  fam_ID_5_1  o     -Inf     3.10
#>  6 fam_ID_6  fam_ID_6_1  o     -Inf     3.31
#>  7 fam_ID_7  fam_ID_7_1  o     -Inf     2.91
#>  8 fam_ID_8  fam_ID_8_1  o        1.71  1.71
#>  9 fam_ID_9  fam_ID_9_1  o     -Inf     3.06
#> 10 fam_ID_10 fam_ID_10_1 o     -Inf     3.17
#> # ℹ 7,990 more rows
#> 

simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE, 
genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, n_sim = 200)
#> $sim_obs
#> # A tibble: 200 × 10
#>    fam_ID       m        f     s1 m_status f_status s1_status m_aoo f_aoo s1_aoo
#>    <chr>    <dbl>    <dbl>  <dbl> <lgl>    <lgl>    <lgl>     <dbl> <dbl>  <dbl>
#>  1 fam_ID… -0.720  0.237   -0.875 FALSE    FALSE    FALSE        43    42     16
#>  2 fam_ID…  0.615 -0.588   -1.67  FALSE    FALSE    FALSE        37    33     11
#>  3 fam_ID… -1.09   0.487    1.17  FALSE    FALSE    FALSE        47    55     29
#>  4 fam_ID… -0.432 -1.35    -0.426 FALSE    FALSE    FALSE        52    49     23
#>  5 fam_ID… -0.322 -0.00678 -0.851 FALSE    FALSE    FALSE        54    52     32
#>  6 fam_ID… -0.647  0.157   -0.153 FALSE    FALSE    FALSE        36    46     17
#>  7 fam_ID… -0.131  1.40     1.86  FALSE    TRUE     TRUE         37    71     54
#>  8 fam_ID…  0.767  0.572    0.391 FALSE    FALSE    FALSE        44    46     26
#>  9 fam_ID… -0.557 -0.578   -1.14  FALSE    FALSE    FALSE        59    62     38
#> 10 fam_ID…  0.720 -0.390   -2.12  FALSE    FALSE    FALSE        53    51     30
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 600 × 5
#>    fam_ID    indiv_ID    role  lower upper
#>    <chr>     <chr>       <chr> <dbl> <dbl>
#>  1 fam_ID_1  fam_ID_1_1  m      -Inf  2.30
#>  2 fam_ID_2  fam_ID_2_1  m      -Inf  2.55
#>  3 fam_ID_3  fam_ID_3_1  m      -Inf  2.13
#>  4 fam_ID_4  fam_ID_4_1  m      -Inf  1.93
#>  5 fam_ID_5  fam_ID_5_1  m      -Inf  1.85
#>  6 fam_ID_6  fam_ID_6_1  m      -Inf  2.59
#>  7 fam_ID_7  fam_ID_7_1  m      -Inf  2.55
#>  8 fam_ID_8  fam_ID_8_1  m      -Inf  2.26
#>  9 fam_ID_9  fam_ID_9_1  m      -Inf  1.68
#> 10 fam_ID_10 fam_ID_10_1 m      -Inf  1.89
#> # ℹ 590 more rows
#> 

simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5, 
n_sim = 200, pop_prev = 0.05)
#> Warning: Neither fam_vec nor n_fam is specified...
#> $sim_obs
#> # A tibble: 200 × 5
#>    fam_ID         g       o o_status o_aoo
#>    <chr>      <dbl>   <dbl> <lgl>    <dbl>
#>  1 fam_ID_1  -0.396 -1.13   FALSE       15
#>  2 fam_ID_2   0.213  0.0171 FALSE       24
#>  3 fam_ID_3  -1.04  -1.49   FALSE       23
#>  4 fam_ID_4  -0.856 -0.870  FALSE       21
#>  5 fam_ID_5  -0.382 -0.635  FALSE       20
#>  6 fam_ID_6  -0.242  0.398  FALSE       29
#>  7 fam_ID_7  -0.219  0.939  FALSE       26
#>  8 fam_ID_8   0.419  0.731  FALSE       34
#>  9 fam_ID_9  -0.681 -0.551  FALSE       16
#> 10 fam_ID_10  0.260 -0.0593 FALSE       34
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 200 × 5
#>    fam_ID    indiv_ID    role  lower upper
#>    <chr>     <chr>       <chr> <dbl> <dbl>
#>  1 fam_ID_1  fam_ID_1_1  o      -Inf  3.57
#>  2 fam_ID_2  fam_ID_2_1  o      -Inf  3.26
#>  3 fam_ID_3  fam_ID_3_1  o      -Inf  3.30
#>  4 fam_ID_4  fam_ID_4_1  o      -Inf  3.37
#>  5 fam_ID_5  fam_ID_5_1  o      -Inf  3.40
#>  6 fam_ID_6  fam_ID_6_1  o      -Inf  3.09
#>  7 fam_ID_7  fam_ID_7_1  o      -Inf  3.19
#>  8 fam_ID_8  fam_ID_8_1  o      -Inf  2.90
#>  9 fam_ID_9  fam_ID_9_1  o      -Inf  3.54
#> 10 fam_ID_10 fam_ID_10_1 o      -Inf  2.90
#> # ℹ 190 more rows
#> 
```
