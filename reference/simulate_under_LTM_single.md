# Simulate under the liability threshold model (single phenotype).

`simulate_under_LTM_single` simulates families and thresholds under the
liability threshold model for a given family structure and a single
phenotype. Please note that it is not possible to simulate different
family structures.

## Usage

``` r
simulate_under_LTM_single(
  fam_vec = c("m", "f", "s1", "mgm", "mgf", "pgm", "pgf"),
  n_fam = NULL,
  add_ind = TRUE,
  h2 = 0.5,
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

  A number representing the liability-scale heritability for a single
  phenotype. Must be non-negative. Note that under the liability
  threshold model, the heritability must also be at most 1. Defaults to
  0.5.

- n_sim:

  A positive number representing the number of simulations. Defaults to
  1000.

- pop_prev:

  A positive number representing the population prevalence, i.e. the
  overall prevalence in the population. Must be smaller than 1. Defaults
  to 0.1.

## Value

If either `fam_vec` or `n_fam` is used as the argument, if it is of the
required format, if the liability-scale heritability `h2` is a number
satisfying \\0 \leq h^2\\, `n_sim` is a strictly positive number, and
`pop_prev` is a positive number that is at most one, then the output
will be a list holding two tibbles. The first tibble, `sim_obs`, holds
the simulated liabilities, the disease status and the current
age/age-of-onset for all family members in each of the `n_sim` families.
The second tibble, `thresholds`, holds the family identifier, the
personal identifier, the role (specified in fam_vec or n_fam) as well as
the lower and upper thresholds for all individuals in all families. Note
that this tibble has the format required in
[`estimate_liability`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md).
In addition, note that if neither `fam_vec` nor `n_fam` are specified,
the function returns the disease status, the current age/age-of-onset,
the lower and upper thresholds, as well as the personal identifier for a
single individual, namely the individual under consideration (called
`o`). If both `fam_vec` and `n_fam` are defined, the user is asked to '
decide on which of the two vectors to use.

## See also

[`construct_covmat`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat.md),
[`simulate_under_LTM_multi`](https://emilmip.github.io/LTFHPlus/reference/simulate_under_LTM_multi.md),
[`simulate_under_LTM`](https://emilmip.github.io/LTFHPlus/reference/simulate_under_LTM.md)

## Examples

``` r
simulate_under_LTM_single()
#> $sim_obs
#> # A tibble: 1,000 × 26
#>    fam_ID          g      o       m      f     s1     mgm    mgf     pgm     pgf
#>    <chr>       <dbl>  <dbl>   <dbl>  <dbl>  <dbl>   <dbl>  <dbl>   <dbl>   <dbl>
#>  1 fam_ID_1  0.540    1.05   0.0659 -0.226  2.14  -1.08   -1.11  -0.354   0.854 
#>  2 fam_ID_2 -0.265    1.64   0.844   0.293 -0.777 -0.179  -0.372  1.99   -1.02  
#>  3 fam_ID_3 -1.10    -2.00  -0.205  -0.820  1.02  -0.944   1.74   0.0201 -1.19  
#>  4 fam_ID_4 -0.344    0.342  1.46    0.313  0.476 -0.0829 -0.107 -0.734   1.44  
#>  5 fam_ID_5  0.152    1.78  -1.21   -0.510  0.949  0.204   0.130 -0.141  -0.187 
#>  6 fam_ID_6 -0.132    0.659 -0.282   0.160 -0.772 -1.52   -1.43  -1.56   -0.911 
#>  7 fam_ID_7  1.43     1.50  -0.489   0.287  1.45  -0.844   0.255  1.53   -0.283 
#>  8 fam_ID_8 -0.107    0.597 -0.277   0.192 -0.997 -1.22    1.29   0.989  -0.477 
#>  9 fam_ID_9  0.00789 -0.364  0.482   1.14   1.73   0.158  -1.03  -0.618  -0.0198
#> 10 fam_ID_…  0.289    0.360 -0.154   0.359  1.42   0.979   0.106  1.57    0.146 
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
#>  1 fam_ID_1  fam_ID_1_1  o     -Inf     2.47
#>  2 fam_ID_2  fam_ID_2_1  o        1.64  1.64
#>  3 fam_ID_3  fam_ID_3_1  o     -Inf     2.91
#>  4 fam_ID_4  fam_ID_4_1  o     -Inf     3.24
#>  5 fam_ID_5  fam_ID_5_1  o        1.78  1.78
#>  6 fam_ID_6  fam_ID_6_1  o     -Inf     3.21
#>  7 fam_ID_7  fam_ID_7_1  o        1.51  1.51
#>  8 fam_ID_8  fam_ID_8_1  o     -Inf     3.48
#>  9 fam_ID_9  fam_ID_9_1  o     -Inf     2.83
#> 10 fam_ID_10 fam_ID_10_1 o     -Inf     2.72
#> # ℹ 7,990 more rows
#> 
simulate_under_LTM_single(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2), 
c("m","mgm","mgf","mhs")))
#> $sim_obs
#> # A tibble: 1,000 × 20
#>    fam_ID       g       o      m    mgm    mgf    mhs1    mhs2 o_status m_status
#>    <chr>    <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl> <lgl>    <lgl>   
#>  1 fam_I… -0.0389  1.00   -0.173  0.122  1.26  -2.11   -0.657  FALSE    FALSE   
#>  2 fam_I…  0.172  -0.145  -0.274  0.642 -0.477  0.680  -0.115  FALSE    FALSE   
#>  3 fam_I… -0.602  -0.985  -2.19  -0.558  1.18  -0.698  -1.54   FALSE    FALSE   
#>  4 fam_I…  0.512   0.971  -0.159  1.08  -1.25  -1.18   -0.685  FALSE    FALSE   
#>  5 fam_I…  0.818   0.879   0.173 -0.874  0.564 -0.424   0.0713 FALSE    FALSE   
#>  6 fam_I… -0.190   1.26   -0.511 -0.366  0.162 -0.335  -0.0346 FALSE    FALSE   
#>  7 fam_I…  0.0386  0.162   2.22   0.251  1.88   0.495   0.134  FALSE    TRUE    
#>  8 fam_I…  0.0758 -0.0225  0.205 -0.129  0.565 -0.671   2.25   FALSE    FALSE   
#>  9 fam_I…  0.170  -0.896  -0.253  2.70   0.778  0.209  -0.183  FALSE    FALSE   
#> 10 fam_I… -0.507  -0.800   0.338 -0.527  0.167 -0.0751  0.660  FALSE    FALSE   
#> # ℹ 990 more rows
#> # ℹ 10 more variables: mgm_status <lgl>, mgf_status <lgl>, mhs1_status <lgl>,
#> #   mhs2_status <lgl>, o_aoo <dbl>, m_aoo <dbl>, mgm_aoo <dbl>, mgf_aoo <dbl>,
#> #   mhs1_aoo <dbl>, mhs2_aoo <dbl>
#> 
#> $thresholds
#> # A tibble: 6,000 × 5
#>    fam_ID    indiv_ID    role  lower upper
#>    <chr>     <chr>       <chr> <dbl> <dbl>
#>  1 fam_ID_1  fam_ID_1_1  o      -Inf  2.83
#>  2 fam_ID_2  fam_ID_2_1  o      -Inf  3.35
#>  3 fam_ID_3  fam_ID_3_1  o      -Inf  3.35
#>  4 fam_ID_4  fam_ID_4_1  o      -Inf  3.48
#>  5 fam_ID_5  fam_ID_5_1  o      -Inf  2.51
#>  6 fam_ID_6  fam_ID_6_1  o      -Inf  2.72
#>  7 fam_ID_7  fam_ID_7_1  o      -Inf  2.79
#>  8 fam_ID_8  fam_ID_8_1  o      -Inf  2.63
#>  9 fam_ID_9  fam_ID_9_1  o      -Inf  3.55
#> 10 fam_ID_10 fam_ID_10_1 o      -Inf  2.55
#> # ℹ 5,990 more rows
#> 
simulate_under_LTM_single(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE, 
h2 = 0.5, n_sim = 500, pop_prev = .05)
#> $sim_obs
#> # A tibble: 500 × 10
#>    fam_ID        m      f      s1 m_status f_status s1_status m_aoo f_aoo s1_aoo
#>    <chr>     <dbl>  <dbl>   <dbl> <lgl>    <lgl>    <lgl>     <dbl> <dbl>  <dbl>
#>  1 fam_ID… -0.558  -1.37  -0.654  FALSE    FALSE    FALSE        45    48     19
#>  2 fam_ID…  0.0812  0.259  0.274  FALSE    FALSE    FALSE        35    40     16
#>  3 fam_ID… -1.55   -0.416 -1.83   FALSE    FALSE    FALSE        37    34     10
#>  4 fam_ID… -0.706   1.83  -0.0682 FALSE    TRUE     FALSE        61    66     33
#>  5 fam_ID… -0.234  -1.20   0.413  FALSE    FALSE    FALSE        54    59     30
#>  6 fam_ID…  1.68    1.30   2.08   TRUE     FALSE    TRUE         80    51     56
#>  7 fam_ID… -1.33    0.592 -1.59   FALSE    FALSE    FALSE        54    64     36
#>  8 fam_ID…  1.26    1.65   0.329  FALSE    TRUE     FALSE        67    92     39
#>  9 fam_ID…  0.341  -0.831  0.0608 FALSE    FALSE    FALSE        51    40     21
#> 10 fam_ID… -1.06   -0.768 -1.92   FALSE    FALSE    FALSE        57    48     27
#> # ℹ 490 more rows
#> 
#> $thresholds
#> # A tibble: 1,500 × 5
#>    fam_ID    indiv_ID    role    lower upper
#>    <chr>     <chr>       <chr>   <dbl> <dbl>
#>  1 fam_ID_1  fam_ID_1_1  m     -Inf     2.48
#>  2 fam_ID_2  fam_ID_2_1  m     -Inf     2.86
#>  3 fam_ID_3  fam_ID_3_1  m     -Inf     2.79
#>  4 fam_ID_4  fam_ID_4_1  m     -Inf     1.93
#>  5 fam_ID_5  fam_ID_5_1  m     -Inf     2.14
#>  6 fam_ID_6  fam_ID_6_1  m        1.68  1.68
#>  7 fam_ID_7  fam_ID_7_1  m     -Inf     2.14
#>  8 fam_ID_8  fam_ID_8_1  m     -Inf     1.81
#>  9 fam_ID_9  fam_ID_9_1  m     -Inf     2.25
#> 10 fam_ID_10 fam_ID_10_1 m     -Inf     2.05
#> # ℹ 1,490 more rows
#> 
simulate_under_LTM_single(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5, 
n_sim = 200, pop_prev = 0.05)
#> Warning: Neither fam_vec nor n_fam is specified...
#> $sim_obs
#> # A tibble: 200 × 5
#>    fam_ID         g       o o_status o_aoo
#>    <chr>      <dbl>   <dbl> <lgl>    <dbl>
#>  1 fam_ID_1  -1.27  -0.439  FALSE       25
#>  2 fam_ID_2  -0.553 -0.0114 FALSE       39
#>  3 fam_ID_3  -0.687 -0.850  FALSE       40
#>  4 fam_ID_4   0.535 -1.09   FALSE       20
#>  5 fam_ID_5   1.70   3.09   TRUE        29
#>  6 fam_ID_6  -0.783 -1.56   FALSE       10
#>  7 fam_ID_7  -0.243  0.0254 FALSE       30
#>  8 fam_ID_8   0.609  0.733  FALSE       13
#>  9 fam_ID_9  -0.178 -0.0371 FALSE       16
#> 10 fam_ID_10 -0.953 -1.09   FALSE       37
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 200 × 5
#>    fam_ID    indiv_ID    role    lower upper
#>    <chr>     <chr>       <chr>   <dbl> <dbl>
#>  1 fam_ID_1  fam_ID_1_1  o     -Inf     3.23
#>  2 fam_ID_2  fam_ID_2_1  o     -Inf     2.71
#>  3 fam_ID_3  fam_ID_3_1  o     -Inf     2.67
#>  4 fam_ID_4  fam_ID_4_1  o     -Inf     3.40
#>  5 fam_ID_5  fam_ID_5_1  o        3.09  3.09
#>  6 fam_ID_6  fam_ID_6_1  o     -Inf     3.73
#>  7 fam_ID_7  fam_ID_7_1  o     -Inf     3.05
#>  8 fam_ID_8  fam_ID_8_1  o     -Inf     3.63
#>  9 fam_ID_9  fam_ID_9_1  o     -Inf     3.54
#> 10 fam_ID_10 fam_ID_10_1 o     -Inf     2.79
#> # ℹ 190 more rows
#> 
```
