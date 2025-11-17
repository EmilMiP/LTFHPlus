# How the covariance is constructed

### What is the input format for LTFHPlus functions?

With version 2.0, updates have been made to the input and the functions
available to estimate the (genetic) liability. Previously, a list entry
format with a set order was expected, where the proband was first, then
followed by father, mother, and any siblings. This limited analysis to
only the immediate family, but if information on, e.g. half-siblings,
grandparents etc, was available, it could not be readily used. Now the
input does not require a set ordering, but instead the user is expected
to provide information on the familial relation to the proband,
e.g. mother, paternal half-sibling, etc. This allows for far more
flexibility for the user to include the familial information that is
available.

### Family Input

The function used to estimate the genetic (or full liability) of an
individual is `estimate_liability`. The family input is input through
`.tbl`, which is a long format where each row is an individual. A role
must accompany each individual. The family relationship to the proband
has its own column.

### Simulate data

From `simulate_under_LTM` an example of the full input data can be seen.
It returns a list, first entry is `sim_obs`, and contains all the
underlying liabilities, status, and age of onset or age for controls.
The second entry is called `thresholds` and it contains a family ID,
individual ID, family relationship to the proband, and a lower and upper
threshold for each individual. The following example simulates a family
with the index person, a mother, a father, and a single sibling. Other
family members can also be used. See the documentation for
[`simulate_under_LTM()`](https://emilmip.github.io/LTFHPlus/reference/simulate_under_LTM.md)
for more information.

``` r
 sims <- simulate_under_LTM(fam_vec = c("m","f","s1"),
                            n_fam = NULL, 
                            add_ind = TRUE, 
                            h2 = 0.5, 
                            n_sim = 10, 
                            pop_prev = .05)
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `max_age = purrr::invoke(pmax, select(., ends_with("_age")))`.
    ## Caused by warning:
    ## ! `invoke()` was deprecated in purrr 1.0.0.
    ## ℹ Please use `exec()` instead.
    ## ℹ The deprecated feature was likely used in the LTFHPlus package.
    ##   Please report the issue at <https://github.com/EmilMiP/LTFHPlus/issues>.

``` r
sims$sim_obs
```

    ## # A tibble: 10 × 14
    ##    fam_ID          g      o       m       f     s1 o_status m_status f_status
    ##    <chr>       <dbl>  <dbl>   <dbl>   <dbl>  <dbl> <lgl>    <lgl>    <lgl>   
    ##  1 fam_ID_1  -1.01   -0.375 -2.47   -0.0384  0.179 FALSE    FALSE    FALSE   
    ##  2 fam_ID_2   0.0466 -1.44  -0.310  -0.307  -0.375 FALSE    FALSE    FALSE   
    ##  3 fam_ID_3  -0.0259  0.520  2.11   -1.59    0.545 FALSE    TRUE     FALSE   
    ##  4 fam_ID_4  -1.32   -1.09  -0.468   0.127  -1.11  FALSE    FALSE    FALSE   
    ##  5 fam_ID_5   0.543   0.608 -0.974   1.07    1.87  FALSE    FALSE    FALSE   
    ##  6 fam_ID_6  -0.632  -1.14  -0.284  -1.10   -1.68  FALSE    FALSE    FALSE   
    ##  7 fam_ID_7   0.874   0.653  0.358   1.74    0.456 FALSE    FALSE    TRUE    
    ##  8 fam_ID_8  -0.601  -1.79  -0.361  -0.396   0.765 FALSE    FALSE    FALSE   
    ##  9 fam_ID_9  -0.135  -0.566 -0.0500 -0.252   0.343 FALSE    FALSE    FALSE   
    ## 10 fam_ID_10  1.53    0.752  0.708   0.248  -1.44  FALSE    FALSE    FALSE   
    ## # ℹ 5 more variables: s1_status <lgl>, o_aoo <dbl>, m_aoo <dbl>, f_aoo <dbl>,
    ## #   s1_aoo <dbl>

``` r
sims$thresholds
```

    ## # A tibble: 40 × 5
    ##    fam_ID    indiv_ID    role  lower upper
    ##    <chr>     <chr>       <chr> <dbl> <dbl>
    ##  1 fam_ID_1  fam_ID_1_1  o      -Inf  3.70
    ##  2 fam_ID_2  fam_ID_2_1  o      -Inf  3.57
    ##  3 fam_ID_3  fam_ID_3_1  o      -Inf  3.54
    ##  4 fam_ID_4  fam_ID_4_1  o      -Inf  3.37
    ##  5 fam_ID_5  fam_ID_5_1  o      -Inf  2.86
    ##  6 fam_ID_6  fam_ID_6_1  o      -Inf  2.79
    ##  7 fam_ID_7  fam_ID_7_1  o      -Inf  3.01
    ##  8 fam_ID_8  fam_ID_8_1  o      -Inf  3.30
    ##  9 fam_ID_9  fam_ID_9_1  o      -Inf  3.60
    ## 10 fam_ID_10 fam_ID_10_1 o      -Inf  3.54
    ## # ℹ 30 more rows

### Covariance Function and examples

We construct the covariance matrix for each family being analysed during
run-time. The covariance function that is used internally in
`estimate_liability` has been updated to allow for a higher degree of
flexibility. This means it is up to the user to provide the familial
relationship, and `construct_covmat` creates the corresponding
covariance matrix based on the heritability and expected genetic overlap
between two individuals.

`construct_covmat` defaults to a family structure with both parents, one
sibling, and the paternal and maternal grandparents. The input format
for `construct_covmat` can be specified in two different ways, either
`fam_vec` (the method used internally in `estimate_liability`) or with
`n_fam`. For `fam_vec` a vector of strings from the list of possible
familial relationships must be provided For the full list, please see
documentation for `construct_covmat`. Family members will then appear in
the covariance matrix in the same order as they appear in `fam_vec`. For
`n_fam` a named vector is provided, where the *names* of the named
vector corresponding to the familial relationship and the values of the
vector corresponds to how often that particular familial role appears.

In order to illustrate the different possible families, we will provide
some examples. If no family information is available, but the age of
onset information is still available, we can use the simplest
covariance, which only contains the genetic and full liability of the
index person:

``` r
# no family members
construct_covmat(fam_vec = NULL, n_fam = NULL, h2 = .5)
```

    ## Warning in construct_covmat_single(fam_vec = fam_vec, n_fam = n_fam, add_ind = add_ind, : Neither fam_vec nor n_fam is specified...

    ##     g   o
    ## g 0.5 0.5
    ## o 0.5 1.0
    ## attr(,"fam_vec")
    ## [1] "g" "o"
    ## attr(,"n_fam")
    ## g o 
    ## 1 1 
    ## attr(,"add_ind")
    ## [1] TRUE
    ## attr(,"h2")
    ## [1] 0.5

The default family contains the index person as well as a father,
mother, one sibling, both maternal and paternal grandparents.

``` r
construct_covmat()
```

    ##         g     o    m    f    s1   mgm   mgf   pgm   pgf
    ## g   0.500 0.500 0.25 0.25 0.250 0.125 0.125 0.125 0.125
    ## o   0.500 1.000 0.25 0.25 0.250 0.125 0.125 0.125 0.125
    ## m   0.250 0.250 1.00 0.00 0.250 0.250 0.250 0.000 0.000
    ## f   0.250 0.250 0.00 1.00 0.250 0.000 0.000 0.250 0.250
    ## s1  0.250 0.250 0.25 0.25 1.000 0.125 0.125 0.125 0.125
    ## mgm 0.125 0.125 0.25 0.00 0.125 1.000 0.000 0.000 0.000
    ## mgf 0.125 0.125 0.25 0.00 0.125 0.000 1.000 0.000 0.000
    ## pgm 0.125 0.125 0.00 0.25 0.125 0.000 0.000 1.000 0.000
    ## pgf 0.125 0.125 0.00 0.25 0.125 0.000 0.000 0.000 1.000
    ## attr(,"fam_vec")
    ## [1] "g"   "o"   "m"   "f"   "s1"  "mgm" "mgf" "pgm" "pgf"
    ## attr(,"n_fam")
    ## 
    ##   f   g   m mgf mgm   o pgf pgm   s 
    ##   1   1   1   1   1   1   1   1   1 
    ## attr(,"add_ind")
    ## [1] TRUE
    ## attr(,"h2")
    ## [1] 0.5

With only a mother and a father

``` r
construct_covmat(fam_vec = c("m", "f"), h2 = .5)
```

    ##      g    o    m    f
    ## g 0.50 0.50 0.25 0.25
    ## o 0.50 1.00 0.25 0.25
    ## m 0.25 0.25 1.00 0.00
    ## f 0.25 0.25 0.00 1.00
    ## attr(,"fam_vec")
    ## [1] "g" "o" "m" "f"
    ## attr(,"n_fam")
    ## 
    ## f g m o 
    ## 1 1 1 1 
    ## attr(,"add_ind")
    ## [1] TRUE
    ## attr(,"h2")
    ## [1] 0.5

In this example, we illustrate the covariance accounting for family
members on either the mother’s or father’s side. Assuming there is no
genetic overlap between the two sides of the family.

``` r
construct_covmat(fam_vec = c("f", "m", "mgm", "pgm", "mhs1", "phs1", "mau", "pau"), h2 = .5)
```

    ##          g     o    f    m   mgm   pgm  mhs1  phs1   mau   pau
    ## g    0.500 0.500 0.25 0.25 0.125 0.125 0.125 0.125 0.125 0.125
    ## o    0.500 1.000 0.25 0.25 0.125 0.125 0.125 0.125 0.125 0.125
    ## f    0.250 0.250 1.00 0.00 0.000 0.250 0.000 0.250 0.000 0.250
    ## m    0.250 0.250 0.00 1.00 0.250 0.000 0.250 0.000 0.250 0.000
    ## mgm  0.125 0.125 0.00 0.25 1.000 0.000 0.125 0.000 0.250 0.000
    ## pgm  0.125 0.125 0.25 0.00 0.000 1.000 0.000 0.125 0.000 0.250
    ## mhs1 0.125 0.125 0.00 0.25 0.125 0.000 1.000 0.000 0.125 0.000
    ## phs1 0.125 0.125 0.25 0.00 0.000 0.125 0.000 1.000 0.000 0.125
    ## mau  0.125 0.125 0.00 0.25 0.250 0.000 0.125 0.000 1.000 0.000
    ## pau  0.125 0.125 0.25 0.00 0.000 0.250 0.000 0.125 0.000 1.000
    ## attr(,"fam_vec")
    ##  [1] "g"    "o"    "f"    "m"    "mgm"  "pgm"  "mhs1" "phs1" "mau"  "pau" 
    ## attr(,"n_fam")
    ## 
    ##   f   g   m mau mgm mhs   o pau pgm phs 
    ##   1   1   1   1   1   1   1   1   1   1 
    ## attr(,"add_ind")
    ## [1] TRUE
    ## attr(,"h2")
    ## [1] 0.5
