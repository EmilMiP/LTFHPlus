# LT-FH++ Example

First, we will load the packages that we will need for the example:

``` r
library(LTFHPlus)
library(dplyr)
library(MASS)
library(future.apply)
library(tidyr)
```

For simplicity’s sake, we will use the same prevalence for both sexes.
However, it is worth noting that in real-world applications, the
cumulative incidence proportions (CIPs) needed for the thresholds will
rarely be identical for both sexes. A common way to estimate the CIP is
with the Aalen-Johansen estimator that accounts for competing events.

### Setting up input

Single trait LT-FH++ requires 3 types of input, namely the
liability-scale heritability, the family relationships and the status,
age or age of onset, sex, and birth year for each family member, and
finally cumulative incidence proportions (CIPs). The more detailed
cumulative incidence proportions, the better. In the LT-FH++ paper, we
used CIPs that were stratified by birth year and sex and fixed the upper
and lower liability threshold for cases to $`\Phi(1 - CIP)`$, since it
provided a very accurate estimate for the full liability. Each CIP curve
was a function of age, where age represented the age of onset or the
current age of a control. If less detailed information, such as CIPs
stratified by sex, but not birth year, we recommend setting the lower
threshold to be $`\Phi(1 - CIP)`$ and the upper limit to infinity for
cases. With this information, the thresholds needed for the
age-dependent liability threshold model can be assigned. We have
provided a function called `prepare_LTFHPlus_input` that can help users
convert the input into a suitable format. An example of how the input
may look and how to use the function can be found in *From CIP and
family to LT-FH++ input*.

### Generating phenotypes

We have implemented a function that allows users to simulate under the
liability threshold model for a given family structure, i.e. mother,
father, sibling, or sibling1, sibling2, child1, child2, etc. The
function simply needs a family structure, heritability, and population
prevalence.For simplicity’s sake, we will use the same prevalence for
both sexes. However, it is worth noting that in real-world applications,
the cumulative incidence proportions needed for the thresholds, will
rarely be identical for both sexes. The function `simulate_under_LTM`
outputs a list with two tibbles. The first one `sim_obs` contains the
underlying liabilities, status, and age or age of onset for each family
and family members. The second tibble `thresholds` contains the
formatted input that is ready to be input to estimate the genetic
liability.

``` r
sims <- simulate_under_LTM(fam_vec = c("m","f","s1"), 
                           n_fam = NULL, 
                           add_ind = TRUE, 
                           h2 = h2, 
                           n_sim = 10, 
                           pop_prev = .1)
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `max_age = purrr::invoke(pmax, select(., ends_with("_age")))`.
    ## Caused by warning:
    ## ! `invoke()` was deprecated in purrr 1.0.0.
    ## ℹ Please use `exec()` instead.
    ## ℹ The deprecated feature was likely used in the LTFHPlus package.
    ##   Please report the issue at <https://github.com/EmilMiP/LTFHPlus/issues>.

``` r
sims
```

    ## $sim_obs
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
    ## 
    ## $thresholds
    ## # A tibble: 40 × 5
    ##    fam_ID    indiv_ID    role  lower upper
    ##    <chr>     <chr>       <chr> <dbl> <dbl>
    ##  1 fam_ID_1  fam_ID_1_1  o      -Inf  3.52
    ##  2 fam_ID_2  fam_ID_2_1  o      -Inf  3.38
    ##  3 fam_ID_3  fam_ID_3_1  o      -Inf  3.35
    ##  4 fam_ID_4  fam_ID_4_1  o      -Inf  3.17
    ##  5 fam_ID_5  fam_ID_5_1  o      -Inf  2.63
    ##  6 fam_ID_6  fam_ID_6_1  o      -Inf  2.55
    ##  7 fam_ID_7  fam_ID_7_1  o      -Inf  2.79
    ##  8 fam_ID_8  fam_ID_8_1  o      -Inf  3.10
    ##  9 fam_ID_9  fam_ID_9_1  o      -Inf  3.42
    ## 10 fam_ID_10 fam_ID_10_1 o      -Inf  3.35
    ## # ℹ 30 more rows

`thresholds` keeps track of families and family members with `fam_id`
and `indiv_id`. Here we simply use dummy variables for both. In
real-world data, the `indiv_id` is usually some pseudonymized identifier
unique to each individual. There is total freedom to set `fam_id` to
whatever, as long as it is unique to each family. One choice could be
the pseudonymized identifier for the index person in a family, or simply
numerate the families from $`1`$ to the total number of families. Next,
`role` identifies each family member’s relationship to the index person.
They follow a system shown in the documentation for `construct_covmat`.
If a family member does not fit one of the shown, a suitable covariance
matrix cannot be constructed. If more detailed family structures are
needed, please contact one of the maintainers of `LTFHPlus`. Finally,
the lower and upper liability thresholds are provided for each family
member. Previously, we utilised list of lists in `R` to link information
to the correct individuals. Going forward, we will still support this
format, but in order to ease the use of non-`R` users, we have opted for
a long format, where each row is an indivdual instead of a family, such
that the input can be generated with the user’s preferred software.

#### Generating your own input data

If you would like to use LTFHPlus, then the above can provide a template
for the input format. The object `sims` is a list of two entries, namely
`sim_obs` and `thresholds`. The values in `sim_obs` are values such as
the true underlying liabilities, simulated ages and onset ages. etc..
The values in `thresholds` are ready to be used in
[`estimate_liability()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md).
Generating this input from ones own register data can be tedious, as it
requires identifying parents (role for parents are $`m`$ and $`f`$),
then one can identify siblings (roles for siblings are
$`s1, s2, etc.`$), and so on for all available roles. The available
roles can be seen in the documentation for
[`estimate_liability()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md).

### Running LT-FH++

With all input parameters generated and a chosen heritability, we have
everything we need to run LT-FH++. We highly recommend the package
**Future** to parallelize the computations for LT-FH++. Since all modern
computers have more than 2 cores available to them, it should be
available to all users. If a user is not able to parallelize, then
`LTFHPlus` can still be used, but it will be far slower. The future
framework also allows for easy use on high performance computing (HPC)
clusters. In this example, we will utilize 4 threads.

``` r
#Setting up parallelization backend
plan(tweak(multisession, workers = nthreads))
#performs LT-FH++ analysis
simu_liab = estimate_liability(.tbl = sims$thresholds,
                               h2 = h2,
                               pid = "indiv_ID",
                               fam_id = "fam_ID",
                               role = "role")
```

    ## The number of workers is 4

``` r
simu_liab
```

    ## # A tibble: 10 × 4
    ##    fam_ID    indiv_ID  genetic_est genetic_se
    ##    <chr>     <chr>           <dbl>      <dbl>
    ##  1 fam_ID_1  fam_ID_1      -0.0327    0.00400
    ##  2 fam_ID_2  fam_ID_2      -0.0271    0.00423
    ##  3 fam_ID_3  fam_ID_3       0.504     0.00357
    ##  4 fam_ID_4  fam_ID_4      -0.0262    0.00449
    ##  5 fam_ID_5  fam_ID_5       0.377     0.00341
    ##  6 fam_ID_6  fam_ID_6      -0.0694    0.00393
    ##  7 fam_ID_7  fam_ID_7       0.409     0.00340
    ##  8 fam_ID_8  fam_ID_8      -0.0198    0.00422
    ##  9 fam_ID_9  fam_ID_9      -0.0188    0.00415
    ## 10 fam_ID_10 fam_ID_10     -0.0156    0.00418

From here the object `simu_liab` can be saved and used in your GWAS
method of choice.
