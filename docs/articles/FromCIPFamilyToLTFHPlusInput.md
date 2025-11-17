# From CIP and family to LT-FH++ input

``` r
library(LTFHPlus)
library(dplyr)
```

### Dummy input

Here we will simply simulate a potential input format. We create `tbæ`,
which contains information on each person to attach thresholds to. It
should contain each family member along with the needed information.
Here, we simply use the proband with no family members. Next, we create
`CIP`, which contains the cumulative incidence proportions. The CIP is
stratified by birth year and sex for illustrative purposes. If users
only have CIPs stratified by sex, it would simply have one fewer
columns. **Please note that all values shown here are only for
illustrative purposes.**

``` r
n_sim = 10
tbl = tibble(
  fam_id = paste0("fam", 1:n_sim),
  pid = 1:n_sim,
  role = rep("o", n_sim),
  sex = sample(x = 0:1, size = n_sim, replace = T),
  status = sample(size = n_sim, x = 0:1, replace = T),
  age = sample(size = n_sim, x = 1:90, replace = T),
  birth_year = 2023 - age,
  aoo = purrr::map2_dbl(.x = status, .y = age, .f = ~ ifelse(.x == 1, sample(size = 1, x = 1:.y), NA))
) %>% 
  print()
```

    ## # A tibble: 10 × 8
    ##    fam_id   pid role    sex status   age birth_year   aoo
    ##    <chr>  <int> <chr> <int>  <int> <int>      <dbl> <dbl>
    ##  1 fam1       1 o         0      1    55       1968     4
    ##  2 fam2       2 o         0      0    43       1980    NA
    ##  3 fam3       3 o         1      0    62       1961    NA
    ##  4 fam4       4 o         0      0    43       1980    NA
    ##  5 fam5       5 o         0      0     5       2018    NA
    ##  6 fam6       6 o         1      1    85       1938    34
    ##  7 fam7       7 o         0      0    44       1979    NA
    ##  8 fam8       8 o         0      0    61       1962    NA
    ##  9 fam9       9 o         0      1    34       1989    25
    ## 10 fam10     10 o         1      1    70       1953    43

``` r
#### THIS IS DUMMY CIP. DO NOT USE FOR REAL-WORLD DATA USE ####
CIP = expand.grid(list(age = 1:100,
                       birth_year = 1900:2024,
                       sex = 0:1)) %>%
  group_by(sex, birth_year) %>%
  mutate(cip = (1:n() - 1)/n() * .1) %>%
  ungroup() %>% 
  print()
```

    ## # A tibble: 25,000 × 4
    ##      age birth_year   sex   cip
    ##    <int>      <int> <int> <dbl>
    ##  1     1       1900     0 0    
    ##  2     2       1900     0 0.001
    ##  3     3       1900     0 0.002
    ##  4     4       1900     0 0.003
    ##  5     5       1900     0 0.004
    ##  6     6       1900     0 0.005
    ##  7     7       1900     0 0.006
    ##  8     8       1900     0 0.007
    ##  9     9       1900     0 0.008
    ## 10    10       1900     0 0.009
    ## # ℹ 24,990 more rows

``` r
#### THIS IS DUMMY CIP. DO NOT USE FOR REAL-WORLD DATA USE ####
```

### Preparing input

Assigning thresholds to each person in `tbl` can now be done with the
function `prepare_LTFHPlus_input`. The thresholds can be assigned in two
ways. The first is matching directly on the combinations of birth year,
sex, and age of each person to the combinations that are present in the
`CIP` object. The second uses interpolation to predict the CIP value
between the observed combinations of birth year, sex, and age that is
present in the CIP object. Currently, only interpolation with
***xgboost*** package is supported. The interpolation can be useful,
since real-world data often lead to ages or age of onsets that can be
expressed as decimals and rounding may lead to large jumps in CIP
values. The outputs below can be subset such that only the required
information is left. For direct input into
[`estimate_liability()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md),
only the family and personal id columns are needed as well as role (if
the graph input is not used) and the lower and upper columns.

#### No interpolation

Without using interpolation, meaning we match on the combinations of
birth year, sex and age that are present in both the `tbl` and `CIP`
objects. The thresholds can then be assigned in the following way:

``` r
tbl2 = prepare_LTFHPlus_input(.tbl = tbl,
                              CIP = CIP, 
                              age_col = "age",
                              aoo_col = "aoo",
                              CIP_merge_columns = c("age","birth_year", "sex"),
                              CIP_cip_col = "cip",
                              status_col = "status",
                              use_fixed_case_thr = F,
                              fam_id_col = "fam_id",
                              personal_id_col = "pid",
                              interpolation = NA,
                              min_CIP_value = 1e-4)
tbl2
```

    ## # A tibble: 10 × 12
    ##    fam_id   pid role    sex status   age birth_year   aoo   cip   thr   lower
    ##    <chr>  <int> <chr> <int>  <int> <dbl>      <dbl> <dbl> <dbl> <dbl>   <dbl>
    ##  1 fam1       1 o         0      1     4       1968     4 0.003  2.75    2.75
    ##  2 fam2       2 o         0      0    43       1980    NA 0.042  1.73 -Inf   
    ##  3 fam3       3 o         1      0    62       1961    NA 0.061  1.55 -Inf   
    ##  4 fam4       4 o         0      0    43       1980    NA 0.042  1.73 -Inf   
    ##  5 fam5       5 o         0      0     5       2018    NA 0.004  2.65 -Inf   
    ##  6 fam6       6 o         1      1    34       1938    34 0.033  1.84    1.84
    ##  7 fam7       7 o         0      0    44       1979    NA 0.043  1.72 -Inf   
    ##  8 fam8       8 o         0      0    61       1962    NA 0.06   1.55 -Inf   
    ##  9 fam9       9 o         0      1    25       1989    25 0.024  1.98    1.98
    ## 10 fam10     10 o         1      1    43       1953    43 0.042  1.73    1.73
    ## # ℹ 1 more variable: upper <dbl>

#### Interpolation with Xgboost

If decimal ages and ages of onset are present, then we can interpolate
the cip values with the xgboost package. The input does not change,
except for `interpolate = "xgboost"`. The current input data does not
contain decimal values, but for illustrative purposes, we will use it
as-is. Parameters can be passed to xgboost through the `bst.params`
variable.

``` r
tbl2_xgb = prepare_LTFHPlus_input(.tbl = tbl,
                                 CIP = CIP, 
                                 age_col = "age",
                                 aoo_col = "aoo",
                                 CIP_merge_columns = c("age","birth_year", "sex"),
                                 CIP_cip_col = "cip",
                                 status_col = "status",
                                 use_fixed_case_thr = F,
                                 fam_id_col = "fam_id",
                                 personal_id_col = "pid",
                                 interpolation = "xgboost", 
                                 xgboost_itr = 30,
                                 min_CIP_value = 1e-4)
```

    ## [1]  train-rmse:0.040135 
    ## [2]  train-rmse:0.028111 
    ## [3]  train-rmse:0.019688 
    ## [4]  train-rmse:0.013790 
    ## [5]  train-rmse:0.009659 
    ## [6]  train-rmse:0.006766 
    ## [7]  train-rmse:0.004740 
    ## [8]  train-rmse:0.003321 
    ## [9]  train-rmse:0.002328 
    ## [10] train-rmse:0.001632 
    ## [11] train-rmse:0.001144 
    ## [12] train-rmse:0.000802 
    ## [13] train-rmse:0.000563 
    ## [14] train-rmse:0.000395 
    ## [15] train-rmse:0.000278 
    ## [16] train-rmse:0.000196 
    ## [17] train-rmse:0.000139 
    ## [18] train-rmse:0.000101 
    ## [19] train-rmse:0.000075 
    ## [20] train-rmse:0.000058 
    ## [21] train-rmse:0.000047 
    ## [22] train-rmse:0.000041 
    ## [23] train-rmse:0.000037 
    ## [24] train-rmse:0.000035 
    ## [25] train-rmse:0.000034 
    ## [26] train-rmse:0.000034 
    ## [27] train-rmse:0.000034 
    ## [28] train-rmse:0.000034 
    ## [29] train-rmse:0.000034 
    ## [30] train-rmse:0.000034

``` r
tbl2_xgb
```

    ## # A tibble: 10 × 13
    ##    fam_id   pid role    sex status   age birth_year   aoo cip_pred event_age
    ##    <chr>  <int> <chr> <int>  <int> <int>      <dbl> <dbl>    <dbl>     <dbl>
    ##  1 fam1       1 o         0      1    55       1968     4  0.00302         4
    ##  2 fam5       5 o         0      0     5       2018    NA  0.00401         5
    ##  3 fam9       9 o         0      1    34       1989    25  0.0240         25
    ##  4 fam6       6 o         1      1    85       1938    34  0.0330         34
    ##  5 fam2       2 o         0      0    43       1980    NA  0.0420         43
    ##  6 fam4       4 o         0      0    43       1980    NA  0.0420         43
    ##  7 fam10     10 o         1      1    70       1953    43  0.0420         43
    ##  8 fam7       7 o         0      0    44       1979    NA  0.0430         44
    ##  9 fam8       8 o         0      0    61       1962    NA  0.0600         61
    ## 10 fam3       3 o         1      0    62       1961    NA  0.0609         62
    ## # ℹ 3 more variables: thr <dbl>, lower <dbl>, upper <dbl>

### Estimating liabilities

From here, the above objects `tbl2` and `tbl2_xgb` can be subset to the
relevant columns and used in
[`estimate_liability()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md).
See ***LT-FH++ Example*** for an example of this.

The objects can also be subset to contain just the family and personal
ID columns, as well as the lower and upper columns, and then used as
input in
[`prepare_graph()`](https://emilmip.github.io/LTFHPlus/reference/prepare_graph.md)
to assign each individual with the threshold information as attributes.
See ***LT-FH++ Graph Example*** for details.
