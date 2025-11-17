# Prepares input for `estimate_liability`

Prepares input for `estimate_liability`

## Usage

``` r
prepare_LTFHPlus_input(
  .tbl,
  CIP,
  age_col,
  aoo_col,
  CIP_merge_columns = c("sex", "birth_year", "age"),
  CIP_cip_col = "cip",
  status_col = "status",
  use_fixed_case_thr = FALSE,
  fam_id_col = "fam_id",
  personal_id_col = "pid",
  interpolation = NULL,
  bst.params = list(max_depth = 10, base_score = 0, nthread = 4, min_child_weight = 10),
  min_CIP_value = 1e-05,
  xgboost_itr = 50
)
```

## Arguments

- .tbl:

  contains family and personal ids and role with a family.

- CIP:

  tibble with population representative cumulative incidence
  proportions. CIP values should be merged by `CIP_columns`.

- age_col:

  name of column with age

- aoo_col:

  name of column with age of onset

- CIP_merge_columns:

  The columns the CIPs are subset by, e.g. CIPs by birth_year, sex.

- CIP_cip_col:

  name of column with CIP values

- status_col:

  Column that contains the status of each family member

- use_fixed_case_thr:

  Should the threshold be fixed for cases? Can be used if CIPs are
  detailed, e.g. stratified by birth_year and sex.

- fam_id_col:

  Column that contains the family ID

- personal_id_col:

  Column that contains the personal ID

- interpolation:

  type of interpolation, defaults to NULL.

- bst.params:

  list of parameters to pass on to xgboost

- min_CIP_value:

  minimum cip value to allow, too low values may lead to numerical
  instabilities.

- xgboost_itr:

  Number of iterations to run xgboost for.

## Value

tibble formatted for `estimate_liability`

## Examples

``` r
tbl = data.frame(
  fam_id = c(1, 1, 1, 1),
  pid = c(1, 2, 3, 4),
  role = c("o", "m", "f", "pgf"),
  sex = c(1, 0, 1, 1),
  status = c(0, 0, 1, 1),
  age = c(22, 42, 48, 78),
  birth_year = 2023 - c(22, 42, 48, 78),
  aoo = c(NA, NA, 43, 45))

cip = data.frame(
  age = c(22, 42, 43, 45, 48, 78),
  birth_year = c(2001, 1981, 1975, 1945, 1975, 1945),
  sex = c(1, 0, 1, 1, 1, 1),
  cip = c(0.1, 0.2, 0.3, 0.3, 0.3, 0.4))

prepare_LTFHPlus_input(.tbl = tbl,
                       CIP = cip, 
                       age_col = "age",
                       aoo_col = "aoo",
                       interpolation = NA)
#>   fam_id pid role sex status age birth_year aoo cip       thr     lower
#> 1      1   1    o   1      0  22       2001  NA 0.1 1.2815516      -Inf
#> 2      1   2    m   0      0  42       1981  NA 0.2 0.8416212      -Inf
#> 3      1   3    f   1      1  43       1975  43 0.3 0.5244005 0.5244005
#> 4      1   4  pgf   1      1  45       1945  45 0.3 0.5244005 0.5244005
#>       upper
#> 1 1.2815516
#> 2 0.8416212
#> 3       Inf
#> 4       Inf
```
