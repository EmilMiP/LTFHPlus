# Attempts to convert the list entry input format to a long format

Attempts to convert the list entry input format to a long format

## Usage

``` r
convert_format(family, threshs, personal_id_col = "pid", role_col = NULL)
```

## Arguments

- family:

  a tibble with two entries, family id and personal id. personal id
  should end in "\_role", if a role column is not present.

- threshs:

  thresholds, with a personal id (without role) as well as the lower and
  upper thresholds

- personal_id_col:

  column name that holds the personal id

- role_col:

  column name that holds the role

## Value

returns a format similar to `prepare_LTFHPlus_input`, which is used by
`estimate_liability`

## Examples

``` r
family <- data.frame(
fam_id = c(1, 1, 1, 1),
pid = c(1, 2, 3, 4),
role = c("o", "m", "f", "pgf")
)

threshs <- data.frame(
  pid = c(1, 2, 3, 4),
  lower = c(-Inf, -Inf, 0.8, 0.7),
  upper = c(0.8, 0.8, 0.8, 0.7)
)

convert_format(family, threshs)
#>   fam_id pid role lower upper
#> 1      1   1    o  -Inf   0.8
#> 2      1   2    m  -Inf   0.8
#> 3      1   3    f   0.8   0.8
#> 4      1   4  pgf   0.7   0.7
```
