# Constructing a covariance matrix for a single phenotype

`construct_covmatc_single` returns the covariance matrix for an
underlying target individual and a variable number of its family members

## Usage

``` r
construct_covmat_single(
  fam_vec = c("m", "f", "s1", "mgm", "mgf", "pgm", "pgf"),
  n_fam = NULL,
  add_ind = TRUE,
  h2 = 0.5
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

  - `pau[0-9]*` (Aunts/Uncles - paternal side).

- n_fam:

  A named vector holding the desired number of family members. See
  [`setNames`](https://rdrr.io/r/stats/setNames.html). All names must be
  picked from the list mentioned above. Defaults to NULL.

- add_ind:

  A logical scalar indicating whether the genetic component of the full
  liability as well as the full liability for the underlying individual
  should be included in the covariance matrix. Defaults to TRUE.

- h2:

  A number representing the squared heritability on liability scale for
  a single phenotype. Must be non-negative and at most 1. Defaults to
  0.5.

## Value

If either `fam_vec` or `n_fam` is used as the argument, if it is of the
required format and `h2` is a number satisfying \\0 \leq h2 \leq 1\\,
then the output will be a named covariance matrix. The number of rows
and columns corresponds to the length of `fam_vec` or `n_fam` (+ 2 if
`add_ind=TRUE`). If both `fam_vec = c()/NULL` and `n_fam = c()/NULL`,
the function returns a \\2 \times 2\\ matrix holding only the
correlation between the genetic component of the full liability and the
full liability for the individual. If both `fam_vec` and `n_fam` are
given, the user is asked to decide on which of the two vectors to use.
Note that the returned object has different attributes, such as
`fam_vec`, `n_fam`, `add_ind` and `h2`.

## Details

This function can be used to construct a covariance matrix for a given
number of family members. Each entry in this covariance matrix equals
the percentage of shared DNA between the corresponding individuals times
the liability-scale heritability \\h^2\\. The family members can be
specified using one of two possible formats.

## See also

[`get_relatedness`](https://emilmip.github.io/LTFHPlus/reference/get_relatedness.md),
[`construct_covmat_multi`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_multi.md),
[`construct_covmat`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat.md)

## Examples

``` r
construct_covmat_single()
#>         g     o    m    f    s1   mgm   mgf   pgm   pgf
#> g   0.500 0.500 0.25 0.25 0.250 0.125 0.125 0.125 0.125
#> o   0.500 1.000 0.25 0.25 0.250 0.125 0.125 0.125 0.125
#> m   0.250 0.250 1.00 0.00 0.250 0.250 0.250 0.000 0.000
#> f   0.250 0.250 0.00 1.00 0.250 0.000 0.000 0.250 0.250
#> s1  0.250 0.250 0.25 0.25 1.000 0.125 0.125 0.125 0.125
#> mgm 0.125 0.125 0.25 0.00 0.125 1.000 0.000 0.000 0.000
#> mgf 0.125 0.125 0.25 0.00 0.125 0.000 1.000 0.000 0.000
#> pgm 0.125 0.125 0.00 0.25 0.125 0.000 0.000 1.000 0.000
#> pgf 0.125 0.125 0.00 0.25 0.125 0.000 0.000 0.000 1.000
#> attr(,"fam_vec")
#> [1] "g"   "o"   "m"   "f"   "s1"  "mgm" "mgf" "pgm" "pgf"
#> attr(,"n_fam")
#> 
#>   f   g   m mgf mgm   o pgf pgm   s 
#>   1   1   1   1   1   1   1   1   1 
#> attr(,"add_ind")
#> [1] TRUE
#> attr(,"h2")
#> [1] 0.5
construct_covmat_single(fam_vec = c("m","mgm","mgf","mhs1","mhs2","mau1"), 
n_fam = NULL, add_ind = TRUE, h2 = 0.5)
#>          g     o    m   mgm   mgf  mhs1  mhs2  mau1
#> g    0.500 0.500 0.25 0.125 0.125 0.125 0.125 0.125
#> o    0.500 1.000 0.25 0.125 0.125 0.125 0.125 0.125
#> m    0.250 0.250 1.00 0.250 0.250 0.250 0.250 0.250
#> mgm  0.125 0.125 0.25 1.000 0.000 0.125 0.125 0.250
#> mgf  0.125 0.125 0.25 0.000 1.000 0.125 0.125 0.250
#> mhs1 0.125 0.125 0.25 0.125 0.125 1.000 0.250 0.125
#> mhs2 0.125 0.125 0.25 0.125 0.125 0.250 1.000 0.125
#> mau1 0.125 0.125 0.25 0.250 0.250 0.125 0.125 1.000
#> attr(,"fam_vec")
#> [1] "g"    "o"    "m"    "mgm"  "mgf"  "mhs1" "mhs2" "mau1"
#> attr(,"n_fam")
#> 
#>   g   m mau mgf mgm mhs   o 
#>   1   1   1   1   1   2   1 
#> attr(,"add_ind")
#> [1] TRUE
#> attr(,"h2")
#> [1] 0.5
construct_covmat_single(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2,2), 
c("m","mgm","mgf","s","mhs")), add_ind = FALSE, h2 = 0.3)
#>         m   mgm   mgf    s1    s2  mhs1  mhs2
#> m    1.00 0.150 0.150 0.150 0.150 0.150 0.150
#> mgm  0.15 1.000 0.000 0.075 0.075 0.075 0.075
#> mgf  0.15 0.000 1.000 0.075 0.075 0.075 0.075
#> s1   0.15 0.075 0.075 1.000 0.150 0.075 0.075
#> s2   0.15 0.075 0.075 0.150 1.000 0.075 0.075
#> mhs1 0.15 0.075 0.075 0.075 0.075 1.000 0.150
#> mhs2 0.15 0.075 0.075 0.075 0.075 0.150 1.000
#> attr(,"fam_vec")
#> [1] "m"    "mgm"  "mgf"  "s1"   "s2"   "mhs1" "mhs2"
#> attr(,"n_fam")
#>   m mgm mgf   s mhs 
#>   1   1   1   2   2 
#> attr(,"add_ind")
#> [1] FALSE
#> attr(,"h2")
#> [1] 0.3
```
