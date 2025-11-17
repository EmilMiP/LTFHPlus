# Positive definite matrices

`correct_positive_definite` verifies that a given covariance matrix is
indeed positive definite by checking that all eigenvalues are positive.
If the given covariance matrix is not positive definite,
`correct_positive_definite` tries to modify the underlying correlation
matrices genetic_corrmat and full_corrmat in order to obtain a positive
definite covariance matrix.

## Usage

``` r
correct_positive_definite(
  covmat,
  correction_val = 0.99,
  correction_limit = 100
)
```

## Arguments

- covmat:

  A symmetric and numeric matrix. If the covariance matrix should be
  corrected, it must have a number of attributes, such as
  `attr(covmat,"fam_vec")`, `attr(covmat,"n_fam")`,
  `attr(covmat,"add_ind")`, `attr(covmat,"h2")`,
  `attr(covmat,"genetic_corrmat")`, `attr(covmat,"full_corrmat")` and
  `attr(covmat,"phenotype_names")`. Any covariance matrix obtained by
  [`construct_covmat`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat.md),
  [`construct_covmat_single`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_single.md)
  or
  [`construct_covmat_multi`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_multi.md)
  will have these attributes by default.

- correction_val:

  A positive number representing the amount by which `genetic_corrmat`
  and `full_corrmat` will be changed, if some eigenvalues are
  non-positive. That is, correction_val is the number that will be
  multiplied to all off_diagonal entries in `genetic_corrmat` and
  `full_corrmat`. Defaults to 0.99.

- correction_limit:

  A positive integer representing the upper limit for the correction
  procedure. Defaults to 100.

## Value

If `covmat` is a symmetric and numeric matrix and all eigenvalues are
positive, `correct_positive_definite` simply returns `covmat`. If some
eigenvalues are not positive and `correction_val` is a positive number,
`correct_positive_definite` tries to convert `covmat` into a positive
definite matrix. If `covmat` has attributes `add_ind`, `h2`,
`genetic_corrmat`, `full_corrmat` and `phenotype_names`,
`correct_positive_definite` computes a new covariance matrix using
slightly modified correlation matrices `genetic_corrmat` and
`full_corrmat`. If the correction is performed successfully, i.e. if the
new covariance matrix is positive definite,the new covariance matrix is
returned. Otherwise, `correct_positive_definite` returns the original
covariance matrix.

## Details

This function can be used to verify that a given covariance matrix is
positive definite. It calculates all eigenvalues in order to investigate
whether they are all positive. This property is necessary for the
covariance matrix to be used as a Gaussian covariance matrix. It is
especially useful to check whether any covariance matrix obtained by
[`construct_covmat_multi`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_multi.md)
is positive definite. If the given covariance matrix is not positive
definite, `correct_positive_definite` tries to modify the underlying
correlation matrices (called `genetic_corrmat` and `full_corrmat` in
[`construct_covmat`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat.md)
or
[`construct_covmat_multi`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_multi.md))
by multiplying all off-diagonal entries in the correlation matrices by a
given number.

## See also

[`construct_covmat`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat.md),
[`construct_covmat_single`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_single.md)
and
[`construct_covmat_multi`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_multi.md).

## Examples

``` r
ntrait <- 2
genetic_corrmat <- matrix(0.6, ncol = ntrait, nrow = ntrait)
diag(genetic_corrmat) <- 1
full_corrmat <- matrix(-0.25, ncol = ntrait, nrow = ntrait)
diag(full_corrmat) <- 1
h2_vec <- rep(0.6, ntrait)
cov <- construct_covmat(fam_vec = c("m", "f"),
  genetic_corrmat = genetic_corrmat,
  h2 = h2_vec,
  full_corrmat = full_corrmat)
cov
#>              g_phenotype1 o_phenotype1 m_phenotype1 f_phenotype1 g_phenotype2
#> g_phenotype1         0.60         0.60         0.30         0.30         0.36
#> o_phenotype1         0.60         1.00         0.30         0.30         0.36
#> m_phenotype1         0.30         0.30         1.00         0.00         0.18
#> f_phenotype1         0.30         0.30         0.00         1.00         0.18
#> g_phenotype2         0.36         0.36         0.18         0.18         0.60
#> o_phenotype2         0.36        -0.25         0.18         0.18         0.60
#> m_phenotype2         0.18         0.18        -0.25         0.00         0.30
#> f_phenotype2         0.18         0.18         0.00        -0.25         0.30
#>              o_phenotype2 m_phenotype2 f_phenotype2
#> g_phenotype1         0.36         0.18         0.18
#> o_phenotype1        -0.25         0.18         0.18
#> m_phenotype1         0.18        -0.25         0.00
#> f_phenotype1         0.18         0.00        -0.25
#> g_phenotype2         0.60         0.30         0.30
#> o_phenotype2         1.00         0.30         0.30
#> m_phenotype2         0.30         1.00         0.00
#> f_phenotype2         0.30         0.00         1.00
#> attr(,"fam_vec")
#> [1] "g" "o" "m" "f"
#> attr(,"n_fam")
#> 
#> f g m o 
#> 1 1 1 1 
#> attr(,"add_ind")
#> [1] TRUE
#> attr(,"h2")
#> [1] 0.6 0.6
#> attr(,"genetic_corrmat")
#>      [,1] [,2]
#> [1,]  1.0  0.6
#> [2,]  0.6  1.0
#> attr(,"full_corrmat")
#>       [,1]  [,2]
#> [1,]  1.00 -0.25
#> [2,] -0.25  1.00
#> attr(,"phenotype_names")
#> [1] "phenotype1" "phenotype2"
correct_positive_definite(cov)
#> The specified covariance matrix is not positive definite. 
#> Trying to correct the covariance matrix...
#> The correction was performed successfully! All off-diagonal entries are corrected by0.656.
#>              g_phenotype1 o_phenotype1 m_phenotype1 f_phenotype1 g_phenotype2
#> g_phenotype1    0.6000000    0.6000000    0.3000000    0.3000000    0.2360373
#> o_phenotype1    0.6000000    1.0000000    0.3000000    0.3000000    0.2360373
#> m_phenotype1    0.3000000    0.3000000    1.0000000    0.0000000    0.1180187
#> f_phenotype1    0.3000000    0.3000000    0.0000000    1.0000000    0.1180187
#> g_phenotype2    0.2360373    0.2360373    0.1180187    0.1180187    0.6000000
#> o_phenotype2    0.2360373   -0.1639148    0.1180187    0.1180187    0.6000000
#> m_phenotype2    0.1180187    0.1180187   -0.1639148    0.0000000    0.3000000
#> f_phenotype2    0.1180187    0.1180187    0.0000000   -0.1639148    0.3000000
#>              o_phenotype2 m_phenotype2 f_phenotype2
#> g_phenotype1    0.2360373    0.1180187    0.1180187
#> o_phenotype1   -0.1639148    0.1180187    0.1180187
#> m_phenotype1    0.1180187   -0.1639148    0.0000000
#> f_phenotype1    0.1180187    0.0000000   -0.1639148
#> g_phenotype2    0.6000000    0.3000000    0.3000000
#> o_phenotype2    1.0000000    0.3000000    0.3000000
#> m_phenotype2    0.3000000    1.0000000    0.0000000
#> f_phenotype2    0.3000000    0.0000000    1.0000000
#> attr(,"fam_vec")
#> [1] "g" "o" "m" "f"
#> attr(,"n_fam")
#> 
#> f g m o 
#> 1 1 1 1 
#> attr(,"add_ind")
#> [1] TRUE
#> attr(,"h2")
#> [1] 0.6 0.6
#> attr(,"genetic_corrmat")
#>           [,1]      [,2]
#> [1,] 1.0000000 0.3933955
#> [2,] 0.3933955 1.0000000
#> attr(,"full_corrmat")
#>            [,1]       [,2]
#> [1,]  1.0000000 -0.1639148
#> [2,] -0.1639148  1.0000000
#> attr(,"phenotype_names")
#> [1] "phenotype1" "phenotype2"
```
