# Constructing a covariance matrix for a variable number of phenotypes

`construct_covmat` returns the covariance matrix for an underlying
target individual and a variable number of its family members for a
variable number of phenotypes. It is a wrapper around
[`construct_covmat_single`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_single.md)
and
[`construct_covmat_multi`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_multi.md).

## Usage

``` r
construct_covmat(
  fam_vec = c("m", "f", "s1", "mgm", "mgf", "pgm", "pgf"),
  n_fam = NULL,
  add_ind = TRUE,
  h2 = 0.5,
  genetic_corrmat = NULL,
  full_corrmat = NULL,
  phen_names = NULL
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
    c("m","f","s1","mgm","mgf","pgm","pgf").

- n_fam:

  A named vector holding the desired number of family members. See
  [`setNames`](https://rdrr.io/r/stats/setNames.html). All names must be
  picked from the list mentioned above. Defaults to NULL.

- add_ind:

  A logical scalar indicating whether the genetic component of the full
  liability as well as the full liability for the underlying individual
  should be included in the covariance matrix. Defaults to TRUE.

- h2:

  Either a number representing the heritability on liability scale for
  one single phenotype or a numeric vector representing the
  liability-scale heritabilities for a positive number of phenotypes.
  All entries in h2 must be non-negative and at most 1.

- genetic_corrmat:

  Either `NULL` or a numeric matrix holding the genetic correlations
  between the desired phenotypes. All diagonal entries must be equal to
  one, while all off-diagonal entries must be between -1 and 1. In
  addition, the matrix must be symmetric. Defaults to NULL.

- full_corrmat:

  Either `NULL` or a numeric matrix holding the full correlations
  between the desired phenotypes. All diagonal entries must be equal to
  one, while all off-diagonal entries must be between -1 and 1. In
  addition, the matrix must be symmetric. Defaults to NULL.

- phen_names:

  Either `NULL` or a character vector holding the phenotype names. These
  names will be used to create the row and column names for the
  covariance matrix. If it is not specified, the names will default to
  phenotype1, phenotype2, etc. Defaults to NULL.

## Value

If either `fam_vec` or `n_fam` is used as the argument, if it is of the
required format, if `add_ind` is a logical scalar and `h2` is a number
satisfying \$\$0 \leq h2 \leq 1\$\$, then the function
`construct_covmat` will return a named covariance matrix, which row- and
column-number corresponds to the length of `fam_vec` or `n_fam` (+ 2 if
`add_ind=TRUE`). However, if `h2` is a numeric vector satisfying \$\$0
\leq h2_i \leq 1\$\$ for all \$\$i \in \\1,...,n_pheno\\\$\$ and if
`genetic_corrmat` and `full_corrmat` are two numeric and symmetric
matrices satisfying that all diagonal entries are one and that all
off-diagonal entries are between -1 and 1, then `construct_covmat` will
return a named covariance matrix, which number of rows and columns
corresponds to the number of phenotypes times the length of `fam_vec` or
`n_fam` (+ 2 if `add_ind=TRUE`). If both `fam_vec` and `n_fam` are equal
to [`c()`](https://rdrr.io/r/base/c.html) or `NULL`, the function
returns either a \\2 \times 2\\ matrix holding only the correlation
between the genetic component of the full liability and the full
liability for the individual under consideration, or a \$\$(2 \times
n_pheno) \times (2\times n_pheno)\$\$ matrix holding the correlation
between the genetic component of the full liability and the full
liability for the underlying individual for all phenotypes. If both
`fam_vec` and `n_fam` are specified, the user is asked to decide on
which of the two vectors to use. Note that the returned object has
different attributes, such as `fam_vec`, `n_fam`, `add_ind` and `h2`.

## Details

This function can be used to construct a covariance matrix for a given
number of family members. If `h2` is a number, each entry in this
covariance matrix equals the percentage of shared DNA between the
corresponding individuals times the liability-scale heritability
\$\$h^2\$\$. However, if `h2` is a numeric vector, and genetic_corrmat
and full_corrmat are two symmetric correlation matrices, each entry
equals either the percentage of shared DNA between the corresponding
individuals times the liability-scale heritability \$\$h^2\$\$ or the
percentage of shared DNA between the corresponding individuals times the
correlation between the corresponding phenotypes. The family members can
be specified using one of two possible formats.

## See also

[`get_relatedness`](https://emilmip.github.io/LTFHPlus/reference/get_relatedness.md),
[`construct_covmat_single`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_single.md),
[`construct_covmat_multi`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_multi.md)

## Examples

``` r
construct_covmat()
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
construct_covmat(fam_vec = c("m","mgm","mgf","mhs1","mhs2","mau1"), 
                 n_fam = NULL, 
                 add_ind = TRUE, 
                 h2 = 0.5)
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
construct_covmat(fam_vec = NULL, 
                 n_fam = stats::setNames(c(1,1,1,2,2), c("m","mgm","mgf","s","mhs")), 
                 add_ind = FALSE,
                 h2 = 0.3)
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
construct_covmat(h2 = c(0.5,0.5), genetic_corrmat = matrix(c(1,0.4,0.4,1), nrow = 2),
                 full_corrmat = matrix(c(1,0.6,0.6,1), nrow = 2))
#>                g_phenotype1 o_phenotype1 m_phenotype1 f_phenotype1
#> g_phenotype1          0.500        0.500         0.25         0.25
#> o_phenotype1          0.500        1.000         0.25         0.25
#> m_phenotype1          0.250        0.250         1.00         0.00
#> f_phenotype1          0.250        0.250         0.00         1.00
#> s1_phenotype1         0.250        0.250         0.25         0.25
#> mgm_phenotype1        0.125        0.125         0.25         0.00
#> mgf_phenotype1        0.125        0.125         0.25         0.00
#> pgm_phenotype1        0.125        0.125         0.00         0.25
#> pgf_phenotype1        0.125        0.125         0.00         0.25
#> g_phenotype2          0.200        0.200         0.10         0.10
#> o_phenotype2          0.200        0.600         0.10         0.10
#> m_phenotype2          0.100        0.100         0.60         0.00
#> f_phenotype2          0.100        0.100         0.00         0.60
#> s1_phenotype2         0.100        0.100         0.10         0.10
#> mgm_phenotype2        0.050        0.050         0.10         0.00
#> mgf_phenotype2        0.050        0.050         0.10         0.00
#> pgm_phenotype2        0.050        0.050         0.00         0.10
#> pgf_phenotype2        0.050        0.050         0.00         0.10
#>                s1_phenotype1 mgm_phenotype1 mgf_phenotype1 pgm_phenotype1
#> g_phenotype1           0.250          0.125          0.125          0.125
#> o_phenotype1           0.250          0.125          0.125          0.125
#> m_phenotype1           0.250          0.250          0.250          0.000
#> f_phenotype1           0.250          0.000          0.000          0.250
#> s1_phenotype1          1.000          0.125          0.125          0.125
#> mgm_phenotype1         0.125          1.000          0.000          0.000
#> mgf_phenotype1         0.125          0.000          1.000          0.000
#> pgm_phenotype1         0.125          0.000          0.000          1.000
#> pgf_phenotype1         0.125          0.000          0.000          0.000
#> g_phenotype2           0.100          0.050          0.050          0.050
#> o_phenotype2           0.100          0.050          0.050          0.050
#> m_phenotype2           0.100          0.100          0.100          0.000
#> f_phenotype2           0.100          0.000          0.000          0.100
#> s1_phenotype2          0.600          0.050          0.050          0.050
#> mgm_phenotype2         0.050          0.600          0.000          0.000
#> mgf_phenotype2         0.050          0.000          0.600          0.000
#> pgm_phenotype2         0.050          0.000          0.000          0.600
#> pgf_phenotype2         0.050          0.000          0.000          0.000
#>                pgf_phenotype1 g_phenotype2 o_phenotype2 m_phenotype2
#> g_phenotype1            0.125        0.200        0.200         0.10
#> o_phenotype1            0.125        0.200        0.600         0.10
#> m_phenotype1            0.000        0.100        0.100         0.60
#> f_phenotype1            0.250        0.100        0.100         0.00
#> s1_phenotype1           0.125        0.100        0.100         0.10
#> mgm_phenotype1          0.000        0.050        0.050         0.10
#> mgf_phenotype1          0.000        0.050        0.050         0.10
#> pgm_phenotype1          0.000        0.050        0.050         0.00
#> pgf_phenotype1          1.000        0.050        0.050         0.00
#> g_phenotype2            0.050        0.500        0.500         0.25
#> o_phenotype2            0.050        0.500        1.000         0.25
#> m_phenotype2            0.000        0.250        0.250         1.00
#> f_phenotype2            0.100        0.250        0.250         0.00
#> s1_phenotype2           0.050        0.250        0.250         0.25
#> mgm_phenotype2          0.000        0.125        0.125         0.25
#> mgf_phenotype2          0.000        0.125        0.125         0.25
#> pgm_phenotype2          0.000        0.125        0.125         0.00
#> pgf_phenotype2          0.600        0.125        0.125         0.00
#>                f_phenotype2 s1_phenotype2 mgm_phenotype2 mgf_phenotype2
#> g_phenotype1           0.10         0.100          0.050          0.050
#> o_phenotype1           0.10         0.100          0.050          0.050
#> m_phenotype1           0.00         0.100          0.100          0.100
#> f_phenotype1           0.60         0.100          0.000          0.000
#> s1_phenotype1          0.10         0.600          0.050          0.050
#> mgm_phenotype1         0.00         0.050          0.600          0.000
#> mgf_phenotype1         0.00         0.050          0.000          0.600
#> pgm_phenotype1         0.10         0.050          0.000          0.000
#> pgf_phenotype1         0.10         0.050          0.000          0.000
#> g_phenotype2           0.25         0.250          0.125          0.125
#> o_phenotype2           0.25         0.250          0.125          0.125
#> m_phenotype2           0.00         0.250          0.250          0.250
#> f_phenotype2           1.00         0.250          0.000          0.000
#> s1_phenotype2          0.25         1.000          0.125          0.125
#> mgm_phenotype2         0.00         0.125          1.000          0.000
#> mgf_phenotype2         0.00         0.125          0.000          1.000
#> pgm_phenotype2         0.25         0.125          0.000          0.000
#> pgf_phenotype2         0.25         0.125          0.000          0.000
#>                pgm_phenotype2 pgf_phenotype2
#> g_phenotype1            0.050          0.050
#> o_phenotype1            0.050          0.050
#> m_phenotype1            0.000          0.000
#> f_phenotype1            0.100          0.100
#> s1_phenotype1           0.050          0.050
#> mgm_phenotype1          0.000          0.000
#> mgf_phenotype1          0.000          0.000
#> pgm_phenotype1          0.600          0.000
#> pgf_phenotype1          0.000          0.600
#> g_phenotype2            0.125          0.125
#> o_phenotype2            0.125          0.125
#> m_phenotype2            0.000          0.000
#> f_phenotype2            0.250          0.250
#> s1_phenotype2           0.125          0.125
#> mgm_phenotype2          0.000          0.000
#> mgf_phenotype2          0.000          0.000
#> pgm_phenotype2          1.000          0.000
#> pgf_phenotype2          0.000          1.000
#> attr(,"fam_vec")
#> [1] "g"   "o"   "m"   "f"   "s1"  "mgm" "mgf" "pgm" "pgf"
#> attr(,"n_fam")
#> 
#>   f   g   m mgf mgm   o pgf pgm   s 
#>   1   1   1   1   1   1   1   1   1 
#> attr(,"add_ind")
#> [1] TRUE
#> attr(,"h2")
#> [1] 0.5 0.5
#> attr(,"genetic_corrmat")
#>      [,1] [,2]
#> [1,]  1.0  0.4
#> [2,]  0.4  1.0
#> attr(,"full_corrmat")
#>      [,1] [,2]
#> [1,]  1.0  0.6
#> [2,]  0.6  1.0
#> attr(,"phenotype_names")
#> [1] "phenotype1" "phenotype2"
```
