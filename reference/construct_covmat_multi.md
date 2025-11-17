# Constructing a covariance matrix for multiple phenotypes

`construct_covmat_multi` returns the covariance matrix for an underlying
target individual and a variable number of its family members for
multiple phenotypes.

## Usage

``` r
construct_covmat_multi(
  fam_vec = c("m", "f", "s1", "mgm", "mgf", "pgm", "pgf"),
  n_fam = NULL,
  add_ind = TRUE,
  genetic_corrmat,
  full_corrmat,
  h2_vec,
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

- genetic_corrmat:

  A numeric matrix holding the genetic correlations between the desired
  phenotypes. All diagonal entries must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric.

- full_corrmat:

  A numeric matrix holding the full correlations between the desired
  phenotypes. All diagonal entries must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric.

- h2_vec:

  A numeric vector representing the liability-scale heritabilities for
  all phenotypes. All entries in h2_vec must be non-negative and at most
  1.

- phen_names:

  A character vector holding the phenotype names. These names will be
  used to create the row and column names for the covariance matrix. If
  it is not specified, the names will default to phenotype1, phenotype2,
  etc. Defaults to NULL.

## Value

If either `fam_vec` or `n_fam` is used as the argument and if it is of
the required format, if `genetic_corrmat` and `full_corrmat` are two
numeric and symmetric matrices satisfying that all diagonal entries are
one and that all off-diagonal entries are between -1 and 1, and if
`h2_vec` is a numeric vector satisfying \\0 \leq h2_i \leq 1\\ for all
\\i \in \\1,...,n_pheno\\\\, then the output will be a named covariance
matrix. The number of rows and columns corresponds to the number of
phenotypes times the length of `fam_vec` or `n_fam` (+ 2 if
`add_ind=TRUE`). If both `fam_vec` and `n_fam` are equal to
[`c()`](https://rdrr.io/r/base/c.html) or `NULL`, the function returns a
\\(2 \times n_pheno) \times (2\times n_pheno)\\ matrix holding only the
correlation between the genetic component of the full liability and the
full liability for the underlying individual for all phenotypes. If both
`fam_vec` and `n_fam` are specified, the user is asked to decide on
which of the two vectors to use. Note that the returned object has a
number different attributes,namely `fam_vec`, `n_fam`, `add_ind`,
`genetic_corrmat`, `full_corrmat`, `h2` and `phenotype_names`.

## Details

This function can be used to construct a covariance matrix for a given
number of family members. Each entry in this covariance matrix equals
either the percentage of shared DNA between the corresponding
individuals times the liability-scale heritability \\h^2\\ or the
percentage of shared DNA between the corresponding individuals times the
correlation between the corresponding phenotypes. That is, for the same
phenotype, the covariance between all combinations of the genetic
component of the full liability and the full liability is given by
\$\$\text{Cov}\left( l_g, l_g \right) = h^2,\$\$ \$\$\text{Cov}\left(
l_g, l_o \right) = h^2,\$\$ \$\$\text{Cov}\left( l_o, l_g \right) =
h^2\$\$ and \$\$\text{Cov}\left( l_o, l_o \right) = 1.\$\$ For two
different phenotypes, the covariance is given by \$\$\text{Cov}\left(
l_g^1, l_g^2 \right) = \rho_g^{1,2},\$\$ \$\$\text{Cov}\left( l_g^1,
l_o^2 \right) = \rho_g^{1,2},\$\$ \$\$\text{Cov}\left( l_o^1, l_g^2
\right) = \rho_g^{1,2}\$\$ and \$\$\text{Cov}\left( l_o^1, l_o^2 \right)
= \rho_g^{1,2} + \rho_e^{1,2},\$\$ where \\l_g^i\\ and \\l_o^i\\ are the
genetic component of the full liability and the full liability for
phenotype \\i\\, respectively, \\\rho_g^{i,j}\\ is the genetic
correlation between phenotype \\i\\ and \\j\\ and \\\rho_e^{1,2}\\ is
the environmental correlation between phenotype \\i\\ and \\j\\. The
family members can be specified using one of two possible formats.

## See also

[`get_relatedness`](https://emilmip.github.io/LTFHPlus/reference/get_relatedness.md),
[`construct_covmat_single`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_single.md)
and
[`construct_covmat`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat.md).

## Examples

``` r
construct_covmat_multi(fam_vec = NULL, 
                       genetic_corrmat = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
                       full_corrmat = matrix(c(1, 0.55, 0.55, 1), nrow = 2),
                       h2_vec = c(0.37,0.44),
                       phen_names = c("p1","p2"))
#> Warning: 
#>  Neither fam_vec nor n_fam is specified...
#>      g_p1 o_p1 g_p2 o_p2
#> g_p1 0.37 0.37 0.50 0.50
#> o_p1 0.37 1.00 0.50 0.55
#> g_p2 0.50 0.50 0.44 0.44
#> o_p2 0.50 0.55 0.44 1.00
#> attr(,"fam_vec")
#> [1] "g" "o"
#> attr(,"n_fam")
#> g o 
#> 1 1 
#> attr(,"add_ind")
#> [1] TRUE
#> attr(,"h2")
#> [1] 0.37 0.44
#> attr(,"genetic_corrmat")
#>      [,1] [,2]
#> [1,]  1.0  0.5
#> [2,]  0.5  1.0
#> attr(,"full_corrmat")
#>      [,1] [,2]
#> [1,] 1.00 0.55
#> [2,] 0.55 1.00
#> attr(,"phenotype_names")
#> [1] "p1" "p2"
construct_covmat_multi(fam_vec = c("m","mgm","mgf","mhs1","mhs2","mau1"), 
                       n_fam = NULL, 
                       add_ind = TRUE,
                       genetic_corrmat = diag(3),
                       full_corrmat = diag(3),
                       h2_vec = c(0.8, 0.65))
#>                 g_phenotype1 o_phenotype1 m_phenotype1 mgm_phenotype1
#> g_phenotype1             0.8          0.8          0.4            0.2
#> o_phenotype1             0.8          1.0          0.4            0.2
#> m_phenotype1             0.4          0.4          1.0            0.4
#> mgm_phenotype1           0.2          0.2          0.4            1.0
#> mgf_phenotype1           0.2          0.2          0.4            0.0
#> mhs1_phenotype1          0.2          0.2          0.4            0.2
#> mhs2_phenotype1          0.2          0.2          0.4            0.2
#> mau1_phenotype1          0.2          0.2          0.4            0.4
#> g_phenotype2             0.0          0.0          0.0            0.0
#> o_phenotype2             0.0          0.0          0.0            0.0
#> m_phenotype2             0.0          0.0          0.0            0.0
#> mgm_phenotype2           0.0          0.0          0.0            0.0
#> mgf_phenotype2           0.0          0.0          0.0            0.0
#> mhs1_phenotype2          0.0          0.0          0.0            0.0
#> mhs2_phenotype2          0.0          0.0          0.0            0.0
#> mau1_phenotype2          0.0          0.0          0.0            0.0
#>                 mgf_phenotype1 mhs1_phenotype1 mhs2_phenotype1 mau1_phenotype1
#> g_phenotype1               0.2             0.2             0.2             0.2
#> o_phenotype1               0.2             0.2             0.2             0.2
#> m_phenotype1               0.4             0.4             0.4             0.4
#> mgm_phenotype1             0.0             0.2             0.2             0.4
#> mgf_phenotype1             1.0             0.2             0.2             0.4
#> mhs1_phenotype1            0.2             1.0             0.4             0.2
#> mhs2_phenotype1            0.2             0.4             1.0             0.2
#> mau1_phenotype1            0.4             0.2             0.2             1.0
#> g_phenotype2               0.0             0.0             0.0             0.0
#> o_phenotype2               0.0             0.0             0.0             0.0
#> m_phenotype2               0.0             0.0             0.0             0.0
#> mgm_phenotype2             0.0             0.0             0.0             0.0
#> mgf_phenotype2             0.0             0.0             0.0             0.0
#> mhs1_phenotype2            0.0             0.0             0.0             0.0
#> mhs2_phenotype2            0.0             0.0             0.0             0.0
#> mau1_phenotype2            0.0             0.0             0.0             0.0
#>                 g_phenotype2 o_phenotype2 m_phenotype2 mgm_phenotype2
#> g_phenotype1          0.0000       0.0000        0.000         0.0000
#> o_phenotype1          0.0000       0.0000        0.000         0.0000
#> m_phenotype1          0.0000       0.0000        0.000         0.0000
#> mgm_phenotype1        0.0000       0.0000        0.000         0.0000
#> mgf_phenotype1        0.0000       0.0000        0.000         0.0000
#> mhs1_phenotype1       0.0000       0.0000        0.000         0.0000
#> mhs2_phenotype1       0.0000       0.0000        0.000         0.0000
#> mau1_phenotype1       0.0000       0.0000        0.000         0.0000
#> g_phenotype2          0.6500       0.6500        0.325         0.1625
#> o_phenotype2          0.6500       1.0000        0.325         0.1625
#> m_phenotype2          0.3250       0.3250        1.000         0.3250
#> mgm_phenotype2        0.1625       0.1625        0.325         1.0000
#> mgf_phenotype2        0.1625       0.1625        0.325         0.0000
#> mhs1_phenotype2       0.1625       0.1625        0.325         0.1625
#> mhs2_phenotype2       0.1625       0.1625        0.325         0.1625
#> mau1_phenotype2       0.1625       0.1625        0.325         0.3250
#>                 mgf_phenotype2 mhs1_phenotype2 mhs2_phenotype2 mau1_phenotype2
#> g_phenotype1            0.0000          0.0000          0.0000          0.0000
#> o_phenotype1            0.0000          0.0000          0.0000          0.0000
#> m_phenotype1            0.0000          0.0000          0.0000          0.0000
#> mgm_phenotype1          0.0000          0.0000          0.0000          0.0000
#> mgf_phenotype1          0.0000          0.0000          0.0000          0.0000
#> mhs1_phenotype1         0.0000          0.0000          0.0000          0.0000
#> mhs2_phenotype1         0.0000          0.0000          0.0000          0.0000
#> mau1_phenotype1         0.0000          0.0000          0.0000          0.0000
#> g_phenotype2            0.1625          0.1625          0.1625          0.1625
#> o_phenotype2            0.1625          0.1625          0.1625          0.1625
#> m_phenotype2            0.3250          0.3250          0.3250          0.3250
#> mgm_phenotype2          0.0000          0.1625          0.1625          0.3250
#> mgf_phenotype2          1.0000          0.1625          0.1625          0.3250
#> mhs1_phenotype2         0.1625          1.0000          0.3250          0.1625
#> mhs2_phenotype2         0.1625          0.3250          1.0000          0.1625
#> mau1_phenotype2         0.3250          0.1625          0.1625          1.0000
#> attr(,"fam_vec")
#> [1] "g"    "o"    "m"    "mgm"  "mgf"  "mhs1" "mhs2" "mau1"
#> attr(,"n_fam")
#> 
#>   g   m mau mgf mgm mhs   o 
#>   1   1   1   1   1   2   1 
#> attr(,"add_ind")
#> [1] TRUE
#> attr(,"h2")
#> [1] 0.80 0.65
#> attr(,"genetic_corrmat")
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> attr(,"full_corrmat")
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> attr(,"phenotype_names")
#> [1] "phenotype1" "phenotype2"
construct_covmat_multi(fam_vec = NULL, 
                       n_fam = stats::setNames(c(1,1,1,2,2), c("m","mgm","mgf","s","mhs")), 
                       add_ind = FALSE,
                       genetic_corrmat = diag(2),
                       full_corrmat = diag(2),
                       h2_vec = c(0.75,0.85))
#>                 m_phenotype1 mgm_phenotype1 mgf_phenotype1 s1_phenotype1
#> m_phenotype1           1.000         0.3750         0.3750        0.3750
#> mgm_phenotype1         0.375         1.0000         0.0000        0.1875
#> mgf_phenotype1         0.375         0.0000         1.0000        0.1875
#> s1_phenotype1          0.375         0.1875         0.1875        1.0000
#> s2_phenotype1          0.375         0.1875         0.1875        0.3750
#> mhs1_phenotype1        0.375         0.1875         0.1875        0.1875
#> mhs2_phenotype1        0.375         0.1875         0.1875        0.1875
#> m_phenotype2           0.000         0.0000         0.0000        0.0000
#> mgm_phenotype2         0.000         0.0000         0.0000        0.0000
#> mgf_phenotype2         0.000         0.0000         0.0000        0.0000
#> s1_phenotype2          0.000         0.0000         0.0000        0.0000
#> s2_phenotype2          0.000         0.0000         0.0000        0.0000
#> mhs1_phenotype2        0.000         0.0000         0.0000        0.0000
#> mhs2_phenotype2        0.000         0.0000         0.0000        0.0000
#>                 s2_phenotype1 mhs1_phenotype1 mhs2_phenotype1 m_phenotype2
#> m_phenotype1           0.3750          0.3750          0.3750        0.000
#> mgm_phenotype1         0.1875          0.1875          0.1875        0.000
#> mgf_phenotype1         0.1875          0.1875          0.1875        0.000
#> s1_phenotype1          0.3750          0.1875          0.1875        0.000
#> s2_phenotype1          1.0000          0.1875          0.1875        0.000
#> mhs1_phenotype1        0.1875          1.0000          0.3750        0.000
#> mhs2_phenotype1        0.1875          0.3750          1.0000        0.000
#> m_phenotype2           0.0000          0.0000          0.0000        1.000
#> mgm_phenotype2         0.0000          0.0000          0.0000        0.425
#> mgf_phenotype2         0.0000          0.0000          0.0000        0.425
#> s1_phenotype2          0.0000          0.0000          0.0000        0.425
#> s2_phenotype2          0.0000          0.0000          0.0000        0.425
#> mhs1_phenotype2        0.0000          0.0000          0.0000        0.425
#> mhs2_phenotype2        0.0000          0.0000          0.0000        0.425
#>                 mgm_phenotype2 mgf_phenotype2 s1_phenotype2 s2_phenotype2
#> m_phenotype1            0.0000         0.0000        0.0000        0.0000
#> mgm_phenotype1          0.0000         0.0000        0.0000        0.0000
#> mgf_phenotype1          0.0000         0.0000        0.0000        0.0000
#> s1_phenotype1           0.0000         0.0000        0.0000        0.0000
#> s2_phenotype1           0.0000         0.0000        0.0000        0.0000
#> mhs1_phenotype1         0.0000         0.0000        0.0000        0.0000
#> mhs2_phenotype1         0.0000         0.0000        0.0000        0.0000
#> m_phenotype2            0.4250         0.4250        0.4250        0.4250
#> mgm_phenotype2          1.0000         0.0000        0.2125        0.2125
#> mgf_phenotype2          0.0000         1.0000        0.2125        0.2125
#> s1_phenotype2           0.2125         0.2125        1.0000        0.4250
#> s2_phenotype2           0.2125         0.2125        0.4250        1.0000
#> mhs1_phenotype2         0.2125         0.2125        0.2125        0.2125
#> mhs2_phenotype2         0.2125         0.2125        0.2125        0.2125
#>                 mhs1_phenotype2 mhs2_phenotype2
#> m_phenotype1             0.0000          0.0000
#> mgm_phenotype1           0.0000          0.0000
#> mgf_phenotype1           0.0000          0.0000
#> s1_phenotype1            0.0000          0.0000
#> s2_phenotype1            0.0000          0.0000
#> mhs1_phenotype1          0.0000          0.0000
#> mhs2_phenotype1          0.0000          0.0000
#> m_phenotype2             0.4250          0.4250
#> mgm_phenotype2           0.2125          0.2125
#> mgf_phenotype2           0.2125          0.2125
#> s1_phenotype2            0.2125          0.2125
#> s2_phenotype2            0.2125          0.2125
#> mhs1_phenotype2          1.0000          0.4250
#> mhs2_phenotype2          0.4250          1.0000
#> attr(,"fam_vec")
#> [1] "m"    "mgm"  "mgf"  "s1"   "s2"   "mhs1" "mhs2"
#> attr(,"n_fam")
#>   m mgm mgf   s mhs 
#>   1   1   1   2   2 
#> attr(,"add_ind")
#> [1] FALSE
#> attr(,"h2")
#> [1] 0.75 0.85
#> attr(,"genetic_corrmat")
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> attr(,"full_corrmat")
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> attr(,"phenotype_names")
#> [1] "phenotype1" "phenotype2"
```
