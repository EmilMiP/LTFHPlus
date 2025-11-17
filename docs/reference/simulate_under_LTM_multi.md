# Simulate under the liability threshold model (multiple phenotypes).

`simulate_under_LTM_multi` simulates families and thresholds under the
liability threshold model for a given family structure and multiple
phenotypes. Please note that it is not possible to simulate different
family structures.

## Usage

``` r
simulate_under_LTM_multi(
  fam_vec = c("m", "f", "s1", "mgm", "mgf", "pgm", "pgf"),
  n_fam = NULL,
  add_ind = TRUE,
  genetic_corrmat = diag(3),
  full_corrmat = diag(3),
  h2_vec = rep(0.5, 3),
  phen_names = NULL,
  n_sim = 1000,
  pop_prev = rep(0.1, 3)
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

- genetic_corrmat:

  A numeric matrix holding the genetic correlations between the desired
  phenotypes. All diagonal entries must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric. Defaults to `diag(3)`.

- full_corrmat:

  A numeric matrix holding the full correlations between the desired
  phenotypes. All diagonal entries must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric. Defaults to `diag(3)`.

- h2_vec:

  A numeric vector holding the liability-scale heritabilities for a
  number of phenotype. All entries must be non-negative. Note that under
  the liability threshold model, the heritabilities must also be at
  most 1. Defaults to `rep(0.5,3)`.

- phen_names:

  A character vector holding the phenotype names. These names will be
  used to create the row and column names for the covariance matrix. If
  it is not specified, the names will default to phenotype1, phenotype2,
  etc. Defaults to `NULL`.

- n_sim:

  A positive number representing the number of simulations. Defaults to
  1000.

- pop_prev:

  A numeric vector holding the population prevalences, i.e. the overall
  prevalences in the population. All entries in `pop_prev` must be
  positive and smaller than 1. Defaults to `rep(.1,3)`.

## Value

If either `fam_vec` or `n_fam` is used as the argument and if it is of
the required format, if `genetic_corrmat` and `full_corrmat` are two
numeric and symmetric matrices satisfying that all diagonal entries are
one and that all off-diagonal entries are between -1 and 1, if the
liability-scale heritabilities in `h2_vec` are numbers satisfying \\0
\leq h^2_i\\ for all \\i \in \\1,...,n_pheno\\\\, `n_sim` is a strictly
positive number, and `pop_prev` is a positive numeric vector such that
all entries are at most one, then the output will be a list containing
lists for each phenotype. The first outer list, which is named after the
first phenotype in `phen_names`, holds the tibble `sim_obs`, which holds
the simulated liabilities, the disease status and the current
age/age-of-onset for all family members in each of the `n_sim` families
for the first phenotype. As the first outer list, the second outer list,
which is named after the second phenotype in `phen_names`, holds the
tibble `sim_obs`, which holds the simulated liabilities, the disease
status and the current age/age-of-onset for all family members in each
of the `n_sim` families for the second phenotype. There is a list
containing `sim_obs` for each phenotype in `phen_names`. The last list
entry, `thresholds`, holds the family identifier, the personal
identifier, the role (specified in fam_vec or n_fam) as well as the
lower and upper thresholds for all individuals in all families and all
phenotypes. Note that this tibble has the format required in
[`estimate_liability`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md).
Finally, note that if neither `fam_vec` nor `n_fam` are specified, the
function returns the disease status, the current age/age-of-onset, the
lower and upper thresholds, as well as the personal identifier for a
single individual, namely the individual under consideration (called
`o`). If both `fam_vec` and `n_fam` are defined, the user is asked to '
decide on which of the two vectors to use.

## See also

[`construct_covmat`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat.md)

## Examples

``` r
simulate_under_LTM_multi()
#> $phenotype1
#> $phenotype1$sim_obs
#> # A tibble: 1,000 × 26
#>    fam_ID    g_phenotype1 o_phenotype1 m_phenotype1 f_phenotype1 s1_phenotype1
#>    <chr>            <dbl>        <dbl>        <dbl>        <dbl>         <dbl>
#>  1 fam_ID_1        -1.45      -1.26           0.561       -0.109        0.444 
#>  2 fam_ID_2         0.436      0.708          0.339       -1.66         0.966 
#>  3 fam_ID_3        -0.724     -0.991         -0.181        0.115       -0.0842
#>  4 fam_ID_4        -0.687     -1.14           0.674       -0.307        0.323 
#>  5 fam_ID_5         0.379      0.101         -1.11         1.07         1.44  
#>  6 fam_ID_6        -1.03      -1.23          -0.949        0.538        0.333 
#>  7 fam_ID_7         0.820      1.02           1.20        -0.222       -1.43  
#>  8 fam_ID_8        -0.502     -1.62          -0.136       -0.408       -0.415 
#>  9 fam_ID_9         0.509      0.00612        0.270       -0.324       -0.255 
#> 10 fam_ID_10        1.54       2.15           1.70         0.477        0.367 
#> # ℹ 990 more rows
#> # ℹ 20 more variables: mgm_phenotype1 <dbl>, mgf_phenotype1 <dbl>,
#> #   pgm_phenotype1 <dbl>, pgf_phenotype1 <dbl>, o_phenotype1_status <lgl>,
#> #   m_phenotype1_status <lgl>, f_phenotype1_status <lgl>,
#> #   s1_phenotype1_status <lgl>, mgm_phenotype1_status <lgl>,
#> #   mgf_phenotype1_status <lgl>, pgm_phenotype1_status <lgl>,
#> #   pgf_phenotype1_status <lgl>, o_phenotype1_aoo <dbl>, …
#> 
#> 
#> $phenotype2
#> $phenotype2$sim_obs
#> # A tibble: 1,000 × 26
#>    fam_ID    g_phenotype2 o_phenotype2 m_phenotype2 f_phenotype2 s1_phenotype2
#>    <chr>            <dbl>        <dbl>        <dbl>        <dbl>         <dbl>
#>  1 fam_ID_1        0.834         0.488        1.32        1.88           0.673
#>  2 fam_ID_2        0.535         1.14         1.32       -0.434          1.47 
#>  3 fam_ID_3       -0.544        -1.98        -0.404       0.275         -0.556
#>  4 fam_ID_4        0.348         0.921       -0.128      -0.512          0.215
#>  5 fam_ID_5        0.330         0.817        1.09        0.709          0.492
#>  6 fam_ID_6       -1.02         -0.234        0.390      -1.30          -1.60 
#>  7 fam_ID_7       -0.0155        0.110        1.25       -0.0646         0.172
#>  8 fam_ID_8       -1.72         -2.22        -0.817      -1.72          -0.255
#>  9 fam_ID_9       -0.224        -0.613       -0.821       0.559          0.247
#> 10 fam_ID_10       1.10          1.58        -0.279      -2.28           0.212
#> # ℹ 990 more rows
#> # ℹ 20 more variables: mgm_phenotype2 <dbl>, mgf_phenotype2 <dbl>,
#> #   pgm_phenotype2 <dbl>, pgf_phenotype2 <dbl>, o_phenotype2_status <lgl>,
#> #   m_phenotype2_status <lgl>, f_phenotype2_status <lgl>,
#> #   s1_phenotype2_status <lgl>, mgm_phenotype2_status <lgl>,
#> #   mgf_phenotype2_status <lgl>, pgm_phenotype2_status <lgl>,
#> #   pgf_phenotype2_status <lgl>, o_phenotype2_aoo <dbl>, …
#> 
#> 
#> $phenotype3
#> $phenotype3$sim_obs
#> # A tibble: 1,000 × 26
#>    fam_ID    g_phenotype3 o_phenotype3 m_phenotype3 f_phenotype3 s1_phenotype3
#>    <chr>            <dbl>        <dbl>        <dbl>        <dbl>         <dbl>
#>  1 fam_ID_1        -0.436      -0.0744       -0.110       -1.79        -1.33  
#>  2 fam_ID_2         0.450      -0.725         1.06        -0.849        0.0225
#>  3 fam_ID_3        -0.512      -2.12         -1.30        -0.153       -1.49  
#>  4 fam_ID_4         0.270       0.0313        0.511       -1.11         0.0146
#>  5 fam_ID_5        -0.956      -0.952         0.715       -0.546        0.681 
#>  6 fam_ID_6         0.295       0.466         0.255        1.77         1.44  
#>  7 fam_ID_7        -0.300      -0.0237       -0.542       -0.652       -2.17  
#>  8 fam_ID_8         0.124      -0.147        -0.782        0.407        0.232 
#>  9 fam_ID_9         1.06        1.91         -0.261        1.33         0.0386
#> 10 fam_ID_10       -0.745      -0.447         1.15        -1.26        -0.139 
#> # ℹ 990 more rows
#> # ℹ 20 more variables: mgm_phenotype3 <dbl>, mgf_phenotype3 <dbl>,
#> #   pgm_phenotype3 <dbl>, pgf_phenotype3 <dbl>, o_phenotype3_status <lgl>,
#> #   m_phenotype3_status <lgl>, f_phenotype3_status <lgl>,
#> #   s1_phenotype3_status <lgl>, mgm_phenotype3_status <lgl>,
#> #   mgf_phenotype3_status <lgl>, pgm_phenotype3_status <lgl>,
#> #   pgf_phenotype3_status <lgl>, o_phenotype3_aoo <dbl>, …
#> 
#> 
#> $thresholds
#> # A tibble: 8,000 × 9
#>    fam_ID    indiv_ID   role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>     <chr>      <chr>            <dbl>            <dbl>            <dbl>
#>  1 fam_ID_1  fam_ID_1_1 o              -Inf                3.31          -Inf   
#>  2 fam_ID_2  fam_ID_2_1 o              -Inf                2.72          -Inf   
#>  3 fam_ID_3  fam_ID_3_1 o              -Inf                2.95          -Inf   
#>  4 fam_ID_4  fam_ID_4_1 o              -Inf                2.76          -Inf   
#>  5 fam_ID_5  fam_ID_5_1 o              -Inf                2.91          -Inf   
#>  6 fam_ID_6  fam_ID_6_1 o              -Inf                3.14          -Inf   
#>  7 fam_ID_7  fam_ID_7_1 o              -Inf                2.63          -Inf   
#>  8 fam_ID_8  fam_ID_8_1 o              -Inf                3.10          -Inf   
#>  9 fam_ID_9  fam_ID_9_1 o              -Inf                2.83          -Inf   
#> 10 fam_ID_10 fam_ID_10… o                 2.13             2.13             1.59
#> # ℹ 7,990 more rows
#> # ℹ 3 more variables: upper_phenotype2 <dbl>, lower_phenotype3 <dbl>,
#> #   upper_phenotype3 <dbl>
#> 

genetic_corrmat <- matrix(0.4, 3, 3)
diag(genetic_corrmat) <- 1
full_corrmat <- matrix(0.6, 3, 3)
diag(full_corrmat) <- 1

simulate_under_LTM_multi(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2,2), 
c("m","mgm","mgf","s","mhs")))
#> $phenotype1
#> $phenotype1$sim_obs
#> # A tibble: 1,000 × 26
#>    fam_ID   g_phenotype1 o_phenotype1 m_phenotype1 mgm_phenotype1 mgf_phenotype1
#>    <chr>           <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
#>  1 fam_ID_1        0.532        0.749       -1.39         -0.625          -1.22 
#>  2 fam_ID_2        0.129       -0.331        0.279         1.98            0.358
#>  3 fam_ID_3       -1.14        -1.09        -0.572         2.23           -0.133
#>  4 fam_ID_4        1.11         1.48         0.791         0.0840         -0.460
#>  5 fam_ID_5        0.116       -0.646       -1.13         -1.03           -0.391
#>  6 fam_ID_6       -0.410       -0.884        1.80          1.40            0.619
#>  7 fam_ID_7       -0.204       -0.680        0.267         0.302           0.218
#>  8 fam_ID_8        1.64         3.45         0.503         2.01           -0.253
#>  9 fam_ID_9       -0.468       -0.160        0.535         0.898          -0.724
#> 10 fam_ID_…        0.733        2.34        -1.74         -0.721           1.51 
#> # ℹ 990 more rows
#> # ℹ 20 more variables: s1_phenotype1 <dbl>, s2_phenotype1 <dbl>,
#> #   mhs1_phenotype1 <dbl>, mhs2_phenotype1 <dbl>, o_phenotype1_status <lgl>,
#> #   m_phenotype1_status <lgl>, mgm_phenotype1_status <lgl>,
#> #   mgf_phenotype1_status <lgl>, s1_phenotype1_status <lgl>,
#> #   s2_phenotype1_status <lgl>, mhs1_phenotype1_status <lgl>,
#> #   mhs2_phenotype1_status <lgl>, o_phenotype1_aoo <dbl>, …
#> 
#> 
#> $phenotype2
#> $phenotype2$sim_obs
#> # A tibble: 1,000 × 26
#>    fam_ID   g_phenotype2 o_phenotype2 m_phenotype2 mgm_phenotype2 mgf_phenotype2
#>    <chr>           <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
#>  1 fam_ID_1       0.475        1.16         -0.571         -0.297          0.615
#>  2 fam_ID_2      -0.0386      -0.0580       -1.03          -1.82           0.135
#>  3 fam_ID_3      -0.355        0.192        -0.771          1.12          -0.935
#>  4 fam_ID_4      -0.518       -0.750         1.05           1.24          -0.852
#>  5 fam_ID_5      -1.44        -1.48          1.19           0.334         -0.704
#>  6 fam_ID_6      -1.38        -2.03         -0.267         -0.270         -2.38 
#>  7 fam_ID_7       0.807       -0.101         1.56           1.83           0.298
#>  8 fam_ID_8      -0.277       -1.48         -1.79          -1.62          -0.218
#>  9 fam_ID_9      -0.914       -0.862        -0.148          0.926         -1.21 
#> 10 fam_ID_…      -1.10        -0.691        -0.749         -0.242         -2.19 
#> # ℹ 990 more rows
#> # ℹ 20 more variables: s1_phenotype2 <dbl>, s2_phenotype2 <dbl>,
#> #   mhs1_phenotype2 <dbl>, mhs2_phenotype2 <dbl>, o_phenotype2_status <lgl>,
#> #   m_phenotype2_status <lgl>, mgm_phenotype2_status <lgl>,
#> #   mgf_phenotype2_status <lgl>, s1_phenotype2_status <lgl>,
#> #   s2_phenotype2_status <lgl>, mhs1_phenotype2_status <lgl>,
#> #   mhs2_phenotype2_status <lgl>, o_phenotype2_aoo <dbl>, …
#> 
#> 
#> $phenotype3
#> $phenotype3$sim_obs
#> # A tibble: 1,000 × 26
#>    fam_ID   g_phenotype3 o_phenotype3 m_phenotype3 mgm_phenotype3 mgf_phenotype3
#>    <chr>           <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
#>  1 fam_ID_1       1.76        3.05           1.25          1.75          -1.10  
#>  2 fam_ID_2       0.364       0.0260         0.515        -0.328          1.01  
#>  3 fam_ID_3      -0.522      -0.167          1.07          0.908          0.0133
#>  4 fam_ID_4       0.712      -0.469         -0.150         0.611          0.540 
#>  5 fam_ID_5       0.358      -0.311          0.901         0.276         -1.04  
#>  6 fam_ID_6       0.0731      0.680          0.180        -0.965         -0.925 
#>  7 fam_ID_7      -0.0185      0.232          0.408        -0.0170         2.09  
#>  8 fam_ID_8      -0.367       0.00631       -0.974        -0.499         -0.132 
#>  9 fam_ID_9      -0.124      -0.973         -1.05          0.0193         0.334 
#> 10 fam_ID_…      -0.0663      0.972         -0.417        -0.858         -0.191 
#> # ℹ 990 more rows
#> # ℹ 20 more variables: s1_phenotype3 <dbl>, s2_phenotype3 <dbl>,
#> #   mhs1_phenotype3 <dbl>, mhs2_phenotype3 <dbl>, o_phenotype3_status <lgl>,
#> #   m_phenotype3_status <lgl>, mgm_phenotype3_status <lgl>,
#> #   mgf_phenotype3_status <lgl>, s1_phenotype3_status <lgl>,
#> #   s2_phenotype3_status <lgl>, mhs1_phenotype3_status <lgl>,
#> #   mhs2_phenotype3_status <lgl>, o_phenotype3_aoo <dbl>, …
#> 
#> 
#> $thresholds
#> # A tibble: 8,000 × 9
#>    fam_ID    indiv_ID   role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>     <chr>      <chr>            <dbl>            <dbl>            <dbl>
#>  1 fam_ID_1  fam_ID_1_1 o              -Inf                3.52             -Inf
#>  2 fam_ID_2  fam_ID_2_1 o              -Inf                2.59             -Inf
#>  3 fam_ID_3  fam_ID_3_1 o              -Inf                2.99             -Inf
#>  4 fam_ID_4  fam_ID_4_1 o                 1.49             1.49             -Inf
#>  5 fam_ID_5  fam_ID_5_1 o              -Inf                3.03             -Inf
#>  6 fam_ID_6  fam_ID_6_1 o              -Inf                2.95             -Inf
#>  7 fam_ID_7  fam_ID_7_1 o              -Inf                2.55             -Inf
#>  8 fam_ID_8  fam_ID_8_1 o                 3.45             3.45             -Inf
#>  9 fam_ID_9  fam_ID_9_1 o              -Inf                3.31             -Inf
#> 10 fam_ID_10 fam_ID_10… o                 2.34             2.34             -Inf
#> # ℹ 7,990 more rows
#> # ℹ 3 more variables: upper_phenotype2 <dbl>, lower_phenotype3 <dbl>,
#> #   upper_phenotype3 <dbl>
#> 

simulate_under_LTM_multi(fam_vec = c("m","f","s1"), add_ind = FALSE, 
genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, n_sim = 100)
#> $phenotype1
#> $phenotype1$sim_obs
#> # A tibble: 100 × 10
#>    fam_ID    m_phenotype1 f_phenotype1 s1_phenotype1 m_phenotype1_status
#>    <chr>            <dbl>        <dbl>         <dbl> <lgl>              
#>  1 fam_ID_1         0.606        0.586         0.380 FALSE              
#>  2 fam_ID_2         2.30         0.253         1.05  TRUE               
#>  3 fam_ID_3         1.28        -1.41          1.25  FALSE              
#>  4 fam_ID_4         0.537        0.630         0.124 FALSE              
#>  5 fam_ID_5         1.15         0.267         0.651 FALSE              
#>  6 fam_ID_6         1.19        -0.547         1.74  FALSE              
#>  7 fam_ID_7        -0.328       -0.973        -0.455 FALSE              
#>  8 fam_ID_8        -0.104        1.43          1.23  FALSE              
#>  9 fam_ID_9         0.155       -0.383        -0.767 FALSE              
#> 10 fam_ID_10       -1.57        -0.220         2.05  FALSE              
#> # ℹ 90 more rows
#> # ℹ 5 more variables: f_phenotype1_status <lgl>, s1_phenotype1_status <lgl>,
#> #   m_phenotype1_aoo <dbl>, f_phenotype1_aoo <dbl>, s1_phenotype1_aoo <dbl>
#> 
#> 
#> $phenotype2
#> $phenotype2$sim_obs
#> # A tibble: 100 × 10
#>    fam_ID    m_phenotype2 f_phenotype2 s1_phenotype2 m_phenotype2_status
#>    <chr>            <dbl>        <dbl>         <dbl> <lgl>              
#>  1 fam_ID_1       0.923          0.761        0.0222 FALSE              
#>  2 fam_ID_2       2.16           0.911        2.27   TRUE               
#>  3 fam_ID_3       1.17           0.799        0.496  FALSE              
#>  4 fam_ID_4       0.238          1.28         0.160  FALSE              
#>  5 fam_ID_5       0.653         -0.647        0.0403 FALSE              
#>  6 fam_ID_6       0.0722        -2.66         0.986  FALSE              
#>  7 fam_ID_7      -0.00295        0.988        0.751  FALSE              
#>  8 fam_ID_8      -0.797          0.689        1.15   FALSE              
#>  9 fam_ID_9      -0.818         -0.268       -0.698  FALSE              
#> 10 fam_ID_10     -1.41          -1.93         0.945  FALSE              
#> # ℹ 90 more rows
#> # ℹ 5 more variables: f_phenotype2_status <lgl>, s1_phenotype2_status <lgl>,
#> #   m_phenotype2_aoo <dbl>, f_phenotype2_aoo <dbl>, s1_phenotype2_aoo <dbl>
#> 
#> 
#> $phenotype3
#> $phenotype3$sim_obs
#> # A tibble: 100 × 10
#>    fam_ID    m_phenotype3 f_phenotype3 s1_phenotype3 m_phenotype3_status
#>    <chr>            <dbl>        <dbl>         <dbl> <lgl>              
#>  1 fam_ID_1       -0.836         0.361        -0.862 FALSE              
#>  2 fam_ID_2        0.384         1.30          1.13  FALSE              
#>  3 fam_ID_3       -0.219         1.43          0.775 FALSE              
#>  4 fam_ID_4        0.0585        0.936         1.21  FALSE              
#>  5 fam_ID_5       -0.346        -0.698         0.978 FALSE              
#>  6 fam_ID_6        2.90          0.393         2.57  TRUE               
#>  7 fam_ID_7       -2.47          0.904         0.204 FALSE              
#>  8 fam_ID_8       -1.26          0.785         0.120 FALSE              
#>  9 fam_ID_9       -0.652        -0.172        -0.473 FALSE              
#> 10 fam_ID_10       0.869        -0.238         1.94  FALSE              
#> # ℹ 90 more rows
#> # ℹ 5 more variables: f_phenotype3_status <lgl>, s1_phenotype3_status <lgl>,
#> #   m_phenotype3_aoo <dbl>, f_phenotype3_aoo <dbl>, s1_phenotype3_aoo <dbl>
#> 
#> 
#> $thresholds
#> # A tibble: 300 × 9
#>    fam_ID    indiv_ID   role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>     <chr>      <chr>            <dbl>            <dbl>            <dbl>
#>  1 fam_ID_1  fam_ID_1_1 m              -Inf                2.83          -Inf   
#>  2 fam_ID_2  fam_ID_2_1 m                 2.30             2.30             2.18
#>  3 fam_ID_3  fam_ID_3_1 m              -Inf                1.74          -Inf   
#>  4 fam_ID_4  fam_ID_4_1 m              -Inf                2.43          -Inf   
#>  5 fam_ID_5  fam_ID_5_1 m              -Inf                1.62          -Inf   
#>  6 fam_ID_6  fam_ID_6_1 m              -Inf                2.26          -Inf   
#>  7 fam_ID_7  fam_ID_7_1 m              -Inf                2.43          -Inf   
#>  8 fam_ID_8  fam_ID_8_1 m              -Inf                1.51          -Inf   
#>  9 fam_ID_9  fam_ID_9_1 m              -Inf                1.54          -Inf   
#> 10 fam_ID_10 fam_ID_10… m              -Inf                1.71          -Inf   
#> # ℹ 290 more rows
#> # ℹ 3 more variables: upper_phenotype2 <dbl>, lower_phenotype3 <dbl>,
#> #   upper_phenotype3 <dbl>
#> 

simulate_under_LTM_multi(fam_vec = c(), n_fam = NULL, add_ind = TRUE, n_sim = 150)
#> $phenotype1
#> $phenotype1$sim_obs
#> # A tibble: 150 × 5
#>    fam_ID    g_phenotype1 o_phenotype1 o_phenotype1_status o_phenotype1_aoo
#>    <chr>            <dbl>        <dbl> <lgl>                          <dbl>
#>  1 fam_ID_1        0.802       1.78    TRUE                              56
#>  2 fam_ID_2       -0.201      -0.538   FALSE                             25
#>  3 fam_ID_3       -0.452       0.795   FALSE                             33
#>  4 fam_ID_4        0.236       0.130   FALSE                             40
#>  5 fam_ID_5       -1.40       -1.36    FALSE                             15
#>  6 fam_ID_6       -0.113       0.0800  FALSE                             18
#>  7 fam_ID_7        0.0234      0.759   FALSE                             22
#>  8 fam_ID_8        0.298       1.39    TRUE                              73
#>  9 fam_ID_9        0.976       1.37    TRUE                              74
#> 10 fam_ID_10      -0.236      -0.00616 FALSE                             36
#> # ℹ 140 more rows
#> 
#> 
#> $phenotype2
#> $phenotype2$sim_obs
#> # A tibble: 150 × 5
#>    fam_ID    g_phenotype2 o_phenotype2 o_phenotype2_status o_phenotype2_aoo
#>    <chr>            <dbl>        <dbl> <lgl>                          <dbl>
#>  1 fam_ID_1         0.117       -1.22  FALSE                             24
#>  2 fam_ID_2        -1.12        -0.290 FALSE                             25
#>  3 fam_ID_3        -0.599       -0.880 FALSE                             33
#>  4 fam_ID_4         0.133       -1.00  FALSE                             40
#>  5 fam_ID_5         0.936       -0.556 FALSE                             15
#>  6 fam_ID_6         0.322        1.03  FALSE                             18
#>  7 fam_ID_7         0.227       -0.982 FALSE                             22
#>  8 fam_ID_8         1.21         1.55  TRUE                              63
#>  9 fam_ID_9        -0.288       -0.722 FALSE                             36
#> 10 fam_ID_10       -1.47        -2.21  FALSE                             36
#> # ℹ 140 more rows
#> 
#> 
#> $phenotype3
#> $phenotype3$sim_obs
#> # A tibble: 150 × 5
#>    fam_ID    g_phenotype3 o_phenotype3 o_phenotype3_status o_phenotype3_aoo
#>    <chr>            <dbl>        <dbl> <lgl>                          <dbl>
#>  1 fam_ID_1       -0.0210     -0.167   FALSE                             24
#>  2 fam_ID_2        0.0714     -0.639   FALSE                             25
#>  3 fam_ID_3        0.0743      0.575   FALSE                             33
#>  4 fam_ID_4       -0.551      -0.739   FALSE                             40
#>  5 fam_ID_5        0.209       1.99    TRUE                              50
#>  6 fam_ID_6       -0.107      -0.286   FALSE                             18
#>  7 fam_ID_7        0.210      -0.00283 FALSE                             22
#>  8 fam_ID_8       -0.975      -2.34    FALSE                             17
#>  9 fam_ID_9        0.0187     -0.440   FALSE                             36
#> 10 fam_ID_10      -0.297      -0.601   FALSE                             36
#> # ℹ 140 more rows
#> 
#> 
#> $thresholds
#> # A tibble: 150 × 9
#>    fam_ID    indiv_ID   role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>     <chr>      <chr>            <dbl>            <dbl>            <dbl>
#>  1 fam_ID_1  fam_ID_1_1 o                 1.78             1.78          -Inf   
#>  2 fam_ID_2  fam_ID_2_1 o              -Inf                3.03          -Inf   
#>  3 fam_ID_3  fam_ID_3_1 o              -Inf                2.72          -Inf   
#>  4 fam_ID_4  fam_ID_4_1 o              -Inf                2.43          -Inf   
#>  5 fam_ID_5  fam_ID_5_1 o              -Inf                3.38          -Inf   
#>  6 fam_ID_6  fam_ID_6_1 o              -Inf                3.28          -Inf   
#>  7 fam_ID_7  fam_ID_7_1 o              -Inf                3.14          -Inf   
#>  8 fam_ID_8  fam_ID_8_1 o                 1.38             1.38             1.56
#>  9 fam_ID_9  fam_ID_9_1 o                 1.37             1.37          -Inf   
#> 10 fam_ID_10 fam_ID_10… o              -Inf                2.59          -Inf   
#> # ℹ 140 more rows
#> # ℹ 3 more variables: upper_phenotype2 <dbl>, lower_phenotype3 <dbl>,
#> #   upper_phenotype3 <dbl>
#> 
```
