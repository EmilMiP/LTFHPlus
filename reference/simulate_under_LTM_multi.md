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
#>  1 fam_ID_1        -0.428     -1.34         -0.0952      -0.791          1.40 
#>  2 fam_ID_2        -0.301     -0.645         1.71        -1.28           0.258
#>  3 fam_ID_3        -1.11      -2.16         -0.869       -2.03          -1.93 
#>  4 fam_ID_4         0.424      1.10          0.922        0.871         -0.895
#>  5 fam_ID_5         0.114     -1.16          1.10        -1.62          -0.263
#>  6 fam_ID_6         0.993      1.45          1.31         0.0951        -0.787
#>  7 fam_ID_7         0.434      1.04         -0.585        2.83           1.05 
#>  8 fam_ID_8        -0.174     -0.00258      -1.58        -1.02          -0.461
#>  9 fam_ID_9         0.424      0.970         1.06        -0.963          1.69 
#> 10 fam_ID_10       -0.154      1.31         -0.958       -0.102          1.08 
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
#>  1 fam_ID_1        0.313       0.00468        0.990       -1.18         0.570 
#>  2 fam_ID_2       -0.510      -0.928          0.570       -0.349       -0.0342
#>  3 fam_ID_3        0.0908     -0.0394         1.63         0.305        2.16  
#>  4 fam_ID_4       -0.941      -1.83           0.512        0.956       -0.377 
#>  5 fam_ID_5       -0.980      -0.356         -1.78        -0.548       -1.06  
#>  6 fam_ID_6       -0.723      -0.798         -1.35        -0.793       -0.999 
#>  7 fam_ID_7       -0.475      -1.20           0.422       -1.64         0.597 
#>  8 fam_ID_8        0.655      -0.592          1.30        -1.26        -0.836 
#>  9 fam_ID_9       -1.02       -0.0343         0.926       -0.677       -0.893 
#> 10 fam_ID_10       0.0292      0.616         -1.31         1.02         0.434 
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
#>  1 fam_ID_1       -0.573       0.0225       -0.481         1.46         -0.288
#>  2 fam_ID_2       -0.104      -0.343         0.848        -0.355         0.226
#>  3 fam_ID_3        0.603      -0.165         0.287        -0.640        -1.40 
#>  4 fam_ID_4        0.538       0.910         0.830        -1.45          0.693
#>  5 fam_ID_5       -0.0690     -0.425        -0.0800        0.844         0.250
#>  6 fam_ID_6       -0.175       0.131        -0.885        -0.675        -0.207
#>  7 fam_ID_7       -0.718      -0.706        -0.377         0.495        -0.166
#>  8 fam_ID_8        0.379       0.373         0.978        -0.412         0.384
#>  9 fam_ID_9        0.233      -0.260        -0.672         0.148        -0.445
#> 10 fam_ID_10      -0.153       0.00952      -0.532        -1.25         -0.611
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
#>  1 fam_ID_1  fam_ID_1_1 o              -Inf                2.72             -Inf
#>  2 fam_ID_2  fam_ID_2_1 o              -Inf                2.95             -Inf
#>  3 fam_ID_3  fam_ID_3_1 o              -Inf                2.76             -Inf
#>  4 fam_ID_4  fam_ID_4_1 o              -Inf                2.91             -Inf
#>  5 fam_ID_5  fam_ID_5_1 o              -Inf                3.14             -Inf
#>  6 fam_ID_6  fam_ID_6_1 o                 1.45             1.45             -Inf
#>  7 fam_ID_7  fam_ID_7_1 o              -Inf                3.10             -Inf
#>  8 fam_ID_8  fam_ID_8_1 o              -Inf                2.83             -Inf
#>  9 fam_ID_9  fam_ID_9_1 o              -Inf                3.55             -Inf
#> 10 fam_ID_10 fam_ID_10… o                 1.32             1.32             -Inf
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
#>  1 fam_ID_1      0.389        -0.254       -0.667         -1.02           0.301 
#>  2 fam_ID_2     -0.404        -1.10        -0.583         -0.591         -1.01  
#>  3 fam_ID_3      0.0540       -0.0958       0.280         -0.908          1.97  
#>  4 fam_ID_4     -0.186         0.210        0.304          0.655          0.757 
#>  5 fam_ID_5      0.401         0.815       -1.05          -0.924         -0.0408
#>  6 fam_ID_6      1.12          0.325       -0.761         -0.714         -0.487 
#>  7 fam_ID_7     -1.37         -1.37        -0.0452        -1.37           0.482 
#>  8 fam_ID_8      2.01          2.16         1.68           0.971          1.04  
#>  9 fam_ID_9      0.00225      -0.612       -0.857         -0.712         -0.340 
#> 10 fam_ID_…     -0.450        -1.83         0.0627        -0.0935        -0.276 
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
#>  1 fam_ID_1      1.12           1.32       -0.555           0.599         0.403 
#>  2 fam_ID_2     -0.121         -0.326      -0.0192         -0.172        -0.828 
#>  3 fam_ID_3     -0.0501        -0.440       0.367           0.543         1.46  
#>  4 fam_ID_4      0.353          0.353      -0.226          -0.376        -0.0302
#>  5 fam_ID_5      0.663          0.826       0.987          -1.09          1.11  
#>  6 fam_ID_6      0.00256        0.520       0.739          -0.805         0.732 
#>  7 fam_ID_7      1.04           1.35        0.897           0.317         0.149 
#>  8 fam_ID_8     -0.754         -0.310      -0.181          -1.56         -0.684 
#>  9 fam_ID_9     -0.271         -0.773       0.844          -0.286         3.73  
#> 10 fam_ID_…      0.827          1.84       -0.243          -0.588        -0.883 
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
#>  1 fam_ID_1     -1.47         -1.77        -1.88           -1.07         -1.56  
#>  2 fam_ID_2      0.529         0.727       -0.154          -0.311        -0.101 
#>  3 fam_ID_3      1.13          1.66         1.88            0.409         0.740 
#>  4 fam_ID_4      0.446         0.784        0.111          -0.224        -1.96  
#>  5 fam_ID_5      0.00553      -0.361        1.15            1.03         -2.09  
#>  6 fam_ID_6     -0.222         0.756       -2.41           -0.525         0.146 
#>  7 fam_ID_7      0.403        -0.0793       0.0359         -1.42         -0.833 
#>  8 fam_ID_8      0.0465        0.284       -0.289          -0.299         0.209 
#>  9 fam_ID_9     -0.330         0.0340      -0.0332          1.05         -0.375 
#> 10 fam_ID_…      0.216         0.286       -0.124          -0.442        -0.0778
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
#>  1 fam_ID_1  fam_ID_1_1 o              -Inf                2.59             1.32
#>  2 fam_ID_2  fam_ID_2_1 o              -Inf                2.99          -Inf   
#>  3 fam_ID_3  fam_ID_3_1 o              -Inf                3.10          -Inf   
#>  4 fam_ID_4  fam_ID_4_1 o              -Inf                3.03          -Inf   
#>  5 fam_ID_5  fam_ID_5_1 o              -Inf                2.95          -Inf   
#>  6 fam_ID_6  fam_ID_6_1 o              -Inf                2.55          -Inf   
#>  7 fam_ID_7  fam_ID_7_1 o              -Inf                2.99             1.35
#>  8 fam_ID_8  fam_ID_8_1 o                 2.18             2.18          -Inf   
#>  9 fam_ID_9  fam_ID_9_1 o              -Inf                3.24          -Inf   
#> 10 fam_ID_10 fam_ID_10… o              -Inf                2.68             1.85
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
#>  1 fam_ID_1        1.04        -1.28         -0.761  FALSE              
#>  2 fam_ID_2       -0.233        0.581         0.365  FALSE              
#>  3 fam_ID_3       -2.32         0.0143       -0.768  FALSE              
#>  4 fam_ID_4       -0.659       -0.0169       -0.0988 FALSE              
#>  5 fam_ID_5       -0.765       -0.0794       -1.73   FALSE              
#>  6 fam_ID_6        1.38         0.0513       -0.288  TRUE               
#>  7 fam_ID_7        2.29        -0.950         0.468  TRUE               
#>  8 fam_ID_8       -0.0787       0.114         1.19   FALSE              
#>  9 fam_ID_9       -0.597        1.08          0.319  FALSE              
#> 10 fam_ID_10       0.812       -0.372         0.463  FALSE              
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
#>  1 fam_ID_1        -0.665      -0.928         -0.751 FALSE              
#>  2 fam_ID_2        -0.587       0.732         -0.109 FALSE              
#>  3 fam_ID_3         1.08       -0.769         -0.565 FALSE              
#>  4 fam_ID_4        -1.52        0.369          0.487 FALSE              
#>  5 fam_ID_5        -0.788       0.773         -0.322 FALSE              
#>  6 fam_ID_6        -1.01       -0.304         -1.47  FALSE              
#>  7 fam_ID_7         0.426      -0.568          0.251 FALSE              
#>  8 fam_ID_8         0.418       0.733          1.17  FALSE              
#>  9 fam_ID_9        -0.419       0.934         -0.139 FALSE              
#> 10 fam_ID_10       -0.314      -0.0139         0.273 FALSE              
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
#>  1 fam_ID_1        0.816        -1.32         -1.25  FALSE              
#>  2 fam_ID_2       -0.326         0.230        -0.565 FALSE              
#>  3 fam_ID_3       -0.900        -0.895        -1.03  FALSE              
#>  4 fam_ID_4       -0.689         0.190        -1.16  FALSE              
#>  5 fam_ID_5       -2.11         -0.136        -1.80  FALSE              
#>  6 fam_ID_6        0.155        -0.248        -1.02  FALSE              
#>  7 fam_ID_7       -0.0697       -1.64         -0.248 FALSE              
#>  8 fam_ID_8        0.168         1.24          0.930 FALSE              
#>  9 fam_ID_9       -0.275         0.337        -0.345 FALSE              
#> 10 fam_ID_10       0.190        -0.236        -0.617 FALSE              
#> # ℹ 90 more rows
#> # ℹ 5 more variables: f_phenotype3_status <lgl>, s1_phenotype3_status <lgl>,
#> #   m_phenotype3_aoo <dbl>, f_phenotype3_aoo <dbl>, s1_phenotype3_aoo <dbl>
#> 
#> 
#> $thresholds
#> # A tibble: 300 × 9
#>    fam_ID    indiv_ID   role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>     <chr>      <chr>            <dbl>            <dbl>            <dbl>
#>  1 fam_ID_1  fam_ID_1_1 m              -Inf                2.59             -Inf
#>  2 fam_ID_2  fam_ID_2_1 m              -Inf                1.47             -Inf
#>  3 fam_ID_3  fam_ID_3_1 m              -Inf                2.43             -Inf
#>  4 fam_ID_4  fam_ID_4_1 m              -Inf                1.74             -Inf
#>  5 fam_ID_5  fam_ID_5_1 m              -Inf                2.09             -Inf
#>  6 fam_ID_6  fam_ID_6_1 m                 1.38             1.38             -Inf
#>  7 fam_ID_7  fam_ID_7_1 m                 2.30             2.30             -Inf
#>  8 fam_ID_8  fam_ID_8_1 m              -Inf                1.56             -Inf
#>  9 fam_ID_9  fam_ID_9_1 m              -Inf                1.47             -Inf
#> 10 fam_ID_10 fam_ID_10… m              -Inf                1.71             -Inf
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
