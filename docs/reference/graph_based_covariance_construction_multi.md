# Constructing covariance matrix from local family graph for multi trait analysis

Function that constructs the genetic covariance matrix given a graph
around a proband and extracts the threshold information from the graph.

## Usage

``` r
graph_based_covariance_construction_multi(
  fam_id,
  pid,
  cur_proband_id,
  cur_family_graph,
  h2_vec,
  genetic_corrmat,
  phen_names,
  add_ind = TRUE
)
```

## Arguments

- fam_id:

  Name of column with the family ID

- pid:

  Name of column of personal ID

- cur_proband_id:

  id of proband

- cur_family_graph:

  local graph of current proband

- h2_vec:

  vector of liability scale heritabilities

- genetic_corrmat:

  matrix with genetic correlations between considered phenotypes. Must
  have same order as h2_vec.

- phen_names:

  Names of the phenotypes, as given in cur_family_graph.

- add_ind:

  whether to add genetic liability of the proband or not. Defaults to
  true.

## Value

list with three elements. The first element is temp_tbl, which contains
the id of the current proband, the family ID and the lower and upper
thresholds for all phenotypes. The second element, cov, is the
covariance matrix of the local graph centred on the current proband. The
third element is newOrder, which is the order of ids from pid and
phen_names pasted together, such that order can be enforced elsewhere
too.

## Examples

``` r
fam <- data.frame(
fam = c(1, 1, 1,1),
id = c("pid", "mom", "dad", "pgf"),
dadcol = c("dad", 0, "pgf", 0),
momcol = c("mom", 0, 0, 0))

thresholds <- data.frame(
  id = c("pid", "mom", "dad", "pgf"),
  lower_1 = c(-Inf, -Inf, 0.8, 0.7),
  upper_1 = c(0.8, 0.8, 0.8, 0.7),
  lower_2 = c(-Inf, 0.3, -Inf, 0.2),
  upper_2 = c(0.3, 0.3, 0.3, 0.2))

graph <- prepare_graph(fam, icol = "id", fcol = "dadcol", mcol = "momcol",
 node_attributes = thresholds)

ntrait <- 2
genetic_corrmat <- matrix(0.2, ncol = ntrait, nrow = ntrait)
diag(genetic_corrmat) <- 1
full_corrmat <- matrix(0.3, ncol = ntrait, nrow = ntrait)
diag(full_corrmat) <- 1
h2_vec <- rep(0.6, ntrait)

graph_based_covariance_construction_multi(fam_id = "fam",
                                          pid = "id",
                                          cur_proband_id = "pid",
                                          cur_family_graph = graph,
                                          h2_vec = h2_vec,
                                          genetic_corrmat = genetic_corrmat,
                                          phen_names = c("1", "2"))
#> $temp_tbl
#> # A tibble: 5 × 6
#>   fam   id    lower_1 lower_2 upper_1 upper_2
#>   <chr> <chr>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1 pid   pid_g  -Inf    -Inf     Inf     Inf  
#> 2 pid   pid    -Inf    -Inf       0.8     0.3
#> 3 pid   mom    -Inf       0.3     0.8     0.3
#> 4 pid   dad       0.8  -Inf       0.8     0.3
#> 5 pid   pgf       0.7     0.2     0.7     0.2
#> 
#> $cov
#>         pid_g_1 pid_1 mom_1 dad_1 pgf_1 pid_g_2 pid_2 mom_2 dad_2 pgf_2
#> pid_g_1    0.60  0.60  0.30  0.30  0.15    0.12  0.12  0.06  0.06  0.03
#> pid_1      0.60  1.00  0.30  0.30  0.15    0.12  0.12  0.06  0.06  0.03
#> mom_1      0.30  0.30  1.00  0.00  0.00    0.06  0.06  0.12  0.00  0.00
#> dad_1      0.30  0.30  0.00  1.00  0.30    0.06  0.06  0.00  0.12  0.06
#> pgf_1      0.15  0.15  0.00  0.30  1.00    0.03  0.03  0.00  0.06  0.12
#> pid_g_2    0.12  0.12  0.06  0.06  0.03    0.60  0.60  0.30  0.30  0.15
#> pid_2      0.12  0.12  0.06  0.06  0.03    0.60  1.00  0.30  0.30  0.15
#> mom_2      0.06  0.06  0.12  0.00  0.00    0.30  0.30  1.00  0.00  0.00
#> dad_2      0.06  0.06  0.00  0.12  0.06    0.30  0.30  0.00  1.00  0.30
#> pgf_2      0.03  0.03  0.00  0.06  0.12    0.15  0.15  0.00  0.30  1.00
#> 
#> $newOrder
#>  [1] "pid_g_1" "pid_1"   "mom_1"   "dad_1"   "pgf_1"   "pid_g_2" "pid_2"  
#>  [8] "mom_2"   "dad_2"   "pgf_2"  
#> 
```
