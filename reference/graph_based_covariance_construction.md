# Constructing covariance matrix from local family graph

Function that constructs the genetic covariance matrix given a graph
around a proband and extracts the threshold information from the graph.

## Usage

``` r
graph_based_covariance_construction(
  pid,
  cur_proband_id,
  cur_family_graph,
  h2,
  add_ind = TRUE
)
```

## Arguments

- pid:

  Name of column of personal ID

- cur_proband_id:

  id of proband

- cur_family_graph:

  local graph of current proband

- h2:

  liability scale heritability

- add_ind:

  whether to add genetic liability of the proband or not. Defaults to
  true.

## Value

list with two elements. The first element is temp_tbl, which contains
the id of the current proband, the family ID and the lower and upper
thresholds. The second element, cov, is the covariance matrix of the
local graph centered on the current proband.

## Examples

``` r
fam <- data.frame(
  id = c("pid", "mom", "dad", "pgf"),
  dadcol = c("dad", 0, "pgf", 0),
  momcol = c("mom", 0, 0, 0))

thresholds <- data.frame(
  id = c("pid", "mom", "dad", "pgf"),
  lower = c(-Inf, -Inf, 0.8, 0.7),
  upper = c(0.8, 0.8, 0.8, 0.7))

graph <- prepare_graph(fam, icol = "id", fcol = "dadcol", mcol = "momcol",
 node_attributes = thresholds)

graph_based_covariance_construction(pid = "id",
                                    cur_proband_id = "pid",
                                    cur_family_graph = graph,
                                    h2 = 0.5)
#> $temp_tbl
#> # A tibble: 5 × 3
#>   id     lower upper
#>   <chr>  <dbl> <dbl>
#> 1 pid_g -Inf   Inf  
#> 2 pid   -Inf     0.8
#> 3 mom   -Inf     0.8
#> 4 dad      0.8   0.8
#> 5 pgf      0.7   0.7
#> 
#> $covmat
#>       pid_g   pid  mom  dad   pgf
#> pid_g 0.500 0.500 0.25 0.25 0.125
#> pid   0.500 1.000 0.25 0.25 0.125
#> mom   0.250 0.250 1.00 0.00 0.000
#> dad   0.250 0.250 0.00 1.00 0.250
#> pgf   0.125 0.125 0.00 0.25 1.000
#> 
```
