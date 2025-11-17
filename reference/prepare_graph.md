# Construct graph from register information

`prepare_graph` constructs a graph based on mother, father, and
offspring links.

## Usage

``` r
prepare_graph(
  .tbl,
  icol,
  fcol,
  mcol,
  node_attributes = NA,
  missingID_patterns = "^0$"
)
```

## Arguments

- .tbl:

  tibble with columns icol, fcol, mcol. Additional columns will be
  attributes in the constructed graph.

- icol:

  column name of column with proband ids.

- fcol:

  column name of column with father ids.

- mcol:

  column name of column with mother ids.

- node_attributes:

  tibble with icol and any additional information, such as sex, lower
  threshold, and upper threshold. Used to assign attributes to each node
  in the graph, e.g. lower and upper thresholds to individuals in the
  graph.

- missingID_patterns:

  string of missing values in the ID columns. Multiple values can be
  used, but must be separated by "\|". Defaults to "^0\$". OBS: "0" is
  NOT enough, since it relies on regex.

## Value

An igraph object. A (directed) graph object based on the links provided
in .tbl, potentially with provided attributes stored for each node.

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

prepare_graph(fam, icol = "id", fcol = "dadcol", mcol = "momcol", node_attributes = thresholds)
#> IGRAPH 40e0e29 DN-- 4 3 -- 
#> + attr: name (v/c), lower (v/n), upper (v/n)
#> + edges from 40e0e29 (vertex names):
#> [1] dad->pid mom->pid pgf->dad
```
