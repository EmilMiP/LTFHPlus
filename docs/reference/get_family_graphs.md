# Automatically identify family members of degree n

This function identifies individuals ndegree-steps away from the proband
in the population graph.

## Usage

``` r
get_family_graphs(
  pop_graph,
  ndegree,
  proband_vec,
  fid = "fid",
  fam_graph_col = "fam_graph",
  mindist = 0,
  mode = "all"
)
```

## Arguments

- pop_graph:

  Population graph from prepare_graph()

- ndegree:

  Number of steps away from proband to include

- proband_vec:

  Vector of proband ids to create family graphs for. Must be strings.

- fid:

  Column name of proband ids in the output.

- fam_graph_col:

  Column name of family graphs in the output.

- mindist:

  Minimum distance from proband to exclude in the graph (experimental,
  untested), defaults to 0, passed directly to make_neighborhood_graph.

- mode:

  Type of distance measure in the graph (experimental, untested),
  defaults to "all", passed directly to make_neighborhood_graph.

## Value

Tibble with two columns, family ids (fid) and family graphs
(fam_graph_col).

## Examples

``` r
# See Vignettes.
```
