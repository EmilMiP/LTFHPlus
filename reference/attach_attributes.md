# Attach attributes to a family graphs

This function attaches attributes to family graphs, such as lower and
upper thresholds, for each family member. This allows for a
user-friendly way to attach personalised thresholds and other per-family
specific attributes to the family graphs.

## Usage

``` r
attach_attributes(
  cur_fam_graph,
  cur_proband,
  pid,
  attr_tbl,
  attr_names,
  proband_cols_to_censor = NA
)
```

## Arguments

- cur_fam_graph:

  An igraph object (neighbourhood graph around a proband) with family
  members up to degree n.

- cur_proband:

  Current proband id (center of the neighbourhood graph).

- pid:

  Column name of personal id (within a family).

- attr_tbl:

  Tibble with family id and attributes for each family member.

- attr_names:

  Names of attributes to be assigned to each node (family member) in the
  graph.

- proband_cols_to_censor:

  Which columns should be made uninformative for the proband? Defaults
  to NA. Used to exclude proband's information for prediction with, e.g.
  c("lower", "upper").

## Value

igraph object (neighbourhood graph around a proband) with updated
attributes for each node in the graph.
