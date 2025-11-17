# Wrapper to attach attributes to family graphs

This function can attach attributes to family graphs, such as lower and
upper thresholds, for each family member. This allows for personalised
thresholds and other per-family specific attributes. This function wraps
around attach_attributes to ease the process of attaching attributes to
family graphs in the standard format.

## Usage

``` r
familywise_attach_attributes(
  family_graphs,
  fam_attr,
  fam_graph_col = "fam_graph",
  attached_fam_graph_col = "masked_fam_graph",
  fid = "fid",
  pid = "pid",
  cols_to_attach = c("lower", "upper"),
  proband_cols_to_censor = NA
)
```

## Arguments

- family_graphs:

  tibble with family ids and family graphs

- fam_attr:

  tibble with attributes for each family member

- fam_graph_col:

  column name of family graphs in family_graphs. defailts to "fam_graph"

- attached_fam_graph_col:

  column name of the updated family graphs with attached attributes.
  defaults to "masked_fam_graph".

- fid:

  column name of family id. Typically contains the name of the proband
  that a family graph is centred on. defaults to "fid".

- pid:

  personal identifier for each individual in a family. Allows for
  multiple instances of the same individual across families. Defaults to
  "pid".

- cols_to_attach:

  columns to attach to the family graphs from fam_attr, typically lower
  and upper thresholds. Mixture input also requires K_i and K_pop.

- proband_cols_to_censor:

  Should proband's upper and lower thresholds be made uninformative?
  Defaults to TRUE. Used to exclude proband's information for
  prediction.

## Value

tibble with family ids and an updated family graph with attached
attributes. If lower and upper thresholds are specified, the input is
ready for estimate_liability().

## Examples

``` r
# See Vignettes.
```
