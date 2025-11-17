# Compute Generational Distances and Kinship Coefficients from a Family Graph

Calculates generational distances and kinship coefficients between all
pairs of individuals represented in a directed family graph. The
function identifies shortest paths between individuals, accounts for
common ancestors, and derives kinship coefficients based on the number
of generations separating each pair.

## Usage

``` r
get_generations(fam_graph)
```

## Arguments

- fam_graph:

  An [igraph](https://r.igraph.org/reference/aaa-igraph-package.html)
  object representing the family structure, where directed edges
  indicate parent–child relationships (from parent to child). See
  get_family_graphs().

## Value

A tibble with one row per unique pair of individuals and the following
columns:

- id1:

  Identifier for the first individual.

- id2:

  Identifier for the second individual.

- k:

  Estimated kinship coefficient between `id1` and `id2`.

- gen.x:

  Mean number of generations separating `id1` from the common ancestor.

- gen.y:

  Mean number of generations separating `id2` from the common ancestor.

The family graph is centered on a proband (typically the same id as
fid), all relations for the proband can be found by selecting only
relations with the proband's id in the column id1.

## Examples

``` r
# see vignette on identifying and labelling relatives
```
