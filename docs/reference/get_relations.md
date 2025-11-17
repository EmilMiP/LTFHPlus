# Compute and Label Pairwise Relationships Across Multiple Family Graphs

Applies
[`get_generations()`](https://emilmip.github.io/LTFHPlus/reference/get_generations.md)
and
[`label_relatives()`](https://emilmip.github.io/LTFHPlus/reference/label_relatives.md)
to a tibble of family graph objects from
[`get_family_graphs()`](https://emilmip.github.io/LTFHPlus/reference/get_family_graphs.md),
returning a unified table of labelled pairwise relationships for all
individuals across the specified families.

## Usage

``` r
get_relations(
  family_graphs,
  fid = "fid",
  family_id_vec = NULL,
  fam_graph_col = "fam_graph"
)
```

## Arguments

- family_graphs:

  A tibble containing family-specific graph objects from
  [`get_family_graphs()`](https://emilmip.github.io/LTFHPlus/reference/get_family_graphs.md)
  (typically of class `igraph`). Each row should correspond to a
  distinct family, with one column containing the graph object and
  another containing the family identifier (typically the proband's id).

- fid:

  A character string specifying the name of the column in
  `family_graphs` that holds the family identifiers. Defaults to
  `"fid"`.

- family_id_vec:

  An optional character or numeric vector specifying which families to
  process. If `NULL` (default), the function will process all families
  in `family_graphs`.

- fam_graph_col:

  A character string specifying the name of the column in
  `family_graphs` that contains the family graph objects. Defaults to
  `"fam_graph"`.

## Value

A tibble containing all labelled pairwise relationships across the
specified families, with columns:

- fid:

  Family identifier (typically the proband's id).

- id1, id2:

  Identifiers for the two individuals being compared.

- gen.x, gen.y:

  Number of generations separating each individual from their most
  recent common ancestor.

- k:

  Kinship coefficient between the pair.

- lab:

  Relationship label (e.g., `"S"`, `"P"`, `"GP"`, `"1C"`, `"HNib"`).

## See also

[`get_generations()`](https://emilmip.github.io/LTFHPlus/reference/get_generations.md)
for computing generational distances and kinship coefficients,
[`label_relatives()`](https://emilmip.github.io/LTFHPlus/reference/label_relatives.md)
for labelling relationships based on generational patterns.

## Examples

``` r
# See vignette
```
