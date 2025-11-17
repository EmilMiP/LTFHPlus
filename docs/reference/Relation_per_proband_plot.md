# Plot the (Average) Number of Relatives per Proband by Relationship Type

Produces a structured visualisation of the (average) number of each
relationship type per proband, based on the labelled pairwise
relationship data from
[`label_relatives()`](https://emilmip.github.io/LTFHPlus/reference/label_relatives.md).
The plot arranges relationship types according to generational distance
and degree of relatedness, providing an intuitive overview of kinship
structure within the study sample.

## Usage

``` r
Relation_per_proband_plot(
  labelled_relations,
  proband_vec,
  reported_info = "both"
)
```

## Arguments

- labelled_relations:

  A tibble or data frame containing pairwise relationship labels and
  their associated metadata. Must include the following columns:

  id1

  :   Identifier for the first individual (typically the proband).

  id2

  :   Identifier for the second individual (the relative).

  gen.x, gen.y

  :   Number of generations separating each individual from their most
      recent common ancestor.

  k

  :   Kinship coefficient between the two individuals.

  lab

  :   Relationship label assigned by
      [`label_relatives()`](https://emilmip.github.io/LTFHPlus/reference/label_relatives.md).

- proband_vec:

  A vector of identifiers for probands. The function restricts the
  analysis to pairs where `id1` is included in this vector.

- reported_info:

  Chose which information is reported on the figure.

  total

  :   shows the total number of relatives of each type across all
      probands.

  average

  :   shows the average number of relatives of each type per proband.

  both

  :   (default) shows both total and average numbers on the plot.

## Value

A ggplot2 object showing the (mean) number of relatives per proband for
each relationship type. The plot can be further modified using standard
ggplot2 functions (e.g., `+ theme()` or `+ labs()`).

## Details

If any relationship types in the input are not recognised in the
predefined mapping (e.g., rare or complex kinships), these are
aggregated and shown as `"Other"`.

## See also

[`label_relatives()`](https://emilmip.github.io/LTFHPlus/reference/label_relatives.md)
for generating the relationship labels used as input.

## Examples

``` r
# see vignette on identifying and labelling relatives
```
