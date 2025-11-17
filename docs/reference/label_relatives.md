# Label Pairwise Relationships Based on Generational Distance and Kinship Coefficient

Assigns standard pedigree relationship labels (e.g., *Parent*, *Child*,
*Sibling*, *Grandparent*, *Cousin*) to all pairs of individuals based on
their generational distances (`gen.x`, `gen.y`) and kinship coefficients
(`k`), typically produced by
[`get_generations()`](https://emilmip.github.io/LTFHPlus/reference/get_generations.md).

## Usage

``` r
label_relatives(tbl)
```

## Arguments

- tbl:

  A tibble or data frame containing at least the following columns:

  fid

  :   Column with family identifier (typically the proband's id).

  id1

  :   Identifier for the first individual.

  id2

  :   Identifier for the second individual.

  gen.x

  :   Number of generations between `id1` and their most recent common
      ancestor with `id2`.

  gen.y

  :   Number of generations between `id2` and their most recent common
      ancestor with `id1`.

  k

  :   Estimated kinship coefficient between the two individuals.

## Value

A tibble with the following columns:

- fid:

  Column with family identifier (typically the proband's id).

- id1:

  Identifier for the first individual.

- id2:

  Identifier for the second individual.

- gen.x:

  Generational distance for `id1`.

- gen.y:

  Generational distance for `id2`.

- k:

  Kinship coefficient between the two individuals.

- lab:

  Assigned relationship label (e.g., `"S"`, `"P"`, `"1C"`, `"H1C"`,
  `"2GP"`, etc.).

## Details

This function derives descriptive relationship labels using generational
differences and kinship patterns. The labels are written in a short-hand
notation, an explaination of a subset is given below:

- *P* - Parent

- *Ch* - Child

- *S* - Sibling

- *GP* - Grandparent

- *Pib* - "Pibling" (parental sibling; aunt/uncle)

- *Nib* - "Nibling" (sibling's child; niece/nephew)

- *GCh* - Grandchild

- *GPib* - Grandpibling (grandparent's sibling)

- *GNib* - Grandnibling (sibling's grandchild)

- *C* - Cousin

- *1C1R* - First Cousin Once Removed

- *2C2R* - Second Cousin Twice Removed

- *H* prefix - Half relationships (e.g., *HS* for Half-Sibling)

## See also

[`get_generations()`](https://emilmip.github.io/LTFHPlus/reference/get_generations.md)
for computing the generational and kinship inputs used by this function.

## Examples

``` r
# see vignette on identifying and labelling relatives
```
