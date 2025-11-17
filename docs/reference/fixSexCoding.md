# Fixing sex coding in trio info

Internal function used to assist in fixing sex coding separately from id
coding type.

## Usage

``` r
fixSexCoding(x, sex_coding = TRUE, dadid, momid)
```

## Arguments

- x:

  current row to check against

- sex_coding:

  logical. Is sex coded as character?

- dadid:

  column name of father ids

- momid:

  column name of mother ids

## Value

appropriate sex coding
