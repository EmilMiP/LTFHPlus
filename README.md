Jette Steinbach
16/02/2022

<!-- README.md is generated from README.Rmd. Please edit that file -->

# LTFHPpp

LTFHpp extends the package
[LTFHplus](https://emilmip.github.io/LTFHPlus/) by Emil M. Pedersen to
allow for more flexible family structures. As LTFHplus, LTFHpp can be
used to estimate an individual’s genetic component of the full
liability, the full liability or both by accounting for the family
history and population prevalences. LTFHpp inherits its efficiency from
LTFHplus, as it uses the same methods, i.e. Gibbs sampler implemented
with [Rcpp](http://www.rcpp.org/). A detailed description of the
liability threshold model conditioned on family history, age of onset
and sex can be found [here](https://doi.org/10.1016/j.ajhg.2022.01.009).

# Installation

You can install LTFHpp by:

``` r
devtools::install_github("JetteS/LTFHpp")
```
