# LTFHPlus

Implementation of LT-FH++. Preprint for LT-FH++ can be found on [bioxriv](https://www.biorxiv.org/content/10.1101/2021.04.20.440585v1)

LT-FH++ can be used to estimate the genetic liability of an individual by accounting for family history and population prevalences. It utilises an efficient Gibbs sampler, which is implemented with the Rcpp package and its highly scaleable. 

# Installation

You can install LTFHPlus by:

```{r}
devtools::install_github("./EmilMiP/LTFHPlus")
```

<!-- badges: start -->
[![R build status](https://github.com/EmilMiP/LTFHPlus/workflows/R-CMD-check/badge.svg)](https://github.com/EmilMiP/LTFHPlus/actions)
<!-- badges: end -->
