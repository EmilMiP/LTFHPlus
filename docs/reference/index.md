# Package index

## Estimating (Genetic) Liability

Functions used to estimate the genetic liability in a range of different
setups with and without family history and prevalence information.

- [`estimate_liability()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md)
  : Estimating the genetic or full liability for a variable number of
  phenotypes

## Helper functions

Collection of functions used as examples of how to include additional
information.

- [`prepare_LTFHPlus_input()`](https://emilmip.github.io/LTFHPlus/reference/prepare_LTFHPlus_input.md)
  :

  Prepares input for `estimate_liability`

- [`prepare_graph()`](https://emilmip.github.io/LTFHPlus/reference/prepare_graph.md)
  : Construct graph from register information

- [`graph_to_trio()`](https://emilmip.github.io/LTFHPlus/reference/graph_to_trio.md)
  : Convert from igraph to trio information

- [`attach_attributes()`](https://emilmip.github.io/LTFHPlus/reference/attach_attributes.md)
  : Attach attributes to a family graphs

- [`familywise_attach_attributes()`](https://emilmip.github.io/LTFHPlus/reference/familywise_attach_attributes.md)
  : Wrapper to attach attributes to family graphs

- [`get_family_graphs()`](https://emilmip.github.io/LTFHPlus/reference/get_family_graphs.md)
  : Automatically identify family members of degree n

- [`convert_age_to_cir()`](https://emilmip.github.io/LTFHPlus/reference/convert_age_to_cir.md)
  : Convert age to cumulative incidence rate

- [`convert_age_to_thresh()`](https://emilmip.github.io/LTFHPlus/reference/convert_age_to_thresh.md)
  : Convert age to threshold

- [`convert_cir_to_age()`](https://emilmip.github.io/LTFHPlus/reference/convert_cir_to_age.md)
  : Convert cumulative incidence rate to age

- [`convert_format()`](https://emilmip.github.io/LTFHPlus/reference/convert_format.md)
  : Attempts to convert the list entry input format to a long format

- [`convert_liability_to_aoo()`](https://emilmip.github.io/LTFHPlus/reference/convert_liability_to_aoo.md)
  : Convert liability to age of onset

- [`convert_observed_to_liability_scale()`](https://emilmip.github.io/LTFHPlus/reference/convert_observed_to_liability_scale.md)
  : Convert the heritability on the observed scale to that on the
  liability scale

- [`simulate_under_LTM()`](https://emilmip.github.io/LTFHPlus/reference/simulate_under_LTM.md)
  : Simulate under the liability threshold model.

- [`truncated_normal_cdf()`](https://emilmip.github.io/LTFHPlus/reference/truncated_normal_cdf.md)
  : CDF for truncated normal distribution.

## Auxiliary functions

Functions used to define covariance matrices.

- [`rtmvnorm.gibbs()`](https://emilmip.github.io/LTFHPlus/reference/rtmvnorm.gibbs.md)
  : Gibbs Sampler for the truncated multivariate normal distribution
- [`get_kinship()`](https://emilmip.github.io/LTFHPlus/reference/get_kinship.md)
  : Construct kinship matrix from graph
- [`get_generations()`](https://emilmip.github.io/LTFHPlus/reference/get_generations.md)
  : Compute Generational Distances and Kinship Coefficients from a
  Family Graph
- [`get_relations()`](https://emilmip.github.io/LTFHPlus/reference/get_relations.md)
  : Compute and Label Pairwise Relationships Across Multiple Family
  Graphs
- [`get_relatedness()`](https://emilmip.github.io/LTFHPlus/reference/get_relatedness.md)
  : Relatedness between a pair of family members
- [`label_relatives()`](https://emilmip.github.io/LTFHPlus/reference/label_relatives.md)
  : Label Pairwise Relationships Based on Generational Distance and
  Kinship Coefficient
- [`Relation_per_proband_plot()`](https://emilmip.github.io/LTFHPlus/reference/Relation_per_proband_plot.md)
  : Plot the (Average) Number of Relatives per Proband by Relationship
  Type
- [`construct_covmat_single()`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_single.md)
  : Constructing a covariance matrix for a single phenotype
- [`construct_covmat_multi()`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat_multi.md)
  : Constructing a covariance matrix for multiple phenotypes
- [`construct_covmat()`](https://emilmip.github.io/LTFHPlus/reference/construct_covmat.md)
  : Constructing a covariance matrix for a variable number of phenotypes
- [`correct_positive_definite()`](https://emilmip.github.io/LTFHPlus/reference/correct_positive_definite.md)
  : Positive definite matrices
- [`graph_based_covariance_construction()`](https://emilmip.github.io/LTFHPlus/reference/graph_based_covariance_construction.md)
  : Constructing covariance matrix from local family graph
- [`graph_based_covariance_construction_multi()`](https://emilmip.github.io/LTFHPlus/reference/graph_based_covariance_construction_multi.md)
  : Constructing covariance matrix from local family graph for multi
  trait analysis

## Legacy & internal functions

functions that are not currently being used, but kept for
compatabilitity and internal functions used by the package.

- [`get_all_combs()`](https://emilmip.github.io/LTFHPlus/reference/get_all_combs.md)
  : construct all combinations of input vector
- [`estimate_gen_liability_ltfh()`](https://emilmip.github.io/LTFHPlus/reference/estimate_gen_liability_ltfh.md)
  : Estimate genetic liability similar to LT-FH
- [`estimate_liability_multi()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability_multi.md)
  : Estimating the genetic or full liability for multiple phenotypes
- [`estimate_liability_single()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability_single.md)
  : Estimating the genetic or full liability
- [`simulate_under_LTM_multi()`](https://emilmip.github.io/LTFHPlus/reference/simulate_under_LTM_multi.md)
  : Simulate under the liability threshold model (multiple phenotypes).
- [`simulate_under_LTM_single()`](https://emilmip.github.io/LTFHPlus/reference/simulate_under_LTM_single.md)
  : Simulate under the liability threshold model (single phenotype).
- [`fixSexCoding()`](https://emilmip.github.io/LTFHPlus/reference/fixSexCoding.md)
  : Fixing sex coding in trio info
