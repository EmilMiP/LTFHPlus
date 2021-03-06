% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariate_functions.R
\name{generate_cov_matrix_noFH}
\alias{generate_cov_matrix_noFH}
\title{Constructs the covariance matrix for LT-FH without Family History, and correlated traits}
\usage{
generate_cov_matrix_noFH(h2_vec, gen_cor_vec)
}
\arguments{
\item{h2_vec}{a vector of heritability for traits to be considered.}

\item{gen_cor_vec}{a vector of genetic correlations between the traits considered. The order is important, and they must fit into a correlation matrix by inserting them rowwise.}
}
\value{
Returns the covariance matrix needed for LT-FH without family history, and correlated traits.
}
\description{
Constructs the covariance matrix for LT-FH without Family History, and correlated traits
}
