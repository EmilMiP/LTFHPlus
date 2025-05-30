% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Helper_functions.R
\name{get_kinship}
\alias{get_kinship}
\title{Construct kinship matrix from graph}
\usage{
get_kinship(fam_graph, h2, index_id = NA, add_ind = TRUE, fix_diag = TRUE)
}
\arguments{
\item{fam_graph}{graph.}

\item{h2}{heritability.}

\item{index_id}{proband id. Only used in conjuction with add_ind = TRUE.}

\item{add_ind}{add genetic liability to the kinship matrix. Defaults to true.}

\item{fix_diag}{Whether to set diagonal to 1 for all entries except for the
genetic liability.}
}
\value{
A kinship matrix.
}
\description{
construct the kinship matrix from a graph representation of a family, centered on an index person (proband).
}
\examples{
fam <- data.frame(
i = c(1, 2, 3, 4),
f = c(3, 0, 4, 0),
m = c(2, 0, 0, 0)
)

thresholds <- data.frame(
  i = c(1, 2, 3, 4),
  lower = c(-Inf, -Inf, 0.8, 0.7),
  upper = c(0.8, 0.8, 0.8, 0.7)
)

graph <- prepare_graph(fam, icol = "i", fcol = "f", mcol = "m", node_attributes = thresholds)

get_kinship(graph, h2 = 0.5, index_id = "1")
get_kinship(graph, h2 = 1, add_ind = FALSE)

}
