# testing construct_covmat_single
# sample random values for input
# assert diagonal
# assert symmetric
# assert positive definite?


# assert all entries between 0 and 1?
# assert row names are correct ?
# assert no relatedness between father and mothers side of family ?

random_symmetric <- function(n, min = -1, max = 1, diag = NULL) {
  x <- matrix(runif(n ^ 2, min = min, max = max), n) 
  ind <- lower.tri(x) 
  x[ind] <- t(x)[ind]
  if (!is.null(diag)) diag(x) <- diag
  return(x)
}

symmetrize <- function(mat) {
  ind <- lower.tri(mat) 
  mat[ind] <- t(mat)[ind]
  diag(mat) <- 1
  return(mat)
}

n_fam <- stats::setNames(c(rbinom(6, 1, 0.6), rbinom(6, 9, 0.15)), c("m", "f", "mgm", "mgf", "pgm", "pgf", "c", "s", "mhs", "phs", "mau", "pau"))
h2_single <- runif(1)
phen_no <- sample(2:4, 1)
h2_multi <- runif(phen_no)

full <- matrix(runif(phen_no ^ 2, min = -1, max = 1), phen_no)
genetic <- full * rbinom(phen_no ^ 2, 1, 0.9) * runif(phen_no ^ 2)

full_corrmat <- symmetrize(full)
genetic_corrmat <- symmetrize(genetic)


