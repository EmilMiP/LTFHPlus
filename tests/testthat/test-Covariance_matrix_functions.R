# testing construct_covmat_single
# sample random values for input
# assert diagonal
# assert symmetric
# assert positive definite?


# assert all entries between 0 and 1?
# assert row names are correct ?
# assert no relatedness between father and mothers side of family ?

SEED <- round(runif(1, 1, 9999))
#SEED <- 6374 #n_phen 2 unique_member_chance 0.3 plural_member_chance 0
set.seed(SEED)

symmetrize <- function(mat, vals) {
  mat[upper.tri(mat, diag = TRUE)] <- vals
  ind <- lower.tri(mat) 
  mat[ind] <- t(mat)[ind]
  diag(mat) <- 1
  return(mat)
}

n_phen_min <- 2
n_phen_max <- 2 #5

unique_member_chance <- 0.3 #0.6
plural_member_chance <- 0 #0.15
family_members <- c("m", "f", "mgm", "mgf", "pgm", "pgf", "c1", "s", "mhs", "phs", "mau", "pau")

n_fam <- stats::setNames(c(rbinom(6, 1, unique_member_chance), rbinom(6, 9, plural_member_chance)), family_members)
h2_single <- runif(1)

n_phen <- 2 #sample(n_phen_min : n_phen_max, 1)
h2_multi <- runif(n_phen)

n <- (n_phen * (n_phen + 1)) / 2
corrs <- rbinom(n, 1, 0.8)
fvals <- runif(n, min = -1, max = 1) * corrs
gvals <- runif(n, min = -1, max = 1) * corrs

full_corrmat <- symmetrize(matrix(nrow = n_phen, ncol = n_phen), fvals)
genetic_corrmat <- symmetrize(matrix(nrow = n_phen, ncol = n_phen), gvals)

cat("\nseed: ", SEED, "\n")
cat(n_fam, "\n", n_phen, "\n", h2_multi, "\n")
print(full_corrmat)
print(genetic_corrmat)

# testing construct_covmat_single

covmat_single <- construct_covmat_single(fam_vec = NULL, n_fam = n_fam, h2 = h2_single)

test_that("construct_covmat_single produces diagonal with h2 in first entry and 1s after", {
  expect_equal(unname(diag(covmat_single)), c(h2_single, rep(1, nrow(covmat_single) - 1)))
})

test_that("construct_covmat_single produces symmetric matrix", {
  expect_true(isSymmetric(covmat_single))
})

test_that("construct_covmat_single produces positive definite matrix", {
  expect_true(!(any(eigen(covmat_single)$values < 0)))
})

# testing construct_covmat_multi

covmat_multi <- construct_covmat_multi(fam_vec = NULL, n_fam = n_fam, h2 = h2_multi, full_corrmat = full_corrmat, genetic_corrmat = genetic_corrmat)

test_that("construct_covmat_multi produces symmetric matrix", {
  expect_true(isSymmetric(covmat_multi))
})

test_that("construct_covmat_multi produces positive definite matrix", {
  cat("\n")
  print(covmat_multi)
  print(eigen(covmat_multi)$values)
  print(det(covmat_multi))
  expect_true(!(any(eigen(covmat_multi)$values < 0)))
})

test_that("construct_covmat_multi produces some diagonal...", {
  
})