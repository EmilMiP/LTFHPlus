
#'
#' Creates the covariance matrix needed for LT-FH++
#'
#' @param h2 Heritability estimate on liability scale to construct the covariance matrix off of.
#' @param n_sib Number of Siblings to include in the covariance matrix.
#'
#' @return None
#'
#' @examples
#' get_cov(.5)
#'
#' @export
#constructs covariance matrix with a baseling of 2 parents and n_sib siblings (with-in disorder):
get_cov = function(h2, n_sib = 0) {
  cov <- matrix(h2/2, 4 + n_sib, 4 + n_sib)
  diag(cov) <- 1
  cov[3,4] <- cov[4,3] <- 0
  cov[1:2, 1] <- cov[1, 1:2] <- h2
  cov
}

#'
#' Helps construct the correlation entries in the covariance matrix needed for LT-FH++ with multiple traits.
#'
#' @param h2_1 Heritability estimate on liability scale for the first phenotype to construct the correlation entries of the covariance matrix off of.
#' @param h2_2 Heritability estimate on liability scale for the second phenotype to construct the correlation entries of the covariance matrix off of.
#' @param rho Correlation between the two phenotypes.
#' @param n_sib Number of Siblings to include in the matrix.
#'
#' @return 
#'
#' @examples
#' get_between_trait_cov(.5, .5, .5)
#'
#' @export
#gets the correlation between two disorders (between disorder)
get_between_trait_cov <- function(h2_1, h2_2, rho, n_sib = 0) {
  corr = get_cov(1, n_sib = n_sib)
  corr * sqrt(h2_1 * h2_2) * rho
}

#'
#' Helps construct the correlation entries in the covariance matrix needed for LT-FH++ with multiple traits.
#'
#' @param corr_mat Matrix containing the heritabilities and correlations for each of the phenotypes provided. Heritabilities should be on the diagonal with the corresponing off-diagonal entry holding the correlation between two phenotypes.
#' @param n_sib Number of Siblings to include in the matrix.
#'
#' @return Returns the full covariance matrix based off of the heritabilities and corralations provided in corr_mat.
#'
#' @examples
#' get_full_cov(matrix(.5, 2, 2))
#'
#' @export
#constructs the full covariate matrix based off of the provided.
get_full_cov = function(corr_mat, n_sib = 0) {

  n_trait = nrow(corr_mat)
  block_size = 4 + n_sib
  cov_size = block_size * n_trait
  full_cov = matrix(NA, ncol = cov_size, nrow = cov_size) 
  
  for (i in 1:n_trait) {
    for (j in 1:n_trait) {
      if (i == j) {
        full_cov[(i - 1) * block_size + 1:block_size, (i - 1) * block_size + 1:block_size] = get_cov(corr_mat[i,i], n_sib = n_sib)
      } else {
        full_cov[(i - 1) * block_size + 1:block_size, (j - 1) * block_size + 1:block_size] = get_between_trait_cov(h2_1 = corr_mat[i,i],
                                                                                                                   h2_2 = corr_mat[j,j],
                                                                                                                   rho = corr_mat[i,j],
                                                                                                                   n_sib = n_sib)
      }
    }
  }
  return(full_cov)
}


#'
#' Verifies that the covariance matrix is positive definite (a properti of gaussian covariance matrices) by checking if all eigen values are positive. 
#'
#' @param full_cov The full covariance matrix for the LT-FH++ method.
#' @param corr_mat Matrix containing the heritabilities and correlations for each of the phenotypes provided. Heritabilities should be on the diagonal with the corresponing off-diagonal entry holding the correlation between two phenotypes.
#' @param correction_val value to multiply to off-diagonal entries in an attempt to force positive definite.
#' @param n_sib Number of Siblings to include in the matrix.
#'
#' @return Returns the full covariance matrix based off of the heritabilities and corralations provided in corr_mat.
#'
#' @examples
#' check_positive_definite(get_cov(.5), matrix(.5, 1, 1))
#'
#' @export
check_positive_definite = function(full_cov, corr_mat, correction_val = .99, n_sib = 0) {
  stop_ctr = 1
  ctr = 1
  cat("Covariance matrix is not positive definite (not all eigen values are positive) - Slightly reducing correlation to mitigate. \n")
  while (any(eigen(full_cov)$values < 0 ) & stop_ctr <= 100) {
    correction_mat = matrix(correction_val, ncol = ncol(corr_mat), nrow = nrow(corr_mat))
    diag(correction_mat) <- 1
    corr_mat = corr_mat * correction_mat
    full_cov = get_full_cov(corr_mat, n_sib = n_sib)
    ctr = ctr + 1
    stop_ctr = stop_ctr + 1
  }
  if (stop_ctr >= 100) {
    cat("something is wrong with the covariance matrix. Please revisit it. \n")
    print(full_cov)
    return(NULL) 
  }
  if (ctr > 1) {
    cat("shrunk off diagonal entries by:", correction_val^ctr, ".\n")
  }
  return(full_cov)
}





#'
#' Constructs the correlation matrix
#'
#' @param names Names of disorders to include
#' @param heri list, data.frame, or tibble with the heritabilities stored under the provded names
#' @param corr similar, but for correlations.
#' @param heritability_col Name of column that contains the heritabilities.
#'
#' @return Returns the correlation matrix needed as input for get_full_cov(), the LT-FH++ method and similar.
#'
#' @examples
#'
#' @export

generate_corr_matrix = function(names, heri, corr, heritability_col = "h2 (Liability Scale)") {
  n_len = length(names)
  corr_matrix = matrix(NA, ncol = n_len, nrow = n_len)
  colnames(corr_matrix) <- rownames(corr_matrix) <- names
  
  for (i in 1:n_len) {
    for (j in i:n_len) {
      if (j > n_len) break
      if (i == j) {
        corr_matrix[i,i] <- heri[[heritability_col]][heri$Phenotype == names[i]]
      } else {
        corr_matrix[i,j] <- corr_matrix[j,i] <- corr$rG[which(paste(names[i], "_", names[j], sep = "") == corr$pheno)]
      }
    }
  }
  return(corr_matrix)
}






if (FALSE) {
  heri = as_tibble(fread("D:/Work/schork_heritability.txt"))
  corr = as_tibble(fread("D:/Work/schork_genetic_correlation.txt"))
  h2_1 = heri$`h2 (Liability Scale)`[7]
  h2_2 = heri$`h2 (Liability Scale)`[3]
  gen_corr  = corr$rG[9]
  corr_mat = diag(c(h2_1, h2_2))
  corr_mat[1,2] <- corr_mat[2,1] <- gen_corr * sqrt(h2_1 * h2_2)
  
  get_full_cov(corr_mat)
  corr_mat = diag(.5, 2, 2)
  corr_mat[1,2] = corr_mat[2, 1] <- 1
  (full_cov = get_full_cov(corr_mat = corr_mat))
  check_positive_definite(full_cov = full_cov, corr_mat = corr_mat )
}
