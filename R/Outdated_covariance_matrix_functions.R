#'
#' Creates the covariance matrix needed for LT-FH++
#' 
#' This function has been retired in favor of \code{\link{construct_covmat}}. This function is kept for legacy reasons.
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
#constructs covariance matrix with a baseline of 2 parents and n_sib siblings (with-in disorder):
get_cov = function(h2, n_sib = 0) {
  
  warning("'get_cov()' was deprecated in LTFHPlus v1.0.0. It is only kept for legacy reasons.\n
          Please consider using 'construct_covmat()' instead.\
          The signature and semantics have changed, see '?construct_covmat'.")
  
  cov <- matrix(h2/2, 4 + n_sib, 4 + n_sib)
  diag(cov) <- 1
  cov[3,4] <- cov[4,3] <- 0
  cov[1:2, 1] <- cov[1, 1:2] <- h2
  cov
}




#'
#' Constructs the covariance matrix for LT-FH without Family History, and correlated traits
#' 
#' This function has been retired in favor of \code{\link{construct_covmat}}. This function is kept for legacy reasons.
#'
#' @param h2_vec a vector of heritability for traits to be considered.
#' @param gen_cor_vec a vector of genetic correlations between the traits considered. The order is important, and they must fit into a correlation matrix by inserting them rowwise.
#'
#' @return Returns the covariance matrix needed for LT-FH without family history, and correlated traits.
#'
#' @export

generate_cov_matrix_noFH = function(h2_vec, gen_cor_vec) {
  
  warning("'generate_cov_matrix_noFH()' was deprecated in LTFHPlus v1.0.0. It is only kept for legacy reasons.\n
          Please consider using 'construct_covmat()' instead.\
          The signature and semantics have changed, see '?construct_covmat'.")
  
  ntraits = length(h2_vec)
  
  #generates initial matrix
  cov_mat = diag(rep(1, ntraits))
  #fills in genetic correlations
  cov_mat[upper.tri(cov_mat)] <- cov_mat[lower.tri(cov_mat)] <- gen_cor_vec
  #use outer product to generate matrix of heritability products
  heri_mat = tcrossprod(sqrt(h2_vec))
  #chaning diagonal to 1, so we can do elementwise multiplication
  diag(heri_mat) = 1
  
  #covariance matrix done for full heritabilities
  cov_mat = cov_mat * heri_mat
  
  #adding the genetic liability for the primary trait
  cov_mat = cbind(0,rbind(0,cov_mat))
  cov_mat[1,-(1:2)] <- cov_mat[-(1:2), 1] <- cov_mat[2,-(1:2)]
  cov_mat[1,1:2] <- cov_mat[2,1] <- h2_vec[1]
  
  return(cov_mat)
}



#'
#' Helps construct the correlation entries in the covariance matrix needed for LT-FH++ with multiple traits.
#' 
#' This function has been retired in favor of \code{\link{construct_covmat}}. This function is kept for legacy reasons.
#'
#' @param h2_1 Heritability estimate on liability scale for the first phenotype to construct the correlation entries of the covariance matrix off of.
#' @param h2_2 Heritability estimate on liability scale for the second phenotype to construct the correlation entries of the covariance matrix off of.
#' @param rho Correlation between the two phenotypes.
#' @param n_sib Number of Siblings to include in the matrix.
#'
#' @examples
#' get_between_trait_cov(.5, .5, .5)
#'
#' @export
#gets the correlation between two disorders (between disorder)
get_between_trait_cov <- function(h2_1, h2_2, rho, n_sib = 0) {
  
  warning("'get_between_trait_cov()' was deprecated in LTFHPlus v1.0.0. It is only kept for legacy reasons.\n
          Please consider using 'construct_covmat()' instead.\
          The signature and semantics have changed, see '?construct_covmat'.")
  
  corr = suppressWarnings(get_cov(1, n_sib = n_sib))
  corr * sqrt(h2_1 * h2_2) * rho
}

#'
#' Helps construct the correlation entries in the covariance matrix needed for LT-FH++ with multiple traits. 
#' 
#' This function has been retired in favor of \code{\link{construct_covmat}}. This function is kept for legacy reasons.
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
  
  warning("'get_full_cov()' was deprecated in LTFHPlus v1.0.0. It is only kept for legacy reasons.\n
          Please consider using 'construct_covmat()' instead.\
          The signature and semantics have changed, see '?construct_covmat'.")
  
  n_trait = nrow(corr_mat)
  block_size = 4 + n_sib
  cov_size = block_size * n_trait
  full_cov = matrix(NA, ncol = cov_size, nrow = cov_size) 
  
  for (i in 1:n_trait) {
    for (j in 1:n_trait) {
      if (i == j) {
        full_cov[(i - 1) * block_size + 1:block_size, (i - 1) * block_size + 1:block_size] = suppressWarnings(get_cov(corr_mat[i,i], n_sib = n_sib))  
      } else {
        full_cov[(i - 1) * block_size + 1:block_size, (j - 1) * block_size + 1:block_size] = suppressWarnings(get_between_trait_cov(h2_1 = corr_mat[i,i],
                                                                                                                                    h2_2 = corr_mat[j,j],
                                                                                                                                    rho = corr_mat[i,j],
                                                                                                                                    n_sib = n_sib))
      }
    }
  }
  return(full_cov)
}