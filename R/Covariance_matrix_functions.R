utils::globalVariables("sq.herit")
#' Relatedness between a pair of family members
#'
#' \code{get_relatedness} returns the relatedness times the
#' squared heritability for a pair of family members
#'
#' This function can be used to get the percentage of shared
#' DNA times the squared heritability \eqn{h^2} for two family members.
#'
#' @param s1,s2 Strings representing the two family members.
#' The strings must be chosen from the following list of strings:
#' - g (Genetic component of full liability)
#' - o (Full liability)
#' - m (Mother)
#' - f (Father)
#' - mgm (Maternal grandmother)
#' - mgf (Maternal grandfather)
#' - pgm (Paternal grandmother)
#' - pgf (Paternal grandfather)
#' - s\[0-9\]* (Full siblings)
#' - mhs\[0-9\]* (Half-siblings - maternal side)
#' - phs\[0-9\]* (Half-siblings - paternal side)
#' - mau\[0-9\]* (Aunts/Uncles - maternal side)
#' - pau\[0-9\]* (Aunts/Uncles - paternal side).
#' @param sq.herit A number representing the squared heritability on liability scale.
#' Must be non-negative and at most 1. Defaults to 0.5
#' 
#' @return If both \code{s1} and \code{s2} are strings chosen from the mentioned list of strings and \code{sq.herit} is a number
#' satisfying \eqn{0 \leq sq.herit \leq 1}, then the output will be a number that equals the percentage of shared 
#' DNA between \code{s1} and \code{s2} times the squared heritability \code{sq.herit}.
#' 
#' @note If you are only interested in the percentage of shared DNA, set \code{sq.herit = 1}.
#' 
#' @examples
#' get_relatedness("g","o")
#' get_relatedness("g","f",sq.herit = 1)
#' get_relatedness("o","s",sq.herit = 0.3)
#'
#'
#' \dontrun{
#' # This will result in errors:
#' get_relatedness("a","b")
#' get_relatedness(m,mhs)
#' }
#' 
#' @importFrom stringr str_detect
#' 
#' @export
get_relatedness <- function(s1,s2, sq.herit=0.5){
  
  # Checking that s1 and s2 are strings
  if(class(s1) != "character") stop("s1 must be a string!")
  if(class(s2) != "character") stop("s2 must be a string!")
  # Convert s1 and s2 to lowercase strings
  s1 <- tolower(s1)
  s2 <- tolower(s2)
  # Checking that s1 and s2 are valid strings
  if(check_valid_relatives(s1)){invisible()}
  if(check_valid_relatives(s2)){invisible()}

  # Checking that the heritability is valid
  if(check_proportion(sq.herit)){invisible()}
  
  # Getting the percentage of shared DNA
  if(str_detect(s1, "^o$")){
    
    if(s1==s2)return(1)
    if(str_detect(s2, "^g$")) return(1*sq.herit)
    if(str_detect(s2, "^m$") | str_detect(s2, "^f$") | str_detect(s2, "^s[0-9]*")) return(0.5*sq.herit)
    if(str_detect(s2, "^[mp]hs[0-9]*") | str_detect(s2, "^[mp]g[mf]") | str_detect(s2, "^[mp]au[0-9]*")) return(0.25*sq.herit)
    
  }else if(str_detect(s1, "^g$")){
    
    if(str_detect(s2, "^[go]$")) return(1*sq.herit)
    if(str_detect(s2, "^m$") | str_detect(s2, "^f$") | str_detect(s2, "^s[0-9]*")) return(0.5*sq.herit)
    if(str_detect(s2, "^[mp]hs[0-9]*") | str_detect(s2, "^[mp]g[mf]") | str_detect(s2, "^[mp]au[0-9]*")) return(0.25*sq.herit)
    
  }else if(str_detect(s1, "^m$")){
    
    if(str_detect(s2, "^m$")) return(1)
    if(str_detect(s2, "^[go]$")|str_detect(s2, "^s[0-9]*")|str_detect(s2, "^mhs[0-9]*")|str_detect(s2, "^mg[mf]$")|str_detect(s2, "^mau[0-9]") ) return(0.5*sq.herit)
    if(str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^phs[0-9]*") | str_detect(s2, "^pau[0-9]*")) return(0)
    
  }else if(str_detect(s1, "^f$")){
    
    if(str_detect(s2, "^f$")) return(1)
    if(str_detect(s2, "^[go]$")|str_detect(s2, "^s[0-9]*")|str_detect(s2, "^phs[0-9]*")|str_detect(s2, "^pg[mf]$")|str_detect(s2, "^pau[0-9]") ) return(0.5*sq.herit)
    if(str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mhs[0-9]*") | str_detect(s2, "^mau[0-9]*")) return(0)
    
  }else if(str_detect(s1, "^s[0-9]*")){
    
    if(s1==s2) return(1)
    if(str_detect(s2, "^[go]$") | str_detect(s2, "^m$") | str_detect(s2, "^f$") | str_detect(s2, "^s[0-9]*")) return(0.5*sq.herit)
    if(str_detect(s2, "^[mp]hs[0-9]*") | str_detect(s2, "^[mp]g[mf]") | str_detect(s2, "^[mp]au[0-9]*")) return(0.25*sq.herit)
    
  }else if(str_detect(s1, "^mg[mf]$")){
    
    if(s1==s2)return(1)
    if(str_detect(s2, "^mg[mf]$")) return(0)
    if(str_detect(s2, "^[go]$")|str_detect(s2, "^s[0-9]*")|str_detect(s2, "^mhs[0-9]*")) return(0.25*sq.herit)
    if(str_detect(s2, "^m$") | str_detect(s2, "^mau[0-9]*") ) return(0.5*sq.herit)
    if(str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$")| str_detect(s2, "^phs[0-9]*") | str_detect(s2, "^pau[0-9]*")) return(0)
    
  }else if(str_detect(s1, "^pg[mf]$")){
    
    if(s1==s2)return(1)
    if(str_detect(s2, "^pg[mf]$")) return(0)
    if(str_detect(s2, "^[go]$")|str_detect(s2, "^s[0-9]*")|str_detect(s2, "^phs[0-9]*")) return(0.25*sq.herit)
    if(str_detect(s2, "^f$") | str_detect(s2, "^pau[0-9]*")) return(0.5*sq.herit)
    if(str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$")| str_detect(s2, "^mhs[0-9]*") | str_detect(s2, "^mau[0-9]*")) return(0)
    
  }else if(str_detect(s1, "^mhs[0-9]*")){
    
    if(s1==s2)return(1)
    if(str_detect(s2, "^[go]$") | str_detect(s2, "^s[0-9]*")|str_detect(s2, "^mg[mf]$")|str_detect(s2, "^mau[0-9]")) return(0.25*sq.herit)
    if(str_detect(s2, "^m$") | str_detect(s2, "^mhs[0-9]*")) return(0.5*sq.herit)
    if(str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^phs[0-9]*") | str_detect(s2, "^pau[0-9]*")) return(0)
    
  }else if(str_detect(s1, "^phs[0-9]*")){
    
    if(s1==s2)return(1)
    if(str_detect(s2, "^g$")|str_detect(s2, "^o$")|str_detect(s2, "^s[0-9]*")|str_detect(s2, "^pg[mf]$")|str_detect(s2, "^pau[0-9]")) return(0.25*sq.herit)
    if(str_detect(s2, "^f$") | str_detect(s2, "^phs[0-9]*")) return(0.5*sq.herit)
    if(str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$")| str_detect(s2, "^mhs[0-9]*") | str_detect(s2, "^mau[0-9]*")) return(0)
    
  }else if(str_detect(s1, "^mau[0-9]*")){
    
    if(s1==s2)return(1)
    if(str_detect(s2, "^g$")|str_detect(s2, "^o$")|str_detect(s2, "^s[0-9]*")| str_detect(s2, "^mhs[0-9]")) return(0.25*sq.herit)
    if(str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mau[0-9]*")) return(0.5*sq.herit)
    if(str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$")| str_detect(s2, "^phs[0-9]*") | str_detect(s2, "^pau[0-9]*")) return(0)
    
  }else if(str_detect(s1, "^pau[0-9]*")){
    
    if(s1==s2)return(1)
    if(str_detect(s2, "^g$")|str_detect(s2, "^o$")|str_detect(s2, "^s[0-9]*")| str_detect(s2, "^phs[0-9]")) return(0.25*sq.herit)
    if(str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pau[0-9]*")) return(0.5*sq.herit)
    if(str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$")| str_detect(s2, "^mhs[0-9]*") | str_detect(s2, "^mau[0-9]*")) return(0)
    
  }else{
    return(NA)
  }
}

#' Constructing a covariance matrix for a single phenotype
#'
#' \code{construct_covmatc_single} returns the covariance matrix for an
#' underlying individual and a variable number of its family members
#'
#' This function can be used to construct a covariance matrix for
#' a given number of family members. Each entry in this covariance 
#' matrix equals the percentage of shared DNA between the corresponding 
#' individuals times the squared heritability \eqn{h^2}. The family members
#' can be specified using one of two possible formats.
#'
#' @param fam_vec A vector of strings holding the different 
#' family members. All family members must be represented by strings from the 
#' following list:
#' 
#' - m (Mother)
#' - f (Father)
#' - mgm (Maternal grandmother)
#' - mgf (Maternal grandfather)
#' - pgm (Paternal grandmother)
#' - pgf (Paternal grandfather)
#' - s\[0-9\]* (Full siblings)
#' - mhs\[0-9\]* (Half-siblings - maternal side)
#' - phs\[0-9\]* (Half-siblings - paternal side)
#' - mau\[0-9\]* (Aunts/Uncles - maternal side)
#' - pau\[0-9\]* (Aunts/Uncles - paternal side).
#'  Defaults to c("m","f","s1","mgm","mgf","pgm","pgf").
#' @param n_fam A named vector holding the desired number of family members.
#' All names must be picked from the list mentioned above. Defaults to NULL.
#' @param add_ind A logical scalar indicating whether the genetic 
#' component of the full liability as well as the full
#' liability for the underlying individual should be included in 
#' the covariance matrix. Defaults to TRUE.
#' @param sq.herit A number representing the squared heritability on liability scale
#' for a single phenotype. Must be non-negative and at most 1.
#' Defaults to 0.5.
#' 
#' @return If either \code{fam_vec} or \code{n_fam} is used as the argument, if it 
#' is of the required format and \code{sq.herit} is a number satisfying 
#' \eqn{0 \leq sq.herit \leq 1}, then the output will be a named covariance matrix. 
#' The number of rows and columns corresponds to the length of \code{fam_vec}
#' or \code{n_fam} (+ 2 if \code{add_ind=TRUE}). 
#' If both \code{fam_vec = c()/NULL} and \code{n_fam = c()/NULL}, the 
#' function returns a \eqn{2 \times 2} matrix holding only the correlation
#' between the genetic component of the full liability and 
#' the full liability for the individual. If both \code{fam_vec} and 
#' \code{n_fam} are given, the user is asked to decide on which 
#' of the two vectors to use.
#' Note that the returned object has different attributes, such as 
#' \code{fam_vec}, \code{n_fam}, \code{add_ind} and \code{sq.herit}.
#' 
#' @examples
#' construct_covmat_single()
#' construct_covmat_single(fam_vec = c("m","mgm","mgf","mhs1","mhs2","mau1"), 
#' n_fam = NULL, add_ind = TRUE, sq.herit = 0.5)
#' construct_covmat_single(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2,2), 
#' c("m","mgm","mgf","s","mhs")), add_ind = FALSE, sq.herit = 0.3)
#' 
#' @seealso \code{\link{get_relatedness}}, \code{\link{construct_covmat_multi}},
#' \code{\link{construct_covmat}}
#' @export
construct_covmat_single <- function(fam_vec = c("m","f","s1","mgm","mgf","pgm","pgf"), n_fam = NULL, 
                                    add_ind = TRUE, sq.herit = 0.5){
  
  # Turn add_ind into a logical 
  add_ind <- as.logical(add_ind)
  
  # Checking that the squared heritability is valid
  if(check_proportion(sq.herit)){invisible()}
  
  # If neither fam_vec nor n_fam are specified, the
  # covariance matrix will simply include the 
  # correlation between the genetic component
  # and the full liability for the individual.
  if(is.null(fam_vec) && is.null(n_fam)){
    
    warning("Neither fam_vec nor n_fam is specified...")
    # Constructing the simple 2x2 matrix
    covmat <- matrix(c(sq.herit,sq.herit,sq.herit,1), nrow = 2)
    # Naming the columns and rows
    colnames(covmat) <- rownames(covmat) <- c("g","o")
    
    # Adding attributes to covmat
    attributes(covmat)$fam_vec <- NULL
    attributes(covmat)$n_fam <- NULL
    attributes(covmat)$add_ind <- add_ind
    attributes(covmat)$sq.herit <- sq.herit
    
    return(covmat)
    
  }else if(!is.null(fam_vec) && !is.null(n_fam)){
    # If both fam_vec and n_fam are specified, the user 
    # needs to decide on which vector to use.
    # The user will be asked to enter y (to use fam_vec) 
    # or n (to use n_fam). If the entered string is not
    # y nor n three times, the function will be aborted.
    count <- 1
    ver <- "a"
    while(!is.null(fam_vec) && !is.null(n_fam)){
      
      if(ver == "y" & (count < 4)){
        n_fam <-  NULL
      }else if(ver == "n" & (count < 4)){
        fam_vec <- NULL
      }else if( count >=4){
        stop("Function aborted...")
      }else{
        count <- count + 1
        ver <- readline(prompt="Both fam_vec and n_fam are specified... \n Do you want to use fam_vec to create the covariance matrix [y/n]?: ")
      }
    }
  }
  
  if(is.null(fam_vec)){
    # If only n_fam is specified, the function uses this
    # vector to create the covariance matrix
    
    # Checking that n_fam is named
    if(is.null(names(n_fam))) stop("n_fam must be a named vector")
    # Checking that all family members are valid strings
    if(check_valid_relatives(names(n_fam))){invisible()}
    
    # Checking that all entries in n_fam are non-negative
    if(any(n_fam<0)) stop("All entries in n_fam must be non-negative!")
    
    # Removing all family members that occur zero times
    n_fam <- n_fam[n_fam > 0]
    
    # Constructing a vector holding all family members
    # (including g and o if add_ind = T).
    if(any(names(n_fam) %in% c("g","o"))){
      n_fam <- n_fam[-which(names(n_fam) %in% c("g","o"))]
    }
    if(add_ind){
      n_fam <- c(stats::setNames(c(1,1), c("g", "o")), n_fam)
    }
    
    # Extracting all family members that can only occur exactly once
    single_indx <- which(stringr::str_detect(names(n_fam), "^[gomf]$") | stringr::str_detect(names(n_fam), "^[mp]g[mf]$"))
    # Extracting all family members that can occur multiple times
    multiple_indx <- setdiff(1:length(n_fam), single_indx)
    
    # Getting the names for all family members that
    # occur only once
    if(length(single_indx)==0){
      fam_vec <- c()
    }else{
      fam_vec <- names(n_fam)[single_indx]
    }
    # And the names for those family members that occur
    # several times
    if(length(multiple_indx)>0){
      for(indx in multiple_indx){
        
        fam_vec <- c(fam_vec, paste0(names(n_fam)[indx],"",1:n_fam[indx]))
      }
    }
    
  }else if(is.null(n_fam)){
    # If only fam_vec is specified, the function uses this
    # vector to create the covariance matrix
    
    # Checking that all family members are represented by valid strings
    if(check_valid_relatives(fam_vec)){invisible()}
    
    # If add_ind = T, the genetic component and the full
    # liability are added to the family members
    if(any(fam_vec %in% c("g","o"))){
      fam_vec <- fam_vec[-which(fam_vec %in% c("g","o"))]
    }
    if(add_ind){
      fam_vec <- c(c("g", "o"), fam_vec)
    }
    
    n_fam <- table(sub("[0-9]*$", "", fam_vec))
  }
    
  # Now that we have a vector holding all desired family members,
  # we can create the covariance matrix
  covmat <- matrix(NA, nrow = length(fam_vec), ncol = length(fam_vec))
  # Changing the row and column names
  rownames(covmat) <- colnames(covmat) <- fam_vec
    
  # Filling in all entries 
  for(mem in fam_vec) {
    covmat[which(rownames(covmat) == mem),] <- sapply(fam_vec, get_relatedness, s1 = mem, sq.herit = sq.herit)
  }
  
  # Adding attributes to covmat
  attributes(covmat)$fam_vec <- fam_vec
  attributes(covmat)$n_fam <- n_fam
  attributes(covmat)$add_ind <- add_ind
  attributes(covmat)$sq.herit <- sq.herit
    
  return(covmat)
}

#' Constructing a covariance matrix for multiple phenotypes
#'
#' \code{construct_covmat_multi} returns the covariance matrix for an 
#' underlying individual and a variable number of its family members
#' for multiple phenotypes.
#'
#' This function can be used to construct a covariance matrix for
#' a given number of family members. Each entry in this covariance 
#' matrix equals either the percentage of shared DNA between the corresponding 
#' individuals times the squared heritability \eqn{h^2} or the 
#' percentage of shared DNA between the corresponding individuals times 
#' the correlation between the corresponding phenotypes. The family members
#' can be specified using one of two possible formats.
#' 
#' @param fam_vec A vector of strings holding the different 
#' family members. All family members must be represented by strings from the 
#' following list:
#' - m (Mother)
#' - f (Father)
#' - mgm (Maternal grandmother)
#' - mgf (Maternal grandfather)
#' - pgm (Paternal grandmother)
#' - pgf (Paternal grandfather)
#' - s\[0-9\]* (Full siblings)
#' - mhs\[0-9\]* (Half-siblings - maternal side)
#' - phs\[0-9\]* (Half-siblings - paternal side)
#' - mau\[0-9\]* (Aunts/Uncles - maternal side)
#' - pau\[0-9\]* (Aunts/Uncles - paternal side).
#'  Defaults to c("m","f","s1","mgm","mgf","pgm","pgf").
#' @param n_fam A named vector holding the desired number of family members.
#' All names must be picked from the list mentioned above. Defaults to NULL.
#' @param add_ind A logical scalar indicating whether the genetic 
#' component of the full liability as well as the full
#' liability for the underlying individual should be included in 
#' the covariance matrix. Defaults to TRUE.
#' @param genetic_corrmat A numeric matrix holding the genetic correlations between the desired 
#' phenotypes. All diagonal entries must be equal to one, while all off-diagonal entries 
#' must be between -1 and 1. In addition, the matrix must be symmetric.
#' @param full_corrmat A numeric matrix holding the full correlations between the desired 
#' phenotypes. All diagonal entries must be equal to one, while all off-diagonal entries 
#' must be between -1 and 1. In addition, the matrix must be symmetric.
#' @param sq.herits A numeric vector representing the squared heritabilities on liability scale
#' for all phenotypes. All entries in sq.herits must be non-negative and at most 1.
#' @param phen_names A character vector holding the phenotype names. These names
#' will be used to create the row and column names for the covariance matrix.
#' If it is not specified, the names will default to phenotype1, phenotype2, etc.
#' Defaults to NULL.
#' 
#' @return If either \code{fam_vec} or \code{n_fam} is used as the argument and if it is of the 
#' required format, if \code{genetic_corrmat} and \code{full_corrmat} are two numeric and symmetric matrices 
#' satisfying that all diagonal entries are one and that all off-diagonal
#' entries are between -1 and 1, and if \code{sq.herits} is a numeric vector satisfying
#' \eqn{0 \leq sq.herits_i \leq 1} for all \eqn{i \in \{1,...,n_pheno\}},
#' then the output will be a named covariance matrix. 
#' The number of rows and columns corresponds to the number of phenotypes times 
#' the length of \code{fam_vec} or \code{n_fam} (+ 2 if \code{add_ind=TRUE}). 
#' If both \code{fam_vec} and \code{n_fam} are equal to \code{c()} or \code{NULL},
#' the function returns a \eqn{(2 \times n_pheno) \times (2\times n_pheno)} 
#' matrix holding only the correlation between the genetic component of the full
#' liability and the full liability for the underlying individual for all
#' phenotypes. If both \code{fam_vec} and \code{n_fam} are specified, the user is asked to 
#' decide on which of the two vectors to use.
#' Note that the returned object has a number different attributes,namely
#' \code{fam_vec}, \code{n_fam}, \code{add_ind}, \code{genetic_corrmat}, \code{full_corrmat},
#' \code{sq.herits} and \code{phenotype_names}.
#' 
#' @examples
#' construct_covmat_multi(fam_vec = NULL, 
#'                        genetic_corrmat = matrix(c(0.9, 0.5, 0.5, 0.8), nrow = 2),
#'                        full_corrmat = matrix(c(0.9, 0.5, 0.5, 0.8), nrow = 2),
#'                        sq.herits = c(0.37,0.44),
#'                        phen_names = c("p1","p2"))
#' construct_covmat_multi(fam_vec = c("m","mgm","mgf","mhs1","mhs2","mau1"), 
#'                        n_fam = NULL, 
#'                        add_ind = TRUE,
#'                        genetic_corrmat = diag(3),
#'                        full_corrmat = diag(3),
#'                        sq.herits = c(0.8, 0.65))
#' construct_covmat_multi(fam_vec = NULL, 
#'                        n_fam = stats::setNames(c(1,1,1,2,2), c("m","mgm","mgf","s","mhs")), 
#'                        add_ind = FALSE,
#'                        genetic_corrmat = diag(2),
#'                        full_corrmat = diag(2),
#'                        sq.herits = c(0.75,0.85))
#' 
#' @seealso \code{\link{get_relatedness}}, \code{\link{construct_covmat_single}} and
#' \code{\link{construct_covmat}}.
#' 
#' @export
construct_covmat_multi <- function(fam_vec = c("m","f","s1","mgm","mgf","pgm","pgf"), n_fam = NULL, add_ind = TRUE, 
                                   genetic_corrmat, full_corrmat, sq.herits, phen_names = NULL){
  
  # Turn add_ind into a logical 
  add_ind <- as.logical(add_ind)
  
  # Checking that the heritabilities are valid
  if(check_proportion(sq.herits)){invisible()}
  
  # Checking that all correlations are valid
  if(check_correlation_matrix(genetic_corrmat)){invisible()}
  if(check_correlation_matrix(full_corrmat)){invisible()}
  # Computing the number of phenotypes
  num_phen <- length(sq.herits)
  
  # Checking that phen_names is either NULL or a valid
  # vector of strings
  if(is.null(phen_names)){
    phen_names <- paste0("phenotype", 1:num_phen)
  }else{
    if(class(phen_names) != "character") phen_names <- as.character(phen_names)
    if(length(phen_names) != num_phen) stop("The number of names in phen_num and the number of phenotypes differ...")
  }

  # If neither fam_vec nor n_fam are specified, the
  # covariance matrix will simply include the 
  # correlation between the genetic component
  # and the full liability for the individual
  # for all phenotypes
  if(is.null(fam_vec) && is.null(n_fam)){
    
    cat("Warning message: \n Neither fam_vec nor n_fam is specified...")
    # Constructing a simple covariance matrix
    covmat <- matrix(NA, nrow = 2*num_phen, ncol = 2*num_phen)
    # Filling in all entries
    for(p1 in 1:num_phen){
      for(p2 in 1:num_phen){
        
        if(p1==p2){
          
          covmat[2*(p1-1) + 1:2, 2*(p2-1) + 1:2] <- matrix(c(sq.herits[p1], sq.herits[p1], sq.herits[p1], 1), nrow = 2)
        }else{
          
          covmat[2*(p1-1) + 1:2, 2*(p2-1) + 1:2] <- matrix(sqrt(sq.herits[p1]*sq.herits[p2])*full_corrmat[p1,p2], nrow = 2, ncol = 2)
        }
      }
    }
    # Changing the row and column names
    colnames(covmat) <- rownames(covmat) <- paste0(c("g_", "o_"), rep(phen_names, each = 2))
    
    # Adding attributes to covmat
    attributes(covmat)$fam_vec <- NULL
    attributes(covmat)$n_fam <- NULL
    attributes(covmat)$add_ind <- add_ind
    attributes(covmat)$sq.herit <- sq.herits
    attributes(covmat)$genetic_corrmat <- genetic_corrmat
    attributes(covmat)$full_corrmat <- full_corrmat
    attributes(covmat)$phenotype_names <- phen_names
    
    return(covmat)
    
  }else if(!is.null(fam_vec) && !is.null(n_fam)){
    # If both fam_vec and n_fam are specified, the user 
    # needs to decide on which vector to use.
    # The user will be asked to enter y (to use fam_vec) 
    # or n (to use n_fam). If the entered string is not
    # y nor n three times, the function will be aborted.
    count <- 1
    ver <- "a"
    while(!is.null(fam_vec) && !is.null(n_fam)){
      
      if(ver == "y" & (count < 4)){
        n_fam <-  NULL
      }else if(ver == "n" & (count < 4)){
        fam_vec <- NULL
      }else if( count >=4){
        stop("Function aborted...")
      }else{
        count <- count + 1
        ver <- readline(prompt="Both fam_vec and n_fam are specified... \n Do you want to use fam_vec to create the covariance matrix [y/n]?: ")
      }
    }
  }
  if(is.null(fam_vec)){
    # If only n_fam is specified, the function uses this
    # vector to create the covariance matrix
    
    # Checking that n_fam is named
    if(is.null(names(n_fam))) stop("n_fam must be a named vector")
    # Checking that all family members are valid strings
    if(check_valid_relatives(names(n_fam))){invisible()}
    
    # Checking that all entries in n_fam are non-negative
    if(any(n_fam<0)) stop("All entries in n_fam must be non-negative!")
    
    # Removing all family members that occur zero times
    n_fam <- n_fam[n_fam > 0]
    
    # Constructing a vector holding all family members
    # (including g and o if add_ind = T).
    if(any(names(n_fam) %in% c("g","o"))){
      n_fam <- n_fam[-which(names(n_fam) %in% c("g","o"))]
    }
    if(add_ind){
      n_fam <- c(stats::setNames(c(1,1), c("g", "o")), n_fam)
    }
    
    # Extracting all family members that can only occur exactly once
    single_indx <- which(stringr::str_detect(names(n_fam), "^[gomf]$") | stringr::str_detect(names(n_fam), "^[mp]g[mf]$"))
    # Extracting all family members that can occur multiple times
    multiple_indx <- setdiff(1:length(n_fam), single_indx)
    
    # Getting the names for all family members that
    # occur only once
    if(length(single_indx)==0){
      fam_vec <- c()
    }else{
      fam_vec <- names(n_fam)[single_indx]
    }
    # And the names for those family members that occur
    # several times
    if(length(multiple_indx)>0){
      for(indx in multiple_indx){
        
        fam_vec <- c(fam_vec, paste0(names(n_fam)[indx],"",1:n_fam[indx]))
      }
    }
    
  }else if(is.null(n_fam)){
    # If only fam_vec is specified, the function uses this
    # vector to create the covariance matrix
    
    # Checking that all family members are represented by valid strings
    if(check_valid_relatives(fam_vec)){invisible()}
    
    # If add_ind = T, the genetic component and the full
    # liability are added to the family members
    if(any(fam_vec %in% c("g","o"))){
      fam_vec <- fam_vec[-which(fam_vec %in% c("g","o"))]
    }
    if(add_ind){
      fam_vec <- c(c("g", "o"), fam_vec)
    }
    
    n_fam <- table(sub("[0-9]*$", "", fam_vec))
  }
  
  # Now that we have a vector holding the desired family
  # members, we can create the covariance matrix.
  covmat <- matrix(NA, nrow = length(fam_vec)*num_phen, ncol = length(fam_vec)*num_phen)
  # Changing the row and column names
  rownames(covmat) <- colnames(covmat) <- paste0(fam_vec,"_", rep(phen_names, each = length(fam_vec)))
  
  # Filling in all entries
  for(p1 in 1:num_phen){
    for(p2 in 1:num_phen){
      for(mem in fam_vec){
        
        if(p1==p2){
          
          covmat[which(rownames(covmat) == paste0(mem,"_", phen_names[p1])), (p2-1)*length(fam_vec) + 1:length(fam_vec)] <- sapply(fam_vec, get_relatedness, s1 = mem, sq.herit = sq.herits[p1])
        }else{
         
          covmat[which(rownames(covmat) == paste0(mem,"_", phen_names[p1])), (p2-1)*length(fam_vec) + 1:length(fam_vec)] <- sapply(fam_vec, get_relatedness, s1 = mem, sq.herit = sqrt(sq.herits[p1]*sq.herits[p2])*genetic_corrmat[p1,p2]) 
        }
       
      }
    }
  }
  
  # Adding attributes to covmat
  attributes(covmat)$fam_vec <- fam_vec
  attributes(covmat)$n_fam <- n_fam
  attributes(covmat)$add_ind <- add_ind
  attributes(covmat)$sq.herit <- sq.herits
  attributes(covmat)$genetic_corrmat <- genetic_corrmat
  attributes(covmat)$full_corrmat <- full_corrmat
  attributes(covmat)$phenotype_names <- phen_names
  
  return(covmat)
}

#' Constructing a covariance matrix for a variable number of
#' phenotypes
#'
#' \code{construct_covmat} returns the covariance matrix for an
#' underlying individual and a variable number of its family members
#' for a variable number of phenotypes. It is a wrapper around 
#' \code{\link{construct_covmat_single}} and \code{\link{construct_covmat_multi}}. 
#'
#' This function can be used to construct a covariance matrix for
#' a given number of family members. If sq.herit is a number,
#' each entry in this covariance matrix equals the percentage 
#' of shared DNA between the corresponding individuals times 
#' the squared heritability \deqn{h^2}. However, if sq.herit is a numeric vector,
#' and genetic_corrmat and full_corrmat are two symmetric correlation matrices,
#' each entry equals either the percentage of shared DNA between the corresponding 
#' individuals times the squared heritability \deqn{h^2} or the 
#' percentage of shared DNA between the corresponding individuals times 
#' the correlation between the corresponding phenotypes. The family members
#' can be specified using one of two possible formats.
#'
#' @param fam_vec A vector of strings holding the different 
#' family members. All family members must be represented by strings from the 
#' following list:
#' - m (Mother)
#' - f (Father)
#' - mgm (Maternal grandmother)
#' - mgf (Maternal grandfather)
#' - pgm (Paternal grandmother)
#' - pgf (Paternal grandfather)
#' - s\[0-9\]* (Full siblings)
#' - mhs\[0-9\]* (Half-siblings - maternal side)
#' - phs\[0-9\]* (Half-siblings - paternal side)
#' - mau\[0-9\]* (Aunts/Uncles - maternal side)
#' - pau\[0-9\]* (Aunts/Uncles - paternal side).
#'  Defaults to c("m","f","s1","mgm","mgf","pgm","pgf").
#' @param n_fam A named vector holding the desired number of family members.
#' All names must be picked from the list mentioned above. Defaults to NULL.
#' @param add_ind A logical scalar indicating whether the genetic 
#' component of the full liability as well as the full
#' liability for the underlying individual should be included in 
#' the covariance matrix. Defaults to TRUE.
#' @param sq.herit Either a number representing the squared heritability 
#' on liability scale for one single phenotype or a numeric vector representing
#' the squared heritabilities on liability scale for a positive number of phenotypes.
#' All entries in sq.herit must be non-negative and at most 1.
#' @param genetic_corrmat Either \code{NULL} or a numeric matrix holding the genetic correlations between the desired 
#' phenotypes. All diagonal entries must be equal to one, while all off-diagonal entries 
#' must be between -1 and 1. In addition, the matrix must be symmetric.
#' Defaults to NULL.
#' @param full_corrmat Either \code{NULL} or a  numeric matrix holding the full correlations between the desired 
#' phenotypes. All diagonal entries must be equal to one, while all off-diagonal entries 
#' must be between -1 and 1. In addition, the matrix must be symmetric.
#' Defaults to NULL.
#' @param phen_names Either \code{NULL} or a character vector holding the phenotype names. These names
#' will be used to create the row and column names for the covariance matrix.
#' If it is not specified, the names will default to phenotype1, phenotype2, etc.
#' Defaults to NULL.
#' 
#' @return If either \code{fam_vec} or \code{n_fam} is used as the argument, if it is of 
#' the required format, if \code{add_ind} is a logical scalar and \code{sq.herit} is a 
#' number satisfying \deqn{0 \leq sq.herit \leq 1}, then the function \code{construct_covmat}
#' will return a named covariance matrix, which row- and column-number 
#' corresponds to the length of \code{fam_vec} or \code{n_fam} (+ 2 if \code{add_ind=TRUE}).
#' However, if \code{sq.herits} is a numeric vector satisfying
#' \deqn{0 \leq sq.herits_i \leq 1} for all \deqn{i \in \{1,...,n_pheno\}} and if 
#' \code{genetic_corrmat} and \code{full_corrmat} are two numeric and symmetric matrices 
#' satisfying that all diagonal entries are one and that all off-diagonal
#' entries are between -1 and 1, then \code{construct_covmat} will return 
#' a named covariance matrix, which number of rows and columns corresponds to the number 
#' of phenotypes times the length of \code{fam_vec} or \code{n_fam} (+ 2 if \code{add_ind=TRUE}). 
#' If both \code{fam_vec} and \code{n_fam} are equal to \code{c()} or \code{NULL},
#' the function returns either a \eqn{2 \times 2} matrix holding only the correlation
#' between the genetic component of the full liability and the full liability for the 
#' individual under consideration, or a \deqn{(2 \times n_pheno) \times (2\times n_pheno)} 
#' matrix holding the correlation between the genetic component of the full
#' liability and the full liability for the underlying individual for all
#' phenotypes. 
#' If both \code{fam_vec} and \code{n_fam} are specified, the user is asked to 
#' decide on which of the two vectors to use.
#' Note that the returned object has different attributes, such as 
#' \code{fam_vec}, \code{n_fam}, \code{add_ind} and \code{sq.herit}.
#' 
#' @examples
#' construct_covmat()
#' construct_covmat(fam_vec = c("m","mgm","mgf","mhs1","mhs2","mau1"), 
#'                  n_fam = NULL, 
#'                  add_ind = TRUE, 
#'                  sq.herit = 0.5)
#' construct_covmat(fam_vec = NULL, 
#'                  n_fam = stats::setNames(c(1,1,1,2,2), c("m","mgm","mgf","s","mhs")), 
#'                  add_ind = FALSE, sq.herit = 0.3)
#' construct_covmat(sq.herit = c(0.5,0.5), genetic_corrmat = matrix(c(1,0.6,0.4,1), nrow = 2),
#' full_corrmat = matrix(c(1,0.7,0.6,1), nrow = 2))
#' 
#' @seealso \code{\link{get_relatedness}}, \code{\link{construct_covmat_single}},
#' \code{\link{construct_covmat_multi}}
#' @export
construct_covmat <- function(fam_vec = c("m","f","s1","mgm","mgf","pgm","pgf"), n_fam = NULL, 
                             add_ind = TRUE, sq.herit = 0.5, genetic_corrmat = NULL, 
                             full_corrmat = NULL, phen_names = NULL){
  
  if(length(sq.herit) == 1){
    
    return(construct_covmat_single(fam_vec = fam_vec, n_fam = n_fam, add_ind = add_ind, sq.herit = sq.herit))
    
  }else{
    
    return(construct_covmat_multi(fam_vec = fam_vec, n_fam = n_fam, add_ind = add_ind, 
                                  genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat,
                                  sq.herits = sq.herit, phen_names = phen_names))
  } 
}
