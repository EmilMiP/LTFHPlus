#' Relatedness between a pair of family members
#'
#' \code{get_relatedness} returns the relatedness times the
#' heritability for a pair of family members
#'
#' This function can be used to get the percentage of shared
#' DNA times the heritability \eqn{h^2} for two family members.
#'
#' @param s1,s2 Strings representing the two family members.
#' The string must be one of the following:
#' - g (Genetic component of full liability)
#' - o (Full liability)
#' - m (Mother)
#' - f (Father)
#' - mgm (Maternal grandmother)
#' - mgf (Maternal grandfather)
#' - mgp\[1-2\]* (Maternal grandparent)
#' - pgm (Paternal grandmother)
#' - pgf (Paternal grandfather)
#' - pgp\[1-2\]* (Paternal grandparent)
#' - s\[0-9\]* (Full siblings)
#' - c\[0-9\]* (children)
#' - mhs\[0-9\]* (Half-siblings - maternal side)
#' - phs\[0-9\]* (Half-siblings - paternal side)
#' - mau\[0-9\]* (Aunts/Uncles - maternal side)
#' - pau\[0-9\]* (Aunts/Uncles - paternal side).
#' @param h2 A number representing the heritability on liability scale.
#' Must be non-negative. Note that under the liability threshold model,
#' the heritability must also be at most 1.
#' 
#' @return If both s1 and s2 are strings chosen from the mentioned list of strings and h2 is a number
#' satisfying 0 <= h2, then the output will be a number that equals the percentage of shared 
#' DNA between \code{s1} and \code{s2} times the heritability \code{h2}.
#' 
#' @note If you are only interested in the percentage of shared DNA, set \code{h2 = 1}.
#' 
#' @examples
#' get_relatedness("g","o")
#' get_relatedness("g","f",h2 = 1)
#' get_relatedness("o","s",h2 = 0.3)
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
get_relatedness <- function(s1,s2, h2=0.5){

  # Checking that s1 and s2 are strings
  if (!is.character(s1)) stop("s1 must be a string!")
  if (!is.character(s2)) stop("s2 must be a string!")
  # Checking that s1 and s2 are valid strings
  check_family_role(s1)
  check_family_role(s2)
  # Checking that the heritability is valid
  if (h2 < 0) stop("The heritability must be non-negative!")
  if (h2 > 1) cat("Warning message: \n Under the liability threshold model, the heritability must be smaller than or equal to 1.")
  
  # Getting the percentage of shared DNA
  if (str_detect(s1, "^o$")) { # offspring
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^g$")) return(1*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^f$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^c[0-9]*")) return(0.5*h2)
    if (str_detect(s2, "^[mp]hs[0-9]*") | str_detect(s2, "^[mp]g[mf]") | str_detect(s2, "^[mp]gp[0-9]*") | str_detect(s2, "^[mp]au[0-9]*")) return(0.25*h2)
    
  } else if (str_detect(s1, "^g$")) { # offspring's genetic liability
    
    if (str_detect(s2, "^[go]$")) return(1*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^f$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^c[0-9]*")) return(0.5*h2)
    if (str_detect(s2, "^[mp]hs[0-9]*") | str_detect(s2, "^[mp]g[mf]") | str_detect(s2, "^[mp]gp[0-9]*") | str_detect(s2, "^[mp]au[0-9]*")) return(0.25*h2)
    
  } else if (str_detect(s1, "^m$")) { # mother
    
    if (str_detect(s2, "^m$")) return(1)
    if (str_detect(s2, "^[go]$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^mhs[0-9]*") | str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mgp[0-9]*") | str_detect(s2, "^mau[0-9]") ) return(0.5*h2)
    if (str_detect(s2, "^c[0-9]*")) return(0.25*h2)
    if (str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pgp[0-9]*") | str_detect(s2, "^phs[0-9]*") | str_detect(s2, "^pau[0-9]*")) return(0)
    
  } else if (str_detect(s1, "^f$")) { # father
    
    if (str_detect(s2, "^f$")) return(1)
    if (str_detect(s2, "^[go]$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^phs[0-9]*") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pgp[0-9]*") | str_detect(s2, "^pau[0-9]") ) return(0.5*h2)
    if (str_detect(s2, "^c[0-9]*")) return(0.25*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mgp[0-9]*") | str_detect(s2, "^mhs[0-9]*") | str_detect(s2, "^mau[0-9]*")) return(0)
    
  } else if (str_detect(s1, "^s[0-9]*")) { # sibling(s)
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^[go]$") | str_detect(s2, "^m$") | str_detect(s2, "^f$") | str_detect(s2, "^s[0-9]*")) return(0.5*h2)
    if (str_detect(s2, "^[mp]hs[0-9]*") | str_detect(s2, "^[mp]g[mf]") | str_detect(s2, "^[mp]gp[0-9]*") | str_detect(s2, "^[mp]au[0-9]*") | str_detect(s2, "^c[0-9]*")) return(0.25*h2)
    
  } else if (str_detect(s1, "^mg[mf]$")) { #maternal grand f/m
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mgp[0-9]*")) return(0)
    if (str_detect(s2, "^[go]$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^mhs[0-9]*")) return(0.25*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^mau[0-9]*") ) return(0.5*h2)
    if (str_detect(s2, "^c[0-9]*")) return(0.125*h2)
    if (str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pgp[0-9]*")| str_detect(s2, "^phs[0-9]*") | str_detect(s2, "^pau[0-9]*")) return(0)
    
  } else if (str_detect(s1, "^pg[mf]$")) { # paternal grand f/m
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pgp[0-9]*")) return(0)
    if (str_detect(s2, "^[go]$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^phs[0-9]*")) return(0.25*h2)
    if (str_detect(s2, "^f$") | str_detect(s2, "^pau[0-9]*")) return(0.5*h2)
    if (str_detect(s2, "^c[0-9]*")) return(0.125*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mgp[0-9]*") | str_detect(s2, "^mhs[0-9]*") | str_detect(s2, "^mau[0-9]*")) return(0)
    
  } else if (str_detect(s1, "^mgp[0-9]*")) { # maternal grandparent
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mgp[0-9]*")) return(0)
    if (str_detect(s2, "^[go]$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^mhs[0-9]*")) return(0.25*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^mau[0-9]*") ) return(0.5*h2)
    if (str_detect(s2, "^c[0-9]*")) return(0.125*h2)
    if (str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pgp[0-9]*") | str_detect(s2, "^phs[0-9]*") | str_detect(s2, "^pau[0-9]*")) return(0)
    
  } else if (str_detect(s1, "^pgp[0-9]*$")) { # paternal grandparent
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pgp[0-9]*")) return(0)
    if (str_detect(s2, "^[go]$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^phs[0-9]*")) return(0.25*h2)
    if (str_detect(s2, "^f$") | str_detect(s2, "^pau[0-9]*")) return(0.5*h2)
    if (str_detect(s2, "^c[0-9]*")) return(0.125*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mgp[0-9]*") | str_detect(s2, "^mhs[0-9]*") | str_detect(s2, "^mau[0-9]*")) return(0)
    
  }else if (str_detect(s1, "^mhs[0-9]*")) { # maternal halfsibling
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^[go]$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mgp[0-9]*") | str_detect(s2, "^mau[0-9]")) return(0.25*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^mhs[0-9]*")) return(0.5*h2)
    if (str_detect(s2, "^c[0-9]*")) return(0.125*h2)
    if (str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pgp[0-9]*") | str_detect(s2, "^phs[0-9]*") | str_detect(s2, "^pau[0-9]*")) return(0)
    
  } else if (str_detect(s1, "^phs[0-9]*")) { # paternal halfsibling
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^g$") | str_detect(s2, "^o$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pgp[0-9]*") | str_detect(s2, "^pau[0-9]")) return(0.25*h2)
    if (str_detect(s2, "^c[0-9]*")) return(0.125*h2)
    if (str_detect(s2, "^f$") | str_detect(s2, "^phs[0-9]*")) return(0.5*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mgp[0-9]*") | str_detect(s2, "^mhs[0-9]*") | str_detect(s2, "^mau[0-9]*")) return(0)
    
  } else if (str_detect(s1, "^mau[0-9]*")) { # maternal aunts and uncles
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^g$") | str_detect(s2, "^o$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^mhs[0-9]")) return(0.25*h2)
    if (str_detect(s2, "^c[0-9]*")) return(0.125*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mgp[0-9]*") | str_detect(s2, "^mau[0-9]*")) return(0.5*h2)
    if (str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pgp[0-9]*") | str_detect(s2, "^phs[0-9]*") | str_detect(s2, "^pau[0-9]*")) return(0)

  } else if (str_detect(s1, "^c[0-9]*")) { # index person's child
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^[go]$") | str_detect(s2, "^c[0-9]*")) return(0.5*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^f$") | str_detect(s2, "^s[0-9]*")) return(0.25*h2)
    if (str_detect(s2, "^[mp]hs[0-9]*") | str_detect(s2, "^[mp]g[mf]") | str_detect(s2, "^[mp]gp[0-9]*") | str_detect(s2, "^[mp]au[0-9]*")) return(0.125*h2) 

  } else if (str_detect(s1, "^pau[0-9]*")) { # paternal aunts and uncles
    
    if (s1 == s2) return(1)
    if (str_detect(s2, "^g$") | str_detect(s2, "^o$") | str_detect(s2, "^s[0-9]*") | str_detect(s2, "^phs[0-9]")) return(0.25*h2)
    if (str_detect(s2, "^c[0-9]*")) return(0.125*h2)
    if (str_detect(s2, "^f$") | str_detect(s2, "^pg[mf]$") | str_detect(s2, "^pgp[0-9]*") | str_detect(s2, "^pau[0-9]*")) return(0.5*h2)
    if (str_detect(s2, "^m$") | str_detect(s2, "^mg[mf]$") | str_detect(s2, "^mgp[0-9]*") |str_detect(s2, "^mhs[0-9]*") | str_detect(s2, "^mau[0-9]*")) return(0)
    
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
#' individuals times the heritability \eqn{h^2}. The family members
#' can be specified using one of two possible formats.
#'
#' @param fam_vec A vector of strings holding the different 
#' family members. All family members must be represented by strings from the 
#' following list:
#' - m (Mother)
#' - f (Father)
#' - mgm (Maternal grandmother)
#' - mgf (Maternal grandfather)
#' - mgp\[1-2\]* (Maternal grandparent)
#' - pgm (Paternal grandmother)
#' - pgf (Paternal grandfather)
#' - pgp\[1-2\]* (Paternal grandparent)
#' - s\[0-9\]* (Full siblings)
#' - c\[0-9\]* (children) 
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
#' @param h2 A number representing the heritability on liability scale
#' for a single phenotype. Must be non-negative. Note that under the liability threshold model,
#' the heritability must also be at most 1.
#' Defaults to 0.5.
#' 
#' @return If either fam_vec or n_fam is used as the argument, if it is of the required format and h2 is a number
#' satisfying 0 <= h2, then the output will be a named covariance matrix. The number of rows and columns
#' correspond to the length of fam_vec or n_fam (+ 2 if add_ind=T). 
#' If both fam_vec = c()/NULL and n_fam = c()/NULL, the 
#' function returns a \eqn{2 \times 2} matrix holding only the correlation
#' between the genetic component of the full liability and 
#' the full liability for the individual. If both fam_vec and 
#' n_fam are given, the user is asked to decide on which 
#' of the two vectors to use.
#' Note that the returned object has different attributes, such as 
#' fam_vec, n_fam, add_ind, h2 and phenotype_names.
#' 
#' @examples
#' construct_covmat_single()
#' construct_covmat_single(fam_vec = c("m","mgm","mgf","mhs1","mhs2","mau1"), 
#' n_fam = NULL, add_ind = TRUE, h2 = 0.5)
#' construct_covmat_single(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2,2), 
#' c("m","mgm","mgf","s","mhs")), add_ind = FALSE, h2 = 0.3)
#' 
#' @seealso \code{\link{get_relatedness}}, \code{\link{construct_covmat_multi}},
#' \code{\link{construct_covmat}}
#' @export
construct_covmat_single <- function(fam_vec = c("m","f","s1","mgm","mgf","pgm","pgf"), n_fam = NULL, add_ind = TRUE, h2 = 0.5){
  # forcing fam_vec to be null if a vector of length 0 is provided. 
  if (length(fam_vec) == 0) fam_vec = NULL
  
  # Turn add_ind into a logical 
  add_ind <- as.logical(add_ind)
  
  # Checking that the heritability is valid
  if(is.null(h2)) stop("The heritability must be specified!")
  if(h2<0) stop("The heritability must be non-negative")
  if(h2>1)cat("Warning message: \n Under the liability threshold model, the heritability must be smaller than or equal to 1.")
  
  # If neither fam_vec nor n_fam are specified, the
  # covariance matrix will simply include the 
  # correlation between the genetic component
  # and the full liability for the individual.
  if(is.null(fam_vec) && is.null(n_fam)){
    
    # Constructing the simple 2x2 matrix
    covmat <- matrix(c(h2,h2,h2,1), nrow = 2)
    # Naming the columns and rows
    colnames(covmat) <- rownames(covmat) <- c("g","o")
    
    # Adding attributes to covmat
    attributes(covmat)$fam_vec <- NULL
    attributes(covmat)$n_fam <- NULL
    attributes(covmat)$add_ind <- add_ind
    attributes(covmat)$h2 <- h2
    attributes(covmat)$phenotype_names <- NULL
    
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
    # Extracting the names of n_fam, which are equal to the family 
    # members
    fam <- names(n_fam)
    # Checking that all family members are valid strings
    check_family_role(fam)

    # Checking that all entries in n_fam are non-negative
    if(any(n_fam<0)) stop("All entries in n_fam must be non-negative!")
    
    # Removing all family members that occur zero times
    n_fam <- n_fam[n_fam > 0]
    
    # Constructing a vector holding all family members
    # (including g and o if add_ind = T).
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
      fam <- c()
    }else{
      fam <- names(n_fam)[single_indx]
    }
    # And the names for those family members that occur
    # several times
    if(length(multiple_indx)>0){
      for(indx in multiple_indx){
        
        fam <- c(fam, paste0(names(n_fam)[indx],"",1:n_fam[indx]))
      }
    }
    
    # Now that we have a vector holding all desired family members,
    # we can create the covariance matrix
    covmat <- matrix(NA, nrow = length(fam), ncol = length(fam))
    # Changing the row and column names
    rownames(covmat) <- colnames(covmat) <- fam
      
    # Filling in all entries 
    for(mem in fam) {
      covmat[which(rownames(covmat) == mem),] <- sapply(fam, get_relatedness, s1 = mem, h2 = h2)
    }
    
    # Adding attributes to covmat
    attributes(covmat)$fam_vec <- NULL
    attributes(covmat)$n_fam <- n_fam
    attributes(covmat)$add_ind <- add_ind
    attributes(covmat)$h2 <- h2
    attributes(covmat)$phenotype_names <- NULL
      
    return(covmat)
    
  }else if(is.null(n_fam)){
    # If only fam_vec is specified, the function uses this
    # vector to create the covariance matrix
    
    # Checking that all family members are represented by valid strings
    check_family_role(fam_vec)

    # If add_ind = T, the genetic component and the full
    # liability are added to the family members
    if(add_ind){
      fam_vec <- c(c("g", "o"), fam_vec)
    }
      
    # Now that we have a vector holding the desired family members,
    # we can build the covariance matrix.
    covmat <- matrix(NA, nrow = length(fam_vec), ncol = length(fam_vec))
    # Changing the row and column names
    rownames(covmat) <- colnames(covmat) <- fam_vec
    
    # Filling in all entries 
    for(mem in fam_vec) {
      covmat[which(rownames(covmat) == mem),] <- sapply(fam_vec, get_relatedness, s1 = mem, h2 = h2)
    }
    
    # Adding attributes to covmat
    attributes(covmat)$fam_vec <- fam_vec
    attributes(covmat)$n_fam <- NULL
    attributes(covmat)$add_ind <- add_ind
    attributes(covmat)$h2 <- h2
    attributes(covmat)$phenotype_names <- NULL
    
    return(covmat)
  }
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
#' individuals times the heritability \eqn{h^2} or the 
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
#' - mgp\[1-2\]* (Maternal grandparent)
#' - pgm (Paternal grandmother)
#' - pgf (Paternal grandfather)
#' - pgp\[1-2\]* (Paternal grandparent)
#' - s\[0-9\]* (Full siblings)
#' - c\[0-9\]* (children)
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
#' @param corrmat A numeric matrix holding the genetic correlation between the desired 
#' number of phenotypes as well as the full correlation. 
#' The full correlations must be given on the diagonal,
#' while the off-diagonal entries must hold the genetic correlations between phenotypes.
#' All correlations must be between -1 and 1.
#' @param h2s A numeric vector representing the heritability on liability scale
#' for all phenotypes. All entries in h2s must be non-negative and at most 1.
#' @param phen_names A character vector holding the phenotype names. These names
#' will be used to create the row and column names for the covariance matrix.
#' If it is not specified, the names will default to phenotype1, phenotype2, etc.
#' Defaults to NULL.
#' 
#' @return If either fam_vec or n_fam is used as the argument, if it is of the 
#' required format and h2 is a numeric matrix satisfying that all 
#' diagonal entries are non-negative and that all off-diagonal
#' entries are between -1 and 1, then the output will be a named covariance matrix. 
#' The number of rows and columns corresponds to the number of phenotypes times 
#' the length of fam_vec or n_fam (+ 2 if add_ind=T). 
#' If both fam_vec and n_fam are equal to c() or NULL, the function returns
#' a \eqn{(2 \times number of phenotypes) \times (2\times number of phenotypes)}{ (2 X number of phenotypes) X (2 X number of phenotypes)} 
#' matrix holding only the correlation between the genetic component of the full
#' liability and the full liability for the underlying individual for all
#' phenotypes. If both fam_vec and n_fam are specified, the user is asked to 
#' decide on which of the two vectors to use.
#' Note that the returned object has different attributes, such as 
#' fam_vec, n_fam, add_ind, h2 and phenotype_names.
#' 
#' @examples
#' construct_covmat_multi(fam_vec = NULL, 
#'                        corrmat = matrix(c(0.9, 0.5, 0.5, 0.8), nrow = 2), 
#'                        h2s = c(0.37,0.44),
#'                        phen_names = c("p1","p2"))
#' construct_covmat_multi(fam_vec = c("m","mgm","mgf","mhs1","mhs2","mau1"), 
#'                        n_fam = NULL, 
#'                        add_ind = TRUE,
#'                        corrmat = diag(3),
#'                        h2s = c(0.8, 0.65))
#' construct_covmat_multi(fam_vec = NULL, 
#'                        n_fam = stats::setNames(c(1,1,1,2,2), c("m","mgm","mgf","s","mhs")), 
#'                        add_ind = FALSE,
#'                        corrmat = diag(2),
#'                        h2s = c(0.75,0.85))
#' 
#' @seealso \code{\link{get_relatedness}}, \code{\link{construct_covmat_single}} and
#' \code{\link{construct_covmat}}.
#' 
#' @export
construct_covmat_multi <- function(fam_vec = c("m","f","s1","mgm","mgf","pgm","pgf"), n_fam = NULL, add_ind = TRUE, 
                                   corrmat, h2s, phen_names = NULL){
  # forcing fam_vec to be null if a vector of length 0 is provided. 
  if (length(fam_vec) == 0) fam_vec = NULL
  
  # Turn add_ind into a logical 
  add_ind <- as.logical(add_ind)
  
  # Checking that the heritabilities are valid
  if(is.null(h2s)) stop("The heritabilities must be specified!")
  if(!is.numeric(h2s) & !is.integer(h2s))stop("The heritabilities must be numeric!")
  if(any(h2s<0)) stop("The heritabilities must be non-negative!")
  if(any(h2s>1)) stop("Under the liability threshold model, the heritabilities must be smaller than or equal to 1!")
  
  # Checking that all correlations are valid
  if(is.null(corrmat)) stop("The correlation matrix corrmat must be specified!")
  if(any(abs(corrmat)>1)) stop("All correlations in cormat must be between -1 and 1!")
  # In addition, corrmat must be symmetric
  if(!isSymmetric.matrix(corrmat)) stop("The matrix corrmat must be symmetric!")
  
  # Computing the number of phenotypes
  num_phen <- length(h2s)
  
  # Checking that phen_names is either NULL or a valid
  # vector of strings
  if(is.null(phen_names)){
    phen_names <- paste0("phenotype", 1:num_phen)
  }else{
    if(!is.character(phen_names)) phen_names <- as.character(phen_names)
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
          
          covmat[2*(p1-1) + 1:2, 2*(p2-1) + 1:2] <- matrix(c(h2s[p1], h2s[p1], h2s[p1], 1), nrow = 2)
        }else{
          
          covmat[2*(p1-1) + 1:2, 2*(p2-1) + 1:2] <- matrix(sqrt(h2s[p1]*h2s[p2])*corrmat[p1,p2], nrow = 2, ncol = 2)
        }
      }
    }
    # Changing the row and column names
    colnames(covmat) <- rownames(covmat) <- paste0(c("g_", "o_"), rep(phen_names, each = 2))
    
    # Adding attributes to covmat
    attributes(covmat)$fam_vec <- NULL
    attributes(covmat)$n_fam <- NULL
    attributes(covmat)$add_ind <- add_ind
    attributes(covmat)$h2 <- h2s
    attributes(covmat)$corrmat <- corrmat
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
    # Extracting the names which are equal to the desired family members
    fam <- names(n_fam)
    # Checking that all family members are valid strings
    check_family_role(fam)

    # Checking that all entries in n_fam are non-negative
    if(any(n_fam<0)) stop("All entries in n_fam must be non-negative!")
    # Removing all family members that occur zero times
    n_fam <- n_fam[n_fam >0]
    # Constructing a vector holding all family members
    # (including g and o if add_ind = T).
    if(add_ind){
      n_fam <- c(stats::setNames(c(1,1), c("g", "o")), n_fam)
    }
    # Extracting all family members that can occur exactly once
    single_indx <- which(stringr::str_detect(names(n_fam), "^[gomf]$") | stringr::str_detect(names(n_fam), "^[mp]g[mf]$"))
    # Extracting all family members that can occur multiple times
    multiple_indx <- setdiff(1:length(n_fam), single_indx)
    # Getting the names for all family members that
    # occur only once
    if(length(single_indx)==0){
      fam <- c()
    }else{
      fam <- names(n_fam)[single_indx]
    }
    # And the names for those family members that occur
    # several times
    if(length(multiple_indx) > 0){
      for(indx in multiple_indx){
        
        fam <- c(fam, paste0(names(n_fam)[indx],"",1:n_fam[indx]))
      }
    }
    
    # Now that we have a vector holding the desired family
    # members, we can create the covariance matrix.
    covmat <- matrix(NA, nrow = length(fam)*num_phen, ncol = length(fam)*num_phen)
    # Changing the row and column names
    rownames(covmat) <- colnames(covmat) <- paste0(fam,"_", rep(phen_names, each = length(fam)))
    
    # Filling in all entries
    for(p1 in 1:num_phen){
      for(p2 in 1:num_phen){
        for(mem in fam){
          
          if(p1==p2){
            
            covmat[which(rownames(covmat) == paste0(mem, "_", phen_names[p1])), (p2-1)*length(fam) + 1:length(fam)] <- sapply(fam, get_relatedness, s1 = mem, h2 = h2s[p1])
          }else{
           
            covmat[which(rownames(covmat) == paste0(mem, "_", phen_names[p1])), (p2-1)*length(fam) + 1:length(fam)] <- sapply(fam, get_relatedness, s1 = mem, h2 = sqrt(h2s[p1]*h2s[p2])*corrmat[p1,p2]) 
          }
         
        }
      }
    }
    
    # Adding attributes to covmat
    attributes(covmat)$fam_vec <- NULL
    attributes(covmat)$n_fam <- n_fam
    attributes(covmat)$add_ind <- add_ind
    attributes(covmat)$h2 <- h2s
    attributes(covmat)$corrmat <- corrmat
    attributes(covmat)$phenotype_names <- phen_names
    
    return(covmat)
    
  }else if(is.null(n_fam)){
    # If only fam_vec is specified, the function uses this
    # vector to create the covariance matrix
    
    # Checking that all family members are represented by valid strings
    check_family_role(fam_vec)

    # If add_ind = T, the genetic component and the full
    # liability are added to the family members
    if(add_ind){
      fam_vec <- c(c("g", "o"), fam_vec)
    }
    
    # Now that we have a vector holding the desired family
    # members, we can build the covariance matrix.
    covmat <- matrix(NA, nrow = length(fam_vec)*num_phen, ncol = length(fam_vec)*num_phen)
    # Changing the row and column names
    rownames(covmat) <- colnames(covmat) <- paste0(fam_vec, "_", rep(phen_names, each = length(fam_vec)))
      
    # Filling in all entries 
    for(p1 in 1:num_phen){
      for(p2 in 1:num_phen){
        
        for(mem in fam_vec) {
          
          if(p1==p2){
            
            covmat[which(rownames(covmat) == paste0(mem,"_", phen_names[p1])), (p2-1)*length(fam_vec) + 1:length(fam_vec)] <- sapply(fam_vec, get_relatedness, s1 = mem, h2 = h2s[p1])
          }else{
            
            covmat[which(rownames(covmat) == paste0(mem,"_", phen_names[p1])), (p2-1)*length(fam_vec) + 1:length(fam_vec)] <- sapply(fam_vec, get_relatedness, s1 = mem, h2 = sqrt(h2s[p1]*h2s[p2])*corrmat[p1,p2])
          }
        }
      }
    }
    
    # Adding attributes to covmat
    attributes(covmat)$fam_vec <- fam_vec
    attributes(covmat)$n_fam <- NULL
    attributes(covmat)$add_ind <- add_ind
    attributes(covmat)$h2 <- h2s
    attributes(covmat)$corrmat <- corrmat
    attributes(covmat)$phenotype_names <- phen_names
    
    return(covmat)
  }
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
#' a given number of family members. If h2 is a number,
#' each entry in this covariance matrix equals the percentage 
#' of shared DNA between the corresponding individuals times 
#' the heritability \eqn{h^2}. However, if h2 is a matrix,
#' each entry equals either the percentage of shared DNA between the corresponding 
#' individuals times the heritability \eqn{h^2} or the 
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
#' - mgp\[1-2\]* (Maternal grandparent)
#' - pgm (Paternal grandmother)
#' - pgf (Paternal grandfather)
#' - pgp\[1-2\]* (Paternal grandparent)
#' - s\[0-9\]* (Full siblings)
#' - c\[0-9\]* (children)
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
#' @param h2 Either a number representing the heritability 
#' on liability scale for one single phenotype or a numeric vector representing
#' the heritabilities on liability scale for a positive number of phenotypes.
#' All entries in h2 must be non-negative and at most 1.
#' @param corrmat Either NULL or a numeric matrix holding the genetic correlation between the desired 
#' number of phenotypes as well as the full correlation. 
#' The full correlations must be given on the diagonal,
#' while the off-diagonal entries must hold the genetic correlations between phenotypes.
#' All correlations must be between -1 and 1.
#' Defaults to NULL.
#' @param phen_names A character vector holding the phenotype names. These names
#' will be used to create the row and column names for the covariance matrix.
#' If it is not specified, the names will default to phenotype1, phenotype2, etc.
#' Defaults to NULL.
#' 
#' @return If either fam_vec or n_fam is used as the argument, if it is of 
#' the required format, if h2 is either a number or a numeric vector satisfying 
#' that all entries are non-negative and corrmat is a numeric matrix satisfying that all 
#' entries are between -1 and 1, then the output will be a named covariance matrix. 
#' The number of rows and columns corresponds to the number of phenotypes times 
#' the length of fam_vec or n_fam (+ 2 if add_ind=T). 
#' If both fam_vec and n_fam are equal to c() or NULL, the function returns
#' a \eqn{(2 \times number of phenotypes) \times (2 \times number of phenotypes)} 
#' matrix holding only the correlation between the genetic component of the full
#' liability and the full liability for the underlying individual for all
#' phenotypes. If both fam_vec and n_fam are specified, the user is asked to 
#' decide on which of the two vectors to use.
#' Note that the returned object has different attributes, such as 
#' fam_vec, n_fam, add_ind, h2 and phenotype_names.
#' 
#' @examples
#' construct_covmat()
#' construct_covmat(fam_vec = c("m","mgm","mgf","mhs1","mhs2","mau1"), 
#'                  n_fam = NULL, 
#'                  add_ind = TRUE, 
#'                  h2 = 0.5)
#' construct_covmat(fam_vec = NULL, 
#'                  n_fam = stats::setNames(c(1,1,1,2,2), c("m","mgm","mgf","s","mhs")), 
#'                  add_ind = FALSE, h2 = 0.3)
#' construct_covmat(h2 = c(0.5,0.5), corrmat = matrix(c(0.6,0.2,0.2,0.4), nrow = 2))
#' 
#' @seealso \code{\link{get_relatedness}}, \code{\link{construct_covmat_single}},
#' \code{\link{construct_covmat_multi}}
#' @export
construct_covmat <- function(fam_vec = c("m","f","s1","mgm","mgf","pgm","pgf"), n_fam = NULL, add_ind = TRUE, h2 = 0.5, corrmat = NULL, phen_names = NULL){
  
  if(length(h2) == 1){
    
    return(construct_covmat_single(fam_vec = fam_vec, n_fam = n_fam, add_ind = add_ind, h2 = h2))
    
  }else{
    
    return(construct_covmat_multi(fam_vec = fam_vec, n_fam = n_fam, add_ind = add_ind, corrmat = corrmat, h2s = h2, phen_names = phen_names))
  } 
}


#' Positive definite matrices
#'
#' \code{correct_positive_definite} verifies that a given covariance matrix
#' is indeed positive definite by checking that all eigenvalues are positive.
#' If the given covariance matrix is not positive definite, 
#' \code{correct_positive_definite} tries to modify the underlying correlation matrix
#' h2 in order to obtain a positive definite covariance matrix
#'
#' This function can be used to verify that a given covariance matrix 
#' is positive definite. It calculates all eigenvalues in order to
#' investigate whether they are all positive. This property is necessary 
#' if the covariance matrix should be used as a Gaussian covariance matrix.
#' It is especially useful to check whether any covariance matrix obtained
#' by \code{\link{construct_covmat_multi}} is positive definite.
#' If the given covariance matrix is not positive definite, \code{correct_positive_definite}
#' tries to modify the underlying correlation matrix (called \code{h2} in
#' \code{\link{construct_covmat}} or \code{\link{construct_covmat_multi}}) by 
#' multiplying all off-diagonal entries in the correlation matrix by a given number.
#'
#' @param covmat A symmetric and numeric matrix. If the covariance matrix 
#' should be corrected, it must have a number of attributes, such as
#' attr(covmat,"fam_vec"), attr(covmat,"n_fam"), attr(covmat,"add_ind"),
#' attr(covmat,"h2") and attr(covmat,"phenotype_names"). Any covariance
#' matrix obtained by \code{\link{construct_covmat}}, \code{\link{construct_covmat_single}} 
#' or \code{\link{construct_covmat_multi}} will have these attributes by default. 
#' 
#' @param correction_val A positive number representing the amount by which
#' h2 will be changed, if not all eigenvalues are positive. That is, correction_val
#' is the number that will be multiplied to all off_diagonal entries in h2.
#' Defaults to 0.99.
#' 
#' @param correction_limit A positive integer representing the upper limit for the correction
#' procedure. Defaults to 100.
#' 
#' @return If covmat is a symmetric and numeric matrix and all eigenvalues are
#' positive, \code{correct_positive_definite} simply returns covmat. If some 
#' eigenvalues are not positive and correction_val is a positive number, 
#' \code{correct_positive_definite} tries to convert covmat into a positive definite
#' matrix. If covmat has attributes "add_ind" and "h2", with h2 being
#' a matrix with at least two rows and columns, \code{correct_positive_definite} 
#' computes a new covariance matrix using a slightly modified h2. If the 
#' correction is performed successfully, i.e. if the new covmat is positive definite,
#' the new covariance matrix is returned. Otherwise, \code{correct_positive_definite}
#' returns the original covariance matrix.
#' 
#' @seealso \code{\link{construct_covmat}}, \code{\link{construct_covmat_single}} and
#' \code{\link{construct_covmat_multi}}.
#' 
#' @export
correct_positive_definite = function(covmat, correction_val = .99, correction_limit = 100) {
  
  # Checking that covmat is symmetric
  if(!isSymmetric.matrix(covmat)) stop("The covariance matrix covmat must be symmetric!")
  # and numeric
  if(!is.numeric(covmat)) stop("The covariance matrix covmat must be numeric!")
  
  # Checking whether all eigenvalues are positive
  if(any(eigen(covmat)$values < 0)){
    
    cat("The specified covariance matrix is not positive definite. \n")
  }else{
    
    return(covmat)
  }
  
  # If some eigenvalues are negative, correction_val must be specified,
  # it must be numeric and positive in order to correct the covariance matrix.
  if(!is.numeric(correction_val)) stop("correction_val must be numeric!")
  if(correction_val <= 0) stop("correction_val must be positive!")
  
  # In addition, covmat must have several attributes holding the 
  # family members (fam_vec or n_fam), a logical add_ind as well as
  # a numeric value or numeric matrix h2, in order to change
  # the correlation matrix.
  if(is.null(attr(covmat,"add_ind")) || is.null(attr(covmat,"h2"))){
    
    warning("The required attributes are missing... The covariance matrix could not be corrected!")
    return(covmat)
  }
  
  # If the covariance matrix is for a single phenotype, it is not
  # possible to correct the covariance matrix.
  if(length(attr(covmat,"h2"))==1){
    warning("The covariance matrix cannot be corrected...")
    return(covmat)
  }
  
  # Furthermore, correction_limit must be a positive integer
  if(!is.numeric(correction_limit)) stop("correction_limit must be numeric!")
  if(correction_limit <= 0) stop("Correction limit must be positive!")
  
  cat("Trying to correct the covariance matrix...\n")
  
  # The covariance matrix will be modified at most correction_limit times.
  n <- 0
  # Storing the old covariance matrix
  old_covmat <- covmat
  # The diagonal of h2 will remain unchanged
  h2s <- attr(covmat,"h2")
  corrmat <- attr(covmat,"corrmat")
  
  # We also extract the vectors holding the family members
  fam_vec <- setdiff(attr(covmat,"fam_vec"), c("g","o"))
  n_fam <- attr(covmat,"n_fam")[stringr::str_detect(names(attr(covmat,"n_fam")), "^[^go]")]
  
  while(any(eigen(covmat)$values < 0) && n <= correction_limit) {
    
    # Changing the corrmat matrix slightly by 
    # multiplying all entries by correction_val.
    corrmat <- corrmat*correction_val
    # Computing a new covariance matrix
    covmat <- construct_covmat(fam_vec = fam_vec, n_fam = n_fam, add_ind = attr(covmat,"add_ind"), corrmat = corrmat, h2 = h2s, phen_names = attr(covmat,"phenotype_names"))
    # Updating 
    n <- n+1
  }
  
  # If the matrix has been modified correction_limit times, 
  # the correction step is aborted
  if(n > correction_limit){
    cat("The covariance matrix could not be corrected. Consider revisiting it.\n")
    return(old_covmat)
  }
  
  cat(paste0("The correction was performed successfully! All off-diagonal entries are corrected by", round(correction_val^n,digits = 3),".\n"))

  return(covmat)
}
