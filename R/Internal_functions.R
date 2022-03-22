#' Checking that a relative is represented by a valid string
#'
#' \code{check_valid_relatives} checks whether a relative is represented
#' by a valid abbreviation.
#'
#' This function can be used to check whether a relative is represented
#' by a valid abbreviation. A valid abbreviation is one of the following:
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
#'
#' @param relative A string representing a relative.
#' The string must be one of the following
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
#' - pau\[0-9\]* (Aunts/Uncles - paternal side)
#' for the function to return TRUE.
#' 
#' @return If \code{relative} is a string chosen from the mentioned list of strings,
#' then the function will return TRUE. Otherwise, the function is aborted.
#' 
#' @examples
#' check_valid_relatives("g")
#' check_valid_relatives("o")
#' check_valid_relatives("mgm")
#'
#'
#' \dontrun{
#' # This will result in errors:
#' get_relatedness("a")
#' get_relatedness(m)
#' }
#' 
#' @importFrom stringr str_detect
#' 
check_valid_relatives <- function(relative){
  
  if(class(relative) != "character"){
    stop(paste0(deparse(substitute(relative)), " must be a string!"))
  }else if(!(str_detect(relative, "^[gomf]$") | str_detect(relative, "^[mp]g[mf]$") | 
       str_detect(relative, "^s[0-9]*") | str_detect(relative, "^[mp]hs[0-9]*")| 
       str_detect(relative, "^[mp]au[0-9]*"))){
       
       stop(paste0(deparse(substitute(relative)), " is not a valid string! Use a string from the following list: \n
  - g (Genetic component of full liability)\n
  - o (Full liability)\n
  - m (Mother)\n
  - f (Father)\n
  - mgm (Maternal grandmother)\n
  - mgf (Maternal grandfather)\n
  - pgm (Paternal grandmother)\n
  - pgf (Paternal grandfather)\n
  - s[0-9]* (Full siblings)\n
  - mhs[0-9]* (Half-siblings - maternal side)\n
  - phs[0-9]* (Half-siblings - paternal side)\n
  - mau[0-9]* (Aunts/Uncles - maternal side)\n
  - pau[0-9]* (Aunts/Uncles - paternal side)."))
  }else{
    return(TRUE)
  }
}

#' Checking that a proportion is valid
#'
#' \code{check_proportion} checks whether a proportion is valid, i.e.
#' whether it is non-negative and at most one
#'
#' This function can be used to check whether a proportion is non-negative
#' and at most one
#'
#' @param prop A number or integer representing the proportion that 
#' needs to be validated.
#' 
#' @return If \code{prop} is a valid proportion that is of class \code{numeric}
#' or \code{integer}, non-negative and at most one,
#' then the function will return TRUE. Otherwise, the function aborts.
#' 
#' @examples
#' check_proportion(0.2)
#' check_proportion(0.04)
#' check_proportion(0)
#' check_proportion(1)
#'
#'
#' \dontrun{
#' # This will result in errors:
#' check_proportion(2)
#' check_proportion(-0.5)
#' }
#' 
check_proportion <- function(prop){
  
  if(is.null(prop)){
    stop(paste0(deparse(substitute(prop)), " must be specified!"))
  }else if(class(prop)!= "numeric" && class(prop)!= "integer"){
    stop(paste0(deparse(substitute(prop)), " must be numeric!"))
  }else if(prop<0){
    stop(paste0(deparse(substitute(prop)), " must be non-negative!"))
  }else if(prop>1){
    stop(paste0(deparse(substitute(prop)), " must be smaller than or equal to 1!"))
  }else{
    return(TRUE)
  }
}


#' Checking that a correlation matrix is valid
#'
#' \code{check_correlation_matrix} checks whether a matrix is a valid 
#' correlation matrix, i.e. whether its diagonal entries are equal to one, 
#' while all off-diagonal entries are between -1 and 1, and whether it is
#' symmetric.
#'
#' This function can be used to check whether a correlation matrix has diagonal
#' entries equal to 1 and off-diagonal entries between -1 and 1 as well as whether 
#' it is symmetric.
#'
#' @param corrmat A numeric matrix holding the correlations. 
#' All diagonal entries must be equal to one, while all off-diagonal entries 
#' must be between -1 and 1. In addition, the matrix must be symmetric.
#' 
#' @return If \code{corrmat} is a valid correlation matrix that is symmetric,
#' has one on all diagonal entries and numbers between -1 and 1 on all off-
#' diagonal entries,
#' then the function will return TRUE. Otherwise, the function will be aborted.
#' 
#' @examples
#' check_correlation_matrix(matrix(c(1,0.4,0.4,1), nrow = 2))
#' check_correlation_matrix(diag(3))
#'
#'
#' \dontrun{
#' # This will result in errors:
#' check_correlation_matrix(matrix(c(0.2,0.4,0.4,0.2), nrow = 2))
#' check_correlation_matrix(matrix(nrow=2, ncol = 2))
#' }
#' 
check_correlation_matrix <- function(corrmat){
  
  if(is.null(corrmat)){
    stop(paste0(deparse(substitute(corrmat)), " must be specified!"))
  }else if(any(diag(corrmat)!= 1)){
    stop(paste0("All diagonal entries in ", deparse(substitute(corrmat))," must be 1!"))
  }else if(any(abs(corrmat)>1)){
    stop(paste0("All off-diagonal entries in ", deparse(substitute(corrmat))," must be between -1 and 1!"))
  }else if(!isSymmetric.matrix(corrmat)){
    stop(paste0(deparse(substitute(corrmat)), " must be symmetric!"))
  }else{
    return(TRUE)
  }
}
