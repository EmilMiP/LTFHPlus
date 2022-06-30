#' Checking that relatives are represented by valid strings
#'
#' \code{validate_relatives} checks whether relatives are represented
#' by valid abbreviations.
#'
#' This function can be used to check whether relatives are represented
#' by valid abbreviations. A valid abbreviation is one of the following:
#' - g (Genetic component of full liability)
#' - o (Full liability)
#' - m (Mother)
#' - f (Father)
#' - c\[0-9\]* (Children)
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
#' @param relatives A string or character vector representing 
#' the relatives.
#' All strings must be chosen among the following abbreviations
#' - g (Genetic component of full liability)
#' - o (Full liability)
#' - m (Mother)
#' - f (Father)
#' - c\[0-9\]* (Children)
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
#' @return If \code{relatives} is a string or character vector such that
#' all strings are chosen from the mentioned list of strings,
#' then the function will return TRUE. Otherwise, the function is aborted.
#' 
#' @examples
#' \dontrun{
#' validate_relatives("g")
#' validate_relatives("o")
#' validate_relatives("mgm")
#' 
#' # This will result in errors:
#' validate_relatives("a")
#' validate_relatives(m)
#' }
#' 
#' @importFrom stringr str_detect
#' 
validate_relatives <- function(relatives){
  
  if(!is.character(relatives)){
    
    stop(paste0(deparse(substitute(relatives)), " must be a string or character vector!"))
    
  }else if(any(!(str_detect(relatives, "^[gomf]$") | str_detect(relatives, "^c[0-9]*")|
                 str_detect(relatives, "^[mp]g[mf]$") | str_detect(relatives, "^s[0-9]*") | 
                 str_detect(relatives, "^[mp]hs[0-9]*")| str_detect(relatives, "^[mp]au[0-9]*")))){
    
    stop(paste0(deparse(substitute(relatives)), " contains invalid abbreviations! Use strings from the following list: \n
  - g (Genetic component of full liability)\n
  - o (Full liability)\n
  - m (Mother)\n
  - f (Father)\n
  - c[0-9]* (Children)\n
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

#' Checking that proportions are valid
#'
#' \code{validate_proportion} checks whether proportions are valid, i.e.
#' whether they are non-negative and at most one.
#'
#' This function can be used to check whether proportions are non-negative
#' and at most one.
#'
#' @param prop A number, integer or numeric vector representing the proportions that 
#' need to be validated.
#' 
#' @return If \code{prop} is a vector holding valid proportions of class \code{numeric}
#' or \code{integer} that are non-negative and at most one,
#' then the function will return TRUE. Otherwise, the function aborts.
#' 
#' @examples
#' \dontrun{
#' validate_proportion(0.2)
#' validate_proportion(0.04)
#' validate_proportion(0)
#' validate_proportion(1)
#' 
#' # This will result in errors:
#' validate_proportion(2)
#' validate_proportion(-0.5)
#' }
#' 
validate_proportion <- function(prop){
  
  if(is.null(prop)){
    
    stop(paste0(deparse(substitute(prop)), " must be specified!"))
    
  }else if(!is.numeric(prop) && !is.integer(prop)){
    
    stop(paste0(deparse(substitute(prop)), " must be numeric!"))
    
  }else if(any(prop<0)){
    
    stop(paste0(deparse(substitute(prop)), " must be non-negative!"))
    
  }else if(any(prop>1)){
    
    stop(paste0(deparse(substitute(prop)), " must be smaller than or equal to 1!"))
    
  }else{
    return(TRUE)
  }
}


#' Checking that a correlation matrix is valid
#'
#' \code{validate_correlation_matrix} checks whether a matrix is a valid 
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
#' \dontrun{
#' validate_correlation_matrix(matrix(c(1,0.4,0.4,1), nrow = 2))
#' validate_correlation_matrix(diag(3))
#' 
#' # This will result in errors:
#' validate_correlation_matrix(matrix(c(0.2,0.4,0.4,0.2), nrow = 2))
#' validate_correlation_matrix(matrix(nrow=2, ncol = 2))
#' }
#' 
validate_correlation_matrix <- function(corrmat){
  
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


#' Constructing age of onset (aoo)
#'
#' \code{construct_aoo} constructs the age of onset (aoo)
#' for a variable number of family members based on their 
#' liability, disease status and current age.
#'
#' @param fam_mem A character vector holding all family members.
#' @param .tbl A tibble holding the liability as well as age and
#' disease status for the set of individuals in \code{fam_mem}.
#' @param pop_prev A positive number representing the population prevalence, i.e. the 
#' overall prevalence in the population.
#' 
#' @return A tibble holding all columns present in .tbl as well
#' as the age of onset or the current age
#' (depending on the disease status) for all individuals 
#' given in \code{fam_mem}. 
#' 
#' @importFrom dplyr %>% rowwise select mutate bind_cols
#' @importFrom rlang :=
construct_aoo <- function(fam_mem,.tbl, pop_prev){
  
  # Removing the genetic component from the 
  # set of family members, if it is present
  i_ind <- setdiff(fam_mem, c("g"))
  
  # Looping over all family members ind i_ind
  lapply(i_ind, function(i){
    
    # Selecting the liability, disease status and age for 
    # individual i, in order to compute the age of onset.
    select(.tbl, c(tidyselect::matches(paste0("^",i,"$")), tidyselect::matches(paste0("^",i,"_[as].*$")))) %>%
      rowwise() %>% 
      mutate(., !!as.symbol(paste0(i,"_aoo")) := ifelse(!!as.symbol(paste0(i,"_status")), 
                                                        round(convert_liability_to_aoo(!!as.symbol(i), dist = "logistic", pop_prev = pop_prev, mid_point = 60, slope = 1/8)),
                                                        !!as.symbol(paste0(i,"_age")))) %>%
      select(., !!as.symbol(paste0(i,"_aoo")))
  }
  ) %>% do.call("bind_cols",.) %>% bind_cols(.tbl,.)
}


#' Computing thresholds
#'
#' \code{construct_thresholds} computes the upper and lower 
#' thresholds for a variable number of family members based on their 
#' disease status and current age or age of onset (depending on 
#' the disease status).
#'
#' @param fam_mem A character vector holding all family members.
#' @param .tbl A tibble holding the family ID, disease status as well
#' as the age of onset or the current age
#' (depending on the disease status).
#' @param pop_prev A positive number representing the population prevalence, i.e. the 
#' overall prevalence in the population.
#' 
#' @return A tibble holding the personal identifier (PID) as well as 
#' the lower and the upper threshold for all individuals
#' present in \code{fam_mem}.
#' 
#' @importFrom dplyr %>% rowwise select mutate bind_rows ungroup
construct_thresholds <- function(fam_mem, .tbl, pop_prev){
  
  # Removing the genetic component from the 
  # set of family members, if it is present
  i_ind <- setdiff(fam_mem, c("g"))
  
  # Looping over all family members ind i_ind
  lapply(seq_along(i_ind), function(nbr){
    i = i_ind[nbr]
    
    # Selecting the family ID, disease status and age/aoo for 
    # individual i, in order to compute the thresholds.
    select(.tbl, c(tidyselect::matches(paste0("^indiv_ID$")), tidyselect::matches(paste0("^",i,"_status$")), tidyselect::matches(paste0("^",i,"_aoo$")))) %>%
      rowwise() %>% 
      mutate(., fam_ID = paste0(indiv_ID,"_member", nbr), 
             upper = convert_age_to_thresh(!!as.symbol(paste0(i,"_aoo")), dist = "logistic", pop_prev = pop_prev, mid_point = 60, slope = 1/8), 
             lower = ifelse(!!as.symbol(paste0(i,"_status")), 
                            convert_age_to_thresh(!!as.symbol(paste0(i,"_aoo")), dist = "logistic", pop_prev = pop_prev, mid_point = 60, slope = 1/8),
                            -Inf)) %>%
      select(., fam_ID, lower, upper) %>% 
      ungroup()
    
  }) %>% do.call("bind_rows",.)
}