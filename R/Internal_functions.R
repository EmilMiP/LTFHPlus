utils::globalVariables("role")
#' Checking that relatives are represented by valid strings
#'
#' \code{validate_relatives} checks whether relatives are represented
#' by valid abbreviations.
#'
#' This function can be used to check whether relatives are represented
#' by valid abbreviations. A valid abbreviation is one of the following:
#' - \code{g} (Genetic component of full liability)
#' - \code{o} (Full liability)
#' - \code{m} (Mother)
#' - \code{f} (Father)
#' - \code{c[0-9]*.[0-9]*} (Children)
#' - \code{mgm} (Maternal grandmother)
#' - \code{mgf} (Maternal grandfather)
#' - \code{pgm} (Paternal grandmother)
#' - \code{pgf} (Paternal grandfather)
#' - \code{s[0-9]*} (Full siblings)
#' - \code{mhs[0-9]*} (Half-siblings - maternal side)
#' - \code{phs[0-9]*} (Half-siblings - paternal side)
#' - \code{mau[0-9]*} (Aunts/Uncles - maternal side)
#' - \code{pau[0-9]*} (Aunts/Uncles - paternal side).
#'
#' @param relatives A string or character vector representing 
#' the relatives.
#' All strings must be chosen among the following abbreviations
#' - \code{g} (Genetic component of full liability)
#' - \code{o} (Full liability)
#' - \code{m} (Mother)
#' - \code{f} (Father)
#' - \code{c[0-9]*.[0-9]*} (Children)
#' - \code{mgm} (Maternal grandmother)
#' - \code{mgf} (Maternal grandfather)
#' - \code{pgm} (Paternal grandmother)
#' - \code{pgf} (Paternal grandfather)
#' - \code{s[0-9]*} (Full siblings)
#' - \code{mhs[0-9]*} (Half-siblings - maternal side)
#' - \code{phs[0-9]*} (Half-siblings - paternal side)
#' - \code{mau[0-9]*} (Aunts/Uncles - maternal side)
#' - \code{pau[0-9]*} (Aunts/Uncles - paternal side)
#' for the function to return TRUE.
#' 
#' @return If \code{relatives} is a string or character vector such that
#' all strings are chosen from the mentioned list of strings,
#' then the function will return TRUE. Otherwise, the function is aborted.
#' 
#' @examples
#' LTFHPlus:::validate_relatives("g")
#' LTFHPlus:::validate_relatives("o")
#' LTFHPlus:::validate_relatives("mgm")
#' 
#' # This will result in errors:
#' try(LTFHPlus:::validate_relatives("a"))
#' try(LTFHPlus:::validate_relatives(m))
#' 
#' @importFrom stringr str_detect
#' @noRd
validate_relatives <- function(relatives){
  
  if(!is.character(relatives)){
    
    stop(paste0(deparse(substitute(relatives)), " must be a string or character vector!"))
    
  }else if(any(!(str_detect(relatives, "^[gomf]$") | str_detect(relatives, "^c[0-9]*.[0-9]*")|
                 str_detect(relatives, "^[mp]g[mf]$") | str_detect(relatives, "^s[0-9]*") | 
                 str_detect(relatives, "^[mp]hs[0-9]*")| str_detect(relatives, "^[mp]au[0-9]*")))){
    
    stop(paste0(deparse(substitute(relatives)), " contains invalid abbreviations! Use strings from the following list: \n
  - g (Genetic component of full liability)\n
  - o (Full liability)\n
  - m (Mother)\n
  - f (Father)\n
  - c[0-9]*.[0-9]* (Children)\n
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
#' @param from_covmat logical variable. Only used internally. allows for skip of negative check.
#' 
#' @return If \code{prop} is a vector holding valid proportions of class \code{numeric}
#' or \code{integer} that are non-negative and at most one,
#' then the function will return TRUE. Otherwise, the function aborts.
#' 
#' @examples
#' LTFHPlus:::validate_proportion(0.2)
#' LTFHPlus:::validate_proportion(0.04)
#' LTFHPlus:::validate_proportion(0)
#' LTFHPlus:::validate_proportion(1)
#' 
#' # This will result in errors:
#' try(LTFHPlus:::validate_proportion(2))
#' try(LTFHPlus:::validate_proportion(-0.5))
#' @noRd
validate_proportion <- function(prop, from_covmat = FALSE){
  if(is.null(prop)){
    
    stop(paste0(deparse(substitute(prop)), " must be specified!"))
    
  }else if(!is.numeric(prop) && !is.integer(prop)){
    
    stop(paste0(deparse(substitute(prop)), " must be numeric!"))
    
  }else if(any(prop<0) & !from_covmat){
    
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
#' LTFHPlus:::validate_correlation_matrix(matrix(c(1,0.4,0.4,1), nrow = 2))
#' LTFHPlus:::validate_correlation_matrix(diag(3))
#' 
#' # This will result in errors:
#' try(LTFHPlus:::validate_correlation_matrix(matrix(c(0.2,0.4,0.4,0.2), nrow = 2)))
#' try(LTFHPlus:::validate_correlation_matrix(matrix(nrow=2, ncol = 2)))
#' @noRd
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
#' @param phen_name Either \code{NULL} or character vector holding the 
#' phenotype name. Must be specified in the multi-trait case.
#' Defaults to \code{NULL}.
#' 
#' @return A tibble holding all columns present in .tbl as well
#' as the age of onset or the current age
#' (depending on the disease status) for all individuals 
#' given in \code{fam_mem}. 
#' 
#' @importFrom dplyr %>% rowwise select mutate bind_cols
#' @importFrom rlang :=
#' @noRd
construct_aoo <- function(fam_mem,.tbl, pop_prev, phen_name = NULL){
  
  # Removing the genetic component from the 
  # set of family members, if it is present
  i_ind <- setdiff(fam_mem, c("g"))
  
  if(is.null(phen_name)){
    
    # Looping over all family members ind i_ind
    lapply(i_ind, function(j){
      
      # Selecting the liability, disease status and age for 
      # individual j, in order to compute the age of onset.
      select(.tbl, c(tidyselect::matches(paste0("^",j,"$")), tidyselect::matches(paste0("^",j,"_[as].*$")))) %>%
        rowwise() %>% 
        mutate(., !!as.symbol(paste0(j,"_aoo")) := ifelse(!!as.symbol(paste0(j,"_status")), 
                                                          round(convert_liability_to_aoo(!!as.symbol(j), dist = "logistic", pop_prev = pop_prev, mid_point = 60, slope = 1/8)),
                                                          !!as.symbol(paste0(j,"_age")))) %>%
        select(., !!as.symbol(paste0(j,"_aoo")))
    }
    ) %>% do.call("bind_cols",.) %>% bind_cols(.tbl,.)
  
  }else{
    
    # Looping over all family members ind i_ind
    lapply(i_ind, function(j){
      
      # Selecting the liability, disease status and age for 
      # individual j, in order to compute the age of onset.
      select(.tbl, tidyselect::starts_with(paste0(j, "_"))) %>%
        rowwise() %>% 
        mutate(., !!as.symbol(paste0(j,"_", phen_name ,"_aoo")) := ifelse(!!as.symbol(paste0(j, "_", phen_name, "_status")), 
                                                                          round(convert_liability_to_aoo(!!as.symbol(paste0(j, "_", phen_name)), dist = "logistic", pop_prev = pop_prev, mid_point = 60, slope = 1/8)),
                                                                          !!as.symbol(paste0(j,"_age")))) %>%
        select(., !!as.symbol(paste0(j,"_", phen_name ,"_aoo")))
    }
    ) %>% do.call("bind_cols",.) %>% bind_cols(.tbl,.)
  }
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
#' @param phen_name Either \code{NULL} or character vector holding the 
#' phenotype name. Must be specified in the multi-trait case.
#' Defaults to \code{NULL}.
#' 
#' @return A tibble holding the personal identifier (PID) as well as 
#' the lower and the upper threshold for all individuals
#' present in \code{fam_mem}.
#' 
#' @importFrom dplyr %>% rowwise select mutate bind_rows ungroup
#' @noRd
construct_thresholds <- function(fam_mem, .tbl, pop_prev, phen_name = NULL){
  
  # Removing the genetic component from the 
  # set of family members, if it is present
  i_ind <- setdiff(fam_mem, c("g"))
  
  if(!is.null(phen_name)){
    
    # Looping over all family members ind i_ind
    lapply(i_ind, function(j){
      
      nbr <- which(i_ind == j)
      
      # Selecting the family ID, disease status and age/aoo for 
      # individual j, in order to compute the thresholds.
      select(.tbl, c(fam_ID, 
                     tidyselect::matches(paste0(j, "_", phen_name, "_status")), 
                     tidyselect::matches(paste0(j, "_", phen_name, "_aoo")))) %>%
        rowwise() %>% 
        mutate(., indiv_ID = paste0(fam_ID,"_", nbr),
                  role = paste0(j),
                  upper = convert_age_to_thresh(!!as.symbol(paste0(j, "_", phen_name, "_aoo")), dist = "logistic", pop_prev = pop_prev, mid_point = 60, slope = 1/8), 
                  lower = ifelse(!!as.symbol(paste0(j, "_", phen_name, "_status")), 
                                  convert_age_to_thresh(!!as.symbol(paste0(j, "_", phen_name, "_aoo")), dist = "logistic", pop_prev = pop_prev, mid_point = 60, slope = 1/8),
                                  -Inf)) %>%
        rename(., !!as.symbol(paste0("lower_", phen_name)) := lower, !!as.symbol(paste0("upper_", phen_name)) := upper) %>% 
        select(., fam_ID, indiv_ID, role, starts_with("lower"), starts_with("upper")) %>% 
        ungroup()
      
    }) %>% do.call("bind_rows",.)
    
  }else{
    
    # Looping over all family members ind i_ind
    lapply(i_ind, function(j){
      
      nbr <- which(i_ind == j)
      
      # Selecting the family ID, disease status and age/aoo for 
      # individual i, in order to compute the thresholds.
      select(.tbl, c(fam_ID, 
                     tidyselect::matches(paste0("^",j,"_status$")), 
                     tidyselect::matches(paste0("^",j,"_aoo$")))) %>%
        rowwise() %>% 
        mutate(., indiv_ID = paste0(fam_ID,"_", nbr),
                  role = paste0(j),
                  upper = convert_age_to_thresh(!!as.symbol(paste0(j,"_aoo")), dist = "logistic", pop_prev = pop_prev, mid_point = 60, slope = 1/8), 
                  lower = ifelse(!!as.symbol(paste0(j,"_status")), 
                                  convert_age_to_thresh(!!as.symbol(paste0(j,"_aoo")), dist = "logistic", pop_prev = pop_prev, mid_point = 60, slope = 1/8),
                                  -Inf)) %>%
        select(., fam_ID, indiv_ID, role, starts_with("lower"), starts_with("upper")) %>% 
        ungroup()
      
    }) %>% do.call("bind_rows",.)
  }
}

#' 
#' Add in missing roles for proband
#' 
#' This function adds missing roles for a proband and fills lower threshold values with -Inf and upper threshold values with Inf.
#' 
#' @param temp_tbl tibble to add missing proband roles to; originates from .tbl of estimate_liability.
#' @param role name of role column
#' @param cur_roles values of the role column
#' @param cur_fam_id current family ID being worked on
#' @param pid name of column with personal IDs
#' @param fam_id name of column with family IDs
#' @param phen_names vector of phenotype names as given in .tbl of estimate_liability. Defaults to NULL (which is single trait).
#' 
#' @return The provided temp_tbl object is returned, but with the missing "g" and/or "o" roles added, where -Inf and Inf values
#' have been used to fill the lower and upper threshold values. If phen_names is provided, a pair of upper and lower values is 
#' provided for each entry in phen_names.
#' 
#' @importFrom dplyr filter pull tibble %>% bind_rows
#' @noRd

add_missing_roles_for_proband = function(temp_tbl, role, cur_roles, cur_fam_id, pid, fam_id, phen_names = NULL) {
  # role types to check for, centered on proband
  to_check_for = c("g", "o")
  
  # roles is already calculated; are they present?
  to_be_added = setdiff(to_check_for,   cur_roles)
  present     = intersect(to_check_for, cur_roles)
  
  # if some present, extract individual ID, if not, get family ID
  if (length(present) > 0 ) {
    i_pid = (temp_tbl %>% filter(!!as.symbol(role) == present) %>% pull(!!as.symbol(pid)))[1]
  } else {
    i_pid = pull(temp_tbl, !!as.symbol(fam_id))[1]
  }
  # suffixes of roles to be added
  id_suffixes = paste0("_",to_be_added) %>% stringr::str_replace_all(., "_o", "")
  
  if ( is.null(phen_names) ) { # single trait
    # construct tibble with desired roles
    tibble(
      !!as.symbol(fam_id) := pull(temp_tbl, !!as.symbol(fam_id))[1],
      !!as.symbol(pid)    := paste0(i_pid, id_suffixes),
      !!as.symbol(role)   := to_be_added,
      lower = rep(-Inf, length(to_be_added)),
      upper = rep( Inf, length(to_be_added))
    ) %>%
      bind_rows(., temp_tbl)
    
  } else { # multi trait
    # constructs id rows, then adds lower and upper thresholds from phen_names provided
    tibble(
      !!as.symbol(fam_id) := pull(temp_tbl, !!as.symbol(fam_id))[1],
      !!as.symbol(pid)    := paste0(i_pid, id_suffixes),
      !!as.symbol(role)   := to_be_added
    ) %>%
      bind_cols(
        tibble(!!!c(stats::setNames(rep(-Inf, length(phen_names)), paste0("lower_", phen_names)),
                    stats::setNames(rep( Inf, length(phen_names)), paste0("upper_", phen_names))))) %>%
      bind_rows(
        .,
        temp_tbl
      )
  }
}