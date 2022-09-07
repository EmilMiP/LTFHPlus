utils::globalVariables("cip")

#' Attempts to convert the list entry input format to a long format
#' 
#' @param family a tibble with two entries, family id and personal id. personal id should end in "_role", if a role column is not present.
#' @param threshs thresholds, with a personal id (without role) as well as the lower and upper thresholds
#' @param personal_id_col column name that holds the personal id
#' @param role_col column name that holds the role
#' 
#' @return returns a format similar to \code{prepare_LTFHPlus_input}, which is used by \code{estimate_liability}
#'
#' @export 

convert_format = function(family, threshs, personal_id_col = "pid", role_col = NULL) {
  # standardising input -----------------------------------------------------
  
  #are there any list columns in family (list entry format)?
  which_list_columns = which(sapply(family, is.list))
  if (length(which_list_columns) > 0) { 
    #updating fam to get rid of list columns
    family = tidyr::unnest(family, cols = names(which_list_columns))
    
    ###  checking if role is present in ID or separate column
    # if "_" is present, a role will be there too.
    if (any(stringr::str_detect(family[[personal_id_col]], "_"))) { #if true, extract role
      #split pid_role in two:
      family[[role_col]] = strsplit(family[[personal_id_col]], "_(?=[^_]+$)", perl=TRUE) %>% sapply(., function(x) x[2])
      family[[personal_id_col]] = strsplit(family[[personal_id_col]], "_(?=[^_]+$)", perl=TRUE) %>% sapply(., function(x) x[1])
      
      cat("We've tried converting from list entries to a long format internally. If you see this print, please run prepare_LTFHPlus_input and use the .tbl input going forward! \n")
    } else {
      if (is.null(role_col)) stop("Please provide family roles for each family member. e.g. father(f), mother(m), siblings (s1-s9), etc.") 
      stop("We weren't able to convert data input automatically. Please use prepare_LTFHPlus_input and use the .tbl input!\n")
    }
  }
  .tbl = left_join(family, threshs, by = personal_id_col)
  return(.tbl)
}

#' Prepares input for \code{estimate_liability}
#' 
#' @param family contains family and personal ids and role with a family.
#' @param CIP tibble with population representative cumulative incidence proportions has the interpretation of "proportion of people that has experienced the trait subset by \code{CIP_columns}. 
#' @param CIP_columns the columns the CIPs are subset by, e.g. CIPs by birth_year, sex. 
#' @param status_col Column that contains the status of each family member 
#' @param use_fixed_case_thr Should the threshold be fixed for cases? Can be used if CIPs are detailed, e.g. stratified by birth_year and sex.
#' @param fam_id_col Column that contains the family ID
#' @param personal_id_col Column that contains the personal ID
#' @param role_col Column that cotnains the role of each individual
#' 
#' @export

prepare_LTFHPlus_input = function(family, CIP, 
                                  CIP_columns = c("sex", "birth_year", "age"), 
                                  status_col = "status", 
                                  use_fixed_case_thr = F, 
                                  fam_id_col = "fam_id", 
                                  personal_id_col = "pid",
                                  role_col = "role") {

  

# Merging CIPs and assigning thresholds -----------------------------------

  family = dplyr::left_join(family, CIP, by = CIP_columns) %>% 
    dplyr::mutate(
      lower = ifelse(!!as.symbol(status_col) == 1, stats::qnorm(cip, lower.tail = F), -Inf),
      upper = ifelse(!!as.symbol(status_col) == 1, ifelse(use_fixed_case_thr, stats::qnorm(cip, lower.tail = F), Inf), stats::qnorm(cip, lower.tail = F))
    )
  

# returning formatted input -----------------------------------------------

  dplyr::select(family, !!as.symbol(fam_id_col), !!as.symbol(personal_id_col), !!as.symbol(role_col), lower, upper)
}

