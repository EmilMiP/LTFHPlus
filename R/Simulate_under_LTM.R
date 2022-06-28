utils::globalVariables("fam_ID")
utils::globalVariables("PID")
utils::globalVariables("lower")
utils::globalVariables("upper")
utils::globalVariables("max_age")
utils::globalVariables("m_max_age")
utils::globalVariables("p_max_age")
utils::globalVariables("indiv_ID")
utils::globalVariables("tmp_names")

#' Simulate under the liability threshold model.
#'
#' \code{simulate_under_LTM} simulates families and thresholds under
#' the liability threshold model for a given family
#' structure.
#'
#' @param fam_vec A vector of strings holding the different 
#' family members. All family members must be represented by strings from the 
#' following list:
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
#'  Defaults to c("m","f","s1","mgm","mgf","pgm","pgf").
#' @param n_fam A named vector holding the desired number of family members.
#' All names must be picked from the list mentioned above. Defaults to \code{NULL}.
#' @param add_ind A logical scalar indicating whether the genetic 
#' component of the full liability as well as the full
#' liability for the underlying target individual should be included in 
#' the covariance matrix. Defaults to \code{TRUE}.
#' @param h2 A number representing the liability-scale heritability 
#' for a single phenotype. Must be non-negative. Note that under 
#' the liability threshold model, the heritability must also be at most 1.
#' Defaults to 0.5.
#' @param n_sim A positive number representing the number of simulations. Defaults to 1000.
#' @param pop_prev A positive number representing the population prevalence, i.e. the 
#' overall prevalence in the population. Must be smaller than 1. Defaults to 0.1.
#' 
#' @return If either \code{fam_vec} or \code{n_fam} is used as the argument, 
#' if it is of the required format, if the liability-scale heritability \code{h2} 
#' is a number satisfying \eqn{0 \leq h^2}, \CODE{n_sim} is a strictly positive number,
#' and \code{pop_prev} is a positive number that is at most one, 
#' then the output will be a list containing three tibbles. 
#' The first tibble, \code{sim_obs}, holds the disease status and the current 
#' age/age-of-onset for all family members in each of the \code{n_sim} families. 
#' The second tibble, \code{thresholds}, holds the lower and upper thresholds 
#' for all individuals in all families. Note that this tibble has the format required in 
#' \code{\link{estimate_liability}}. 
#' The last tibble, \code{fam_ID}, connects the personal identifiers
#' for all individuals to the family identifiers. As \code{thresholds}, \code{fam_ID} 
#' has the format required in \code{\link{estimate_liability}}.
#' Note that if neither \code{fam_vec} nor \code{n_fam} are specified, the function 
#' returns the disease status, the current age/age-of-onset, the lower and upper 
#' thresholds, as well as the personal identifier for a single individual, namely 
#' the individual under consideration (called \code{o}).
#' If both \code{fam_vec} and \code{n_fam} are defined, the user is asked to '
#' decide on which of the two vectors to use.
#' 
#' @examples
#' simulate_under_LTM()
#' simulate_under_LTM(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2,2), 
#' c("m","mgm","mgf","s","mhs")))
#' simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE, 
#' h2 = 0.5, n_sim = 500, pop_prev = .05)
#' simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5, 
#' n_sim = 200, pop_prev = 0.05)
#' 
#' @seealso \code{\link{construct_covmat}}
#' 
#' @importFrom dplyr %>% bind_cols select relocate mutate rowwise n across
#' @importFrom tmvtnorm rtmvnorm
#' 
#' @export
simulate_under_LTM <- function(fam_vec = c("m","f","s1","mgm","mgf","pgm","pgf"), 
                               n_fam = NULL, 
                               add_ind = TRUE, 
                               h2 = 0.5, 
                               n_sim=1000, 
                               pop_prev = .1){
  
  # Turning add_ind into class logical
  add_ind <- as.logical(add_ind)
  
  # Checking that the heritability is valid
  if(validate_proportion(h2)){invisible()}
  
  # Checking that n_sim is a number
  if(!is.numeric(n_sim) && !is.integer(n_sim)) stop("The number of simulations n_sim must be numeric!")
  
  # Checking that n_sim is strictly positive
  if(n_sim <=0)stop("n_sim must be a positive number!")
  
  # Checking that pop_prev is valid
  if(validate_proportion(pop_prev)){invisible()}

  # Computing the covariance matrix.
  # If both fam_vec and n_fam are empty, construct_covmat
  # would usually return a warning. We supress this warning here.
  if(is.null(fam_vec) && is.null(n_fam)){
    
    covmat <- suppressWarnings(construct_covmat_single(fam_vec = NULL, n_fam = NULL, add_ind = add_ind, h2 = h2))
  }else{
    covmat <- construct_covmat_single(fam_vec = fam_vec, n_fam = n_fam, add_ind = add_ind, h2 = h2)
  }
  
  # Simulating n_sim liabilities for the each family member
  liabs <- tmvtnorm::rtmvnorm(n = n_sim, mean = replicate(ncol(covmat), 0), sigma = covmat)
  # Adding the column names
  colnames(liabs) <- colnames(covmat)
  
  # Turning the matrix into a tibble and adding the individual ID
  liabs <- tibble::as_tibble(liabs) %>%
    mutate(indiv_ID = paste0("ID", 1:n())) %>%
    relocate(., indiv_ID)
  
  # Adding the disease status for all individuals.
  # RemarK: across() can be used to apply a function (.fns)
  # to a subset of columns (.cols) and storing the resulting
  # columns under pre-specified names (.names).
  # .cols uses the same syntax as select().
  liabs <- mutate(liabs, across(.cols = -c(tidyselect::matches("^g$"), tidyselect::matches("^indiv_ID$")), 
                                .fns = ~ .x > qnorm(pop_prev, lower.tail = FALSE),
                                .names = "{.col}_status" ))
  
  # Adding the age for all individuals.
  # We begin by adding the age for the individual as well as
  # its full- and half-siblings, if available.
  liabs <- mutate(liabs, across(.cols = c(tidyselect::matches("^o$"), tidyselect::matches("^s[0-9]*$"), tidyselect::matches("^[mp]hs[0-9]*$")), 
                                .fns = ~sample(10:40, size = n(), replace = TRUE),
                                .names = "{.col}_age"))
  
  liabs <- mutate(liabs, across(.cols = c(tidyselect::matches("^c[0-9]*$")), 
                                .fns = ~ pmax(o_age - sample(10:20, size = n(), replace = TRUE), 0), # pmax to avoid negative ages.
                                .names = "{.col}_age"))
  
  
  # In order for the parents to have a reasonable age,
  # their age depends on the age of their oldest child,
  # if children are available.
  if(any(stringr::str_detect(colnames(liabs), ".*_age$"))){
    
    liabs <- liabs %>% 
      mutate(., max_age = purrr::invoke(pmax, dplyr::select(., tidyselect::ends_with("_age"))))
  }else{
    
    liabs <- mutate(liabs, max_age = 0)
  }
  
  liabs <- liabs %>%
    mutate(., across(.cols = c(tidyselect::matches("^[mf]$")), 
                     .fns = ~sample(18:30, size = n(), replace = TRUE) + max_age,
                     .names = "{.col}_age")) %>%
    mutate(., across(.cols = c(tidyselect::matches("^[mp]au[0-9]*$")), 
                     .fns = ~sample(12:40, size = n(), replace = TRUE) + max_age,
                     .names = "{.col}_age"))
  
  # In order for the grandparents to have a reasonable age,
  # their age depends on the age of their oldest child,
  # if children are available.
  if(any(stringr::str_detect(colnames(liabs), "^m_age$") | stringr::str_detect(colnames(liabs), "^mau[0-9]*_age$"))){
    
    liabs <- liabs %>% 
      mutate(., m_max_age = purrr::invoke(pmax, dplyr::select(., tidyselect::matches("^m_age$"), tidyselect::matches("^mau[0-9]*_age$"))))
  }else{
    
    liabs <- mutate(liabs, m_max_age = 25)
  }
  
  if(any(stringr::str_detect(colnames(liabs), "^p_age$") | stringr::str_detect(colnames(liabs), "^pau[0-9]*_age$"))){
    
    liabs <- liabs %>% 
      mutate(., p_max_age = purrr::invoke(pmax, dplyr::select(., tidyselect::matches("^f_age$"), tidyselect::matches("^pau[0-9]*_age$"))))
  }else{
    
    liabs <- mutate(liabs, p_max_age = 25)
  }
   
  liabs <- liabs %>%
    dplyr::mutate(., across(.cols = c(tidyselect::matches("^[pm]g[mf]$")), 
                     .fns = ~sample(15:30, size = n(), replace = TRUE) + m_max_age,
                     .names = "{.col}_age")) %>%
    dplyr::mutate(., across(.cols = c(tidyselect::matches("^[pm]gp[0-9]*")), 
                            .fns = ~sample(15:30, size = n(), replace = TRUE) + p_max_age,
                            .names = "{.col}_age")) %>%
    dplyr::select(., -c(max_age, m_max_age, p_max_age))
  
  # Adding age of onset for all individuals having the disease
  liabs <- liabs %>% mutate(., construct_aoo(fam_mem = colnames(covmat), .tbl = ., pop_prev = pop_prev))
  
  # Constructing thresholds
  threshs <- construct_thresholds(fam_mem = colnames(covmat), .tbl = liabs, pop_prev = pop_prev)
  
  # Constructing the personal identifiers for all
  # family members in each family.
  fam_ID <- dplyr::select(liabs, indiv_ID) %>% 
    mutate(., fam_ID = lapply(indiv_ID, function(i){
      tmp_names = setdiff(colnames(covmat),c("g"))
      paste0(i,"_member", seq_along(tmp_names),"_", tmp_names)
    }))
  
  return(list(sim_obs = liabs, 
              thresholds = threshs, 
              fam_ID = fam_ID))
}
