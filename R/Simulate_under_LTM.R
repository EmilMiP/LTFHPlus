utils::globalVariables("fam_ID")
utils::globalVariables("PID")
utils::globalVariables("lower")
utils::globalVariables("upper")
utils::globalVariables("max_age")
utils::globalVariables("m_max_age")
utils::globalVariables("p_max_age")
utils::globalVariables("indiv_ID")
utils::globalVariables("tmp_names")

#' Simulate under the liability threshold model (single phenotype).
#'
#' \code{simulate_under_LTM_single} simulates families and thresholds under
#' the liability threshold model for a given family structure and a single 
#' phenotype. Please note that it is not possible to simulate different 
#' family structures. 
#'
#' @param fam_vec A vector of strings holding the different 
#' family members. All family members must be represented by strings from the 
#' following list:
#' - \code{m} (Mother)
#' - \code{f} (Father)
#' - \code{c[0-9]*.[0-9]\*} (Children)
#' - \code{mgm} (Maternal grandmother)
#' - \code{mgf} (Maternal grandfather)
#' - \code{pgm} (Paternal grandmother)
#' - \code{pgf} (Paternal grandfather)
#' - \code{s[0-9]*} (Full siblings)
#' - \code{mhs[0-9]*} (Half-siblings - maternal side)
#' - \code{phs[0-9]*} (Half-siblings - paternal side)
#' - \code{mau[0-9]*} (Aunts/Uncles - maternal side)
#' - \code{pau[0-9]*} (Aunts/Uncles - paternal side).
#'  Defaults to \code{c("m","f","s1","mgm","mgf","pgm","pgf")}.
#' @param n_fam A named vector holding the desired number of family members.
#' See  \code{\link[stats]{setNames}}.
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
#' is a number satisfying \eqn{0 \leq h^2}, \code{n_sim} is a strictly positive number,
#' and \code{pop_prev} is a positive number that is at most one, 
#' then the output will be a list holding two tibbles. 
#' The first tibble, \code{sim_obs}, holds the simulated liabilities, the disease
#' status and the current age/age-of-onset for all family members in each of the 
#' \code{n_sim} families. 
#' The second tibble, \code{thresholds}, holds the family identifier, the personal
#' identifier and the lower and upper thresholds for all individuals in all families. 
#' Note that this tibble has the format required in \code{\link{estimate_liability}}. 
#' In addition, note that if neither \code{fam_vec} nor \code{n_fam} are specified, the function 
#' returns the disease status, the current age/age-of-onset, the lower and upper 
#' thresholds, as well as the personal identifier for a single individual, namely 
#' the individual under consideration (called \code{o}).
#' If both \code{fam_vec} and \code{n_fam} are defined, the user is asked to '
#' decide on which of the two vectors to use.
#' 
#' @examples
#' simulate_under_LTM_single()
#' simulate_under_LTM_single(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2), 
#' c("m","mgm","mgf","mhs")))
#' simulate_under_LTM_single(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE, 
#' h2 = 0.5, n_sim = 500, pop_prev = .05)
#' simulate_under_LTM_single(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5, 
#' n_sim = 200, pop_prev = 0.05)
#' 
#' @seealso \code{\link{construct_covmat}}, \code{\link{simulate_under_LTM_multi}}, \code{\link{simulate_under_LTM}}
#' 
#' @importFrom dplyr %>% bind_cols select relocate mutate rowwise n across
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom tidyselect matches starts_with ends_with
#' @importFrom stringr str_detect
#' 
#' @export
simulate_under_LTM_single <- function(fam_vec = c("m","f","s1","mgm","mgf","pgm","pgf"), 
                                      n_fam = NULL, 
                                      add_ind = TRUE, 
                                      h2 = 0.5, 
                                      n_sim=1000, 
                                      pop_prev = .1){
  # Making sure that the input is valid ----------------------------------------------
  
  # If fam_vec or n_fam is a vector of length zero, it is set to
  # NULL instead
  if(length(fam_vec) == 0) fam_vec <- NULL
  if(length(n_fam) == 0) n_fam <- NULL
  
  # Turning add_ind into class logical
  add_ind <- as.logical(add_ind)
  
  # Checking that the heritability is valid
  if(validate_proportion(h2)){invisible()}
  
  # Checking that n_sim is a number
  if(!is.numeric(n_sim) && !is.integer(n_sim)) stop("The number of simulations n_sim must be numeric!")
  
  # Checking that n_sim is strictly positive
  if(n_sim <=0) stop("n_sim must be a positive number!")
  
  # Checking that pop_prev is valid
  if(validate_proportion(pop_prev)){invisible()}

  # Computing the covariance matrix.
  covmat <- construct_covmat_single(fam_vec = fam_vec, n_fam =n_fam, add_ind = add_ind, h2 = h2)
  
  # Simulating n_sim liabilities for the each family member.
  # The resulting tibble has n_sim rows and the same number 
  # of columns as covmat.
  liabs <- tmvtnorm::rtmvnorm(n = n_sim, mean = replicate(ncol(covmat), 0), sigma = covmat)

  # Adding the column names
  colnames(liabs) <- attributes(covmat)$fam_vec
  
  # Turning the matrix into a tibble and adding the family ID
  liabs <- tibble::as_tibble(liabs) %>%
    mutate(fam_ID = paste0("fam_ID_", 1:n())) %>%
    relocate(., fam_ID)
  
  # Adding the disease status for all individuals.
  # Remark: across() can be used to apply a function (.fns)
  # to a subset of columns (.cols) and storing the resulting
  # columns under pre-specified names (.names).
  # .cols uses the same syntax as select().
  liabs <- mutate(liabs, across(.cols = -c(matches("^g$"), matches("^fam_ID$")), 
                                .fns = ~ .x > qnorm(pop_prev, lower.tail = FALSE),
                                .names = "{.col}_status" ))
  
  # Adding the age for all individuals.
  # We begin by adding the age for all children,
  # if children are available.
  liabs <- mutate(liabs, across(.cols = matches(paste0("^c[0-9]*.[0-9]*$")), 
                                .fns = ~ sample(0:20, size = n(), replace = TRUE),
                                .names = "{.col}_age"))
  
  # Next, we add the age for the target individual as well as
  # its full- and half-siblings, if available.
  # If the target individual is supposed to have children,
  # we set its minimum age to 18 and the age depend on 
  # the children's age.
  min_age <- ifelse(any(str_detect(attributes(covmat)$fam_vec, "^c[0-9]*.[0-9]*$")), 18,10)
  
  if(any(str_detect(colnames(liabs), ".*_age$"))){
    
    liabs <- liabs %>% mutate(., max_age = purrr::invoke(pmax, select(., ends_with("_age"))))
  }else{
    
    liabs <- mutate(liabs, max_age = 0)
  }
  
  liabs <- mutate(liabs, across(.cols = matches("^o$"), 
                                .fns = ~ sample(min_age:40, size = n(), replace = TRUE) + max_age,
                                .names = "{.col}_age"))

  liabs <- mutate(liabs, across(.cols = c(matches("^s[0-9]*$"), matches("^[mp]hs[0-9]*$")), 
                                .fns = ~ sample(min_age:40, size = n(), replace = TRUE),
                                .names = "{.col}_age"))
  
  # In order for the parents to have a reasonable age,
  # their age depends on the age of their oldest child,
  # if children are available.
  if(any(str_detect(colnames(liabs), ".*_age$"))){
    
    liabs <- liabs %>% mutate(., max_age = purrr::invoke(pmax, select(., ends_with("_age"))))
  }else{
    
    liabs <- mutate(liabs, max_age = 0)
  }
  
  liabs <- liabs %>%
    mutate(., across(.cols = c(matches("^[mf]$")), 
                     .fns = ~sample(18:30, size = n(), replace = TRUE) + max_age,
                     .names = "{.col}_age")) %>%
    mutate(., across(.cols = c(matches("^[mp]au[0-9]*$")), 
                     .fns = ~sample(12:40, size = n(), replace = TRUE) + max_age,
                     .names = "{.col}_age"))
  
  # In order for the grandparents to have a reasonable age,
  # their age depends on the age of their oldest child,
  # if children are available.
  if(any(str_detect(colnames(liabs), "^m_age$") | str_detect(colnames(liabs), "^mau[0-9]*_age$"))){
    
    liabs <- liabs %>% mutate(., m_max_age = purrr::invoke(pmax, select(., matches("^m_age$"), matches("^mau[0-9]*_age$"))))
  }else{
    
    liabs <- mutate(liabs, m_max_age = 25)
  }
  
  if(any(str_detect(colnames(liabs), "^f_age$") | str_detect(colnames(liabs), "^pau[0-9]*_age$"))){
    
    liabs <- liabs %>% mutate(., p_max_age = purrr::invoke(pmax, select(., matches("^f_age$"), matches("^pau[0-9]*_age$"))))
  }else{
    
    liabs <- mutate(liabs, p_max_age = 25)
  }
   
  liabs <- liabs %>%
    mutate(., across(.cols = matches("^mg[mf]$"), 
                     .fns = ~sample(15:30, size = n(), replace = TRUE) + m_max_age,
                     .names = "{.col}_age")) %>%
    mutate(., across(.cols = matches("^pg[mf]$"), 
                            .fns = ~sample(15:30, size = n(), replace = TRUE) + p_max_age,
                            .names = "{.col}_age")) %>%
    select(., -c(max_age, m_max_age, p_max_age))
  
  # Finally, we can add the age of onset for all individuals 
  # having the disease
  liabs <- liabs %>% mutate(., construct_aoo(fam_mem = attributes(covmat)$fam_vec, .tbl = ., pop_prev = pop_prev))
  
  # Constructing thresholds
  threshs <- construct_thresholds(fam_mem = attributes(covmat)$fam_vec, .tbl = liabs, pop_prev = pop_prev)
  
  # Returning the simulated thresholds
  return(list(sim_obs = select(liabs, -c(ends_with("_age"))), 
              thresholds = threshs))
}


#' Simulate under the liability threshold model (multiple phenotypes).
#'
#' \code{simulate_under_LTM_multi} simulates families and thresholds under
#' the liability threshold model for a given family structure and multiple
#' phenotypes. Please note that it is not possible to simulate different 
#' family structures. 
#'
#' @param fam_vec A vector of strings holding the different 
#' family members. All family members must be represented by strings from the 
#' following list:
#' - \code{m} (Mother)
#' - \code{f} (Father)
#' - \code{c[0-9]*.[0-9]\*} (Children)
#' - \code{mgm} (Maternal grandmother)
#' - \code{mgf} (Maternal grandfather)
#' - \code{pgm} (Paternal grandmother)
#' - \code{pgf} (Paternal grandfather)
#' - \code{s[0-9]*} (Full siblings)
#' - \code{mhs[0-9]*} (Half-siblings - maternal side)
#' - \code{phs[0-9]*} (Half-siblings - paternal side)
#' - \code{mau[0-9]*} (Aunts/Uncles - maternal side)
#' - \code{pau[0-9]*} (Aunts/Uncles - paternal side).
#'  Defaults to \code{c("m","f","s1","mgm","mgf","pgm","pgf")}.
#' @param n_fam A named vector holding the desired number of family members.
#' See  \code{\link[stats]{setNames}}.
#' All names must be picked from the list mentioned above. Defaults to \code{NULL}.
#' @param add_ind A logical scalar indicating whether the genetic 
#' component of the full liability as well as the full
#' liability for the underlying target individual should be included in 
#' the covariance matrix. Defaults to \code{TRUE}.
#' @param genetic_corrmat A numeric matrix holding the genetic correlations 
#' between the desired phenotypes. All diagonal entries must be equal to one, 
#' while all off-diagonal entries must be between -1 and 1. In addition, 
#' the matrix must be symmetric.
#' Defaults to \code{diag(3)}.
#' @param full_corrmat A numeric matrix holding the full correlations 
#' between the desired phenotypes. All diagonal entries must be equal to 
#' one, while all off-diagonal entries must be between -1 and 1. In addition, 
#' the matrix must be symmetric.
#' Defaults to \code{diag(3)}.
#' @param h2_vec A numeric vector holding the liability-scale heritabilities 
#' for a number of phenotype. All entries must be non-negative. Note that under 
#' the liability threshold model, the heritabilities must also be at most 1.
#' Defaults to \code{rep(0.5,3)}.
#' @param phen_names A character vector holding the phenotype names. These names
#' will be used to create the row and column names for the covariance matrix.
#' If it is not specified, the names will default to phenotype1, phenotype2, etc.
#' Defaults to \code{NULL}.
#' @param n_sim A positive number representing the number of simulations. Defaults to 1000.
#' @param pop_prev A numeric vector holding the population prevalences, i.e. the 
#' overall prevalences in the population. All entries in \code{pop_prev} must be positive
#' and smaller than 1. Defaults to \code{rep(.1,3)}.
#' 
#' @return If either \code{fam_vec} or \code{n_fam} is used as the argument and if it is of the 
#' required format, if \code{genetic_corrmat} and \code{full_corrmat} are two numeric 
#' and symmetric matrices satisfying that all diagonal entries are one and that all 
#' off-diagonal entries are between -1 and 1, if the liability-scale heritabilities in 
#' \code{h2_vec} are numbers satisfying \eqn{0 \leq h^2_i} for all \eqn{i \in \{1,...,n_pheno\}}, 
#' \code{n_sim} is a strictly positive number, and \code{pop_prev} is a positive numeric 
#' vector such that all entries are at most one,
#' then the output will be a list containing lists for each phenotype.
#' The first outer list, which is named after the first phenotype in \code{phen_names}, 
#' holds two tibbles, namely \code{sim_obs}, which holds the simulated liabilities, the
#' disease status and the current age/age-of-onset for all family members in each of the \code{n_sim} 
#' families for the first phenotype. The second tibble, \code{thresholds}, holds 
#' the family identifier, the personal identifier and the lower and upper thresholds 
#' for all individuals in all families, again for the first phenotype. Note that this 
#' tibble has the format required in \code{\link{estimate_liability}}. 
#' As the first outer list, the second outer list, which is named after the second 
#' phenotype in \code{phen_names}, holds two tibbles. \code{sim_obs}, which holds 
#' the  simulated liabilities, the disease status and the current age/age-of-onset 
#' for all family members in each of the \code{n_sim} families for the second phenotype, 
#' and \code{thresholds}, which holds the family identifier, the personal
#' identifier and the lower and upper thresholds for all individuals in all families
#' for the second phenotype.
#' Once more, note that this tibble has the format required in 
#' \code{\link{estimate_liability}}.
#' There is a list containing \code{sim_obs} and \code{thresholds} for each 
#' phenotype in \code{phen_names}. 
#' Finally, note that if neither \code{fam_vec} nor \code{n_fam} are specified, the function 
#' returns the disease status, the current age/age-of-onset, the lower and upper 
#' thresholds, as well as the personal identifier for a single individual, namely 
#' the individual under consideration (called \code{o}).
#' If both \code{fam_vec} and \code{n_fam} are defined, the user is asked to '
#' decide on which of the two vectors to use.
#' 
#' @examples
#' simulate_under_LTM_multi()
#' 
#' genetic_corrmat <- matrix(0.4, 3, 3)
#' diag(genetic_corrmat) <- 1
#' full_corrmat <- matrix(0.6, 3, 3)
#' diag(full_corrmat) <- 1
#' 
#' simulate_under_LTM_multi(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2,2), 
#' c("m","mgm","mgf","s","mhs")))
#' 
#' simulate_under_LTM_multi(fam_vec = c("m","f","s1"), add_ind = FALSE, 
#' genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, n_sim = 100)
#' 
#' simulate_under_LTM_multi(fam_vec = c(), n_fam = NULL, add_ind = TRUE, n_sim = 150)
#' 
#' @seealso \code{\link{construct_covmat}}
#' 
#' @importFrom dplyr %>% bind_cols select relocate mutate rowwise n across
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom tidyselect matches contains ends_with
#' @importFrom stringr str_detect
#' 
#' @export
simulate_under_LTM_multi <- function(fam_vec = c("m","f","s1","mgm","mgf","pgm","pgf"), 
                                      n_fam = NULL, 
                                      add_ind = TRUE,
                                      genetic_corrmat = diag(3),
                                      full_corrmat = diag(3),
                                      h2_vec = rep(.5,3), 
                                      phen_names = NULL,
                                      n_sim = 1000, 
                                      pop_prev = rep(.1,3)){
  
  # If fam_vec or n_fam is a vector of length zero, it is set to
  # NULL instead
  if(length(fam_vec) == 0) fam_vec <- NULL
  if(length(n_fam) == 0) n_fam <- NULL
  
  # The same holds for the vector holding phenotype names
  if(length(phen_names) == 0) phen_names <- NULL
  
  # Turning add_ind into class logical
  add_ind <- as.logical(add_ind)
  
  # Checking that the heritability is valid
  if(validate_proportion(h2_vec)){invisible()}
  
  # Checking that all correlations are valid
  if(validate_correlation_matrix(genetic_corrmat)){invisible()}
  if(validate_correlation_matrix(full_corrmat)){invisible()}
  
  # Computing the number of phenotypes
  num_phen <- length(h2_vec)
  
  # Checking that phen_names is either NULL or a valid
  # vector of strings
  if(is.null(phen_names)){
    phen_names <- paste0("phenotype", 1:num_phen)
  }else{
    if(!is.character(phen_names)) phen_names <- as.character(phen_names)
    if(length(phen_names) != num_phen) stop("The number of names in phen_num and the number of phenotypes differ...")
  }
  
  # Checking that n_sim is a number
  if(!is.numeric(n_sim) && !is.integer(n_sim)) stop("The number of simulations n_sim must be numeric!")
  
  # Checking that n_sim is strictly positive
  if(n_sim <=0)stop("n_sim must be a positive number!")
  
  # Checking that pop_prev is valid
  if(validate_proportion(pop_prev)){invisible()}
  
  # Computing the covariance matrix.
  covmat <- suppressWarnings(construct_covmat_multi(fam_vec = fam_vec, n_fam = n_fam, add_ind = add_ind, 
                                                    genetic_corrmat = genetic_corrmat, 
                                                    full_corrmat = full_corrmat, h2_vec = h2_vec, 
                                                    phen_names = phen_names))
  
  # Simulating n_sim liabilities for the each family member and each
  # phenotype. The resulting tibble has n_sim rows and the same number
  # of columns as covmat.
  liabs <- tmvtnorm::rtmvnorm(n = n_sim, mean = replicate(ncol(covmat), 0), sigma = covmat)
  # Adding the column names
  colnames(liabs) <- colnames(covmat)
  
  # Turning the matrix into a tibble and adding the individual ID
  liabs <- tibble::as_tibble(liabs) %>%
    mutate(fam_ID = paste0("fam_ID_", 1:n())) %>%
    relocate(., fam_ID)
  
  # Adding the disease status for all individuals.
  # RemarK: across() can be used to apply a function (.fns)
  # to a subset of columns (.cols) and storing the resulting
  # columns under pre-specified names (.names).
  # .cols uses the same syntax as select().
  for(i in 1:num_phen){
    
    liabs <- mutate(liabs, across(.cols = c(contains(phen_names[i]), -matches("^g_.*")), 
                                  .fns = ~ .x > qnorm(pop_prev[i], lower.tail = FALSE),
                                  .names = "{.col}_status" ))
  }
  
  # Adding the age for all individuals.
  # We begin by adding the age for all children,
  # if children are available.
  liabs <- mutate(liabs, across(.cols = matches(paste0("^c[0-9]*.[0-9]*_", phen_names[1],"$")), 
                                .fns = ~ sample(0:20, size = n(), replace = TRUE),
                                .names = "{gsub(phen_names[1],'', {col}, fixed = TRUE)}age"))
  
  
  # Next, we add the age for the target individual as well as
  # its full- and half-siblings, if available.
  # If the target individual is supposed to have children,
  # we set its minimum age to 18 and the age depend on 
  # the children's age.
  min_age <- ifelse(any(str_detect(attributes(covmat)$fam_vec, "^c[0-9]*.[0-9]*$")), 18,10)
  
  if(any(str_detect(colnames(liabs), ".*_age$"))){
    
    liabs <- liabs %>% 
      mutate(., max_age = purrr::invoke(pmax, select(., ends_with("_age"))))
  }else{
    
    liabs <- mutate(liabs, max_age = 0)
  }
  
  liabs <- mutate(liabs, across(.cols = matches(paste0("^o_", phen_names[1],"$")), 
                                .fns = ~ sample(min_age:40, size = n(), replace = TRUE) + max_age,
                                .names = "{gsub(phen_names[1],'', {col}, fixed = TRUE)}age"))
  
  liabs <- mutate(liabs, across(.cols = c(matches(paste0("^s[0-9]*_", phen_names[1],"$")), matches(paste0("^[mp]hs[0-9]*_", phen_names[1],"$"))), 
                                .fns = ~ sample(min_age:40, size = n(), replace = TRUE),
                                .names = "{gsub(phen_names[1],'', {col}, fixed = TRUE)}age"))
  
  
  # In order for the parents to have a reasonable age,
  # their age depends on the age of their oldest child.
  if(any(str_detect(colnames(liabs), ".*_age$"))){
    
    liabs <- liabs %>% 
      mutate(., max_age = purrr::invoke(pmax, select(., ends_with("_age"))))
  }else{
    
    liabs <- mutate(liabs, max_age = 0)
  }
  
  liabs <- liabs %>%
    mutate(., across(.cols = matches(paste0("^[mf]_", phen_names[1],"$")), 
                     .fns = ~sample(18:30, size = n(), replace = TRUE) + max_age,
                     .names = "{gsub(phen_names[1],'', {col}, fixed = TRUE)}age")) %>%
    mutate(., across(.cols = matches(paste0("^[mp]au[0-9]*_", phen_names[1],"$")), 
                     .fns = ~sample(12:40, size = n(), replace = TRUE) + max_age,
                     .names = "{gsub(phen_names[1],'', {col}, fixed = TRUE)}age"))
  
  # In order for the grandparents to have a reasonable age,
  # their age depends on the age of their oldest child,
  # if children are available.
  if(any(str_detect(colnames(liabs), "^m_age$") | str_detect(colnames(liabs), "^mau[0-9]*_age$"))){
    
    liabs <- liabs %>% 
      mutate(., m_max_age = purrr::invoke(pmax, select(., matches("^m_age$"), matches("^mau[0-9]*_age$"))))
  }else{
    
    liabs <- mutate(liabs, m_max_age = 25)
  }
  
  if(any(str_detect(colnames(liabs), "^f_age$") | str_detect(colnames(liabs), "^pau[0-9]*_age$"))){
    
    liabs <- liabs %>% 
      mutate(., p_max_age = purrr::invoke(pmax, select(., matches("^f_age$"), matches("^pau[0-9]*_age$"))))
  }else{
    
    liabs <- mutate(liabs, p_max_age = 25)
  }
  
  liabs <- liabs %>%
    dplyr::mutate(., across(.cols = matches(paste0("^mg[mf]_", phen_names[1],"$")), 
                            .fns = ~sample(15:30, size = n(), replace = TRUE) + m_max_age,
                            .names = "{gsub(phen_names[1],'', {col}, fixed = TRUE)}age")) %>%
    dplyr::mutate(., across(.cols = matches(paste0("^pg[mf]_", phen_names[1],"$")), 
                            .fns = ~sample(15:30, size = n(), replace = TRUE) + p_max_age,
                            .names = "{gsub(phen_names[1],'', {col}, fixed = TRUE)}age")) %>%
    dplyr::select(., -c(max_age, m_max_age, p_max_age))
  
  
  # Finally, we can add the age of onset for all individuals 
  # having a disease and construct the thresholds for
  # all family members in each family. 
  fam_vec = unique(gsub("_.*","",colnames(covmat)))
    
  res <- lapply(seq_along(phen_names), function(i){
    
    i_liabs <- select(liabs, fam_ID, contains(phen_names[i]), ends_with("_age")) %>% 
      mutate(., construct_aoo(fam_mem = fam_vec, .tbl = ., pop_prev = pop_prev[i], phen_name = phen_names[i]))
  
    i_threshs <- construct_thresholds(fam_mem = fam_vec, .tbl = i_liabs, pop_prev = pop_prev[i], phen_name = phen_names[i])
    
    return(list(sim_obs = select(i_liabs, -c(ends_with("_age"))),
                thresholds = i_threshs))
  })
  
  # Renaming the list components
  names(res) <- phen_names
  # Returning the simulated thresholds
  return(res)
}

#' Simulate under the liability threshold model.
#'
#' \code{simulate_under_LTM} simulates families and thresholds under
#' the liability threshold model for a given family structure and a 
#' variable number of phenotypes.Please note that it is not possible 
#' to simulate different family structures. 
#' 
#' This function can be used to simulate the case-control status, the current 
#' age and age-of-onset as well as the lower and upper thresholds for
#' a variable number of phenotypes for all family members in each of 
#' the \code{n_sim} families. 
#' If \code{h2} is a number, \code{simulate_under_LTM} simulates the case-
#' control status, the current age and age-of-onset as well as thresholds
#' for a single phenotype.
#' However, if \code{h2} is a numeric vector, if \code{genetic_corrmat} and 
#' \code{full_corrmat} are two symmetric correlation matrices, and if 
#' \code{phen_names} and \code{pop_prev} are to numeric vectors holding 
#' the phenotype names and the population prevalences, respectively, then 
#' \code{simulate_under_LTM} simulates the case-control status, the current 
#' age and age-of-onset as well as thresholds for two or more (correlated) 
#' phenotypes.
#' The family members can be specified using one of two possible formats.
#'
#' @param fam_vec A vector of strings holding the different 
#' family members. All family members must be represented by strings from the 
#' following list:
#' - \code{m} (Mother)
#' - \code{f} (Father)
#' - \code{c[0-9]*.[0-9]\*} (Children)
#' - \code{mgm} (Maternal grandmother)
#' - \code{mgf} (Maternal grandfather)
#' - \code{pgm} (Paternal grandmother)
#' - \code{pgf} (Paternal grandfather)
#' - \code{s[0-9]*} (Full siblings)
#' - \code{mhs[0-9]*} (Half-siblings - maternal side)
#' - \code{phs[0-9]*} (Half-siblings - paternal side)
#' - \code{mau[0-9]*} (Aunts/Uncles - maternal side)
#' - \code{pau[0-9]*} (Aunts/Uncles - paternal side).
#'  Defaults to \code{c("m","f","s1","mgm","mgf","pgm","pgf")}.
#' @param n_fam A named vector holding the desired number of family members.
#' See  \code{\link[stats]{setNames}}.
#' All names must be picked from the list mentioned above. Defaults to \code{NULL}.
#' @param add_ind A logical scalar indicating whether the genetic 
#' component of the full liability as well as the full
#' liability for the underlying target individual should be included in 
#' the covariance matrix. Defaults to \code{TRUE}.
#' @param h2 Either a number or a numeric vector holding the liability-scale 
#' heritability(ies) for one or more phenotypes. All entries in \code{h2} must 
#' be non-negative. Note that under the liability threshold model, the 
#' heritabilities must also be at most 1. Defaults to 0.5.
#' @param genetic_corrmat Either \code{NULL} or a numeric matrix holding the 
#' genetic correlations between the desired phenotypes. Must be specified, if
#' \code{length(h2)}\eqn{>0}, and will be ignored if \code{h2} is a number.
#' All diagonal entries in \code{genetic_corrmat} must be equal to one, 
#' while all off-diagonal entries must be between -1 and 1. In addition, 
#' the matrix must be symmetric.
#' Defaults to \code{NULL}.
#' @param full_corrmat Either \code{NULL} or a numeric matrix holding the 
#' full correlations between the desired phenotypes. Must be specified, if
#' \code{length(h2)}\eqn{>0}, and will be ignored if \code{h2} is a number.
#' All diagonal entries in \code{full_corrmat} must be equal to one, while 
#' all off-diagonal entries must be between -1 and 1. In addition, the 
#' matrix must be symmetric.
#' Defaults to \code{NULL}.
#' @param phen_names Either \code{NULL} or character vector holding the 
#' phenotype names. These names will be used to create the row and column 
#' names for the covariance matrix. Must be specified, if \code{length(h2)}
#' \eqn{> 0}, and will be ignored if \code{h2} is a number.
#' If it is not specified, the names will default to phenotype1, phenotype2, etc.
#' Defaults to \code{NULL}.
#' @param n_sim A positive number representing the number of simulations. Defaults to 1000.
#' @param pop_prev Either a number or a numeric vector holding the population 
#' prevalence(s), i.e. the overall prevalence(s) in the population. 
#' All entries in \code{pop_prev} must be positive
#' and smaller than 1. Defaults to 0.1.
#' 
#' @return If either \code{fam_vec} or \code{n_fam} is used as the argument, 
#' if it is of the required format, if the liability-scale heritability \code{h2} 
#' is a number satisfying \eqn{0 \leq h^2}, \code{n_sim} is a strictly positive number,
#' and \code{pop_prev} is a positive number that is at most one, 
#' then the output will be a list containing two tibbles. 
#' The first tibble, \code{sim_obs}, holds the simulated liabilities, the disease
#' status and the current age/age-of-onset for all family members in each of the 
#' \code{n_sim} families. 
#' The second tibble, \code{thresholds}, holds the family identifier, the personal
#' identifier and the lower and upper thresholds for all individuals in all families. 
#' Note that this tibble has the format required in \code{\link{estimate_liability}}.
#' If either \code{fam_vec} or \code{n_fam} is used as the argument and if it is of the 
#' required format, if \code{genetic_corrmat} and \code{full_corrmat} are two numeric 
#' and symmetric matrices satisfying that all diagonal entries are one and that all 
#' off-diagonal entries are between -1 and 1, if the liability-scale heritabilities in 
#' \code{h2_vec} are numbers satisfying \eqn{0 \leq h^2_i} for all \eqn{i \in \{1,...,n_pheno\}}, 
#' \code{n_sim} is a strictly positive number, and \code{pop_prev} is a positive numeric 
#' vector such that all entries are at most one, then the output will be a list containing 
#' the following lists.
#' The first outer list, which is named after the first phenotype in \code{phen_names}, 
#' holds two tibbles, namely \code{sim_obs}, which holds the simulated liabilities, the
#' disease status and the current age/age-of-onset for all family members in each of the \code{n_sim} 
#' families for the first phenotype. The second tibble, \code{thresholds}, holds 
#' the family identifier, the personal identifier and the lower and upper thresholds 
#' for all individuals in all families, again for the first phenotype. Note that this 
#' tibble has the format required in \code{\link{estimate_liability}}. 
#' As the first outer list, the second outer list, which is named after the second 
#' phenotype in \code{phen_names}, holds two tibbles. \code{sim_obs}, which holds 
#' the  simulated liabilities, the disease status and the current age/age-of-onset 
#' for all family members in each of the \code{n_sim} families for the second phenotype, 
#' and \code{thresholds}, which holds the family identifier, the personal
#' identifier and the lower and upper thresholds for all individuals in all families
#' for the second phenotype.
#' Once more, note that this tibble has the format required in 
#' \code{\link{estimate_liability}}.
#' There is a list containing \code{sim_obs} and \code{thresholds} for each 
#' phenotype in \code{phen_names}. 
#' Finally, note that if neither \code{fam_vec} nor \code{n_fam} are specified, the function 
#' returns the disease status, the current age/age-of-onset, the lower and upper 
#' thresholds, as well as the personal identifier for a single individual, namely 
#' the individual under consideration (called \code{o}).
#' If both \code{fam_vec} and \code{n_fam} are defined, the user is asked to '
#' decide on which of the two vectors to use.
#' 
#' @examples
#' simulate_under_LTM()
#' 
#' genetic_corrmat <- matrix(0.4, 3, 3)
#' diag(genetic_corrmat) <- 1
#' full_corrmat <- matrix(0.6, 3, 3)
#' diag(full_corrmat) <- 1
#' 
#' simulate_under_LTM(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2,2), 
#' c("m","mgm","mgf","s","mhs")))
#' 
#' simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE, 
#' genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, n_sim = 200)
#' 
#' simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5, 
#' n_sim = 200, pop_prev = 0.05)
#' 
#' @seealso \code{\link{construct_covmat}} \code{\link{simulate_under_LTM_single}}
#' \code{\link{simulate_under_LTM_multi}}
#' 
#' @export
simulate_under_LTM <- function(fam_vec = c("m","f","s1","mgm","mgf","pgm","pgf"), 
                               n_fam = NULL, 
                               add_ind = TRUE,
                               h2 = 0.5, 
                               genetic_corrmat = NULL,
                               full_corrmat = NULL,
                               phen_names = NULL,
                               n_sim = 1000, 
                               pop_prev = 0.1){
  
  if(length(h2) == 1){
    
    return(simulate_under_LTM_single(fam_vec = fam_vec, 
                                      n_fam = n_fam, 
                                      add_ind = add_ind, 
                                      h2 = h2, 
                                      n_sim = n_sim, 
                                      pop_prev = pop_prev))
    
  }else{
    
    return(simulate_under_LTM_multi(fam_vec = fam_vec, 
                                    n_fam = n_fam, 
                                    add_ind = add_ind,
                                    genetic_corrmat = genetic_corrmat,
                                    full_corrmat = full_corrmat,
                                    h2_vec = h2, 
                                    phen_names = phen_names,
                                    n_sim = n_sim, 
                                    pop_prev = pop_prev))
  }
}
