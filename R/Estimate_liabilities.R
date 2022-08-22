utils::globalVariables("lower")
utils::globalVariables("upper")
utils::globalVariables("pb")
utils::globalVariables("phenotype")
utils::globalVariables("est")
utils::globalVariables("se")
utils::globalVariables("name")

#' Estimating the genetic or full liability 
#'
#' \code{estimate_liability_single} estimates the genetic component of the full
#' liability and/or the full liability for a number of individuals based
#' on their family history.
#'
#' This function can be used to estimate either the genetic component of the 
#' full liability, the full liability or both. It is possible to input either 
#' 
#' @param .tbl A matrix, list or data frame that can be converted into a tibble.
#' Must have at least five columns that hold the family identifier, the personal 
#' identifier, the role and the lower and upper thresholds. Note that the 
#' role must be one of the following abbreviations 
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
#' Defaults to \code{NULL}.
#' @param family A matrix, list or data frame that can be converted into a tibble.
#' Must have at least two columns that hold the family identifier and the corresponding
#' personal identifiers, respectively. That is, for each family in \code{fam_id} there should
#' be a list holding all individuals belonging to that family in \code{pid}. Note that the 
#' personal identifiers for all individuals must have a special format. It must end
#' on \code{_?}, where \code{?} is one of the following abbreviations 
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
#' See also \code{\link{construct_covmat}}.
#' Defaults to \code{NULL}.
#' @param threshs A matrix, list or data frame that can be converted into a tibble.
#' Must have at least three columns, one holding the personal identifier for all individuals,
#' and the remaining two holding the lower and upper thresholds, respectively. The latter
#' should be called "lower" and "upper", respectively. Defaults to \code{NULL}.
#' @param h2 A number representing the heritability on liability scale
#' for a single phenotype. Must be non-negative. Note that under the liability threshold model,
#' the heritability must also be at most 1.
#' Defaults to 0.5.
#' @param  pid A string holding the name of the column in \code{.tbl} (or \code{family} and 
#' \code{threshs}) that hold the personal identifier(s). Defaults to "PID".
#' @param fam_id A string holding the name of the column in \code{.tbl} or \code{family} that
#' holds the family identifier. Defaults to "fam_ID".
#' @param role A string holding the name of the column in \code{.tbl} that
#' holds the role. Each role must be chosen from the following list of abbreviations 
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
#' Defaults to "role".
#' @param out A character or numeric vector indicating whether the genetic component
#' of the full liability, the full liability or both should be returned. If \code{out = c(1)} or 
#' \code{out = c("genetic")}, the genetic liability is estimated and returned. If \code{out = c(2)} or 
#' \code{out = c("full")}, the full liability is estimated and returned. If \code{out = c(1,2)} or 
#' \code{out = c("genetic", "full")}, both components are estimated and returned. 
#' Defaults to \code{c(1)}.
#' @param tol A number that is used as the convergence criterion for the Gibbs sampler.
#' Equals the standard error of the mean. That is, a tolerance of 0.2 means that the 
#' standard error of the mean is below 0.2. Defaults to 0.01.
#' 
#' @return If \code{family} and \code{threshs} are two matrices, lists or 
#' data frames that can be converted into tibbles, if \code{family} has two 
#' columns named like the strings represented in \code{pid} and \code{fam_id}, if 
#' \code{threshs} has a column named like the string given in \code{pid} as 
#' well as a column named "lower" and a column named "upper" and if the 
#' liability-scale heritability \code{h2}, \code{out}, \code{tol} and 
#' \code{always_add} are of the required form, then the function returns a 
#' tibble with either four or six columns (depending on the length of out).
#' The first two columns correspond to the columns \code{fam_id} and \code{pid} '
#' present in \code{family}. 
#' If \code{out} is equal to \code{c(1)} or \code{c("genetic")}, the third 
#' and fourth column hold the estimated genetic liability as well as the 
#' corresponding standard error, respectively. 
#' If \code{out} equals \code{c(2)} or \code{c("full")}, the third and 
#' fourth column hold the estimated full liability as well as the 
#' corresponding standard error, respectively. 
#' If \code{out} is equal to \code{c(1,2)} or \code{c("genetic","full")},
#' the third and fourth column hold the estimated genetic liability as 
#' well as the corresponding standard error, respectively, while the fifth and
#' sixth column hold the estimated full liability as well as the corresponding 
#' standard error, respectively.
#' 
#' @examples
#' sims <- simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, 
#' add_ind = TRUE, h2 = 0.5, n_sim=10, pop_prev = .05)
#' #
#' estimate_liability_single(.tbl = sims$thresholds, 
#' h2 = 0.5, pid = "indiv_ID", fam_id = "fam_ID", role = "role", out = c(1), 
#' tol = 0.01)
#' # 
#' sims <- simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, 
#' h2 = 0.5, n_sim=10, pop_prev = .05)
#' #
#' estimate_liability_single(.tbl = sims$thresholds, 
#' h2 = 0.5, pid = "indiv_ID", fam_id = "fam_ID", role = "role",
#' out = c("genetic"), tol = 0.01)
#' 
#' @seealso \code{\link[future.apply]{future_apply}}, \code{\link{estimate_liability_multi}},
#' \code{\link{estimate_liability}}
#' 
#' @importFrom dplyr %>% pull bind_rows bind_cols select filter
#' @importFrom stringr str_ends str_split str_subset str_detect
#' @importFrom rlang :=
#' 
#' @export
estimate_liability_single <- function(.tbl = NULL, 
                                      family = NULL, 
                                      threshs = NULL, 
                                      h2 = 0.5, 
                                      pid = "PID", 
                                      fam_id = "fam_ID", 
                                      role = "role",
                                      out = c(1), 
                                      tol = 0.01){
  
# Making sure input is valid ----------------------------------------------
  
  # Turning .tbl into a tibble
  # if it is not of class tbl
  if (!is.null(.tbl) && !tibble::is_tibble(.tbl))  .tbl <- tibble::as_tibble(.tbl)
  
  # If .tbl is NULL and family and threshs are given as the input,
  # the data must be converted to the correct format
  if (is.null(.tbl) || nrow(.tbl) == 0) .tbl <- convert_format(family = family, threshs = threshs, personal_id_col = pid, role_col = role)
  
  # # Turning progress into class logical
  # progress <- as.logical(progress)
  
  # Turning pid, fam_id and role into strings
  pid <- as.character(pid)
  fam_id <- as.character(fam_id)
  role <- as.character(role)
  
  # Checking that the heritability is valid
  if(validate_proportion(h2)){invisible()}
  
  # Checking that .tbl has three columns named pid, fam_id and role
  if(!(pid %in% colnames(.tbl))) stop(paste0("The column ", pid," does not exist in the tibble .tbl..."))
  if(!(fam_id %in% colnames(.tbl))) stop(paste0("The column ", fam_id," does not exist in the tibble .tbl..."))
  if(!(role %in% colnames(.tbl))) stop(paste0("The column ", role," does not exist in the tibble .tbl..."))
  
  # In addition, we check that two columns named lower and upper are present
  if(any(!c("lower","upper") %in% colnames(.tbl))) stop("The tibble .tbl must include two columns named 'lower' and 'upper'!")
  
  # Checking that tol is valid
  if(!is.numeric(tol) && !is.integer(tol)) stop("The tolerance must be numeric!")
  if(tol <= 0) stop("The tolerance must be strictly positive!")
  
  # Checking that out is either a character vector or a
  # numeric vector 
  if(is.numeric(out)){
    
    out <- intersect(out, c(1,2))
    
  }else if(is.character(out)){
    
    out <- c("genetic", "full")[rowSums(sapply(out, grepl, x = c("genetic", "full"))) > 0]
    out[out == "genetic"] <- 1
    out[out == "full"] <- 2
    out <- as.numeric(out)
  }else{
    stop("out must be a numeric or character vector!")
  }
  
  # Checking whether out is empty
  if(length(out) == 0){
    
    cat("Warning message: \n out is not of the required format! \n The function will return the estimated genetic liability!")
    out <- c(1)
  }
  
  # Sorting out
  out <- sort(out)
  
  # If the tibble consists of more than the required columns, 
  # we select only the relevant ones.
  .tbl <- select(.tbl, !!as.symbol(fam_id), !!as.symbol(pid), !!as.symbol(role), tidyselect::starts_with("lower"), tidyselect::starts_with("upper"))
  
  # Finally, we also check whether all lower thresholds are 
  # smaller than or equal to the upper thresholds
  if(any(pull(.tbl, lower) > pull(.tbl, upper))){
    cat("Warning message: \n Some lower thresholds are larger than the corresponding upper thresholds! \n
  The lower and upper thresholds will be swapped...")
    
    swapping_indx <- which(pull(.tbl, lower) > pull(.tbl, upper))
    
    .tbl$lower[swapping_indx] <- .tbl$lower[swapping_indx] + .tbl$upper[swapping_indx]
    .tbl$upper[swapping_indx] <- .tbl$lower[swapping_indx] - .tbl$upper[swapping_indx]
    .tbl$lower[swapping_indx] <- .tbl$lower[swapping_indx] - .tbl$upper[swapping_indx]
  }
  

  # Extracting the family identifiers
  fam_list <- unique(pull(.tbl, !!as.symbol(fam_id)))

  
  cat(paste0("The number of workers is ", future::nbrOfWorkers(), "\n"))
  
  # # If progress = TRUE, a progress bar will be displayed
  # if(progress){
  #   pb <- utils::txtProgressBar(min = 0, max = nrow(family), style = 3, char = "=")
  #   j <- 0
  # }
  #   
  gibbs_res <- future.apply::future_lapply(X = 1:length(fam_list), FUN = function(i){
    
    # Extract the personal numbers and roles for all
    # family members
    pids <- pull(filter(.tbl, !!as.symbol(fam_id) == fam_list[i]), !!as.symbol(pid))
    roles <- pull(filter(.tbl, !!as.symbol(fam_id) == fam_list[i]), !!as.symbol(role))
    
    # Constructing the covariance matrix.
    cov <- construct_covmat(fam_vec = roles, n_fam = NULL, add_ind = TRUE, h2 = h2)
    
    # Extracting the thresholds for all family members
    # including the thresholds for o and/or g
    temp_tbl = filter(.tbl, !!as.symbol(fam_id) == fam_list[i])
    
    # If the thresholds for o and/or g are missing, they will
    # be set to -Inf (lower threshold) and Inf (upper 
    # threshold)
    
    if(setdiff(c("g","o"), pull(temp_tbl, !!as.symbol(role))) == "g") {
      
      if("o" %in% pull(temp_tbl, !!as.symbol(role))) i_pid <- pull(filter(temp_tbl, !!as.symbol(role)== "o"), !!as.symbol(pid))
      if(!"o" %in% pull(temp_tbl, !!as.symbol(role))) i_pid <- paste0(pull(temp_tbl, !!as.symbol(fam_id))[1], "_g")
      
      temp_tbl <- bind_rows(tibble::tibble(!!as.symbol(fam_id) := pull(temp_tbl, !!as.symbol(fam_id))[1], 
                                           !!as.symbol(pid) := i_pid, 
                                           !!as.symbol(role) := "g",
                                           lower = -Inf, 
                                           upper = Inf), 
                            temp_tbl)
    }else if(setdiff(c("g","o"), pull(temp_tbl, !!as.symbol(role))) == "o"){
      
      if("g" %in% pull(temp_tbl, !!as.symbol(role))) i_pid <- pull(filter(temp_tbl, !!as.symbol(role)== "g"), !!as.symbol(pid))
      if(!"g" %in% pull(temp_tbl, !!as.symbol(role))) i_pid <- paste0(pull(temp_tbl, !!as.symbol(fam_id))[1], "_o")
      
      temp_tbl <- bind_rows(tibble::tibble(!!as.symbol(fam_id) := pull(temp_tbl, !!as.symbol(fam_id))[1], 
                                           !!as.symbol(pid) := i_pid, 
                                           !!as.symbol(role) := "o",
                                           lower = -Inf, 
                                           upper = Inf), 
                            temp_tbl)
    }else if(sort(setdiff(c("g","o"), pull(temp_tbl, !!as.symbol(role)))) == c("g","o")){
      
      i_pid <- pull(temp_tbl, !!as.symbol(fam_id))[1]
      
      temp_tbl <- bind_rows(tibble::tibble(!!as.symbol(fam_id) := rep(pull(temp_tbl, !!as.symbol(fam_id))[1],2), 
                                           !!as.symbol(pid) := c(paste0(i_pid, "_g"), paste0(i_pid, "_o")), 
                                           !!as.symbol(role) := c("g","o"),
                                           lower = c(-Inf, -Inf), 
                                           upper = c(Inf, Inf)), 
                            temp_tbl)
    }else{}
    
    # Now that we have extracted all the relevant information, we 
    # only need to order the observations before we can run
    # Gibbs sampler, as g and o need to be the first two 
    # observations.
    first_indx <- match(c("g","o"), pull(temp_tbl, !!as.symbol(role)))
    other_indx <- setdiff(1:length(pull(temp_tbl, !!as.symbol(role))), first_indx)
    temp_tbl <- temp_tbl[c(first_indx, other_indx),]
    
    
    # Setting the variables needed for Gibbs sampler
    fixed <- (pull(temp_tbl,upper) - pull(temp_tbl,lower)) < 1e-04
    std_err <- rep(Inf, length(out))
    names(std_err) <- c("genetic", "full")[out]
    n_gibbs <- 1
    
    # Running Gibbs sampler
    while(any(std_err > tol)){
      
      if(n_gibbs == 1){
        
        est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, covmat = cov, lower = pull(temp_tbl, lower), upper = pull(temp_tbl, upper),
                                    fixed = fixed, out = out, burn_in = 1000) %>% 
          `colnames<-`(c("genetic", "full")[out]) %>% 
          tibble::as_tibble()
        
      }else{
        
        est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, covmat = cov, lower = pull(temp_tbl, lower), upper = pull(temp_tbl, upper),
                                    fixed = fixed, out = out, burn_in = 1000) %>%
          `colnames<-`(c("genetic", "full")[out]) %>%
          tibble::as_tibble() %>% 
          bind_rows(est_liabs)
      }
      
      # Computing the standard error
      std_err <- batchmeans::bmmat(est_liabs)[,2]
      # Adding one to the counter
      n_gibbs <- n_gibbs +1
    }
    
    # If all standard errors are below the tolerance, 
    # the estimated liabilities as well as the corresponding 
    # standard error can be returned
    res <- tibble::as_tibble(batchmeans::bmmat(est_liabs), rownames = "out") %>% 
      tidyr::pivot_longer(., cols = c(est,se)) %>% 
      mutate(., name = paste0(out, "_", name), .keep = "unused") %>% 
      tibble::deframe(.)
    
    res <- c(unlist(select(filter(temp_tbl, !!as.symbol(role)=="o"), !!as.symbol(fam_id), !!as.symbol(pid))),
             res)
    
    return(res)
    
    # If progress = TRUE, a progress bar will be displayed
    # if(progress){
    #   j <- j+1
    #   utils::setTxtProgressBar(pb, j)
    # }
    
    
  }, future.seed = TRUE) %>% 
    do.call("bind_rows", .)
    
  # # Close the connection 
  # if(progress) close(pb)  

  return(gibbs_res)
}


#' Estimating the genetic or full liability for multiple phenotypes
#' 
#' \code{estimate_liability_multi} estimates the genetic component of the full
#' liability and/or the full liability for a number of individuals based
#' on their family history for a variable number of phenotypes.
#'
#' This function can be used to estimate either the genetic component of the 
#' full liability, the full liability or both for a variable number of traits.
#'
#' @param .tbl A matrix, list or data frame that can be converted into a tibble.
#' Must have at least seven columns that hold the family identifier, the personal 
#' identifier, the role and the lower and upper thresholds for all phenotypes
#' of interest. Note that the role must be one of the following abbreviations 
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
#' Defaults to \code{NULL}.
#' @param family A matrix, list or data frame that can be converted into a tibble.
#' Must have at least two columns that hold the family identifier and the corresponding
#' personal identifiers, respectively. That is, for each family in fam_id there should
#' be a list holding all individuals belonging to that family in pid. Note that the 
#' personal identifiers for all individuals must have a special format. It must be end
#' with _?, where ? is one of the following abbreviations 
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
#' See also \code{\link{construct_covmat}}.
#' Defaults to \code{NULL}.
#' @param threshs A matrix, list or data frame that can be converted into a tibble.
#' Must have at least five columns; one holding the personal identifier for all individuals,
#' and the remaining four holding the lower and upper thresholds for the first and second
#' phenotype, respectively. It must be possible to tie each pair of lower and upper thresholds
#' to a specific phenotype uniquely. This is done easily by adding _{name_of_phenotype} to
#' the column names lower and upper, e.g. lower_p1 and upper_p1 for the lower and upper
#' thresholds corresponding to the first phenotype. Defaults to \code{NULL}.
#' @param genetic_corrmat A numeric matrix holding the genetic correlations between the desired 
#' phenotypes. All diagonal entries must be equal to one, while all off-diagonal entries 
#' must be between -1 and 1. In addition, the matrix must be symmetric.
#' @param full_corrmat A numeric matrix holding the full correlations between the desired 
#' phenotypes. All diagonal entries must be equal to one, while all off-diagonal entries 
#' must be between -1 and 1. In addition, the matrix must be symmetric.
#' @param h2_vec A numeric vector representing the liability-scale heritabilities
#' for all phenotypes. All entries in h2_vec must be non-negative and at most 1.
#' @param phen_names A character vector holding the phenotype names. These names
#' will be used to create the row and column names for the covariance matrix.
#' If it is not specified, the names will default to phenotype1, phenotype2, etc.
#' Defaults to NULL.
#' @param  pid A string holding the name of the column in \code{family} and 
#' \code{threshs} that hold the personal identifier(s). Defaults to "PID".
#' @param fam_id A string holding the name of the column in \code{family} that
#' holds the family identifier. Defaults to "fam_ID".
#' @param role A string holding the name of the column in \code{.tbl} that
#' holds the role.Each role must be chosen from the following list of abbreviations 
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
#' Defaults to "role".
#' @param out A character or numeric vector indicating whether the genetic component
#' of the full liability, the full liability or both should be returned. If \code{out = c(1)} or 
#' \code{out = c("genetic")}, the genetic liability is estimated and returned. If \code{out = c(2)} or 
#' \code{out = c("full")}, the full liability is estimated and returned. If \code{out = c(1,2)} or 
#' \code{out = c("genetic", "full")}, both components are estimated and returned. 
#' Defaults to \code{c(1)}.
#' @param tol A number that is used as the convergence criterion for the Gibbs sampler.
#' Equals the standard error of the mean. That is, a tolerance of 0.2 means that the 
#' standard error of the mean is below 0.2. Defaults to 0.01.
#' 
#' @return If \code{family} and \code{threshs} are two matrices, lists or data frames 
#' that can be converted into tibbles, if \code{family} has two columns named like 
#' the strings represented in \code{pid} and \code{fam_id}, if \code{threshs} has a 
#' column named like the string given in \code{pid} as well as a column named \code{"lower"} 
#' and a column named \code{"upper"} and if the liability-scale heritabilities in \code{h2_vec}, 
#' \code{genetic_corrmat}, \code{full_corrmat}, \code{out} and \code{tol} are of the 
#' required form, then the function returns a tibble with at least six columns (depending 
#' on the length of out).
#' The first two columns correspond to the columns \code{fam_id} and \code{pid} present in 
#' the tibble \code{family}. 
#' If \code{out} is equal to \code{c(1)} or \code{c("genetic")}, the third and fourth columns 
#' hold the estimated genetic liability as well as the corresponding standard error for the 
#' first phenotype, respectively. 
#' If \code{out} equals \code{c(2)} or \code{c("full")}, the third and fourth columns hold 
#' the estimated full liability as well as the corresponding standard error for the first 
#' phenotype, respectively. 
#' If \code{out} is equal to \code{c(1,2)} or \code{c("genetic","full")}, the third and 
#' fourth columns hold the estimated genetic liability as well as the corresponding standard 
#' error for the first phenotype, respectively, while the fifth and sixth columns hold the 
#' estimated full liability as well as the corresponding standard error for the first 
#' phenotype, respectively.
#' The remaining columns hold the estimated genetic liabilities and/or the estimated full 
#' liabilities as well as the corresponding standard errors for the remaining phenotypes.
#' 
#' @examples 
#' genetic_corrmat <- matrix(0.4, 3, 3)
#' diag(genetic_corrmat) <- 1
#' full_corrmat <- matrix(0.6, 3, 3)
#' diag(full_corrmat) <- 1
#' #
#' sims <- simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = TRUE, 
#' genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, h2 = rep(.5,3), 
#' n_sim = 100, pop_prev = rep(.1,3))
#' # 
#' estimate_liability_multi(.tbl = sims$thresholds, h2_vec = rep(.5,3), 
#' genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat,
#' pid = "indiv_ID", fam_id = "fam_ID", role = "role", out = c(1), 
#' tol = 0.01)
#' 
#' 
#' @seealso \code{\link[future.apply]{future_apply}}, \code{\link{estimate_liability_single}},
#' \code{\link{estimate_liability}}
#' 
#' @importFrom dplyr %>% pull bind_rows bind_cols select row_number rename arrange left_join
#' @importFrom rlang :=
#' @importFrom stringr str_detect
#' @importFrom tidyselect starts_with
#' 
#' @export
estimate_liability_multi <- function(.tbl = NULL, 
                                     family = NULL, 
                                     threshs = NULL, 
                                     h2_vec, 
                                     genetic_corrmat, 
                                     full_corrmat,
                                     phen_names = NULL, 
                                     pid = "PID",
                                     fam_id = "fam_ID",
                                     role = "role",
                                     out = c(1), 
                                     tol = 0.01){
  
  # Making sure input is valid--------------------------------------------
  
  # Turning .tbl into a tibble
  # if it is not of class tbl
  if(!tibble::is_tibble(.tbl))  .tbl <- tibble::as_tibble(.tbl)
  
  # If .tbl is NULL and family and threshs are given as the input,
  # the data must be converted to the correct format
  if (is.null(.tbl) || nrow(.tbl) == 0) .tbl <- convert_format(family = family, threshs = threshs, personal_id_col = pid, role_col = role)
  
  # Turning progress into class logical
  # progress <- as.logical(progress)
  
  # Turning pid, fam_id and role into strings
  pid <- as.character(pid)
  fam_id <- as.character(fam_id)
  role <- as.character(role)
  
  # Checking that the heritability is valid
  if(validate_proportion(h2_vec)){invisible()}
  
  # Checking that all correlations are valid
  if(validate_correlation_matrix(genetic_corrmat)){invisible()}
  if(validate_correlation_matrix(full_corrmat)){invisible()}
  
  # Checking that .tbl has three columns named pid_col, fam_id and role
  if(!(pid %in% colnames(.tbl))) stop(paste0("The column ", pid," does not exist in the tibble .tbl..."))
  if(!(fam_id %in% colnames(.tbl))) stop(paste0("The column ", fam_id," does not exist in the tibble .tbl..."))
  if(!(role %in% colnames(.tbl))) stop(paste0("The column ", role," does not exist in the tibble .tbl..."))
  
  # In addition, we check that the two columns lower and upper are present
  if(any(!c("lower","upper") %in% gsub("_.*", "", colnames(.tbl), ignore.case = TRUE))) stop("The tibble .tbl must include two columns named 'lower' and 'upper'!")
  
  # Checking that tol is valid
  if(!is.numeric(tol) && !is.integer(tol)) stop("The tolerance must be numeric!")
  if(tol <= 0) stop("The tolerance must be strictly positive!")
  
  # Checking that out is either a character vector or a
  # numeric vector 
  if(is.numeric(out)){
    
    out <- intersect(out, c(1,2))
    
  }else if(is.character(out)){
    
    out <- c("genetic", "full")[rowSums(sapply(out, grepl, x = c("genetic", "full"))) > 0]
    out[out == "genetic"] <- 1
    out[out == "full"] <- 2
    out <- as.numeric(out)
  }else{
    stop("out must be a numeric or character vector!")
  }
  
  # Checking whether out is empty
  if(length(out) == 0){
    
    cat("Warning message: \n out is not of the required format! \n The function will return the estimated genetic liability!")
    out <- c(1)
  }
  
  # Sorting out
  out <- sort(out)
  
  # If the tibble consists of more than the required columns, 
  # we select only the relevant ones.
  .tbl <- select(.tbl, !!as.symbol(fam_id), !!as.symbol(pid), !!as.symbol(role), tidyselect::starts_with("lower"), tidyselect::starts_with("upper"))
  
  # Now we can extract the number of phenotypes
  n_pheno <- length(h2_vec)
  
  if(ncol(.tbl) != (2*n_pheno + 3)) stop("Something is wrong with the number of phenotypes... \n 
The number of pairs of lower and upper thresholds is not equal to the number of phenotypes specified in h2_vec...\
Does all columns have the required names?")
  
  # As well as the phenotype names
  phen_names <- gsub("lower_", "", colnames(.tbl)[str_detect(colnames(.tbl), "^lower_")], ignore.case = TRUE)
  
  # Finally, we also check whether all lower thresholds are 
  # smaller than or equal to the upper thresholds
  for ( pheno in phen_names ) {
    
    
    if(any(pull(.tbl, !!as.symbol(paste0("lower_", pheno))) > pull(.tbl, !!as.symbol(paste0("upper_", pheno))))){
      cat("Warning message: \n Some lower thresholds are larger than the corresponding upper thresholds! \n
The lower and upper thresholds will be swapped...")
      
      swapping_indx <- which(pull(.tbl, !!as.symbol(paste0("lower_", pheno))) > pull(.tbl, !!as.symbol(paste0("upper_", pheno))))
      
      .tbl <- mutate(.tbl, !!as.symbol(paste0("lower_", pheno)) := ifelse(row_number() %in% swapping_indx, !!as.symbol(paste0("lower_", pheno)) + !!as.symbol(paste0("upper_", pheno)), !!as.symbol(paste0("lower_", pheno)))) %>% 
        mutate(., !!as.symbol(paste0("upper_", pheno)) := ifelse(row_number() %in% swapping_indx, !!as.symbol(paste0("lower_", pheno)) - !!as.symbol(paste0("upper_", pheno)), !!as.symbol(paste0("upper_", pheno)))) %>% 
        mutate(., !!as.symbol(paste0("lower_", pheno)) := ifelse(row_number() %in% swapping_indx, !!as.symbol(paste0("lower_", pheno)) - !!as.symbol(paste0("upper_", pheno)), !!as.symbol(paste0("lower_", pheno))))
    }
  }
  
  #  Extracting the family identifiers
  fam_list <- unique(pull(.tbl, !!as.symbol(fam_id)))
  
  cat(paste0("The number of workers is ", future::nbrOfWorkers(), "\n"))
  
  # If progress = TRUE, a progress bar will be displayed
  # if(progress){
  #   pb <- utils::txtProgressBar(min = 0, max = nrow(family), style = 3, char = "=")
  #   j <- 0
  # }
  
  gibbs_res <- future.apply::future_lapply(X= 1:length(fam_list), FUN = function(i){
    
    # Extract the personal numbers and roles for all
    # family members
    pids <- pull(filter(.tbl, !!as.symbol(fam_id) == fam_list[i]), !!as.symbol(pid))
    roles <- pull(filter(.tbl, !!as.symbol(fam_id) == fam_list[i]), !!as.symbol(role))
    
    # Constructing the covariance matrix.
    cov <- construct_covmat(fam_vec = roles, n_fam = NULL, add_ind = TRUE, 
                            genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat,
                            h2 = h2_vec, phen_names = phen_names)
    
    # Extracting the thresholds for all family members,
    # including the thresholds for o and/or g,
    # and all phenotypes
    temp_tbl = filter(.tbl, !!as.symbol(fam_id) == fam_list[i])
    
    # If the thresholds for o and/or g are missing, they will
    # be set to -Inf (lower threshold) and Inf (upper 
    # threshold)
    
    if(setdiff(c("g","o"), pull(temp_tbl, !!as.symbol(role))) == "g") {
      
      if("o" %in% pull(temp_tbl, !!as.symbol(role))) i_pid <- pull(filter(temp_tbl, !!as.symbol(role)== "o"), !!as.symbol(pid))
      if(!"o" %in% pull(temp_tbl, !!as.symbol(role))) i_pid <- paste0(pull(temp_tbl, !!as.symbol(fam_id))[1], "_g")
      
      temp_tbl <- bind_rows(bind_cols(tibble::tibble(!!as.symbol(fam_id) := pull(temp_tbl, !!as.symbol(fam_id))[1], 
                                                     !!as.symbol(pid) := i_pid, 
                                                     !!as.symbol(role) := "g"),
                                      tibble::tibble(!!!c(stats::setNames(rep(-Inf, n_pheno), paste0("lower_", phen_names)), stats::setNames(rep(Inf, n_pheno), paste0("upper_", phen_names))))), 
                            temp_tbl)
    }else if(setdiff(c("g","o"), pull(temp_tbl, !!as.symbol(role))) == "o"){
      
      if("g" %in% pull(temp_tbl, !!as.symbol(role))) i_pid <- pull(filter(temp_tbl, !!as.symbol(role)== "g"), !!as.symbol(pid))
      if(!"g" %in% pull(temp_tbl, !!as.symbol(role))) i_pid <- paste0(pull(temp_tbl, !!as.symbol(fam_id))[1], "_o")
      
      temp_tbl <- bind_rows(bind_cols(tibble::tibble(!!as.symbol(fam_id) := pull(temp_tbl, !!as.symbol(fam_id))[1], 
                                                     !!as.symbol(pid) := i_pid, 
                                                     !!as.symbol(role) := "o"),
                                      tibble::tibble(!!!c(stats::setNames(rep(-Inf, n_pheno), paste0("lower_", phen_names)), stats::setNames(rep(Inf, n_pheno), paste0("upper_", phen_names))))), 
                            temp_tbl)
    }else if(sort(setdiff(c("g","o"), pull(temp_tbl, !!as.symbol(role)))) == c("g","o")){
      
      i_pid <- pull(temp_tbl, !!as.symbol(fam_id))[1]
      
      temp_tbl <- bind_rows(bind_cols(tibble::tibble(!!as.symbol(fam_id) := pull(temp_tbl, !!as.symbol(fam_id))[1], 
                                                     !!as.symbol(pid) := i_pid, 
                                                     !!as.symbol(role) := "g"),
                                      tibble::tibble(!!!c(stats::setNames(rep(-Inf, n_pheno), paste0("lower_", phen_names)), stats::setNames(rep(Inf, n_pheno), paste0("upper_", phen_names))))), 
                            bind_cols(tibble::tibble(!!as.symbol(fam_id) := pull(temp_tbl, !!as.symbol(fam_id))[1], 
                                                     !!as.symbol(pid) := i_pid, 
                                                     !!as.symbol(role) := "o"),
                                      tibble::tibble(!!!c(stats::setNames(rep(-Inf, n_pheno), paste0("lower_", phen_names)), stats::setNames(rep(Inf, n_pheno), paste0("upper_", phen_names))))), 
                            temp_tbl)
    }else{}
    
    
    
    # Now that we have extracted all the relevant information, we 
    # only need to order the observations before we can turn the 
    # format from wide to long and run Gibbs sampler.
    # This is due to the fact that g and o need to be the first two 
    # observations.
    first_indx <- match(c("g","o"), pull(temp_tbl, !!as.symbol(role)))
    other_indx <- setdiff(1:length(pull(temp_tbl, !!as.symbol(role))), first_indx)
    
    temp_tbl <- mutate(temp_tbl, !!as.symbol(role) := factor(!!as.symbol(role), levels = pull(temp_tbl,!!as.symbol(role))[c(first_indx, other_indx)])) %>% 
      arrange(., !!as.symbol(role))
    
    # Now, we need to lengthen the data
    fam_threshs <- select(temp_tbl, -c(starts_with("upper"))) %>% 
      tidyr::pivot_longer(., cols = starts_with("lower"), names_to = "phenotype", values_to = "lower") %>% 
      mutate(., phenotype = gsub("lower_","",phenotype))
    
    fam_threshs <- select(temp_tbl, -c(starts_with("lower"))) %>% 
      tidyr::pivot_longer(., cols = starts_with("upper"), names_to = "phenotype", values_to = "upper") %>% 
      mutate(., phenotype = gsub("upper_","",phenotype)) %>% 
      left_join(fam_threshs,., by = stats::setNames(c(fam_id,pid,role,"phenotype"), c(fam_id,pid,role,"phenotype"))) %>% 
      mutate(., phenotype = factor(phenotype, levels = phen_names)) %>% 
      arrange(., phenotype, !!as.symbol(role))
    
    
    # Setting the variables needed for Gibbs sampler
    fixed <- (pull(fam_threshs,upper) - pull(fam_threshs,lower)) < 1e-04
    std_err <- matrix(Inf, ncol =  length(out), nrow = n_pheno)
    colnames(std_err) <- c("genetic", "full")[out]
    n_gibbs <- 1
    
    # And change the variable out, as we need to extract the
    # genetic and/or full liability for each phenotype.
    # And as all thresholds in fam_threshs are ordered according
    # to the phenotype and the role, we need to extract every 
    # (number of individuals)'th entry starting from the entries
    # specified in out.
    updated_out <- sort(sapply(out, function(k) k + length(levels(pull(fam_threshs, !!as.symbol(role))))*(0:(n_pheno-1))))
    
    # Running Gibbs sampler
    while(any(std_err > tol)){
      
      if(n_gibbs == 1){
        
        est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, covmat = cov, lower = pull(fam_threshs, lower), upper = pull(fam_threshs, upper),
                                    fixed = fixed, out = updated_out, burn_in = 1000) %>% 
          `colnames<-`(paste0(rep(c("genetic", "full")[out], n_pheno),"_", rep(phen_names, each = length(out)))) %>% 
          tibble::as_tibble()
        
      }else{
        
        est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, covmat = cov, lower = pull(fam_threshs, lower), upper = pull(fam_threshs, upper),
                                    fixed = fixed, out = updated_out, burn_in = 1000) %>% 
          `colnames<-`(paste0(rep(c("genetic", "full")[out], n_pheno),"_", rep(phen_names, each = length(out)))) %>% 
          tibble::as_tibble() %>% 
          bind_rows(est_liabs,.)
      }
      
      # Computing the standard error
      std_err[1:n_pheno,1:length(out)] <- matrix(batchmeans::bmmat(est_liabs)[,2], ncol =  length(out), nrow = n_pheno, byrow = TRUE)
      # Adding one to the counter
      n_gibbs <- n_gibbs +1
    }
    
    # # If progress = TRUE, a progress bar will be displayed
    # if(progress){
    #   j <- j+1
    #   utils::setTxtProgressBar(pb, j)
    # }
    
    # If all standard errors are below the tolerance, 
    # the estimated liabilities as well as the corresponding 
    # standard error can be returned
    res <- tibble::as_tibble(batchmeans::bmmat(est_liabs), rownames = "out") %>% 
      tidyr::pivot_longer(., cols = c(est,se)) %>% 
      mutate(., name = paste0(out, "_", name), .keep = "unused") %>% 
      tibble::deframe(.)
    
    res <- c(unlist(select(filter(temp_tbl, !!as.symbol(role)=="o"), !!as.symbol(fam_id), !!as.symbol(pid))),
             res)
    
    return(res)
    
  }, future.seed = TRUE) %>% 
    do.call("bind_rows",.)
  
  # # Close the connection 
  # if(progress) close(pb)  
  
  return(gibbs_res)
}


#' Estimating the genetic or full liability for a variable number of
#' phenotypes
#' 
#' \code{estimate_liability} estimates the genetic component of the full
#' liability and/or the full liability for a number of individuals based
#' on their family history for one or more phenotypes.  It is a wrapper around 
#' \code{\link{estimate_liability_single}} and \code{\link{estimate_liability_multi}}.
#'
#' This function can be used to estimate either the genetic component of the 
#' full liability, the full liability or both for a variable number of traits.
#'
#' @param .tbl A matrix, list or data frame that can be converted into a tibble.
#' Must have at least five columns that hold the family identifier, the personal 
#' identifier, the role and the lower and upper thresholds for all phenotypes
#' of interest. Note that the role must be one of the following abbreviations 
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
#' Defaults to \code{NULL}.
#' @param family A matrix, list or data frame that can be converted into a tibble.
#' Must have at least two columns that hold the family identifier and the corresponding
#' personal identifiers, respectively. That is, for each family in fam_id there should
#' be a list holding all individuals belonging to that family in pid. Note that the 
#' personal identifiers for all individuals must have a special format. It must be end
#' with _?, where ? is one of the following abbreviations 
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
#' See also \code{\link{construct_covmat}}.
#' Defaults to \code{NULL}.
#' @param threshs A matrix, list or data frame that can be converted into a tibble.
#' Must have at least five columns; one holding the personal identifier for all individuals,
#' and the remaining four holding the lower and upper thresholds for the first and second
#' phenotype, respectively. It must be possible to tie each pair of lower and upper thresholds
#' to a specific phenotype uniquely. This is done easily by adding _{name_of_phenotype} to
#' the column names lower and upper, e.g. lower_p1 and upper_p1 for the lower and upper
#' thresholds corresponding to the first phenotype. Defaults to \code{NULL}.
#' @param h2 Either a number representing the heritability on liability scale for a 
#' single phenotype, or a numeric vector representing the liability-scale heritabilities
#' for all phenotypes. All entries in h2 must be non-negative and at most 1.
#' @param  pid A string holding the name of the column in \code{family} and 
#' \code{threshs} that hold the personal identifier(s). Defaults to \code{"PID"}.
#' @param fam_id A string holding the name of the column in \code{family} that
#' holds the family identifier. Defaults to \code{"fam_ID"}.
#' @param role A string holding the name of the column in \code{.tbl} that
#' holds the role.Each role must be chosen from the following list of abbreviations 
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
#' Defaults to "role".
#' @param out A character or numeric vector indicating whether the genetic component
#' of the full liability, the full liability or both should be returned. If \code{out = c(1)} or 
#' \code{out = c("genetic")}, the genetic liability is estimated and returned. If \code{out = c(2)} or 
#' \code{out = c("full")}, the full liability is estimated and returned. If code{out = c(1,2)} or 
#' \code{out = c("genetic", "full")}, both components are estimated and returned. 
#' Defaults to \code{c(1)}.
#' @param tol A number that is used as the convergence criterion for the Gibbs sampler.
#' Equals the standard error of the mean. That is, a tolerance of 0.2 means that the 
#' standard error of the mean is below 0.2. Defaults to 0.01.
#' @param genetic_corrmat Either \code{NULL} (if \code{h2} is a number) or a numeric 
#' matrix (if \code{h2} is a vector of length > 1) holding the genetic correlations 
#' between the desired phenotypes. All diagonal entries must be equal to one, while 
#' all off-diagonal entries must be between -1 and 1. In addition, the matrix must 
#' be symmetric. Defaults to \code{NULL}.
#' @param full_corrmat Either \code{NULL} (if \code{h2} is a number) or a numeric 
#' matrix (if \code{h2} is a vector of length > 1) holding the full correlations 
#' between the desired phenotypes. All diagonal entries must be equal to one, while 
#' all off-diagonal entries must be between -1 and 1. In addition, the matrix must 
#' be symmetric. Defaults to \code{NULL}.
#' @param phen_names Either \code{NULL} or a character vector holding the phenotype 
#' names. These names will be used to create the row and column names for the 
#' covariance matrix. If it is not specified, the names will default to 
#' phenotype1, phenotype2, etc. Defaults to NULL.
#' 
#' @return If \code{family} and \code{threshs} are two matrices, lists or 
#' data frames that can be converted into tibbles, if \code{family} has two 
#' columns named like the strings represented in \code{pid} and \code{fam_id}, if 
#' \code{threshs} has a column named like the string given in \code{pid} as 
#' well as a column named "lower" and a column named "upper" and if the 
#' liability-scale heritability \code{h2} is a number (\code{length(h2)=1}), 
#' and \code{out}, \code{tol} and 
#' \code{always_add} are of the required form, then the function returns a 
#' tibble with either four or six columns (depending on the length of out).
#' The first two columns correspond to the columns \code{fam_id} and \code{pid} '
#' present in \code{family}. 
#' If \code{out} is equal to \code{c(1)} or \code{c("genetic")}, the third 
#' and fourth column hold the estimated genetic liability as well as the 
#' corresponding standard error, respectively. 
#' If \code{out} equals \code{c(2)} or \code{c("full")}, the third and 
#' fourth column hold the estimated full liability as well as the 
#' corresponding standard error, respectively. 
#' If \code{out} is equal to \code{c(1,2)} or \code{c("genetic","full")},
#' the third and fourth column hold the estimated genetic liability as 
#' well as the corresponding standard error, respectively, while the fifth and
#' sixth column hold the estimated full liability as well as the corresponding 
#' standard error, respectively.
#' If \code{h2} is a numeric vector of length greater than 1 and if 
#' \code{genetic_corrmat}, \code{full_corrmat}, \code{out} and \code{tol} are of the 
#' required form, then the function returns a tibble with at least six columns (depending 
#' on the length of out).
#' The first two columns correspond to the columns \code{fam_id} and \code{pid} present in 
#' the tibble \code{family}. 
#' If \code{out} is equal to \code{c(1)} or \code{c("genetic")}, the third and fourth columns 
#' hold the estimated genetic liability as well as the corresponding standard error for the 
#' first phenotype, respectively. 
#' If \code{out} equals \code{c(2)} or \code{c("full")}, the third and fourth columns hold 
#' the estimated full liability as well as the corresponding standard error for the first 
#' phenotype, respectively. 
#' If \code{out} is equal to \code{c(1,2)} or \code{c("genetic","full")}, the third and 
#' fourth columns hold the estimated genetic liability as well as the corresponding standard 
#' error for the first phenotype, respectively, while the fifth and sixth columns hold the 
#' estimated full liability as well as the corresponding standard error for the first 
#' phenotype, respectively.
#' The remaining columns hold the estimated genetic liabilities and/or the estimated full 
#' liabilities as well as the corresponding standard errors for the remaining phenotypes.
#' 
#' @seealso \code{\link[future.apply]{future_apply}}, \code{\link{estimate_liability_single}},
#' \code{\link{estimate_liability_multi}}
#' 
#' @importFrom dplyr %>% pull bind_rows bind_cols select row_number rename
#' @importFrom rlang :=
#' 
#' @export
estimate_liability <- function(.tbl = NULL,
                               family = NULL, 
                               threshs = NULL, 
                               h2 = 0.5, 
                               pid = "PID", 
                               fam_id = "fam_ID", 
                               role = "role",
                               out = c(1), 
                               tol = 0.01,
                               genetic_corrmat = NULL, 
                               full_corrmat = NULL, 
                               phen_names = NULL){
  
  if(length(h2) == 1){
    
    return(estimate_liability_single(.tbl = .tbl, 
                                     family = family, 
                                     threshs = threshs, 
                                     h2 = h2, 
                                     pid = pid, 
                                     fam_id = fam_id, 
                                     role = role,
                                     out = out, 
                                     tol = tol))

  }else{
    
    return(estimate_liability_multi(.tbl = .tbl, 
                                    family = family, 
                                    threshs = threshs, 
                                    h2_vec = h2, 
                                    genetic_corrmat = genetic_corrmat, 
                                    full_corrmat = full_corrmat,
                                    phen_names = phen_names, 
                                    pid = pid,
                                    fam_id = fam_id,
                                    role = role,
                                    out = out, 
                                    tol = tol))
  } 
}
