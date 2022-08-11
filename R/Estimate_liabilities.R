utils::globalVariables("lower")
utils::globalVariables("upper")

#' Estimating the genetic or full liability 
#'
#' \code{estimate_liability_single} estimates the genetic component of the full
#' liability and/or the full liability for a number of individuals based
#' on their family history.
#'
#' This function can be used to estimate either the genetic component of the 
#' full liability, the full liability or both. 
#'
#' @param family A matrix, list or data frame that can be converted into a tibble.
#' Must have at least two columns that hold the family identifier and the corresponding
#' personal identifiers, respectively. That is, for each family in \code{fam_id} there should
#' be a list holding all individuals belonging to that family in \code{pid}. Note that the 
#' personal identifiers for all individuals must have a special format. It must end
#' on \code{_?}, where \code{?} is one of the following abbreviations 
#' - g (Genetic component of the full liability)
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
#' See also \code{\link{construct_covmat}}.
#' @param threshs A matrix, list or data frame that can be converted into a tibble.
#' Must have at least three columns, one holding the personal identifier for all individuals,
#' and the remaining two holding the lower and upper thresholds, respectively. The latter
#' should be called "lower" and "upper", respectively. 
#' @param h2 A number representing the heritability on liability scale
#' for a single phenotype. Must be non-negative. Note that under the liability threshold model,
#' the heritability must also be at most 1.
#' Defaults to 0.5.
#' @param  pid A string holding the name of the column in \code{family} and 
#' \code{threshs} that hold the personal identifier(s). Defaults to "PID".
#' @param fam_id A string holding the name of the column in \code{family} that
#' holds the family identifier. Defaults to "fam_ID".
#' @param out A character or numeric vector indicating whether the genetic component
#' of the full liability, the full liability or both should be returned. If \code{out = c(1)} or 
#' \code{out = c("genetic")}, the genetic liability is estimated and returned. If \code{out = c(2)} or 
#' \code{out = c("full")}, the full liability is estimated and returned. If \code{out = c(1,2)} or 
#' \code{out = c("genetic", "full")}, both components are estimated and returned. 
#' Defaults to \code{c(1)}.
#' @param tol A number that is used as the convergence criterion for the Gibbs sampler.
#' Equals the standard error of the mean. That is, a tolerance of 0.2 means that the 
#' standard error of the mean is below 0.2. Defaults to 0.01.
#' @param always_add A character vector or \code{NULL}. If \code{always_add = c("g","o")}, both the genetic component 
#' of the full liability as well as the full liability will be added to the list of family members. 
#' If always_add equals \code{"g"} or \code{"o"}, the genetic component of the full liability or the full liability
#' will be added, respectively. If \code{always_add = NULL}, no component will be added.
#' Defaults to \code{c("g","o")}.
#' @param progress A logical scalar indicating whether the function should display
#' a progress bar. Defaults to \code{FALSE}.
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
#' 
#' estimate_liability_single(family = sims$fam_ID, threshs = sims$thresholds, 
#' h2 = 0.5, pid = "indiv_ID", fam_id = "fam_ID", out = c(1), tol = 0.01,
#' always_add = c("g","o"), progress = TRUE)
#' # 
#' sims <- simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, 
#' h2 = 0.5, n_sim=10, pop_prev = .05)
#' 
#' estimate_liability_single(family = sims$fam_ID, threshs = sims$thresholds, 
#' h2 = 0.5, pid = "indiv_ID", fam_id = "fam_ID", out = c("genetic"), 
#' tol = 0.01, always_add = c("g","o"))
#' 
#' @seealso \code{\link[future.apply]{future_apply}}
#' 
#' @importFrom dplyr %>% pull bind_rows bind_cols select filter
#' @importFrom stringr str_ends str_split str_subset str_detect
#' @importFrom rlang :=
#' 
#' @export
estimate_liability_single <- function(family, threshs, h2 = 0.5, pid = "PID", fam_id = "fam_ID", out = c(1), 
                               tol = 0.01, always_add = c("g","o"), progress = FALSE){
  
# Making sure input is valid ----------------------------------------------

  # If always_add is a vector of length zero, it is set to
  # NULL instead
  if(length(always_add) == 0) always_add <- NULL
  
  # Turning progress into class logical
  progress <- as.logical(progress)
  
  # Turning pid and fam_id into strings
  pid <- as.character(pid)
  fam_id <- as.character(fam_id)
  
  #Turning family and threshs into tibbles
  # if they are not of class tbl
  if(!tibble::is_tibble(family))  family <- tibble::as_tibble(family)
  if(!tibble::is_tibble(threshs)) threshs <- tibble::as_tibble(threshs)
  
  # Checking that the heritability is valid
  if(validate_proportion(h2)){invisible()}
  
  # Checking that family has two columns named pid_col and fam_id
  if(!(pid %in% colnames(family))) stop(paste0("The column ", pid," does not exist in the tibble family..."))
  if(!(fam_id %in% colnames(family))) stop(paste0("The column ", fam_id," does not exist in the tibble family..."))
  
  # And that pid is also present in the tibble threshs
  if(!(pid %in% colnames(threshs))) stop(paste0("The column ", pid," does not exist in the tibble threshs..."))
  
  # In addition, we check that threshs has two columns named lower and upper
  if(any(!c("lower","upper") %in% colnames(threshs))) stop("The tibble threshs must include two columns named 'lower' and 'upper'!")
  
  # Checking that tol is valid
  if(!is.numeric(tol) && !is.integer(tol)) stop("The tolerance must be numeric!")
  if(tol <= 0) stop("The tolerance must be strictly positive!")
  
  # Checking that always_add is a vector of strings
  if(!is.character(always_add) && !is.null(always_add)) stop("always_add must be of class character!") 
  always_add <- intersect(always_add, c("g","o"))
  
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
  family <- select(family, !!as.symbol(fam_id), !!as.symbol(pid))
  threshs <- select(threshs, !!as.symbol(pid), tidyselect::starts_with("lower"), tidyselect::starts_with("upper"))
  
  # Finally, we also check whether all lower thresholds are 
  # smaller than or equal to the upper thresholds
  if(any(pull(threshs, lower) > pull(threshs, upper))){
    cat("Warning message: \n Some lower thresholds are larger than the corresponding upper thresholds! \n
        The lower and upper thresholds will be swapped...")
    
    swapping_indx <- which(pull(threshs, lower) > pull(threshs, upper))
    
    threshs$lower[swapping_indx] <- threshs$lower[swapping_indx] + threshs$upper[swapping_indx]
    threshs$upper[swapping_indx] <- threshs$lower[swapping_indx] - threshs$upper[swapping_indx]
    threshs$lower[swapping_indx] <- threshs$lower[swapping_indx] - threshs$upper[swapping_indx]
  }
  

  # Extracting the families
  fam_list <- pull(family, !!as.symbol(pid))
  
  cat(paste0("The number of workers is ", future::nbrOfWorkers(), "\n"))
  
  # # If progress = TRUE, a progress bar will be displayed
  # if(progress){
  #   pb <- utils::txtProgressBar(min = 0, max = nrow(family), style = 3, char = "=")
  #   j <- 0
  # }
  #   
  gibbs_res <- future.apply::future_lapply(X = 1:nrow(family), FUN = function(i){
    
    # Extract family members
    fam <- fam_list[[i]]
    # Remove individual o and/or g from the set (if present)
    fam <- setdiff(gsub("^.*_", "", fam), c("g","o"))
    
    # Constructing the covariance matrix
    # If always_add holds "g" or "o", add_ind must be TRUE
    cov <- construct_covmat(fam_vec = fam, n_fam = NULL, add_ind = length(always_add), h2 = h2)
    
    # Extracting the thresholds for all family members 
    # and all phenotypes
    fam_threshs = threshs[match(fam_list[[i]][!gsub("^.*_", "", fam_list[[i]]) %in% always_add], pull(threshs,!!as.symbol(pid))), ]
    
    ##################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(setdiff(always_add, intersect(gsub(paste0("^.*_"), "", fam_list[[i]]), always_add)) == "g"){
      
      cov <- cov[-which(str_detect(colnames(cov), "^g")),-which(str_detect(colnames(cov), "^g"))]
      
    }else if(setdiff(c("g","o"), intersect(gsub(paste0("^.*_"), "", fam_list[[i]]), c("g","o"))) == "o"){
      
      cov <- cov[-which(str_detect(colnames(cov), "^o")),-which(str_detect(colnames(cov), "^o"))]
    }
    
    # Adding the individuals present in always_add
    thr <- threshs %>% filter(!!as.symbol(fam_id) %in% fam_list[[i]][gsub("^.*_", "", fam_list[[i]]) %in% always_add])
    
    if (nrow(thr) == 2) {
      fam_threshs <- bind_rows(thr, fam_threshs)
    } else if (nrow(thr) == 0) {
      fam_threshs <- bind_rows(tibble::tibble(!!as.symbol(fam_id) := c("g","o"), lower = c(-Inf,-Inf), upper = c(Inf, Inf)), 
                               fam_threshs)
    } else if (any(str_ends(pull(thr,!!as.symbol(fam_id)) %>% str_subset(fam, .),"_g"))) {
      
      fam_threshs <- bind_rows(tibble::add_row(thr, !!as.symbol(fam_id) := "o", lower = -Inf, upper = Inf),
                               fam_threshs)
    } else if (any(str_ends(pull(thr,!!as.symbol(fam_id)) %>% str_subset(fam, .),"_o"))){
      
      fam_threshs <- bind_rows(tibble::add_row(thr, !!as.symbol(fam_id) := "g", lower = -Inf, upper = Inf, .before = 1),
                               fam_threshs)
    }
    
    # Setting the variables needed for Gibbs sampler
    fixed <- (pull(fam_threshs,upper) - pull(fam_threshs,lower)) < 1e-04
    std_err <- rep(Inf, length(out))
    names(std_err) <- c("genetic", "full")[out]
    n_gibbs <- 1
    
    # Running Gibbs sampler
    while(any(std_err > tol)){
      
      if(n_gibbs == 1){
        
        est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, covmat = cov, lower = pull(fam_threshs, lower), upper = pull(fam_threshs, upper),
                                    fixed = fixed, out = out, burn_in = 1000) %>% 
          `colnames<-`(c("genetic", "full")[out]) %>% 
          tibble::as_tibble()
        
      }else{
        
        est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, covmat = cov, lower = pull(fam_threshs, lower), upper = pull(fam_threshs, upper),
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
    return(stats::setNames(c(t(batchmeans::bmmat(est_liabs))), paste0(rep(c("Posterior_genetic", "Posterior_full")[out], each = 2), "_", c("liab", "std_err"))))
    
    # If progress = TRUE, a progress bar will be displayed
    # if(progress){
    #   j <- j+1
    #   utils::setTxtProgressBar(pb, j)
    # }
    
    
  }, future.seed = TRUE) %>% 
    do.call("bind_rows", .)
    
  # # Close the connection 
  # if(progress) close(pb)  
  
  # Finally, we can add all estimated liabilities as well
  # as their estimated standard errors to the tibble holding
  # the family information
  family <- bind_cols(family, gibbs_res)

  return(family)
}



#' Estimating the genetic or full liability for a variable number of
#' phenotypes
#' FUNCTIONS FOR MULTIPLE TRAITS IS STILL NOT COMPLETE. PLEASE DO NOT USE.
#' 
#' \code{estimate_liability} estimates the genetic component of the full
#' liability and/or the full liability for a number of individuals based
#' on their family history for one or more phenotypes.  It is a wrapper around 
#' \code{\link{estimate_liability_single}} and \code{\link{estimate_liability_multi}}.
#'
#' This function can be used to estimate either the genetic component of the 
#' full liability, the full liability or both for a variable number of traits.
#'
#' @param family A matrix, list or data frame that can be converted into a tibble.
#' Must have at least two columns that hold the family identifier and the corresponding
#' personal identifiers, respectively. That is, for each family in fam_id there should
#' be a list holding all individuals belonging to that family in pid. Note that the 
#' personal identifiers for all individuals must have a special format. It must be end
#' with _?, where ? is one of the following abbreviations 
#' - g (Genetic component of the full liability)
#' - o (full liability)
#' - m (Mother)
#' - f (Father)
#' -c\[0-9\]* (Children)
#' - mgm (Maternal grandmother)
#' - mgf (Maternal grandfather)
#' - pgm (Paternal grandmother)
#' - pgf (Paternal grandfather)
#' - s\[0-9\]* (Full siblings)
#' - mhs\[0-9\]* (Half-siblings - maternal side)
#' - phs\[0-9\]* (Half-siblings - paternal side)
#' - mau\[0-9\]* (Aunts/Uncles - maternal side)
#' - pau\[0-9\]* (Aunts/Uncles - paternal side).
#' See also \code{\link{construct_covmat}}.
#' @param threshs A matrix, list or data frame that can be converted into a tibble.
#' Must have at least five columns; one holding the personal identifier for all individuals,
#' and the remaining four holding the lower and upper thresholds for the first and second
#' phenotype, respectively. It must be possible to tie each pair of lower and upper thresholds
#' to a specific phenotype uniquely. This is done easily by adding _{name_of_phenotype} to
#' the column names lower and upper, e.g. lower_p1 and upper_p1 for the lower and upper
#' thresholds corresponding to the first phenotype. 
#' @param h2 Either a number representing the heritability on liability scale for a 
#' single phenotype, or a numeric vector representing the liability-scale heritabilities
#' for all phenotypes. All entries in h2 must be non-negative and at most 1.
#' @param  pid A string holding the name of the column in \code{family} and 
#' \code{threshs} that hold the personal identifier(s). Defaults to \code{"PID"}.
#' @param fam_id A string holding the name of the column in \code{family} that
#' holds the family identifier. Defaults to \code{"fam_ID"}.
#' @param out A character or numeric vector indicating whether the genetic component
#' of the full liability, the full liability or both should be returned. If \code{out = c(1)} or 
#' \code{out = c("genetic")}, the genetic liability is estimated and returned. If \code{out = c(2)} or 
#' \code{out = c("full")}, the full liability is estimated and returned. If code{out = c(1,2)} or 
#' \code{out = c("genetic", "full")}, both components are estimated and returned. 
#' Defaults to \code{c(1)}.
#' @param tol A number that is used as the convergence criterion for the Gibbs sampler.
#' Equals the standard error of the mean. That is, a tolerance of 0.2 means that the 
#' standard error of the mean is below 0.2. Defaults to 0.01.
#' @param always_add A character vector or \code{NULL}. If \code{always_add = c("g","o")}, both the genetic component 
#' of the full liability as well as the full liability will be added to the list of family members. 
#' If always_add equals \code{"g"} or \code{"o"}, the genetic component of the full liability or the full liability
#' will be added, respectively. If \code{always_add = NULL}, no component will be added.
#' Defaults to \code{c("g","o")}.
#' @param progress A logical scalar indicating whether the function should display
#' a progress bar. Defaults to \code{FALSE}.
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
estimate_liability <- function(family, threshs, h2 = 0.5, pid = "PID", fam_id = "fam_ID", 
                               out = c(1), tol = 0.01, always_add = c("g","o"), progress = FALSE,
                               genetic_corrmat = NULL, full_corrmat = NULL, phen_names = NULL){
  
  if(length(h2) == 1){
    
    return(estimate_liability_single(family = family, threshs=threshs, h2 = h2, pid = pid, fam_id = fam_id, 
                              out = out, tol = tol, always_add = always_add, progress = progress))
    
  }else{
    
    return(estimate_liability_multi(family = family, threshs = threshs, h2_vec = h2, 
                                    genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat,
                                    phen_names = phen_names, always_add = always_add, pid = pid,
                                    fam_id = fam_id, out = out, tol = tol, progress = progress))
  } 
}
