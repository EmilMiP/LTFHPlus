utils::globalVariables("lower")
utils::globalVariables("upper")

#' Estimating the genetic or full liability 
#'
#' \code{estimate_liability} estimates the genetic component of the full
#' liability and/or the full liability for a number of individuals based
#' on their family history.
#'
#' This function can be used to estimate either the genetic component of the 
#' full liability, the full liability or both. 
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
#' and the remaining two holding the lower and upper thresholds, respectively.
#' @param h2 A number representing the heritability on liability scale
#' for a single phenotype. Must be non-negative. Note that under the liability threshold model,
#' the heritability must also be at most 1.
#' Defaults to 0.5.
#' @param  pid A string holding the name of the column in \code{family} and 
#' \code{threshs} that hold the personal identifier(s). Defaults to "PID".
#' @param fam_id A string holding the name of the column in \code{family} that
#' holds the family identifier. Defaults to "fam_ID".
#' @param out A character or numeric vector indicating whether the genetic component
#' of the full liability, the full liability or both should be returned. If out = c(1) or 
#' out = c("genetic"), the genetic liability is estimated and returned. If out = c(2) or 
#' out = c("full"), the full liability is estimated and returned. If out = c(1,2) or 
#' out = c("genetic", "full"), both components are estimated and returned. 
#' Defaults to c(1).
#' @param tol A number that is used as the convergence criterion for the Gibbs sampler.
#' Equals the standard error of the mean. That is, a tolerance of 0.2 means that the 
#' standard error of the mean is below 0.2. Defaults to 0.01.
#' @param always_add A character vector or NULL. If always_add = c("g","o"), both the genetic component 
#' of the full liability as well as the full liability will be added to the list of family members. 
#' If always_add equals "g" or "o", the genetic component of the full liability or the full liability
#' will be added, respectively. If always_add = NULL, no component will be added.
#' Defaults to c("g","o").
#' 
#' @return If family and threshs are two matrices, lists or data frames that can be converted into
#' tibbles, if family has two columns named like the strings represented in pid and fam_id, if 
#' threshs has a column named like the string given in pid as well as a column named "lower" and 
#' a column named "upper" and if the heritability, out, tol and always_add are of the required form,
#' then the function returns a tibble with either four or six columns (depending on the length of out).
#' The first two columns correspond to the columns fam_id and pid from family. 
#' If out is equal to c(1) or c("genetic"), the third and fourth column hold the estimated genetic 
#' liability as well as the corresponding standard error, respectively. 
#' If out equals c(2) or c("full"), the third and fourth column hold the estimated full liability 
#' as well as the corresponding standard error, respectively. 
#' If out is equal to c(1,2) or c("genetic","full"), the third and fourth column hold the estimated 
#' genetic liability as well as the corresponding standard error, respectively, while the fifth and
#' sixth column hold the estimated full liability as well as the corresponding standard error, respectively.
#' 
#' @examples
#' sims <- simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, 
#' add_ind = TRUE, h2 = 0.5, n_sim=10, pop_prev = .05)
#' 
#' estimate_liability(family = sims$fam_ID, threshs = sims$thresholds, 
#' h2 = 0.5, pid = "indiv_ID", fam_id = "fam_ID", out = c(1), tol = 0.01,
#' always_add = c("g","o"))
#' # 
#' sims <- simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, 
#' h2 = 0.5, n_sim=10, pop_prev = .05)
#' estimate_liability(family = sims$fam_ID, threshs = sims$thresholds, 
#' h2 = 0.5, pid = "indiv_ID", fam_id = "fam_ID", out = c("genetic"), 
#' tol = 0.01, always_add = c("g","o"))
#' 
#' @seealso \code{\link[future.apply]{future_apply}}
#' 
#' @importFrom dplyr %>% pull bind_rows bind_cols select
#' @importFrom rlang :=
#' 
#' @export
estimate_liability <- function(family, threshs, h2 = 0.5, pid = "PID", fam_id = "fam_ID", out = c(1), 
                               tol = 0.01, always_add = c("g","o")){
  

# Making sure input is valid ----------------------------------------------

  # Turning pid and fam_id into strings
  pid <- as.character(pid)
  fam_id <- as.character(fam_id)
  
  #Turning family and threshs into tibbles
  # if they are not of class tbl
  if(!tibble::is_tibble(family))  family <- tibble::as_tibble(family)
  if(!tibble::is_tibble(threshs)) threshs <- tibble::as_tibble(threshs)
  
  
  # Checking that the heritability is valid
  if(!is.numeric(h2) && !is.integer(h2) )stop("The heritability must be numeric!")
  if(h2<0)stop("The heritability must be non-negative!")
  if(h2>1)stop("Under the liability threshold model, the heritability must be smaller than or equal to 1!")
  
  # Checking that family has two columns named pid_col and fam_id
  if(!(pid %in% colnames(family))) stop(paste0("The column ", pid," does not exist in the tibble family..."))
  if(!(fam_id %in% colnames(family))) stop(paste0("The column ", fam_id," does not exist in the tibble family..."))
  
  # And that pid is also present in the tibble threshs
  if(!(fam_id %in% colnames(threshs))) stop(paste0("The column ", fam_id," does not exist in the tibble threshs..."))
  
  # In addition, we check that threshs has two columns named lower and upper
  if(any(!c("lower","upper") %in% colnames(threshs))) stop("The tibble threshs must include two columns named 'lower' and 'upper'!")
  
  # Checking that tol is valid
  if(!is.numeric(tol)) stop("The tolerance must be numeric!")
  if(tol <= 0) stop("The tolerance must be strictly positive!")
  
  # Checking that always_add is a vector of strings
  if(!is.character(always_add) & !is.null(always_add)) stop("always_add must be of class character!") 
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
  threshs <- select(threshs, !!as.symbol(fam_id), lower, upper)
  
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
  fam_list <- pull(family, !!as.symbol(fam_id))
  
  cat(paste0("The number of workers is ", future::nbrOfWorkers(), "\n"))
    
  gibbs_res <- future.apply::future_lapply(X= 1:nrow(family), FUN = function(i){
    # Extract family members
    fam <- unlist(fam_list[i])
    # Remove individual o and/or g from the set (if present)
    full_fam <- setdiff(gsub(paste0("^.*_"), "", fam), c("g","o"))
    # Constructing the covariance matrix
    # If always_add holds "g" or "o", add_ind must be TRUE
    cov <- construct_covmat(fam_vec = full_fam, n_fam = NULL, add_ind = length(always_add), h2 = h2)
    
    # Extracting the thresholds for all family members 
    fam_threshs = threshs[match(fam[!(str_detect(fam, "^.*_g$") | str_detect(fam, "^.*_o$"))], pull(threshs,!!as.symbol(fam_id))), ]
    # Adding the individuals present in always_add
    thr <- threshs[match(fam[str_detect(fam, paste0("^.*_[", paste(always_add, collapse = "") , "]"))], pull(threshs,!!as.symbol(fam_id))), ]
    
    if(nrow(thr) == 2){
      
      fam_threshs <- bind_rows(thr, fam_threshs)
    }else if(nrow(thr) == 0){
      
      fam_threshs <- bind_rows(tibble::tibble(!!as.symbol(fam_id) := c("g","o"), lower = c(-Inf,-Inf), upper = c(Inf, Inf)), 
                               fam_threshs)
    }else if(str_detect(pull(thr,!!as.symbol(fam_id)),"^.*_g$")){
      
      fam_threshs <- bind_rows(tibble::add_row(thr, !!as.symbol(fam_id) := "o", lower = -Inf, upper =Inf),
                               fam_threshs)
    }else if(str_detect(pull(thr,!!as.symbol(fam_id)),"^.*_o$")){
      
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
    
  }, future.seed = TRUE) %>% 
    do.call("bind_rows", .)
    
    
  # Finally, we can add all estimated liabilities as well
  # as their estimated standard errors to the tibble holding
  # the family information
  family <- bind_cols(family, gibbs_res)

  return(family)
}