#' Estimating the genetic or full liability for multiple phenotypes
#' using prevalence information
#' 
#' FUNCTIONS FOR MULTIPLE TRAITS IS STILL NOT COMPLETE. PLEASE DO NOT USE.
#'
#' \code{estimate_liability_prevalence} estimates the genetic component of 
#' the full liability and/or the full liability for a number of individuals 
#' based solely on prevalence information for a variable number of 
#' phenotypes.
#'
#' This function can be used to estimate either the genetic component of the 
#' full liability, the full liability or both for a variable number of traits,
#' when no family history is available, but prevalence information for each 
#' phenotype can be obtained.
#'
#' @param status A matrix, list or data frame that can be converted into a tibble.
#' Must have at least three columns; one holding the personal identifier for all individuals,
#' and the remaining holding the phenotype status for the first and second
#' phenotype, respectively. It must be possible to tie each status variable
#' to a specific phenotype uniquely. The function will use the column names to create
#' phenotype names. 
#' @param h2_vec A numeric vector representing the heritabilites on liability scale
#' for all phenotypes. All entries in \code{h2_vec} must be non-negative and at most 1.
#' @param genetic_corrmat A numeric matrix holding the genetic correlations between the desired 
#' phenotypes. All diagonal entries must be equal to one, while all off-diagonal entries 
#' must be between -1 and 1. In addition, the matrix must be symmetric.
#' Defaults to NULL.
#' @param full_corrmat A  numeric matrix holding the full correlations between the desired 
#' phenotypes. All diagonal entries must be equal to one, while all off-diagonal entries 
#' must be between -1 and 1. In addition, the matrix must be symmetric.
#' Defaults to NULL.
#' @param prevalences A numeric, non-negative vector holding the prevalences. 
#' All prevalences must be at most one. 
#' @param  pid A string holding the name of the column in \code{status} that hold 
#' the personal identifier. Defaults to "PID".
#' @param out A character or numeric vector indicating whether the genetic component
#' of the full liability, the full liability or both should be returned. If out = c(1) or 
#' out = c("genetic"), the genetic liability is estimated and returned. If out = c(2) or 
#' out = c("full"), the full liability is estimated and returned. If out = c(1,2) or 
#' out = c("genetic", "full"), both components are estimated and returned. 
#' Defaults to c(1).
#' @param tol A number that is used as the convergence criterion for the Gibbs sampler.
#' Equals the standard error of the mean. That is, a tolerance of 0.2 means that the 
#' standard error of the mean is below 0.2. Defaults to 0.01.
#' @param parallel A logical scalar indicating whether computations should be performed parallel.
#' In order for this to be possible, the user must install the library "future.apply" and create a plan
#' (see \code{\link[future.apply]{future_apply}}). Defaults to FALSE.
#' @param progress A logical scalar indicating whether the function should display
#' a progress bar. Defaults to \code{FALSE}.
#' 
#' 
#' @return If \code{status} is a matrix, list or data frame that can be converted into a
#' tibble and that has a column named \code{PID} and if the heritabilities, corrmat, 
#' out and tol are of the required form, then the function returns a tibble with at least
#' five columns (depending on the length of out). 
#' The first column corresponds to the columns pid. 
#' If out is equal to c(1) or c("genetic"), the second and third columns hold the estimated genetic 
#' liability as well as the corresponding standard error for the first phenotype, respectively. 
#' If out equals c(2) or c("full"), the second and third columns hold the estimated full liability 
#' as well as the corresponding standard error for the first phenotype, respectively. 
#' If out is equal to c(1,2) or c("genetic","full"), the second and third columns hold the estimated 
#' genetic liability as well as the corresponding standard error for the first phenotype, while the fourth and
#' fifth columns hold the estimated full liability as well as the corresponding standard error for the
#' same phenotype. 
#' The remaining columns hold the estimated genetic liabilities and/or the estimated full liabilities
#' as well as the corresponding standard errors for the remaining phenotypes.
#' 
#' @seealso \code{\link[future.apply]{future_apply}}
#' 
#' @importFrom dplyr %>% pull bind_rows bind_cols select row_number rename
#' @importFrom rlang :=
#' 
#' @noRd
estimate_liability_prevalence = function(status, h2_vec, genetic_corrmat, full_corrmat,
                                         prevalences, pid = "PID", out = c(1), tol = 0.01, 
                                         parallel = FALSE, progress = FALSE){
  # Turning parallel into class logical
  parallel <- as.logical(parallel)

  # Turning pid into a string
  pid <- as.character(pid)
  
  # Turning threshs into a 
  if(!tibble::is_tibble(status)) status <- tibble::as_tibble(status)
  
  # Checking that the heritabilities are valid
  if(validate_proportion(h2_vec)){invisible()}

  # Checking that all correlations are valid
  if(validate_correlation_matrix(genetic_corrmat)){invisible()}
  if(validate_correlation_matrix(full_corrmat)){invisible()}
  # Checking that all prevalences are valid
  if(validate_proportion(prevalences)){invisible()}

  # And that pid is also present in the tibble threshs
  if(!(pid %in% colnames(status))) stop(paste0("The column ", pid," does not exist in the tibble status"))
  
  # Checking that tol is valid
  if(!is.numeric(tol)) stop("The tolerance must be numeric!")
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
    
    warning("out is not of the required format! \n The function will return the estimated genetic liability!")
    out <- c(1)
  }
  
  # Sorting out
  out <- sort(out)

  
  # Now we can extract the number of phenotypes
  n_pheno <- length(h2_vec)
  
  if(ncol(status) != (n_pheno + 1)) stop("Something is wrong with the number of phenotypes... \n 
The number of columns in status is not equal to the number of phenotypes specified in h2_vec...\
Does all columns have the required names?")
  
  # As well as the phenotype names
  pheno_names <- colnames(status)[-1]
  
  
  # If parallel = TRUE, future_lapply needs to be installed 
  if(parallel){
    if(!("future.apply" %in% library()$results[,1])){
      stop("In order to use a parallelized version of this function, the package future.apply must be installed.")
    }
  }

  
  if(parallel){
    
    message(paste0("The number of workers is ", future::nbrOfWorkers(), "\n"))
    
    # As all families have the same structure (each family solely includes
    # the individual itself), we can use the same covariance matrix for 
    # all families
    # Constructing the covariance matrix
    cov <- construct_covmat(fam_vec = c(), n_fam = NULL, add_ind = TRUE, 
                            genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat,
                            h2 = h2_vec, phen_names = pheno_names)
    
    gibbs_res <- future.apply::future_sapply(X= 1:nrow(status), FUN = function(i){
      
      cur_config <- status[i,]
      
      # Setting the variables needed for Gibbs sampler
      lower <- rep(-Inf, 2*n_pheno)
      upper <- rep(Inf, 2*n_pheno)
      fixed <- rep(FALSE, 2*n_pheno)
      std_err <- rep(Inf, length(out))
      names(std_err) <- c("genetic", "full")[out]
      n_gibbs <- 1
      
      # Running Gibbs sampler
      while(any(std_err > tol)){
        
        if(n_gibbs == 1){
          
          est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, covmat = cov, lower = lower, upper = upper,
                                      fixed = fixed, out = out, burn_in = 1000) %>% 
            `colnames<-`(c("genetic", "full")[out]) %>% 
            tibble::as_tibble()
          
        }else{
          
          est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, covmat = cov, lower =lower, upper = upper,
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
      return(stats::setNames(c(t(batchmeans::bmmat(est_liabs))), paste0(rep(c("Posterior_genetic", "Posterior_full")[out], each = 2), ".", c("liab", "std_err"))))
      
      
    }, future.seed = TRUE)
    
    
  }else{
    
    # As all families have the same structure (each family solely includes
    # the individual itself), we can use the same covariance matrix for 
    # all families
    # Constructing the covariance matrix
    cov <- construct_covmat(fam_vec = c(), n_fam = NULL, add_ind = TRUE, 
                            genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat,
                            h2 = h2_vec, phen_names = pheno_names)
    
    gibbs_res <- sapply(X = 1:nrow(family), FUN = function(i){
      
      
      cur_config <- status[i,]
      
      # Setting the variables needed for Gibbs sampler
      lower <- rep(-Inf, 2*n_pheno)
      upper <- rep(Inf, 2*n_pheno)
      fixed <- rep(FALSE, 2*n_pheno)
      std_err <- rep(Inf, length(out))
      names(std_err) <- c("genetic", "full")[out]
      n_gibbs <- 1
      
      # Running Gibbs sampler
      while(any(std_err > tol)){
        
        if(n_gibbs == 1){
          
          est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, covmat = cov, lower = lower, upper = upper,
                                      fixed = fixed, out = out, burn_in = 1000) %>% 
            `colnames<-`(c("genetic", "full")[out]) %>% 
            tibble::as_tibble()
          
        }else{
          
          est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, covmat = cov, lower =lower, upper = upper,
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
      
    })
  }
  

  # Finally, we can add all estimated liabilities as well
  # as their estimated standard errors to the tibble holding
  # the family information
  family <- bind_cols(status, tibble::as_tibble(t(gibbs_res)))
  
  return(family)
}
