#' Estimating the genetic or full liability 
#'
#' \code{estimate_liability} estimate the genetic component of the full
#' liability and/or the full liability for a number of individuals based
#' on their family history.
#'
#' This function can be used to estimate either the genetic component of the 
#' full liability, the full liability or both. 
#'
#' @param family A matrix, list or data frame that can be converted into a tibble.
#' Must have at least two columns that hold the family identifier and the corresponding
#' personal identifiers, respectively. That is, for each family in fam_id there should
#' be a list holding all individuals belonging to that family in pid.
#' @param threshs A matrix, list or data frame that can be converted into a tibble.
#' Must have at least three columns, one holding the personal identifier for all individuals,
#' and the remaining two holding the lower and upper thresholds, respectively.
#' @param sq.herit A number representing the squared heritability on liability scale
#' for a single phenotype. Must be non-negative and at most 1.
#' Defaults to 0.5.
#' @param  pid A string holding the name of the column in \code{family} and 
#' \code{threshs} that hold the personal identifier(s). Defaults to "PID".
#' @param fam_id A string holding the name of the column in \code{family} that
#' holds the family identifier.
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
#' (see \code{\link{future.apply::future_apply()}}).
#' @param always_add A character vector or NULL. If always_ad = c("g","o"), both the genetic component 
#' of the full liability as well as the full liability will be added to the list of family members. 
#' If always_add equals "g" or "o", the genetic component of the full liability or the full liability
#' will be added, respectively. If always_add = NULL, no component will be added.
#' 
#' @return If family and threshs are two matrices, lists or data frames that can be converted into
#' tibbles, if family has two columns named like the strings represented in pid and fam_id, if 
#' threshs has a column named like the string given in pid as well as a colimn named "lower" and 
#' a column named "upper" and if the squared heritability, out, tol and always_add are of the required form,
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
#' 
#' estimate_liability()
#' 
#' @seealso \code{\link{future.apply::future_apply()}}
#' 
#' @export
estimate_liability <- function(family, threshs, sq.herit = 0.5, pid = "PID", fam_id = "fam_ID", out = c(1), tol = 0.01, parallel = TRUE, always_add = c("g","o")){
  
  # Turning parallel into class logical
  parallel <- as.logical(parallel)
  
  # Turning pid and fam_id into strings
  pid <- as.character(pid)
  fam_id <- as.character(fam_id)
  
  #Turning family and threshs into tibbles
  family <- as_tibble(family)
  threshs <- as_tibble(threshs)
  
  # Checking that the heritability is valid
  if(sq.herit<0)stop("The squared heritability must be non-negative")
  if(sq.herit>1)stop("The squared heritability must be smaller than or equal to 1.")
  # Checking that family has two columns named pid_col and fam_id
  if(!(pid %in% colnames(family))) stop(paste0("The colum ", pid," does not exist in the tibble family..."))
  if(!(fam_id %in% colnames(family))) stop(paste0("The colum ", fam_id," does not exist in the tibble family..."))
  # And that pid is also present in the tibble threshs
  if(!(pid %in% colnames(threshs))) stop(paste0("The colum ", pid," does not exist in the tibble threshs..."))
  # In addition, we check that threshs has two columns named lower and upper
  if(any(!c("lower","upper") %in% colnames(threshs))) stop("The tibble threshs must include two columns named 'lower' and 'upper'!")
  # Checking that tol is valid
  if(class(tol) != "numeric") stop("The tolerance must be numeric!")
  if(tol <= 0) stop("The tolerance must be strictly positive!")
  # Checking that always_add is a vector of strings
  if(class(always_add) != "character" & class(always_add) != "NULL" ) stop("always_add must be of class character!") 
  always_add <- intersect(always_add, c("g","o"))
  # Checking that out is either a character vector or a
  # numeric vector 
  if(class(out) == "numeric"){
    
    out <- intersect(out, c(1,2))
    
  }else if(class(out) == "character"){
    
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
  threshs <- select(threshs, !!as.symbol(pid), lower, upper)
  
  # If parallel = T, future_lapply needs to be installed 
  if(parallel){
    if(!("future.apply" %in% library()$results[,1])){
      
      count <- 1
      ver <- "a"
      while(!(ver %in% c("y","n"))){
        
        if(count < 4){
          
          count <- count + 1 
          ver <- readline(prompt="In order to use a parallelized version of this function, the package future.apply must be installed. \n 
                          Do you want to install future.apply now [y/n]?: ")
        }else{
          
          stop("Function aborted...")
        }
      }
      
      if(ver == "y"){
        
        install.packages("future.apply")
      }else{
        
        parallel <- FALSE
      }
    }
  }
  
  # Extracting the families
  fam_list <- pull(family, !!as.symbol(pid))
  
  if(parallel){
    
    cat(paste0("The number of workers is ", nbrOfWorkers()))
    
    gibbs_res <- future.apply::future_sapply(X= 1:nrow(family), FUN = function(i){
      # Extract family members
      fam <- unlist(fam_list[i])
      # Remove individual o and/or g from the set (if present)
      full_fam <- setdiff(gsub(paste0("^[0-9]*_"), "", fam), c("g","o"))
      # And add the individuals that should be added
      full_fam <- c(always_add, full_fam)
      # Constructing the covariance matrix
      cov <- construct_covmat(fam_vec = full_fam, n_fam = NULL, add_ind = FALSE, sq.herit = sq.herit)
      
      # Extracting the thresholds for all family members 
      fam_threshs = threshs[match(fam[!(str_detect(fam, "^[0-9]*_g") | str_detect(fam, "^[0-9]*_o"))], pull(threshs,!!as.symbol(pid))), ]
      # Adding the individuals present in always_add
      thr <- threshs[match(fam[str_detect(fam, paste0("^[0-9]*_[", paste(always_add, collapse = "") , "]"))], pull(threshs,!!as.symbol(pid))), ]
      
      if(nrow(thr) == 2){
        
        fam_threshs <- bind_rows(thr, fam_threshs)
      }else if(nrow(thr) == 0){
        
        fam_threshs <- bind_rows(tibble(!!as.symbol(pid) := c("g","o"), lower = c(-Inf,-Inf), upper = c(Inf, Inf)), 
                                 fam_threshs)
      }else if(str_detect(pull(thr,!!as.symbol(pid)),"^[0-9]*_g")){
        
        fam_threshs <- bind_rows(add_row(thr, !!as.symbol(pid) := "o", lower = -Inf, upper =Inf),
                                 fam_threshs)
      }else if(str_detect(pull(thr,!!as.symbol(pid)),"^[0-9]*_o")){
        
        fam_threshs <- bind_rows(add_row(thr, !!as.symbol(pid) := "g", lower = -Inf, upper =Inf, .before = 1),
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
          
          est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, sigma = cov, lower = pull(fam_threshs, lower), upper = pull(fam_threshs, upper),
                                      fixed = fixed, ind = out, burn_in = 1000) %>% 
            `colnames<-`(c("genetic", "full")[out]) %>% 
            as_tibble()
          
        }else{
          
          est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, sigma = cov, lower = pull(fam_threshs, lower), upper = pull(fam_threshs, upper),
                                      fixed = fixed, ind = out, burn_in = 1000) %>%
            `colnames<-`(c("genetic", "full")[out]) %>%
            as_tibble() %>% 
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
      return(setNames(c(t(batchmeans::bmmat(est_liabs))), paste0(rep(c("Posterior_genetic", "Posterior_full")[out], each = 2), ".", c("liab", "std_err"))))
      
    }, future.seed = TRUE)
    
    
  }else{
    
    gibbs_res <- sapply(X = 1:nrow(family), FUN = function(i){
      
      # Extract family members
      fam <- fam_list[[i]]
      # Remove individual o and/or g from the set (if present)
      full_fam <- setdiff(gsub(paste0("^[0-9]*_"), "", fam), c("g","o"))
      # And add the individuals that should be added
      full_fam <- c(always_add, full_fam)
      # Constructing the covariance matrix
      cov <- construct_covmat(fam_vec = full_fam, n_fam = NULL, add_ind = FALSE, sq.herit = sq.herit)
      
      # Extracting the thresholds for all family members 
      fam_threshs = threshs[match(fam[!(str_detect(fam, "^[0-9]*_g") | str_detect(fam, "^[0-9]*_o"))], pull(threshs,!!as.symbol(pid))), ]
      # Adding the individuals present in always_add
      thr <- threshs[match(fam[str_detect(fam, paste0("^[0-9]*_[", paste(always_add, collapse = "") , "]"))], pull(threshs,!!as.symbol(pid))), ]
      
      if(nrow(thr) == 2){
        
        fam_threshs <- bind_rows(thr, fam_threshs)
      }else if(nrow(thr) == 0){
        
        fam_threshs <- bind_rows(tibble(!!as.symbol(pid) := c("g","o"), lower = c(-Inf,-Inf), upper = c(Inf, Inf)), 
                                 fam_threshs)
      }else if(str_detect(pull(thr,!!as.symbol(pid)),"^[0-9]*_g")){
        
        fam_threshs <- bind_rows(add_row(thr, !!as.symbol(pid) := "o", lower = -Inf, upper =Inf),
                                 fam_threshs)
      }else if(str_detect(pull(thr,!!as.symbol(pid)),"^[0-9]*_o")){
        
        fam_threshs <- bind_rows(add_row(thr, !!as.symbol(pid) := "g", lower = -Inf, upper =Inf, .before = 1),
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
          
          est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, sigma = cov, lower = pull(fam_threshs, lower), upper = pull(fam_threshs, upper),
                                      fixed = fixed, ind = out, burn_in = 1000) %>% 
            `colnames<-`(c("genetic", "full")[out]) %>% 
            as_tibble()
          
        }else{
          
          est_liabs <- rtmvnorm.gibbs(n_sim = 1e+05, sigma = cov, lower = pull(fam_threshs, lower), upper = pull(fam_threshs, upper),
                                      fixed = fixed, ind = out, burn_in = 1000) %>%
            `colnames<-`(c("genetic", "full")[out]) %>%
            as_tibble() %>% 
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
      return(setNames(c(t(batchmeans::bmmat(est_liabs))), paste0(rep(c("Posterior_genetic", "Posterior_full")[out], each = 2), "_", c("liab", "std_err"))))
    })
  }
  
  # Finally, we can add all estimated liabilities as well
  # as their estimated standard errors to the tibble holding
  # the family information
  family <- bind_cols(family, as_tibble(t(gibbs_res)))

  return(family)
}