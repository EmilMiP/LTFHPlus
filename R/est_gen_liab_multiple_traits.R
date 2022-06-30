#' Estimating the genetic or full liability for multiple phenotypes
#' FUNCTIONS FOR MULTIPLE TRAITS IS STILL NOT COMPLETE. PLEASE DO NOT USE.
#' 
#' \code{estimate_liability_multi} estimates the genetic component of the full
#' liability and/or the full liability for a number of individuals based
#' on their family history for a variable number of phenotypes.
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
#' @seealso \code{\link[future.apply]{future_apply}}
#' 
#' @importFrom dplyr %>% pull bind_rows bind_cols select row_number rename
#' @importFrom rlang :=
#' 
#' @export
estimate_liability_multi <- function(family, threshs, h2_vec, genetic_corrmat, full_corrmat,
                                     phen_names = NULL, always_add, pid = "PID",
                                     fam_id = "fam_ID", out = c(1), tol = 0.01, 
                                     progress = FALSE){
  
# Making sure input is valid--------------------------------------------
  
  # If always_add is a vector of length zero, it is set to
  # NULL instead
  if(length(always_add) == 0) always_add <- NULL
  
  # Turning parallel into class logical
  progress <- as.logical(progress)

  # Turning pid and fam_id into strings
  pid <- as.character(pid)
  fam_id <- as.character(fam_id)
  
  #Turning family and threshs into tibbles
  if(!tibble::is_tibble(family)) family <- tibble::as_tibble(family)
  if(!tibble::is_tibble(threshs)) threshs <- tibble::as_tibble(threshs)
  
  # Checking that the heritability is valid
  if(validate_proportion(h2_vec)){invisible()}
  
  # Checking that all correlations are valid
  if(validate_correlation_matrix(genetic_corrmat)){invisible()}
  if(validate_correlation_matrix(full_corrmat)){invisible()}
  
  # Checking that family has two columns named pid_col and fam_id
  if(!(pid %in% colnames(family))) stop(paste0("The column ", pid," does not exist in the tibble family..."))
  if(!(fam_id %in% colnames(family))) stop(paste0("The column ", fam_id," does not exist in the tibble family..."))
  
  # And that pid is also present in the tibble threshs
  if(!(pid %in% colnames(threshs))) stop(paste0("The column ", pid," does not exist in the tibble threshs..."))
  
  # In addition, we check that threshs has columns named lower and upper
  if(any(!c("lower","upper") %in% sub("_.*$","",colnames(threshs)))) stop("The tibble threshs must include two columns named 'lower' and 'upper'!")
  
  # Checking that tol is valid
  if(!is.numeric(tol)) stop("The tolerance must be numeric!")
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
  
  # Now we can extract the number of phenotypes
  n_pheno <- length(h2_vec)
  if(ncol(threshs) != (2*n_pheno + 1)) stop("Something is wrong with the number of phenotypes... \n 
The number of pairs of lower and upper thresholds is not equal to the number of phenotypes specified in h2_vec...\
Does all columns have the required names?")
  
  # As well as the phenotype names
  pheno_names <- sub("lower_", "", colnames(threshs)[stringr::str_detect(colnames(threshs), "^lower_")], ignore.case = TRUE)
  
  # Finally, we also check whether all lower thresholds are 
  # smaller than or equal to the upper thresholds
  for ( pheno in pheno_names ) {
  
    
    if(any(pull(threshs, !!as.symbol(paste0("lower_", pheno))) > pull(threshs, !!as.symbol(paste0("upper_", pheno))))){
      cat("Warning message: \n Some lower thresholds are larger than the corresponding upper thresholds! \n
The lower and upper thresholds will be swapped...")
      
      swapping_indx <- which(pull(threshs, !!as.symbol(paste0("lower_", pheno))) > pull(threshs, !!as.symbol(paste0("upper_", pheno))))
      
      threshs <- mutate(threshs, !!as.symbol(paste0("lower_", pheno)) := ifelse(row_number() %in% swapping_indx, !!as.symbol(paste0("lower_", pheno)) + !!as.symbol(paste0("upper_", pheno)), !!as.symbol(paste0("lower_", pheno)))) %>% 
        mutate(., !!as.symbol(paste0("upper_", pheno)) := ifelse(row_number() %in% swapping_indx, !!as.symbol(paste0("lower_", pheno)) - !!as.symbol(paste0("upper_", pheno)), !!as.symbol(paste0("upper_", pheno)))) %>% 
        mutate(., !!as.symbol(paste0("lower_", pheno)) := ifelse(row_number() %in% swapping_indx, !!as.symbol(paste0("lower_", pheno)) - !!as.symbol(paste0("upper_", pheno)), !!as.symbol(paste0("lower_", pheno))))
    }
  }
  
  # Extracting the families
  fam_list <- pull(family, !!as.symbol(fam_id))
  
  cat(paste0("The number of workers is ", future::nbrOfWorkers(), "\n"))
  
  # If progress = TRUE, a progress bar will be displayed
  if(progress){
    pb <- utils::txtProgressBar(min = 0, max = nrow(family), style = 3, char = "=")
    j <- 0
  }
  
  gibbs_res <- future.apply::future_sapply(X= 1:nrow(family), FUN = function(i){
    
    # Extract family members
    fam <- unlist(fam_list[i])

    # Constructing the covariance matrix
    cov <- construct_covmat(fam_vec = fam, n_fam = NULL, add_ind = length(intersect(gsub(paste0("^.*_"), "", fam), c("g","o"))), 
                            genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat,
                            h2 = h2_vec, phen_names = pheno_names)

    if(setdiff(c("g","o"), intersect(gsub(paste0("^.*_"), "", fam), c("g","o"))) == "g"){
      
      cov <- cov[-which(stringr::str_detect(colnames(cov), "^g_")),-which(stringr::str_detect(colnames(cov), "^g_"))]
    
    }else if(setdiff(c("g","o"), intersect(gsub(paste0("^.*_"), "", fam), c("g","o"))) == "o"){
        
      cov <- cov[-which(stringr::str_detect(colnames(cov), "^o_")),-which(stringr::str_detect(colnames(cov), "^o_"))]
    }
    
    # Extracting the thresholds for all family members 
    # and all phenotypes
    intermed_res = threshs[match(fam, pull(threshs,!!as.symbol(pid))), ]
    
    for(pheno in pheno_names){
      
      if(which(pheno_names == pheno) == 1){
        
        fam_threshs <- select(intermed_res, !!as.symbol(pid), !!as.symbol(paste0("lower_", pheno)), !!as.symbol(paste0("upper_", pheno))) %>% 
          rename(lower = !!as.symbol(paste0("lower_", pheno)), upper = !!as.symbol(paste0("upper_", pheno)))
      
      }else{
        
        fam_threshs <- select(intermed_res, !!as.symbol(pid), !!as.symbol(paste0("lower_", pheno)), !!as.symbol(paste0("upper_", pheno))) %>% 
          rename(lower = !!as.symbol(paste0("lower_", pheno)), upper = !!as.symbol(paste0("upper_", pheno))) %>% 
          bind_rows(fam_threshs,.)
      }

    }
    rm(intermed_res)
      
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
    
    # If progress = TRUE, a progress bar will be displayed
    if(progress){
      j <- j+1
      utils::setTxtProgressBar(pb, j)
    }
    
    # If all standard errors are below the tolerance, 
    # the estimated liabilities as well as the corresponding 
    # standard error can be returned
    return(stats::setNames(c(t(batchmeans::bmmat(est_liabs))), paste0(rep(c("Posterior_genetic", "Posterior_full")[out], each = 2), ".", c("liab", "std_err"))))
    
  }, future.seed = TRUE) %>% 
    do.call("bind_rows",.)
  
  # Close the connection 
  if(progress) close(pb)  
  
  # Finally, we can add all estimated liabilities as well
  # as their estimated standard errors to the tibble holding
  # the family information
  family <- bind_cols(family, gibbs_res)

  return(family)
}
