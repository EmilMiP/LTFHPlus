utils::globalVariables("cir")
utils::globalVariables("thresh")

#' Positive definite matrices
#'
#' \code{correct_positive_definite} verifies that a given covariance matrix
#' is indeed positive definite by checking that all eigenvalues are positive.
#' If the given covariance matrix is not positive definite, 
#' \code{correct_positive_definite} tries to modify the underlying correlation matrices
#' genetic_corrmat and full_corrmat in order to obtain a positive definite 
#' covariance matrix.
#'
#' This function can be used to verify that a given covariance matrix 
#' is positive definite. It calculates all eigenvalues in order to
#' investigate whether they are all positive. This property is necessary 
#' for the covariance matrix to be used as a Gaussian covariance matrix.
#' It is especially useful to check whether any covariance matrix obtained
#' by \code{\link{construct_covmat_multi}} is positive definite.
#' If the given covariance matrix is not positive definite, \code{correct_positive_definite}
#' tries to modify the underlying correlation matrices (called \code{genetic_corrmat} and 
#' \code{full_corrmat} in \code{\link{construct_covmat}} or \code{\link{construct_covmat_multi}}) by 
#' multiplying all off-diagonal entries in the correlation matrices by a given number.
#'
#' @param covmat A symmetric and numeric matrix. If the covariance matrix 
#' should be corrected, it must have a number of attributes, such as
#' \code{attr(covmat,"fam_vec")}, \code{attr(covmat,"n_fam")}, 
#' \code{attr(covmat,"add_ind")}, \code{attr(covmat,"h2")},
#' \code{attr(covmat,"genetic_corrmat")}, \code{attr(covmat,"full_corrmat")}
#' and \code{attr(covmat,"phenotype_names")}. Any covariance matrix 
#' obtained by \code{\link{construct_covmat}}, \code{\link{construct_covmat_single}} 
#' or \code{\link{construct_covmat_multi}} will have these attributes by default. 
#' @param correction_val A positive number representing the amount by which
#' \code{genetic_corrmat} and \code{full_corrmat} will be changed, if some
#' eigenvalues are non-positive. That is, correction_val is the number that will be
#' multiplied to all off_diagonal entries in \code{genetic_corrmat} and \code{full_corrmat}.
#' Defaults to 0.99.
#' @param correction_limit A positive integer representing the upper limit for the correction
#' procedure. Defaults to 100.
#' 
#' @return If \code{covmat} is a symmetric and numeric matrix and all eigenvalues are
#' positive, \code{correct_positive_definite} simply returns \code{covmat}. If some 
#' eigenvalues are not positive and \code{correction_val} is a positive number, 
#' \code{correct_positive_definite} tries to convert \code{covmat} into a positive definite
#' matrix. If \code{covmat} has attributes \code{add_ind}, \code{h2},
#' \code{genetic_corrmat}, \code{full_corrmat} and \code{phenotype_names}, 
#' \code{correct_positive_definite} computes a new covariance matrix using slightly
#' modified correlation matrices \code{genetic_corrmat} and \code{full_corrmat}.
#' If the correction is performed successfully, i.e. if the new covariance matrix 
#' is positive definite,the new covariance matrix is returned. 
#' Otherwise, \code{correct_positive_definite} returns the original covariance matrix.
#' 
#' @seealso \code{\link{construct_covmat}}, \code{\link{construct_covmat_single}} and
#' \code{\link{construct_covmat_multi}}.
#' 
#' @export
correct_positive_definite = function(covmat, correction_val = .99, correction_limit = 100) {
  
  # Checking that covmat is symmetric
  if(!isSymmetric.matrix(covmat)) stop("The covariance matrix covmat must be symmetric!")
  # and numeric
  if(!is.numeric(covmat)) stop("The covariance matrix covmat must be numeric!")
  
  # Checking whether all eigenvalues are positive
  if(any(eigen(covmat)$values < 0)){
    
    cat("The specified covariance matrix is not positive definite. \n")
  }else{
    
    return(covmat)
  }
  
  # If some eigenvalues are negative, correction_val must be specified,
  # it must be numeric and positive in order to correct the covariance matrix.
  if(!is.numeric(correction_val) && !is.integer(correction_val)) stop("correction_val must be numeric!")
  if(correction_val <= 0) stop("correction_val must be positive!")
  
  # In addition, covmat must have several attributes holding the 
  # family members (fam_vec or n_fam), a logical add_ind as well as
  # a numeric value or numeric matrix h2 and numeric matrices
  # genetic_corrmat and full_corrmat in order to change
  # the correlation matrix.
  if(is.null(attr(covmat,"add_ind")) || is.null(attr(covmat,"h2")) || 
     is.null(attr(covmat,"genetic_corrmat")) || is.null(attr(covmat,"full_corrmat"))){
    
    warning("The required attributes are missing... The covariance matrix could not be corrected!")
    return(covmat)
  }
  
  # If the covariance matrix is for a single phenotype, it is not
  # possible to correct the covariance matrix.
  if(length(attr(covmat,"h2"))==1){
    warning("The covariance matrix cannot be corrected...")
    return(covmat)
  }
  
  # Furthermore, correction_limit must be a positive number
  if(!is.numeric(correction_limit)) stop("correction_limit must be numeric!")
  if(correction_limit <= 0) stop("Correction limit must be positive!")
  
  cat("Trying to correct the covariance matrix...\n")
  
  # The covariance matrix will be modified at most correction_limit times.
  n <- 0
  # Storing the old covariance matrix
  old_covmat <- covmat
  
  # We also extract the vectors holding the family members
  fam_vec <- setdiff(attr(covmat,"fam_vec"), c("g","o"))
  #n_fam <- attr(covmat,"n_fam")[stringr::str_detect(names(attr(covmat,"n_fam")), "^[^go]")]
  
  while(any(eigen(covmat)$values < 0) && n <= correction_limit) {
    
    # Changing the correlation matrices slightly by 
    # multiplying all entries by correction_val.
    genetic_corrmat <- attr(covmat,"genetic_corrmat")*correction_val
    diag(genetic_corrmat) <- 1
    full_corrmat <- attr(covmat,"full_corrmat")*correction_val
    diag(full_corrmat) <- 1
    # Computing a new covariance matrix
    covmat <- construct_covmat(fam_vec = fam_vec, 
                               n_fam = NULL,
                               add_ind = attr(covmat,"add_ind"), 
                               genetic_corrmat = genetic_corrmat,
                               full_corrmat = full_corrmat,
                               h2 = attr(covmat,"h2"), 
                               phen_names = attr(covmat,"phenotype_names"))
    # Updating 
    n <- n+1
  }
  
  # If the matrix has been modified correction_limit times, 
  # the correction step is aborted
  if(n > correction_limit){
    cat("The covariance matrix could not be corrected. Consider revisiting it.\n")
    return(old_covmat)
  }
  
  cat(paste0("The correction was performed successfully! All off-diagonal entries are corrected by", round(correction_val^n,digits = 3),".\n"))
  
  return(covmat)
}

#' A simplied version of correct_positive_definite
#' 
#' multiplies off-diagonal elements of a matrix with a correction value until all eigen values are positive or correction limit is reached.
#' 
#' @param covmat Covariance matrix that needs to be positive definite
#' @param correction_limit The maximum number of times to correct off-diagonal elements
#' @param correction_val The value to multiply on the off-diagnoal elements.
#' 
#' @return Returns list with the corrected covmat and number of iterations used to achieve positive definite.
#' 

correct_positive_definite_simplified = function(covmat, correction_limit = 100, correction_val = 0.99) {
  eigen_val = eigen(covmat)

  if ( any(eigen_val$values < 0) ) {
    og_diag = diag(covmat)
    n = 0

    while (any(eigen(covmat)$values < 0) & n <= correction_limit ) {
      covmat = covmat*correction_val
      diag(covmat) = og_diag
      n = n + 1
    }

    if (eigen(eigen_val$values < 0)) stop("Unable to enforce a positive definite covariance matrix.")
    return(list(covmat = covmat, nitr = n))
  } else {
    return(list(covmat = covmat, nitr = 0))
  }
}



#' CDF for truncated normal distribution.
#' 
#' \code{truncated_normal_cdf} computes the cumulative density
#' function for a truncated normal distribution. 
#' 
#' This function can be used to compute the value of the cumulative
#' density function for a truncated normal distribution given an
#' individual's true underlying liability. 
#'
#' @param liability  A number representing the individual's 
#' true underlying liability. 
#' @param lower A number representing the lower cutoff point for the 
#' truncated normal distribution. Defaults to 1.645 
#' (stats::qnorm(0.05, lower.tail = F)).
#' @param upper A number representing the upper cutoff point of the 
#' truncated normal distribution. Must be greater or equal to lower.
#' Defaults to Inf.
#' 
#' @return If liability is a number and the lower and upper cutoff points
#' are numbers satisfying lower <= upper, then \code{truncated_normal_cdf}
#' returns the probability that the liability will take on a value less than
#' or equal to \code{liability}.
#'
#' @examples
#' \dontrun{
#'  curve(truncated_normal_cdf, from = qnorm(0.05, lower.tail = F), to = 3.5, 
#'        xname = "Underlying liability")
#' }
#' @export
truncated_normal_cdf = function(liability, lower = stats::qnorm(0.05, lower.tail = F), upper = Inf) {
  
  # Checking that the liability is valid
  if(!is.numeric(liability) && !is.integer(liability)) stop("The liability must be numeric!")
  # Checking that the lower and upper cutoff points are valid
  if(!is.numeric(lower)&& !is.integer(lower)) stop("The lower cutoff point must be numeric!")
  if(!is.numeric(upper) && !is.integer(upper)) stop("The upper cutoff point must be numeric!")
  if(upper < lower){
    cat("The upper cutoff point is below the lower cutoff point! \n 
The upper and lower cutoff points will be swapped...")
    
    temp <- lower
    lower <- upper
    upper <- temp
  }
  
  return(stats::pnorm(liability) - stats::pnorm(lower)) / (stats::pnorm(upper) - stats::pnorm(lower))
}


#' Convert age to cumulative incidence rate
#' 
#' \code{convert_age_to_cir} computes the cumulative incidence 
#' rate from a person's age. 
#' 
#' Given a person's age, \code{convert_age_to_cir} can be used
#' to compute the cumulative incidence rate (cir), which is given
#' by the formula 
#' \deqn{pop_prev / (1 + exp((mid_point - age) * slope))}
#'
#' @param age A non-negative number representing the individual's age.
#' @param pop_prev A positive number representing the overall
#' population prevalence. Must be at most 1. Defaults to 0.1.
#' @param mid_point A positive number representing the mid point 
#' logistic function. Defaults to 60. 
#' @param slope A number holding the rate of increase.
#' Defaults to 1/8.
#' 
#' @return If age and mid_point are positive numbers, if pop_prev
#' is a positive number between 0 and 1 and if slope is a valid number,
#' then \code{convert_age_to_cir} returns a number, which is equal to
#' the cumulative incidence rate.
#'
#' @examples 
#' \dontrun{
#' curve(convert_age_to_cir, from = 10, to = 110, xname = "age")
#' }
#' @export
convert_age_to_cir = function(age, pop_prev = .1, mid_point = 60, slope = 1/8) {
  
  # Checking that age is valid
  if(!is.numeric(age) && !is.integer(age)) stop("The age must be numeric!")
  if(age<0) stop("The age must be non-negative!")
  if(age >=150) warning("At this point, it is unrealistic to be of age 150 or older...")
  
  # Checking that pop_prev is valid
  if(!is.numeric(pop_prev) && !is.integer(pop_prev)) stop("The population prevalence pop_prev must be numeric!")
  if(pop_prev<=0)stop("The population prevalence pop_prev must be positive!")
  if(pop_prev>1)stop("The population prevalence pop_prev must be smaller or equal to 1!")
  
  # Checking that mid_point is valid
  if(!is.numeric(mid_point) && !is.integer(mid_point)) stop("The mid point mid_point must be numeric!")
  if(mid_point<=0)stop("The mid point mid_point must be positive!")
  
  # Checking that slope is valid
  if(!is.numeric(slope) && !is.integer(slope)) stop("The slope must be numeric!")
  
  return(pop_prev / (1 + exp((mid_point - age) * slope)))
}

#' Convert age to threshold
#' 
#' \code{convert_age_to_thresh} computes the threshold
#' from a person's age using either the logistic function
#' or the truncated normal distribution
#' 
#' Given a person's age, \code{convert_age_to_thresh} can be used
#' to first compute the cumulative incidence rate (cir), which is 
#' then used to compute the threshold using either the
#' logistic function or the truncated normal distribution.
#' Under the logistic function, the formula used to compute
#' the threshold from an individual's age is given by
#' \deqn{qnorm(pop_prev / (1 + exp((mid_point - age) * slope)), lower.tail = F)},
#' while it is given by 
#' \deqn{qnorm((1 - (age-min_age)/max_age) * (pnorm(upper) - pnorm(lower)) + pnorm(lower))}
#' under the truncated normal distribution.
#'
#' @param age A non-negative number representing the individual's age.
#' @param dist A string indicating which distribution to use. 
#' If dist = "logistic", the logistic function will be used to 
#' compute the age of onset.
#' If dist = "normal", the truncated normal distribution will be used instead. 
#' Defaults to "logistic".
#' @param pop_prev Only necessary if dist = "logistic". A positive number representing the overall
#' population prevalence. Must be at most 1. Defaults to 0.1.
#' @param mid_point Only necessary if dist = "logistic". A positive number representing the mid point 
#' logistic function. Defaults to 60. 
#' @param slope Only necessary if dist = "logistic". A number holding the rate of increase.
#' Defaults to 1/8.
#' @param min_age Only necessary if dist = "normal". A positive number representing the individual's earliest age. 
#' Defaults to 10.
#' @param max_age Only necessary if dist = "normal". A positive number representing the individual's latest age. 
#' Must be greater than min_aoo. Defaults to 90.
#' @param lower Only necessary if dist = "normal". A number representing the lower cutoff point for the 
#' truncated normal distribution. Defaults to 1.645 
#' (stats::qnorm(0.05, lower.tail = F)).
#' @param upper Only necessary if dist = "normal". A number representing the upper cutoff point of the 
#' truncated normal distribution. Must be greater or equal to lower.
#' Defaults to Inf.
#' 
#' @return If age is a positive number and all other necessary arguments are valid,
#' then \code{convert_age_to_thresh} returns a number, which is equal to
#' the threshold.
#' 
#' @examples 
#' # curve(convert_age_to_thresh, from = 10, to = 110, xname = "Age")
#' # curve(convert_age_to_thresh, from = 10, to = 100, xname = "Age")
#' @export
convert_age_to_thresh = function(age, dist = "logistic", pop_prev = .1, mid_point = 60, slope = 1/8,
                                 min_age = 10, max_age = 90, lower = stats::qnorm(0.05, lower.tail = F), upper = Inf) {
  # Checking that age is valid
  if(!is.numeric(age)&& !is.integer(age)) stop("The age must be numeric!")
  if(age<0) stop("The age must be non-negative!")
  
  # Checking that dist is either logistic or normal.
  if(is.character(dist)){
    
    valid_is <- grep("logistic|normal", dist)
    
  }else{
    stop("dist must be a string!")
  }
  
  # Checking whether dist is empty or of length >1.
  if(length(valid_is) == 0){
    
    cat("Warning message: \n dist is not of the required format! \n The function will use the logistic function to compute the age of onset!")
    dist <- "logistic"
  }else {
    if(length(valid_is) > 1) cat("Warning message: \n dist is not of the required format! \n The function will use the first valid entry to compute the age of onset!")
    dist <- dist[valid_is[1]]
  }
  
  # If dist = logistic, the logistic function will be used to compute the age of onset
  if(dist == "logistic"){
    
    # Checking that pop_prev is valid
    if(!is.numeric(pop_prev) && !is.integer(pop_prev)) stop("The population prevalence pop_prev must be numeric!")
    if(pop_prev<=0)stop("The population prevalence pop_prev must be positive!")
    if(pop_prev>1)stop("The population prevalence pop_prev must be smaller or equal to 1!")
    
    # Checking that mid_point is valid
    if(!is.numeric(mid_point) && !is.integer(mid_point)) stop("The mid point mid_point must be numeric!")
    if(mid_point<=0)stop("The mid point mid_point must be positive!")
    
    # Checking that slope is valid
    if(!is.numeric(slope) && !is.integer(slope)) stop("The slope must be numeric!")
    
    # Computing the threshold
    return(stats::qnorm(pop_prev / (1 + exp((mid_point - age) * slope)), lower.tail = F))
    
  }
  
  # if dist = normal, the truncated normal distribution will be used to compute the age of onset
  if(dist == "normal"){
    
    # Checking that min_age and max_age are valid.
    if(!is.numeric(min_age) && !is.integer(min_age)) stop("The earliest age min_age must be numeric!")
    if(!is.numeric(max_age) && !is.integer(max_age)) stop("The latest age max_age must be numeric!")
    if(min_age <= 0) stop("The earliest age min_age must be positive!")
    if(max_age <= 0) stop("The latest age max_age must be positive!")
    if(min_age > max_age){
      cat("The latest age max_age is below the earliest age min_age! \n 
The earliest and latest age will be swapped...")
      
      min_age <- min_age + max_age
      max_age <- min_age - max_age
      min_age <- min_age - max_age
    }
    
    # Checking that the lower and upper cutoff points are valid
    if(!is.numeric(lower) && !is.integer(lower)) stop("The lower cutoff point must be numeric!")
    if(!is.numeric(upper) && !is.integer(upper)) stop("The upper cutoff point must be numeric!")
    if(upper < lower){
      cat("The upper cutoff point is below the lower cutoff point! \n 
The upper and lower cutoff points will be swapped...")
      
      temp <- lower
      lower <- upper
      upper <- temp
    }
    # Computing the threshold
    return(stats::qnorm((1 - (age-min_age)/max_age) * (stats::pnorm(upper) - stats::pnorm(lower)) + stats::pnorm(lower)))
  }
}


#' Convert cumulative incidence rate to age
#' 
#' \code{convert_cir_to_age} computes the age
#' from a person's cumulative incidence rate. 
#' 
#' Given a person's cumulative incidence rate (cir), \code{convert_cir_to_age} 
#' can be used to compute the corresponding age, which is given by
#' \deqn{mid_point - \log(pop_prev/cir - 1) * 1/slope}
#'
#' @param cir A positive number representing the individual's cumulative 
#' incidence rate.
#' @param pop_prev A positive number representing the overall
#' population prevalence. Must be at most 1 and must be larger than 
#' cir. Defaults to 0.1.
#' @param mid_point A positive number representing the mid point 
#' logistic function. Defaults to 60. 
#' @param slope A number holding the rate of increase.
#' Defaults to 1/8.
#' 
#' @return If cir and mid_point are positive numbers, if pop_prev
#' is a positive number between 0 and 1 and if slope is a valid number,
#' then \code{convert_cir_to_age} returns a number, which is equal to
#' the current age.
#' 
#' @examples 
#' \dontrun{
#' curve(convert_cir_to_age, from = 0, to = 0.1, xname = "Cumulative incidence rate")
#' }
#' @export
convert_cir_to_age = function(cir, pop_prev = .1, mid_point = 60, slope = 1/8) {
  
  # Checking that age is valid
  if(!is.numeric(cir) && !is.integer(cir)) stop("The cumulative incidence rate cir must be numeric!")
  if(cir<=0) stop("The cumulative incidence rate cir must be positive!")
  
  # Checking that pop_prev is valid
  if(!is.numeric(pop_prev) && !is.integer(pop_prev)) stop("The population prevalence pop_prev must be numeric!")
  if(pop_prev<=0)stop("The population prevalence pop_prev must be positive!")
  if(pop_prev>1)stop("The population prevalence pop_prev must be smaller or equal to 1!")
  
  # Checking that mid_point is valid
  if(!is.numeric(mid_point) && !is.integer(mid_point)) stop("The mid point mid_point must be numeric!")
  if(mid_point<=0)stop("The mid point mid_point must be positive!")
  
  # Checking that slope is valid
  if(!is.numeric(slope) && !is.integer(slope)) stop("The slope must be numeric!")
  
  if(cir >= pop_prev){
    
    return(NA)
  }else{
    
    res <- mid_point - log(pop_prev/cir - 1)* 1/slope
    
    if(res > 0) return(res)
    if(res <= 0) return(0)
  }
}


#' Convert liability to age of onset
#' 
#' \code{convert_liability_to_aoo} computes the age
#' of onset from an individual's true underlying liability using 
#' either the logistic function or the truncated normal distribution.
#' 
#' Given a person's cumulative incidence rate (cir), \code{convert_liability_to_aoo} 
#' can be used to compute the corresponding age. Under the logistic function,
#' the age is given by
#' \deqn{mid_point - log(pop_prev/cir - 1) * 1/slope},
#' while it is given by 
#' \deqn{(1 - truncated_normal_cdf(liability = liability, lower = lower , upper = upper)) * max_aoo + min_aoo}
#' under the truncated normal distribution.
#' 
#' @param liability A number representing the individual's 
#' true underlying liability. 
#' @param dist A string indicating which distribution to use. 
#' If dist = "logistic", the logistic function will be used to 
#' compute the age of onset.
#' If dist = "normal", the truncated normal distribution will be used instead. 
#' Defaults to "logistic".
#' @param pop_prev Only necessary if dist = "logistic". A positive number representing the overall
#' population prevalence. Must be at most 1. Defaults to 0.1.
#' @param mid_point Only necessary if dist = "logistic". A positive number representing the mid point 
#' logistic function. Defaults to 60. 
#' @param slope Only necessary if dist = "logistic". A number holding the rate of increase.
#' Defaults to 1/8.
#' @param min_aoo Only necessary if dist = "normal". A positive number representing the individual's earliest age of onset. 
#' Defaults to 10.
#' @param max_aoo Only necessary if dist = "normal". A positive number representing the individual's latest age of onset. 
#' Must be greater than min_aoo. Defaults to 90.
#' @param lower Only necessary if dist = "normal". A number representing the lower cutoff point for the 
#' truncated normal distribution. Defaults to 1.645 
#' (stats::qnorm(0.05, lower.tail = F)).
#' @param upper Only necessary if dist = "normal". A number representing the upper cutoff point of the 
#' truncated normal distribution. Must be greater or equal to lower.
#' Defaults to Inf.
#' 
#' @return If liability is a number and all other necessary arguments are valid,
#' then \code{convert_liability_to_aoo} returns a positive number, which is equal to
#' the age of onset.
#' 
#' @examples 
#' \dontrun{
#' curve(convert_liability_to_aoo, from = 1.3, to = 3.5, 
#'       xname = "Underlying liability")
#' curve(convert_liability_to_aoo, from = qnorm(0.05, lower.tail = F), to = 3.5, 
#'       xname = "Underlying liability", dist = "normal")
#' }
#' @export
convert_liability_to_aoo = function(liability, dist = "logistic", pop_prev = .1, mid_point = 60, slope = 1/8,
                                    min_aoo = 10, max_aoo = 90, lower = stats::qnorm(0.05, lower.tail = F), upper = Inf ) {
  
  # Checking that liability is valid
  if(!is.numeric(liability) && !is.integer(liability)) stop("The liability must be numeric!")
  
  # Checking that dist is either logistic or normal.
  if(is.character(dist)){
    
    valid_is <- grep("logistic|normal", dist)
    
  }else{
    stop("dist must be a string!")
  }
  
  # Checking whether dist is empty or of length >1.
  if(length(valid_is) == 0){
    
    cat("Warning message: \n dist is not of the required format! \n The function will use the logistic function to compute the age of onset!")
    dist <- "logistic"
  }else {
    
    if(length(valid_is) > 1) cat("Warning message: \n dist is not of the required format! \n The function will use the first valid entry to compute the age of onset!")
    dist <- dist[valid_is[1]]
  }
  
  # If dist = logistic, the logistic function will be used to compute the age of onset
  if(dist == "logistic"){
    
    # Checking that pop_prev is valid
    if(!is.numeric(pop_prev) && !is.integer(pop_prev)) stop("The population prevalence pop_prev must be numeric!")
    if(pop_prev<=0)stop("The population prevalence pop_prev must be positive!")
    if(pop_prev>1)stop("The population prevalence pop_prev must be smaller or equal to 1!")
    
    # Checking that mid_point is valid
    if(!is.numeric(mid_point) && !is.integer(mid_point)) stop("The mid point mid_point must be numeric!")
    if(mid_point<=0)stop("The mid point mid_point must be positive!")
    
    # Checking that slope is valid
    if(!is.numeric(slope) && !is.integer(slope)) stop("The slope must be numeric!")
    
    # Computing the age of onset
    if(stats::pnorm(liability, lower.tail = F) >= pop_prev){
      return(NA)
    }else{
       
      res <- mid_point - log(pop_prev/stats::pnorm(liability, lower.tail = F) - 1)* 1/slope
      
      if(res > 0) return(res)
      if(res <= 0) return(0)
    }
  }
  
  # if dist = normal, the truncated normal distribution will be used to compute the age of onset
  if(dist == "normal"){
    
    # Checking that min_aoo and max_aoo are valid.
    if(!is.numeric(min_aoo) && !is.integer(min_aoo)) stop("The earliest age of onset min_aoo must be numeric!")
    if(!is.numeric(max_aoo) && !is.integer(max_aoo)) stop("The latest age of onset max_aoo must be numeric!")
    if(min_aoo <= 0) stop("The earliest age of onset min_aoo must be positive!")
    if(max_aoo <= 0) stop("The latest age of onset max_aoo must be positive!")
    if(min_aoo > max_aoo){
      cat("The latest age of onset max_aoo is below the earliest age of onset min_aoo! \n 
The earliest and latest age of onset will be swapped...")
      
      min_aoo <- min_aoo + max_aoo
      max_aoo <- min_aoo - max_aoo
      min_aoo <- min_aoo - max_aoo
    }
    
    # Checking that the lower and upper cutoff points are valid
    if(!is.numeric(lower) && !is.integer(lower)) stop("The lower cutoff point must be numeric!")
    if(!is.numeric(upper) && !is.integer(upper)) stop("The upper cutoff point must be numeric!")
    if(upper < lower){
      cat("The upper cutoff point is below the lower cutoff point! \n 
The upper and lower cutoff points will be swapped...")

      temp <- lower
      lower <- upper
      upper <- temp
    }
    
    # Computing the age of onset
    res <- (1 - truncated_normal_cdf(liability = liability, lower = lower , upper = upper)) * max_aoo + min_aoo
    
    return(res)
  }
}

#' Convert the heritability on the observed scale to that on the liability scale
#' 
#' \code{convert_observed_to_liability_scale} transforms the heritability on the 
#' observed scale to the heritability on the liability scale.
#' 
#' This function can be used to transform the heritability on the observed 
#' scale to that on the liability scale. \code{convert_observed_to_liability_scale}
#' uses either Equation 17 (if prop_cases = NULL) or Equation 23 from
#' Sang Hong Lee, Naomi R. Wray, Michael E. Goddard and Peter M. Visscher, "Estimating
#' Missing Heritability for Diseases from Genome-wide Association Studies",
#' The American Journal of Human Genetics, Volume 88, Issue 3, 2011, pp. 294-305,
#' \doi{10.1016/j.ajhg.2011.02.002} to transform the heritability on the observed
#' scale to the heritability on the liability scale.
#' 
#' @param obs_h2 A number or numeric vector representing the liability-scale 
#' heritability(ies)on the observed scale. Must be non-negative and at most 1. 
#' Defaults to 0.5
#' @param pop_prev A number or numeric vector representing the population prevalence(s). All 
#' entries must be non-negative and at most one.
#' If it is a vector, it must have the same length as obs_h2. Defaults to 0.05. 
#' @param prop_cases Either NULL or a number or a numeric vector representing the proportion
#' of cases in the sample. All entries must be non-negative and at most one. 
#' If it is a vector, it must have the same length as obs_h2. Defaults to 0.5.
#' 
#' If \code{obs_h2}, \code{pop_prev} and \code{prop_cases} are non-negative numbers 
#' that are at most one, the function returns the heritability on the liability
#' scale using Equation 23 from 
#' Sang Hong Lee, Naomi R. Wray, Michael E. Goddard and Peter M. Visscher, "Estimating
#' Missing Heritability for Diseases from Genome-wide Association Studies",
#' The American Journal of Human Genetics, Volume 88, Issue 3, 2011, pp. 294-305,
#' \doi{10.1016/j.ajhg.2011.02.002}.
#' If \code{obs_h2}, \code{pop_prev} and \code{prop_cases} are non-negative numeric
#' vectors where all entries are at most one, the function returns a vector of the same
#' length as obs_h2. Each entry holds to the heritability on the liability
#' scale which was obtained from the corresponding entry in obs_h2 using Equation 23.
#' If \code{obs_h2} and \code{pop_prev} are non-negative numbers that are at most
#' one and \code{prop_cases} is \code{NULL}, the function returns the heritability 
#' on the liability scale using Equation 17 from 
#' Sang Hong Lee, Naomi R. Wray, Michael E. Goddard and Peter M. Visscher, "Estimating
#' Missing Heritability for Diseases from Genome-wide Association Studies",
#' The American Journal of Human Genetics, Volume 88, Issue 3, 2011, pp. 294-305,
#' \doi{10.1016/j.ajhg.2011.02.002}. 
#' If \code{obs_h2} and \code{pop_prev} are non-negative numeric vectors such that
#' all entries are at most one, while \code{prop_cases} is \code{NULL},
#' \code{convert_observed_to_liability_scale} returns a vector of the same
#' length as obq_h2. Each entry holds to the liability-scale heritability that 
#' was obtained from the corresponding entry in obs_h2 using Equation 17.
#' 
#' @examples 
#' convert_observed_to_liability_scale()
#' convert_observed_to_liability_scale(prop_cases=NULL)
#' convert_observed_to_liability_scale(obs_h2 = 0.8, pop_prev = 1/44, 
#'                                     prop_cases = NULL)
#' convert_observed_to_liability_scale(obs_h2 = c(0.5,0.8), 
#'                                     pop_prev = c(0.05, 1/44), 
#'                                     prop_cases = NULL)
#' 
#' @references
#' Sang Hong Lee, Naomi R. Wray, Michael E. Goddard, Peter M. Visscher (2011, March). Estimating
#' Missing Heritability for Diseases from Genome-wide Association Studies. In The American Journal 
#' of Human Genetics (Vol. 88, Issue 3, pp. 294-305). \doi{10.1016/j.ajhg.2011.02.002}
#' 
#' @export
convert_observed_to_liability_scale <- function(obs_h2 = 0.5, pop_prev = 0.05, prop_cases = 0.5){
  
  # Checking that the observed heritabilities are valid
  if(!is.numeric(obs_h2) && !is.integer(obs_h2)) stop("The observed heritability(ies) must be numeric!")
  if(any(obs_h2<0))stop("The observed heritability(ies) must be non-negative!")
  if(any(obs_h2>1))stop("The observed heritability(ies) must be smaller than or equal to one!")
  # Checking that the population prevalences are valid
  if(!is.numeric(pop_prev) && !is.integer(pop_prev))stop("The population prevalence(s) must be numeric!")
  if(any(pop_prev<0))stop("The population prevalence(s) must be non-negative!")
  if(any(pop_prev>1))stop("The population prevalence(s) must be smaller than or equal to one!")
  
  # Defining the variable z, which is the height of the truncated
  # normal curve at the point t, and where t is the truncated point,
  # such that the fraction of observations larger than t is equal to 
  # the population prevalence pop_prev. That is
  # t = qnorm(pop_prev, lower.tail = FALSE)
  z <- stats::dnorm(stats::qnorm(pop_prev, lower.tail = FALSE))
  
  # Using one of the two possible transformations depending on whether 
  # prop_cases is NULL or not.
  if(is.null(prop_cases)){
    
    return(obs_h2 * (pop_prev*(1-pop_prev))/(z^2))
  }else{
    
    # Checking that the proportions of cases are valid
    if(!is.numeric(prop_cases) && !is.integer(prop_cases)) stop("The proportion(s) of cases must be numeric!")
    if(any(prop_cases<0))stop("The proportion(s) of cases must be non-negative!")
    if(any(prop_cases>1))stop("The proportion(s) of cases must be smaller than or equal to one!")
    
    return(obs_h2 * (pop_prev*(1-pop_prev))/(z^2) * (pop_prev*(1-pop_prev))/(prop_cases*(1-prop_cases)))
  }
}


# graph helper functions --------------------------------------------------

#' 
#' Add pseudo point to graph
#' 
#' Internal function. Identifies the proband and makes a duplicated point that will act as genetic liability in kinship matrix
#' 
#' @param fam_graph family graph center on proband
#' @param index_id id of proband.
#' 
#' @return returns family graph with a pseudo point added that will have the same connections as the proband point. For kinship construction, the relationship to the proband must be adjusted to 1 \* h2 (and not 0.5 \* h2).
#' @export
#' 

add_gen_liab_to_graph = function(fam_graph, index_id) {
  # find all index edges
  all_edges = igraph::as_ids(igraph::E(fam_graph))
  # extract only edges that directly match the index_id.
  # edges are of the form str1|str2; to avoid matching the wrong substring,
  # we include string end and start as well as the divider "|" in the match
  
  # find edges that start with index_id
  startWith_edges = stringr::str_detect(all_edges, paste0("^", index_id, "\\|")) 
  # find edges that end with index_id
  endWith_edges = stringr::str_detect(all_edges, paste0("\\|", index_id, "$"))
  
  index_edges = all_edges[startWith_edges | endWith_edges ]
  # subset edges to and from index with new id - here with "_g" added
  
  # we need to handle the no edges case, i.e. one point graphs.
  if (length(index_edges)  > 0) {
    gen_liab_edges = stringr::str_replace(index_edges, index_id, paste0(index_id, "_g")) %>% 
      # format into vector with start and end of edge at positions 1 & 2, 3 & 4, etc
      # as defined by igraph.
      stringr::str_split("\\|") %>% 
      # rbind is much faster than cbind
      do.call("rbind", .) %>%
      # transpose for 2xN
      t() %>% 
      # c for column wise conversion to vector - each edge is a column
      c()
  }
  
  fam_graph %>% 
    igraph::add.vertices(., nv = 1, name = paste0(index_id, "_g")) %>% 
    # adding edges found above and edge between _g and index_id
    igraph::add.edges(., edges = c(if (length(index_edges) > 0) gen_liab_edges,
                                   paste0(rep(index_id, 4), c("", "_g", "_g", ""))))
}


#' Construct kinship matrix from graph
#' 
#' construct the kinship matrix from a graph representation of a family, centered on an index person (proband). 
#' 
#' 
#' @param fam_graph graph.
#' @param h2 heritability.
#' @param index_id proband id. Only used in conjuction with add_ind = TRUE.
#' @param add_ind add genetic liability to the kinship matrix. Defaults to true.
#' @param fix_diag Whether to set diagonal to 1 for all entries except for the 
#' genetic liability.
#' 
#' @return A kinship matrix. 
#' 
#' @export
#' 

get_kinship = function(fam_graph, h2, index_id = NA, add_ind = TRUE, fix_diag = TRUE) {
  if (add_ind) {
    index_id_g = paste0(index_id, "_g")
    # adding in point to act as genetic liability for index person
    # need to adjust individual and genetic liab point's internal relationship later
    fam_graph = add_gen_liab_to_graph(fam_graph, index_id)
    
  }
  
  # in-distances
  ph = igraph::distances(fam_graph, mode = "in")
  rnames = rownames(ph)
  nrows = nrow(ph)
  
  # construct dummy matrix
  new_distances = matrix(NA, ncol = nrows, nrow = nrows)
  rownames(new_distances) = rnames
  colnames(new_distances) = rnames
  
  # distance to self is always 0
  diag(new_distances) = 0
  
  # fill dummy matrix
  for (i in 1:(nrows - 1)) {
    for (j in (i + 1):nrows) {
      new_distances[i,j] <- new_distances[j,i] <- min(ph[i, ] + ph[j, ]) 
    }
  }
  
  # calculating kinship matrix, _g needs adjustment
  kinship = 0.5^new_distances
  if (add_ind) {
    # multipying with h2 and fixing diagonal to 1 - except for gen liab
    kinship[index_id, index_id_g] <- kinship[index_id_g, index_id] <- kinship[index_id_g, index_id_g] <- 1
    kinship = kinship * h2
    if (fix_diag) { # fixing diagonal at 1, except the genetic component
      diag(kinship)[-which(rownames(kinship) == index_id_g)] <- 1  
    }
    
  } else {
    kinship = kinship * h2
    
    if (fix_diag) diag(kinship) <- 1
  }
  
  return(kinship)
}


#' 
#' construct all combinations of input vector
#' 
#' pastes together all combinations of input vector
#' 
#' @param vec vector of strings
#' 
#' @return A vector of strings is returned.
#' 
#' @examples
#' \dontrun{
#' get_all_combs(letters[1:3])
#' }
#' 
#' @export
get_all_combs = function(vec) {
  # outer: creates all combinations of entries in vec
  # then we select only the off diagonal entries
  outer(vec, vec, paste, sep = "_")[!diag(TRUE, nrow = length(vec))]
  # results in length(vec) * (length(vec) - 1) entries
  # since it is the number of entries in the upper *and* lower triangle
  # number of entries in either is given by length(vec) * (length(vec) - 1)/2
}
