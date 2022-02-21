utils::globalVariables("cir")
utils::globalVariables("thresh")
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
#' # curve(truncated_normal_cdf, from = qnorm(0.05, lower.tail = F), to = 3.5, xname = "Underlying liability")
#' @export
truncated_normal_cdf = function(liability, lower = stats::qnorm(0.05, lower.tail = F), upper = Inf) {
  
  # Checking that the liability is valid
  if(class(liability) != "numeric") stop("The liability must be numeric!")
  # Checking that the lower and upper cutoff points are valid
  if(class(lower) != "numeric") stop("The lower cutoff point must be numeric!")
  if(class(upper) != "numeric") stop("The upper cutoff point must be numeric!")
  if(upper < lower){
    cat("The upper cutoff point is below the lower cutoff point! \n 
The upper and lower cutoff points will be swapped...")
    
    lower <- lower + upper
    upper <- lower - upper
    lower <- lower - upper
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
#' by the formula \eqn{ pop_prev / (1 + exp((mid_point - age) * slope))}
#'
#' @param age A positive number representing the individual's age.
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
#' 
#' # curve(convert_age_to_cir, from = 10, to = 110, xname = "age")
#' @export
convert_age_to_cir = function(age, pop_prev = .1, mid_point = 60, slope = 1/8) {
  
  # Checking that age is valid
  if(class(age) != "numeric") stop("The age must be numeric!")
  if(age<=0) stop("The age must be positive!")
  if(age >=150) warning("At this point, it is unrealistic to be of age 150 or older!")
  
  # Checking that pop_prev is valid
  if(class(pop_prev) != "numeric") stop("The population prevalence pop_prev must be numeric!")
  if(pop_prev<=0)stop("The population prevalence pop_prev must be positive!")
  if(pop_prev>1)stop("The population prevalence pop_prev must be smaller or equal to 1!")
  
  # Checking that mid_point is valid
  if(class(mid_point) != "numeric") stop("The mid point mid_point must be numeric!")
  if(mid_point<=0)stop("The mid point mid_point must be positive!")
  
  # Checking that slope is valid
  if(class(slope) != "numeric") stop("The slope must be numeric!")
  
  cir <- pop_prev / (1 + exp((mid_point - age) * slope))
  
  return(cir)
}

#' Convert age to threshold
#' 
#' \code{convert_age_to_thresh} computes the threshold
#' from a person's age using either the logistic function
#' or the truncated normal distribution
#' 
#' Given a person's age, \code{convert_age_to_thresh} can be used
#' to first compute the cumulative incidence rate (cir), which is 
#' then used to compute the threshold using  either 
#' 
#' 
#' 
#'  Overall, the formula used to compute
#' the threshold from an individual's age is given by
#' \eqn{qnorm(pop_prev / (1 + exp((mid_point - age) * slope)), lower.tail = F)}
#'
#' @param age A positive number representing the individual's age.
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
  if(class(age) != "numeric") stop("The age must be numeric!")
  if(age<=0) stop("The age must be positive!")
  if(age >=150) warning("At this point, it is unrealistic to be of age 150 or older!")
  
  # Checking that dist is either logistic or normal.
  if(class(dist) == "character"){
    
    dist <- c("logistic", "normal")[rowSums(sapply(dist, grepl, x = c("logistic", "normal"))) > 0]
    
  }else{
    stop("dist must be a string!")
  }
  
  # Checking whether dist is empty or of length >1.
  if(length(dist) == 0){
    
    cat("Warning message: \n dist is not of the required format! \n The function will use the logistic function to compute the age of onset!")
    dist <- "logistic"
  }else if(length(dist) >1){
    
    cat("Warning message: \n dist is not of the required format! \n The function will use the first valid entry to compute the age of onset!")
    dist <- dist[1]
  }
  
  # If dist = logistic, the logistic function will be used to compute the age of onset
  if(dist == "logistic"){
    
    # Checking that pop_prev is valid
    if(class(pop_prev) != "numeric") stop("The population prevalence pop_prev must be numeric!")
    if(pop_prev<=0)stop("The population prevalence pop_prev must be positive!")
    if(pop_prev>1)stop("The population prevalence pop_prev must be smaller or equal to 1!")
    
    # Checking that mid_point is valid
    if(class(mid_point) != "numeric") stop("The mid point mid_point must be numeric!")
    if(mid_point<=0)stop("The mid point mid_point must be positive!")
    
    # Checking that slope is valid
    if(class(slope) != "numeric") stop("The slope must be numeric!")
    
    # Computing the threshold
    thresh <- stats::qnorm(pop_prev / (1 + exp((mid_point - age) * slope)), lower.tail = F)
    return(thresh)
    
  }
  
  # if dist = normal, the truncated normal distribution will be used to compute the age of onset
  if(dist == "normal"){
    
    # Checking that min_age and max_age are valid.
    if(class(min_age) != "numeric") stop("The earliest age min_age must be numeric!")
    if(class(max_age) != "numeric") stop("The latest age max_age must be numeric!")
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
    if(class(lower) != "numeric") stop("The lower cutoff point must be numeric!")
    if(class(upper) != "numeric") stop("The upper cutoff point must be numeric!")
    if(upper < lower){
      cat("The upper cutoff point is below the lower cutoff point! \n 
The upper and lower cutoff points will be swapped...")
      
      lower <- lower + upper
      upper <- lower - upper
      lower <- lower - upper
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
#' Given a person's cumulative incidence rate (cir), \code{{convert_cir_to_age} 
#' can be used to compute the corresponding age, which is given by
#' \eqn{mid_point - log(pop_prev/cir - 1) * 1/slope}
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
#' # curve(convert_cir_to_age, from = 0, to = 0.1, xname = "Cumulative incidence rate")
#' @export
convert_cir_to_age = function(cir, pop_prev = .1, mid_point = 60, slope = 1/8) {
  
  # Checking that age is valid
  if(class(cir) != "numeric") stop("The cumulative incidence rate cir must be numeric!")
  if(cir<=0) stop("The cumulative incidence rate cir must be positive!")
  
  # Checking that pop_prev is valid
  if(class(pop_prev) != "numeric") stop("The population prevalence pop_prev must be numeric!")
  if(pop_prev<=0)stop("The population prevalence pop_prev must be positive!")
  if(pop_prev>1)stop("The population prevalence pop_prev must be smaller or equal to 1!")
  
  # Checking that mid_point is valid
  if(class(mid_point) != "numeric") stop("The mid point mid_point must be numeric!")
  if(mid_point<=0)stop("The mid point mid_point must be positive!")
  
  # Checking that slope is valid
  if(class(slope) != "numeric") stop("The slope must be numeric!")
  
  if(cir >= pop_prev){
    
    return(NA)
  }else{
    
    return(mid_point - log(pop_prev/cir - 1)* 1/slope)
  }
}


#' Convert liability to age of onset
#' 
#' \code{convert_liability_to_aoo} computes the age
#' of onset from an individual's true underlying liability using 
#' either the logistic function or the truncated normal distribution.
#' 
#' Given a person's cumulative incidence rate (cir), \code{{convert_liability_to_aoo} 
#' can be used to compute the corresponding age. Under the logistic function,
#' the age is given by
#' \eqn{mid_point - log(pop_prev/cir - 1) * 1/slope,}
#' while it is given by 
#' \eqn{(1 - truncated_normal_cdf(liability = liability, lower = lower , upper = upper)) * max_aoo + min_aoo}
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
#' # curve(convert_liability_to_aoo, from = 1.3, to = 3.5, xname = "Underlying liability")
#' # curve(convert_liability_to_aoo, from = qnorm(0.05, lower.tail = F), to = 3.5, xname = "Underlying liability", dist = "normal")
#' @export
convert_liability_to_aoo = function(liability, dist = "logistic", pop_prev = .1, mid_point = 60, slope = 1/8,
                                    min_aoo = 10, max_aoo = 90, lower = stats::qnorm(0.05, lower.tail = F), upper = Inf ) {
  
  # Checking that liability is valid
  if(class(liability) != "numeric") stop("The liability must be numeric!")
  
  # Checking that dist is either logistic or normal.
  if(class(dist) == "character"){
    
    dist <- c("logistic", "normal")[rowSums(sapply(dist, grepl, x = c("logistic", "normal"))) > 0]
    
  }else{
    stop("dist must be a string!")
  }
  
  # Checking whether dist is empty or of length >1.
  if(length(dist) == 0){
    
    cat("Warning message: \n dist is not of the required format! \n The function will use the logistic function to compute the age of onset!")
    dist <- "logistic"
  }else if(length(dist) >1){
    
    cat("Warning message: \n dist is not of the required format! \n The function will use the first valid entry to compute the age of onset!")
    dist <- dist[1]
  }
  
  # If dist = logistic, the logistic function will be used to compute the age of onset
  if(dist == "logistic"){
    
    # Checking that pop_prev is valid
    if(class(pop_prev) != "numeric") stop("The population prevalence pop_prev must be numeric!")
    if(pop_prev<=0)stop("The population prevalence pop_prev must be positive!")
    if(pop_prev>1)stop("The population prevalence pop_prev must be smaller or equal to 1!")
    
    # Checking that mid_point is valid
    if(class(mid_point) != "numeric") stop("The mid point mid_point must be numeric!")
    if(mid_point<=0)stop("The mid point mid_point must be positive!")
    
    # Checking that slope is valid
    if(class(slope) != "numeric") stop("The slope must be numeric!")
    
    # Computing the age of onset
    if(stats::pnorm(liability, lower.tail = F) >= pop_prev){
      return(NA)
    }else{
      return(mid_point - log(pop_prev/stats::pnorm(liability, lower.tail = F) - 1)* 1/slope)
    }
  }
  
  # if dist = normal, the truncated normal distribution will be used to compute the age of onset
  if(dist == "normal"){
    
    # Checking that min_aoo and max_aoo are valid.
    if(class(min_aoo) != "numeric") stop("The earliest age of onset min_aoo must be numeric!")
    if(class(max_aoo) != "numeric") stop("The latest age of onset max_aoo must be numeric!")
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
    if(class(lower) != "numeric") stop("The lower cutoff point must be numeric!")
    if(class(upper) != "numeric") stop("The upper cutoff point must be numeric!")
    if(upper < lower){
      cat("The upper cutoff point is below the lower cutoff point! \n 
The upper and lower cutoff points will be swapped...")
      
      lower <- lower + upper
      upper <- lower - upper
      lower <- lower - upper
    }
    
    # Computing the age of onset
    return((1 - truncated_normal_cdf(liability = liability, lower = lower , upper = upper)) * max_aoo + min_aoo)
  }
  
}