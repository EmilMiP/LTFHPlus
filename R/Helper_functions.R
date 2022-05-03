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
#' \dontrun{
#'  curve(truncated_normal_cdf, from = qnorm(0.05, lower.tail = F), to = 3.5, 
#'        xname = "Underlying liability")
#' }
#' @export
truncated_normal_cdf = function(liability, lower = stats::qnorm(0.05, lower.tail = F), upper = Inf) {
  
  # Checking that the liability is valid
  if(class(liability) != "numeric" && class(liability) != "integer") stop("The liability must be numeric!")
  # Checking that the lower and upper cutoff points are valid
  if(class(lower) != "numeric" && class(lower) != "integer") stop("The lower cutoff point must be numeric!")
  if(class(upper) != "numeric" && class(upper) != "integer") stop("The upper cutoff point must be numeric!")
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
  if(class(age) != "numeric" && class(age) != "integer") stop("The age must be numeric!")
  if(any(age<0)) stop("The age must be non-negative!")
  if(any(age >=150)) warning("At this point, it is unrealistic to be of age 150 or older!")
  
  # Checking that pop_prev is valid
  if(class(pop_prev) != "numeric" && class(pop_prev) != "integer") stop("The population prevalence pop_prev must be numeric!")
  if(any(pop_prev<=0))stop("The population prevalence pop_prev must be positive!")
  if(any(pop_prev>1))stop("The population prevalence pop_prev must be smaller or equal to 1!")
  
  # Checking that mid_point is valid
  if(class(mid_point) != "numeric" && class(mid_point) != "integer") stop("The mid point mid_point must be numeric!")
  if(mid_point<=0)stop("The mid point mid_point must be positive!")
  
  # Checking that slope is valid
  if(class(slope) != "numeric" && class(slope) != "integer") stop("The slope must be numeric!")
  
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
  if(class(age) != "numeric" && class(age) != "integer") stop("The age must be numeric!")
  if(any(age<0)) stop("The age must be non-negative!")
  
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
    if(class(pop_prev) != "numeric" && class(pop_prev) != "integer") stop("The population prevalence pop_prev must be numeric!")
    if(any(pop_prev<=0))stop("The population prevalence pop_prev must be positive!")
    if(any(pop_prev>1))stop("The population prevalence pop_prev must be smaller or equal to 1!")
    
    # Checking that mid_point is valid
    if(class(mid_point) != "numeric" && class(mid_point) != "integer") stop("The mid point mid_point must be numeric!")
    if(any(mid_point<=0))stop("The mid point mid_point must be positive!")
    
    # Checking that slope is valid
    if(class(slope) != "numeric" && class(slope) != "integer") stop("The slope must be numeric!")
    
    # Computing the threshold
    thresh <- stats::qnorm(pop_prev / (1 + exp((mid_point - age) * slope)), lower.tail = F)
    return(thresh)
    
  }
  
  # if dist = normal, the truncated normal distribution will be used to compute the age of onset
  if(dist == "normal"){
    
    # Checking that min_age and max_age are valid.
    if(class(min_age) != "numeric" && class(min_age) != "integer") stop("The earliest age min_age must be numeric!")
    if(class(max_age) != "numeric" && class(max_age) != "integer") stop("The latest age max_age must be numeric!")
    if(any(min_age <= 0)) stop("The earliest age min_age must be positive!")
    if(any(max_age <= 0)) stop("The latest age max_age must be positive!")
    if(min_age > max_age){
      cat("The latest age max_age is below the earliest age min_age! \n 
The earliest and latest age will be swapped...")
      
      min_age <- min_age + max_age
      max_age <- min_age - max_age
      min_age <- min_age - max_age
    }
    
    # Checking that the lower and upper cutoff points are valid
    if(class(lower) != "numeric" && class(lower) != "integer") stop("The lower cutoff point must be numeric!")
    if(class(upper) != "numeric" && class(upper) != "integer") stop("The upper cutoff point must be numeric!")
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
#' Given a person's cumulative incidence rate (cir), \code{convert_cir_to_age} 
#' can be used to compute the corresponding age, which is given by
#' \deqn{mid_point - log(pop_prev/cir - 1) * 1/slope}
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
  if(class(cir) != "numeric" && class(cir) != "integer") stop("The cumulative incidence rate cir must be numeric!")
  if(cir<=0) stop("The cumulative incidence rate cir must be positive!")
  
  # Checking that pop_prev is valid
  if(class(pop_prev) != "numeric" && class(pop_prev) != "integer") stop("The population prevalence pop_prev must be numeric!")
  if(any(pop_prev<=0))stop("The population prevalence pop_prev must be positive!")
  if(any(pop_prev>1))stop("The population prevalence pop_prev must be smaller or equal to 1!")
  
  # Checking that mid_point is valid
  if(class(mid_point) != "numeric" && class(mid_point) != "integer") stop("The mid point mid_point must be numeric!")
  if(any(mid_point<=0))stop("The mid point mid_point must be positive!")
  
  # Checking that slope is valid
  if(class(slope) != "numeric" && class(slope) != "integer") stop("The slope must be numeric!")
  
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
  if(class(liability) != "numeric" && class(liability) != "integer") stop("The liability must be numeric!")
  
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
    if(class(pop_prev) != "numeric" && class(pop_prev) != "integer") stop("The population prevalence pop_prev must be numeric!")
    if(any(pop_prev<=0))stop("The population prevalence pop_prev must be positive!")
    if(any(pop_prev>1))stop("The population prevalence pop_prev must be smaller or equal to 1!")
    
    # Checking that mid_point is valid
    if(class(mid_point) != "numeric" && class(mid_point) != "integer") stop("The mid point mid_point must be numeric!")
    if(mid_point<=0)stop("The mid point mid_point must be positive!")
    
    # Checking that slope is valid
    if(class(slope) != "numeric" && class(slope) != "integer") stop("The slope must be numeric!")
    
    # Computing the age of onset
    if(any(stats::pnorm(liability, lower.tail = F) >= pop_prev)){
      return(NA)
    }else{
      return(mid_point - log(pop_prev/stats::pnorm(liability, lower.tail = F) - 1)* 1/slope)
    }
  }
  
  # if dist = normal, the truncated normal distribution will be used to compute the age of onset
  if(dist == "normal"){
    
    # Checking that min_aoo and max_aoo are valid.
    if(class(min_aoo) != "numeric" && class(min_aoo) != "integer") stop("The earliest age of onset min_aoo must be numeric!")
    if(class(max_aoo) != "numeric" && class(max_aoo) != "integer") stop("The latest age of onset max_aoo must be numeric!")
    if(any(min_aoo <= 0)) stop("The earliest age of onset min_aoo must be positive!")
    if(any(max_aoo <= 0)) stop("The latest age of onset max_aoo must be positive!")
    if(any(min_aoo > max_aoo)){
      cat("The latest age of onset max_aoo is below the earliest age of onset min_aoo! \n 
The earliest and latest age of onset will be swapped...")
      
      min_aoo <- min_aoo + max_aoo
      max_aoo <- min_aoo - max_aoo
      min_aoo <- min_aoo - max_aoo
    }
    
    # Checking that the lower and upper cutoff points are valid
    if(class(lower) != "numeric" && class(lower) != "integer") stop("The lower cutoff point must be numeric!")
    if(class(upper) != "numeric" && class(upper) != "integer") stop("The upper cutoff point must be numeric!")
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
#' @param obs_h2 A number or numeric vector representing the heritability(ies)
#' on the observed scale. Must be non-negative and at most 1. Defaults to 0.5
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
#' length as obq_h2. Each entry holds to the heritability on the liability
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
#' length as obq_h2. Each entry holds to the heritability on the liability
#' scale which was obtained from the corresponding entry in obs_h2 using Equation 17.
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
  if(class(obs_h2) != "numeric" && class(obs_h2) != "integer") stop("The observed heritability(ies) must be numeric!")
  if(any(obs_h2<0))stop("The observed heritability(ies) must be non-negative!")
  if(any(obs_h2>1))stop("The observed heritability(ies) must be smaller than or equal to one!")
  # Checking that the population prevalences are valid
  if(class(pop_prev) != "numeric" && class(pop_prev) != "integer")stop("The population prevalence(s) must be numeric!")
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
    if(class(prop_cases) != "numeric" && class(prop_cases) != "integer") stop("The proportion(s) of cases must be numeric!")
    if(any(prop_cases<0))stop("The proportion(s) of cases must be non-negative!")
    if(any(prop_cases>1))stop("The proportion(s) of cases must be smaller than or equal to one!")
    
    return(obs_h2 * (pop_prev*(1-pop_prev))/(z^2) * (pop_prev*(1-pop_prev))/(prop_cases*(1-prop_cases)))
  }
}
