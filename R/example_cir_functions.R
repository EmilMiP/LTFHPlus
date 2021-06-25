utils::globalVariables("cir")
utils::globalVariables("thr")

# Logistic Function -------------------------------------------------------


#'
#' Calculate cumulative incidence rate based off of a person's age (Used for Examples)
#'
#' @param age age of individual.
#' @param pop_prev Overall prevalence in the population.
#' @param age_mid The mid point of the logistic function
#' @param slope The rate of increase.
#' 
#' @return Returns the estimated genetic liabilities.
#'
#' @examples 
#' # See R/Example/example_nosib.R for an example of use and input.
#' # curve(age_to_cir, from = 10, to = 110, xname = "age")
#' @export

age_to_cir = function(age, pop_prev = .1, age_mid = 60, slope = 1/8) {
  pop_prev / (1 + exp((age_mid - age) * slope))
}

#'
#' Calculate thresholds based off of a person's age (Used for Examples)
#'
#' @param age age of individual.
#' @param pop_prev Overall prevalence in the population.
#' @param age_mid The mid point of the logistic function
#' @param slope The rate of increase.
#'
#' @return Returns the estimated genetic liabilities.
#'
#' @examples 
#' # See R/Example/example_nosib.R for an example of use and input.
#' # curve(age_to_thres, from = 10, to = 110, xname = "age")
#' @export
#' 

age_to_thres = function(age, pop_prev = .1, age_mid = 60, slope = 1/8) {
  stats::qnorm(pop_prev / (1 + exp((age_mid - age) * slope)), lower.tail = F)
}

#'
#' Calculate a person's age based off of the cumulative incidence rate  (Used for Examples)
#'
#' @param cir cumulative incidence rate for an individual
#' @param pop_prev Overall prevalence in the population.
#' @param age_mid The mid point of the logistic function
#' @param slope The rate of increase.
#'
#' @return Returns the estimated genetic liabilities.
#'
#' @examples 
#' # See R/Example/example_nosib.R for an example of use and input.
#' # curve(cir_to_age, from = 0.00, to = .1, xname = "cir")
#' @export
#' 
cir_to_age = function(cir, pop_prev = .1, age_mid = 60, slope = 1/8) {
  age_mid - log(pop_prev/cir - 1) * 1/slope
}




#'
#' Calculate a person's age of onset based off of the true underlying liability (Used for Examples)
#'
#' @param liab true underlying liability of an individual.
#' @param pop_prev Overall prevalence in the population.
#' @param age_mid The mid point of the logistic function
#' @param slope The rate of increase.
#'
#' @return Returns the estimated genetic liabilities.
#'
#' @examples 
#' # See R/Example/example_nosib.R for an example of use and input.
#' # curve(liab_to_aoo, from = 1.3, to = 3.5, xname = "liab")
#' @export

liab_to_aoo = function(liab, pop_prev = .1, age_mid = 60, slope = 1/8) {
  age_mid - log(pop_prev/stats::pnorm(liab, lower.tail = F) - 1)* 1/slope
}



# Truncated Normal Distribution -------------------------------------------

#'
#' Cumulative density function of a truncated normal distribution (Used for Examples)
#'
#' @param liab true underlying liability for an individual
#' @param a the lower cutoff point of the truncated normal distribution
#' @param b the upper cutoff point of the truncated normal distribution
#' 
#' @return Returns the estimated genetic liabilities.
#'
#' @examples 
#' # See R/Example/example_nosib.R for an example of use and input.
#' # curve(trunc_normal_cdf, from = qnorm(0.05, lower.tail = F), to = 3.5, xname = "liab")
#' @export

trunc_normal_cdf = function(liab, a = stats::qnorm(0.05, lower.tail = F), b = Inf) {
  (stats::pnorm(liab) - stats::pnorm(a)) / (stats::pnorm(b) - stats::pnorm(a))
}


#'
#' calculate the age of onset for cases based off of the true underlying liability using the truncated normal distribution (Used for Examples)
#'
#' @param liab true underlying liability for an individual
#' @param min_aoo earliest age of onset to consider
#' @param max_aoo latest age of onset to consider
#' @param a the lower cutoff point of the truncated normal distribution
#' @param b the upper cutoff point of the truncated normal distribution
#' 
#' @return Returns the estimated genetic liabilities.
#'
#' @examples 
#' # See R/Example/example_nosib.R for an example of use and input.
#' # curve(liab_to_aoo_case_trunc_normal, from = qnorm(0.05, lower.tail = F), to = 3.5, xname = "liab")
#' @export

liab_to_aoo_case_trunc_normal = function(liab, min_aoo = 10, max_aoo = 100 - min_aoo, a = stats::qnorm(0.05, lower.tail = F), b = Inf) {
  (1 - trunc_normal_cdf(liab, a = a, b = b)) * max_aoo + min_aoo
}

#'
#' calculate the thresholds of an individual based off of an individual's age using the truncated normal distribution (Used for Examples)
#'
#' @param age true underlying liability for an individual
#' @param min_age minimum age to consider
#' @param max_age maximum age to consider
#' @param a the lower cutoff point of the truncated normal distribution
#' @param b the upper cutoff point of the truncated normal distribution
#' 
#' @return Returns the estimated genetic liabilities.
#'
#' @examples 
#' # See R/Example/example_nosib.R for an example of use and input.
#' # curve(age_to_thres_trunc_normal, from = 10, to = 100, xname = "age")
#' @export
#' 

age_to_thres_trunc_normal = function(age, min_age = 10, max_age = 100 - min_age, a = stats::qnorm(0.05, lower.tail = F), b = Inf) {
  stats::qnorm((1 - (age-min_age)/max_age) * (stats::pnorm(b) - stats::pnorm(a)) + stats::pnorm(a))
}






#'
#' Estimating the cumulative incidence rate in sample
#'
#' @param data tibble containing the status column, age, and id of individuals in "indivs" and "ids.
#' @param indivs vector of prefixes to construct status and age columns from 
#' @param ids unique ids for every individual in indivs. These ids will be the identifier for a given individual in the output.
#' 
#' @return Returns cir for the provided data. Estimated in-sample.
#'
#' @export
#' 

#age of onset to liability. simulated age is age of onset if indiv is a case.
est_cir = function(data,
                   indivs = c("child", "father", "mother"),
                   ids = c("FID", "pid_f", "pid_m")) {
  cat("This function has not been tested properly, use at your own discretion! \n")
  res = dplyr::tibble()
  for(i in seq_along(indivs)) {
    indiv = indivs[i]
    stat_col = paste(indiv, "_stat", sep = "")
    age_col  = paste(indiv, "_age", sep = "")

    if (indivs[i] %in% c("father", "mother")) {
      ph = data %>%
        dplyr::arrange(!!as.name(age_col)) %>%
        dplyr::mutate(cir = (cumsum(!!as.name(stat_col)) + 1)/dplyr::n()) %>%
        dplyr::mutate(thr = stats::qnorm(cir, lower.tail = FALSE)) %>%
        dplyr::select(!!as.name(ids[i]), thr)
    } else {
      sex_col = paste(indiv, "_sex", sep = "")
      ph = data %>%
        dplyr::group_by(!!as.name(sex_col)) %>%
        dplyr::arrange(!!as.name(age_col)) %>%
        dplyr::mutate(cir = (cumsum(!!as.name(stat_col)) + 1)/dplyr::n()) %>%
        dplyr::mutate(thr = stats::qnorm(cir, lower.tail = FALSE)) %>%
        dplyr::ungroup() %>%
        dplyr::select(!!as.name(ids[i]), thr)
    }
    colnames(ph) = c("FID", "thr")
    res = dplyr::bind_rows(res, ph)
  }
  return(res)
}
