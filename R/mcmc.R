#' Calaculate parameter matrix from a set of effect sizes
#'
#' @param simData Simulated data object
#' @param effect Vector of effect sizes - from normal distribution
#'
#' @return
#' @export
#'
#' @examples
#' data <- read_in_pools_rc("./extdata/test.pool_rc", fread("./extdata/test.pool_info"), "./extdata/test.gwas", 10)
#' get_parameters(data, 0.1)
#' sim_data <- produce_sim_data(10, 100, 100)
#' get_parameters(sim_data, 0.1)
get_parameters <- function(simData, effect) {
  p <- rowMeans(simData$mia)

  exptMean	<- get_exp_mean(p, effect)

  exptVar	<- get_exp_var(p, effect, exptMean)

  #calculate the combined effects over loci
  combMean <- sum(exptMean)
  combVar	<- sum(exptVar)

  # calculate the values excluding the focal locus
  exclMean	<- combMean	- exptMean
  exclVar	<- combVar	- exptVar

  return(
    list(
      combMean = combMean,
      combVar = combVar,
      exclMean = exclMean,
      exclVar = exclVar,
      obsFreq = p
    )
  )
}

#' Calculate Expected Mean
#'
#' @param p
#' @param effect
#'
#' @return
#' @export
#'
#' @examples
#' get_exp_mean(0.2, 0.1)
get_exp_mean <- function(p, effect){
  exp_mean <- (2 * p ^ 2 + 2 * p * (1 - p)) * effect
  return(exp_mean)
}

#' Calculate Expected Varience
#'
#' @param p
#' @param exp_mean
#'
#' @return
#' @export
#'
#' @examples
#' get_exp_var(0.2, 0.01, 0.1)
get_exp_var <- function(p, effect, exp_mean){
  exp_var <- (p ^ 2  * (2 * effect - exp_mean) ^ 2 +
     2 * p * (1 - p) * (effect - exp_mean) ^ 2 +
     (1 - p) ^ 2 * exp_mean ^ 2)
  return(exp_var)
}
