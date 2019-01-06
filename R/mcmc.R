#' Calaculate parameter matrix from a set of effect sizes
#'
#' @param simData Simulated data object
#' @param effect Vector of effect sizes - from normal distribution
#'
#' @return
#' @export
#'
#' @examples
get_parameters <- function(simData, effect) {
  p <-
    (simData$mHigh + simData$mLow) / 4 #this uses the mean observed frequency across the two pools
  q <- 1 - p

  exptMean	<- (2 * p ^ 2 + 2 * p * q) * effect

  exptVar	<- (p ^ 2  * (2 * effect - exptMean) ^ 2 +
                2 * p * q * (effect - exptMean) ^ 2 +
                q ^ 2 * exptMean ^ 2)
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
get_exp_var <- function(p, exp_mean){
  exp_var <- (p ^ 2  * (2 * effect - exptMean) ^ 2 +
     2 * p * (1 - p) * (effect - exptMean) ^ 2 +
     q ^ 2 * exptMean ^ 2)
  return(exp_var)
}
