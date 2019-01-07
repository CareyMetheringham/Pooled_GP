#' Calaculate parameter matrix from a set of effect sizes
#'
#' @param simData Simulated data object
#' @param effect Vector of effect sizes - from normal distribution
#'
#' @return
#' @export
#'
#' @examples
#' test_data <- read_in_pools_rc("./extdata/test.pool_rc", fread("./extdata/test.pool_info"), "./extdata/test.gwas", 10)
#' get_parameters(test_data, 0.1)
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
      obsFreq = as.vector(p)
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

run_mcmc <- function(my_data, effect, chain = 10^3){
  for (i in 1:chain){
    params <- get_parameters(my_data, effect)
    l <- sample(1:length(my_data$snp_id), 1)
    e <- -1000:1000 / 200
    f <- params$obsFreq[[l]]
    m <- params$exclMean[[l]]
    var_g <- params$exclVar[[l]]
    var_e <- calc_var_e(var_g)
    sd_combined <- sqrt(var_g + var_e)
    p_high <- calc_p_high(f, m, sd_combined, e)
    obs_high <- get_obs_high(my_data)[l]
    #stopgap use of whole numbers for likelihood function
    num_total <- round(f * 1000)
    num_high <- round(obs_high * 1000)
    likelihood <- dbinom(num_high, num_total, p_high)
    #diffuse prior?
    prior <- rbeta(length(p_high), 1, 1)
    #calculate posterior
    post <- prior * likelihood
    x <- which.max(post)
    effect[l] <- e[x]
  }
  return(effect)
}

#need to be able to vary cutoff value??
calc_p_high <- function(f, m, sd_comb, e){
  p_high <- f * pnorm(m + 2 * e, 0.5, sd_comb) + (1 - f) * pnorm(m + e, 0.5, sd_comb)
  return(p_high)
}

#should this be different??? wrong equation
calc_p_low <- function(f, m, sd_comb, e){
  p_high <- f * pnorm(m + 2 * e, 0.5, sd_comb) + (1 - f) * pnorm(m + e, 0.5, sd_comb)
  return(p_high)
}

get_obs_high <- function(my_data, high = '2'){
  mia_high <- my_data$mia[, my_data$y == high]
  mean_high <- rowMeans(mia_high)
  return(mean_high)
}

#' Calculate Enviromental Varience on Assumed h2
#'
#' @param var_g
#' @param h2
#'
#' @return
#' @export
#'
#' @examples
#' calc_var_e(0.1)
#' calc_var_e(0.1, 1) # should give 0 write test
calc_var_e <- function(var_g, h2 = 0.3){
  var_e <- (var_g * ( 1 - h2 )) / h2
  return(var_e)
}


# run_mcmc <-
#   function(chainLength,
#            simData,
#            effect,
#            nInd,
#            nPop,
#            checkpoint) {
#     checkpoints <- seq(from = 1,
#                        to = chainLength,
#                        length.out = checkpoint)
#     rValues <- list()
#     c <- 1
#     for (i in 1:chainLength) {
#       params <-
#         get_parameters(simData, effect) #generate parameters for the starting vector of random effect sizes
#       #select one locus at random
#       l <-
#         sample(1:simData$numVarSites, 1) #pick one locus at random from the varient loci
#       #e <- -nLoci:nLoci / 2 #a vector of proposed effect sizes - THIS NEEDS CHANGING!!!!!
#       e <-
#         -100:100 / 200
#
#       f <- params$obsFreq[[l]] # observed frequency of the allele
#       q <- 1 - f #observed frequency of the alternate allele
#       m <-
#         params$exclMean[[l]] #estimated average effect of all other alleles (why does this appear fixed?)
#       Vg <-
#         params$exclVar[[l]]  #genetic varience of all other alleles
#
#       Ve <- (Vg*(1-h2))/h2 #calculate Ve based on assumed h2
#
#       sdCombined <-
#         sqrt(Vg + Ve) # need to add enviromental varience
#
#       #probability of the allele being in the upper pool
#       pHi <-
#         f * pnorm(m + 2 * e, 0.5, sdCombined) + q * pnorm(m + e, 0.5, sdCombined)
#       pLo <- f * pnorm(m + 2 * e, cutoff , sdCombined) + q * pnorm(m + e, cutoff, sdCombined)
#
#       #pHi/pLo ??
#
#       #observed frequency of the allele in the upper pool
#       obsHi <- simData$mHigh[[l]] / 2
#       obsLo <- simData$mLow[[l]] / 2
#
#       n <-
#         f * nPop * nInd * 4 #number of times the allele is observed
#       y <-obsHi * nPop * nInd * 2
#
#       #y <- obsHi * nPop * nInd * 2 #number of times the allele is observed in the upper pool
#
#       #liklihood function
#       liklihood <-
#         dbinom(round(y), round(n), pHi) #does not take non int values
#
#       #prior ??? - needs work
#       prior <- rbeta(length(pHi), 1, 1)
#
#       #posterior?
#       post <- prior * liklihood
#       #plot(post~e,pch="*")
#
#       x <- which.max(post)
#       effect[l] <-
#         e[x]
#
#       if (i %in% checkpoints) {
#         r <- get_r_squared(simData, effect)
#         rValues[[c]] <- r
#         info <- paste("Checkpoint",c, "- Iteration", i, "- r2 =", r, sep = " ")
#         print(info)
#         c = c + 1
#       }
#     }
#     return(list(estEffect = effect/nLoci,
#                 rValues = rValues))
#   }
