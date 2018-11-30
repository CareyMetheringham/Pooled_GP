sim_test_pop <- function(es, af, h2 = 0.5, test_ind = 100){
  num_sites <- length(es)
  test_pop_gt <- t(generate_population(test_ind, num_sites, af))
  sim_snp_id <- paste("snp", 1:num_sites, sep = "_")
  sim_ind_id <- paste("ind", 1:test_ind, sep = "_")
  rownames(test_pop_gt) <- sim_snp_id
  colnames(test_pop_gt) <- sim_ind_id
  return(test_pop_gt)
}


#'
#' #simulate a test population of individuals
#' testPop <-
#'   simulate_test_pop(simParams, paramList$nLoci, nInd = 150, paramList$h2, simData$fixedSites, paramList$cutoff)
#'
#' #' Simulate Test Population
#' #' @param simParams
#' #' @param nLoci
#' #' @param nInd
#' #' @param h2
#' #' @param fixed
#' #' @return a single test population of nInd individuals
#' #' @export
#' #' @examples
#' simulate_test_pop <- function(simParams, nLoci, nInd, h2, fixed, cutoff) {
#'   lociP <- matrix(nrow = nInd, ncol = nLoci)
#'   for (j in 1:nInd) {
#'     #each loci can take one of three values (0-2)
#'     lociP[j, ] <- create_loci(nLoci, simParams$allelicF)
#'   }
#'   #calculate breeding values
#'   popBV <- (lociP %*% simParams$effectSizes)
#'   #heritability factor used to calculate enviormental varience
#'   env <-
#'     rnorm(length(popBV), 0, calculate_varience(popBV, h2))
#'
#'   phenotype <- popBV + env
#'
#'   topInd <- phenotype > quantile(phenotype, (1 - cutoff))
#'   lowInd <- phenotype < quantile(phenotype, cutoff)
#'
#'   #calculate the difference in alleles between the two pools
#'   High <- lociP[topInd, ]
#'   Low <- lociP[lowInd, ]
#'
#'   return(
#'     list(
#'       genotype = lociP[, fixed == FALSE ],
#'       phenotype = phenotype,
#'       breedingValue = popBV,
#'       highInd = High,
#'       lowInd = Low
#'     )
#'   )
#' }
