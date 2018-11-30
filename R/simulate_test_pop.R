sim_test_pop <- function(es, af, h2, test_ind){
  num_sites <- length(es)
  test_pop_gt <- generate_population(test_ind, num_sites, af)
}



#simulate a test population of individuals
testPop <-
  simulate_test_pop(simParams, paramList$nLoci, nInd = 150, paramList$h2, simData$fixedSites, paramList$cutoff)

#' Simulate Test Population
#' @param simParams
#' @param nLoci
#' @param nInd
#' @param h2
#' @param fixed
#' @return a single test population of nInd individuals
#' @export
#' @examples
simulate_test_pop <- function(simParams, nLoci, nInd, h2, fixed, cutoff) {
  lociP <- matrix(nrow = nInd, ncol = nLoci)
  for (j in 1:nInd) {
    #each loci can take one of three values (0-2)
    lociP[j, ] <- create_loci(nLoci, simParams$allelicF)
  }
  #calculate breeding values
  popBV <- (lociP %*% simParams$effectSizes)
  #heritability factor used to calculate enviormental varience
  env <-
    rnorm(length(popBV), 0, calculate_varience(popBV, h2))

  phenotype <- popBV + env

  topInd <- phenotype > quantile(phenotype, (1 - cutoff))
  lowInd <- phenotype < quantile(phenotype, cutoff)

  #calculate the difference in alleles between the two pools
  High <- lociP[topInd, ]
  Low <- lociP[lowInd, ]

  return(
    list(
      genotype = lociP[, fixed == FALSE ],
      phenotype = phenotype,
      breedingValue = popBV,
      highInd = High,
      lowInd = Low
    )
  )
}
