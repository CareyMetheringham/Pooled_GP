#' Generate Distribution of Effect Sizes
#' To be used for training pools and test population
#' @param snps total number of sites in simulation
#'
#' @return effect size of each locus
#' @export
#'
#' @examples
generate_effect_size <- function(snps) {
  #the effect size of haplotypes at these loci follow a normal distribution
  true_es <- rnorm(snps, 0, .1) #why 0.1 too small
  return(true_es)
}

#' Generate Distribution of Allelic Frequencies for Minor Allele
#'
#' @param snps total number of sites in simulation
#'
#' @return frequency of minor allele at each locus
#' @export
#'
#' @examples
generate_allelic_freqency <- function(snps) {
  #the frequency of alleles varies between loci and follows a beta distribution
  allelic_freq <- rbeta(snps, .1, .1) #gives ~40% variable sites
  return(allelic_freq)
}

#' Create Loci from Binomial Distribution
#' @param snps number od snps being simulated
#' @param allelic_freq frequency of minor allele at each locus
#' @return state of each locus
#' @export
#' @examples
#' create_loci(100,rbeta(100, .1, .1))
create_loci <- function(snps, allelic_freq){
  loci <- rbinom(snps, 2, allelic_freq)
  return(loci)
}

#' Simulate SNPs for a Population
#'
#' @param ind number of individuals in the population
#' @param snps number od snps being simulated
#' @param allelic_freq frequency of minor allele at each locus
#'
#' @return matrix containg loci states in a single population
#' @export
#'
#' @examples generate_population(10,100, generate_allelic_freqency(100))
generate_population <- function(ind, snps, allelic_freq){
  loci_of_pop <- matrix(nrow = ind, ncol = snps)
    for (i in 1:ind){
      loci_of_pop[i, ] <- create_loci(snps, allelic_freq)
    }
  return(loci_of_pop)
}

#' Calculate Enviromental Variance of BV from h2
#' @param bv vector containg breeding values
#' @param h2 estimate of hertibility
#' @return estimate of varience
#' @export
#' @examples
#' calculate_varience(bv=rnorm(100), h2=0.5)
calculate_varience <- function(bv, h2){
  est_var <- sd(bv) * sqrt( (1 - h2) / h2)
  return(est_var)
}

#' Get True Breeding Value
#'
#' @param loci_of_pop matrix containg loci states in a single population
#' @param true_es true effect size of each locus
#'
#' @return true breeding value of each individual in the population
#' @export
#'
#' @examples get_true_bv(generate_population(10, 10, generate_allelic_freqency(10)), generate_effect_size(10))
get_true_bv <- function(loci_of_pop, true_es){
  true_bv <- loci_of_pop %*% true_es
  scale_bv <- scale(true_bv)
  return(scale_bv)
}

#' Get Phenotypic Value
#'
#' @param bv vector containg breeding values
#' @param est_var estimate of varience
#'
#' @return vector containing phenotypes
#' @export
#'
#' @examples get_phenotype(bv = rnorm(100), calculate_varience(bv = rnorm(100), h2=0.5))
get_phenotype <- function(bv, est_var){
  env_var <- rnorm(length(bv), 0, est_var)
  pheno <- bv + env_var
  return(pheno)
}

simulate_training_pools <- function(){

}
