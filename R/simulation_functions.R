#' Generate the distribution of effect sizes and allelic frequency
#' To be used for training pools and test population
#' @param nLoci total number of sites in simulation
#'
#' @return a list cntaining effectSizes and allelicF
#'
#' @examples
generate_effect_size <- function(nLoci) {
  #the effect size of haplotypes at these loci follow a normal distribution
  es <- rnorm(nLoci, 0, .1) #why 0.1 too small
  return(es)
}

generate_allelic_freqency <- function(nLoci) {
  #the frequency of alleles varies between loci and follows a beta distribution
  allelicF <- rbeta(nLoci, .1, .1) #gives ~40% variable sites
  return(allelicF)
}
