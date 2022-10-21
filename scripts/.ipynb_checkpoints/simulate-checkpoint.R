#Simulate object for testing GP prediction

#Libraries
library("argparse")

######## INPUT ARGS ###########
#Define input arguments and parse
parser <- ArgumentParser(description= 'This simulates population for GP')
parser$add_argument('--n_pop', '-p', help= 'Number of populations')
parser$add_argument('--n_ind', '-i', help= 'Number of individuals')
parser$add_argument('--n_site', '-s', help= 'Number of sites')
parser$add_argument('--h2', '-h', help= 'Heritability')
parser$add_argument('--maf', '-m', help= 'MAF')
parser$add_argument('--threshold', '-t', help= 'Threshold')
parser$add_argument('--output1', '-o1', help= 'Simulated training data')
parser$add_argument('--output2', '-o2', help= 'simulated test data')
xargs<- parser$parse_args()

######## FUNCTIONS ###########

#' Generate Distribution of Effect Sizes
#' To be used for training pools and test population
#' @param num_sites total number of num_sites in simulation
#' @return effect size of each locus
generate_effect_size <- function(num_sites) {
  #the effect size of haplotypes at these loci follow a normal distribution
  true_es <- rnorm(num_sites, 0, .1) #why 0.1 too small
  return(true_es)
}

#' Generate Distribution of Allelic Frequencies for Minor Allele
#' @param num_sites total number of num_sites in simulation
#' @return frequency of minor allele at each locus
generate_allelic_freqency <- function(num_sites) {
  #the frequency of alleles varies between loci and follows a beta distribution
  allelic_freq <- rbeta(num_sites, .1, .1) #gives ~40% variable num_sites
  return(allelic_freq)
}

#' Create Loci from Binomial Distribution
#' @param num_sites number od num_sites being simulated
#' @param allelic_freq frequency of minor allele at each locus
#' @return state of each locus
create_loci <- function(num_sites, allelic_freq){
  loci <- rbinom(num_sites, 2, allelic_freq)
  return(loci)
}

#' Simulate num_sites for a Population
#' @param ind number of individuals in the population
#' @param num_sites number od num_sites being simulated
#' @param allelic_freq frequency of minor allele at each locus
#' @return matrix containg loci states in a single population
#' @examples generate_population(10,100, generate_allelic_freqency(100))
generate_population <- function(ind, num_sites, allelic_freq){
  loci_of_pop <- matrix(nrow = ind, ncol = num_sites)
    for (i in 1:ind){
      loci_of_pop[i, ] <- create_loci(num_sites, allelic_freq)
    }
  return(loci_of_pop)
}

#' Calculate Enviromental Variance of BV from h2
#' @param bv vector containg breeding values
#' @param h2 estimate of hertibility
#' @return estimate of varience
calculate_varience <- function(bv, h2){
  est_var <- sd(bv) * sqrt( (1 - h2) / h2)
  return(est_var)
}

#' Get Phenotypic Value
#' @param bv vector containg breeding values
#' @param est_var estimate of varience
#' @return vector containing phenotypes
#' @examples get_phenotype(bv = rnorm(100), calculate_varience(bv = rnorm(100), h2=0.5))
get_phenotype <- function(bv, est_var){
  env_var <- rnorm(length(bv), 0, est_var)
  pheno <- bv + env_var
  return(pheno)
}

#' Simulate Genotype for a Group of Populations
#' @param num_pop number of populations
#' @param num_ind number of individuals in each population
#' @param num_sites number of sites under examination
#' @param allelic_freq frequency of minor allele at each site
#' @return a list of genotype matrices
#' @examples
#' generate_genotypes(num_pop = 10, num_ind = 10, num_sites = 100, generate_allelic_freqency(100))
generate_genotypes <- function(num_pop, num_ind, num_sites, allelic_freq){
  snp_list <- list()
  for (i in 1:num_pop){
    snp_list[[i]] <- generate_population(num_ind, num_sites, allelic_freq)
  }
  return(snp_list)
}

#' Calculate Breeding Values for a group of Populations
#' @param snp_list list of genotype matrices - one per population
#' @param effect_sizes a vector of effect sizes for each site
#' @return list containing vectors of breeding values for each population
get_breeding_values <- function(snp_list, effect_sizes){
  bv_list <- list()
  for (i in 1:length(snp_list)){
    bv_list[[i]] <- snp_list[[i]] %*% effect_sizes
  }
  return(bv_list)
}

#' Simulate the Training Populations
#' @param num_pop number of populations
#' @param num_ind number of individuals in each population
#' @param num_sites number of sites under examination
#' @return a list containing effect size of each site, its allelic frequency, the genotype matrix and a list with
#' vectors of breeding values and phenotypes(continuous)
sim_training_pops <- function(num_pop, num_ind, num_sites, h2 = 0.3){
  effect_sizes <- generate_effect_size(num_sites)
  allelic_freq <- generate_allelic_freqency(num_sites)
  genotypes <- generate_genotypes(num_pop, num_ind, num_sites, allelic_freq)
  genotype_matrix <- make_gt_matrix(genotypes)
  bv <- get_breeding_values(genotypes, effect_sizes)
  phenotypes <- list()
  for ( i in 1:num_pop){
    phenotypes[[i]] <- get_phenotype(bv[[i]], calculate_varience(bv[[i]], h2))
  }
  return(list(es = effect_sizes,
              af = allelic_freq,
              gt_list = genotypes,
              gt_matrix = genotype_matrix,
              bv = bv,
              pt = phenotypes))
}

#' Create a gentype matrix from the geno_list
#' Could potentially be simplified - needs test function
#' @param geno_list a list of genotype matrices
#' @return one matrix containing all genotypes
#' @examples
#' make_gt_matrix(sim_training_pops(5, 10, 10, 0.5)$gt_list)
make_gt_matrix <- function(geno_list){
  pop <- length(geno_list)
  ind <- nrow(geno_list[[1]])
  sites <- ncol(geno_list[[1]])
  geno_matrix <- matrix(nrow = pop * ind, ncol = sites)
  k <- 1
  for (i in 1:pop) {
    geno_matrix[ k:(ind + k - 1), ] <- geno_list[[i]]
    k <- k + ind
  }
  return(geno_matrix)
}

#' Simulate the test population
#' @param es vector containing estimated effect sizes - eg from training pop
#' @param af vector containing allelic frequency - eg from training pop
#' @param h2 heritability
#' @param test_ind number of individuals in test population
sim_test_pop <- function(es, af, h2 = 0.3, test_ind = 100){
  gt <- sim_test_gt(test_ind, af)
  bv <- es %*% gt
  est_var <- calculate_varience(bv, h2)
  ph <- get_phenotype(bv, est_var)
  return(list(gt = gt,
              bv = bv,
              ph = ph
         ))
}

#' Simulate gt matrix for test population
#' @param num_ind number of individuals
#' @param af vector of allelic frequency
#' @return matrix of genotypes: 0, 1, 2
sim_test_gt <- function(num_ind, af){
  num_sites <- length(af)
  test_pop_gt <- t(generate_population(num_ind, num_sites, af))
  sim_snp_id <- paste("snp", 1:num_sites, sep = "_")
  sim_ind_id <- paste("ind", 1:num_ind, sep = "_")
  rownames(test_pop_gt) <- sim_snp_id
  colnames(test_pop_gt) <- sim_ind_id
  return(test_pop_gt)
}

######### CREATE SIMULATAIONS ##########

# Create training data object
training_data <- produce_sim_data(n_pop, n_ind, n_site, h2, MAF, threshold)
# Save training data object
save(training_data, file = xargs$output1)

# Create test data object
test_data <- sim_test_pop(training_data$es, training_data$af, h2, 200)
# Save test data object
save(test_data, file = xargs$output2)