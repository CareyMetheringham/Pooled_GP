#' Indentify Variant and Invaraiant Sites
#'
#' @param gt matrix of genotypes for all populations
#' @param MAF minor allele frequency cutoff
#'
#' @return #boolean vector where TRUE is a varient site
#' @export
#'
#' @examples
#' find_varient_sites(sim_training_pops(10, 100, 100, 0.5)$gt_matrix, 0.01)
find_varient_sites <- function(gt, MAF = 0.01) {
  sums <- colSums(gt)
  total_ind <- nrow(gt)
  #minumum nInd of times allele can occur
  min_count <- (total_ind * MAF)
  max_count <- (total_ind * (2 - (MAF)))
  variant <- (sums >= min_count & sums <= max_count)
  return(variant)
}

#' Samples from Population
#'
#' @param pheno matrix of dim 1*n where n is the number of individuals
#' in a population
#' @param threshold the proportion of the population taken from each extreme
#' default set to 0.5 (50:50) split
#'
#' @return two boolean vectors identifying which individuals
#' to use in high and low pools
#' @export
#'
#' @examples
get_samples <- function(pheno, threshold = 0.5){
  hi_pheno <- pheno > quantile(pheno, (1 - threshold))
  lo_pheno <- pheno < quantile(pheno, threshold)
  return(list(hi = hi_pheno,
              lo = lo_pheno))
}

#' Get Pooled Allele Frequency
#'
#' @param gt matrix
#' @param sample boolean list
#' @param variant_sites boolean vector
#'
#' @return mean frequency of the minor allele at each varient site
#' @export
#'
#' @examples
pool_freq <- function(gt, sample, variant_sites){
  freq <- colMeans(gt[unlist(sample) == TRUE, ])
  v_freq <- freq[variant_sites == TRUE]
  return(v_freq)
}

get_pools <- function(sim, MAF = 0.01, threshold = 0.5){
  variant_sites <- find_varient_sites(sim$gt_matrix, MAF)
  num_var_sites <- sum(variant_sites)
  print(num_var_sites)
  num_pop <- length(sim$bv)
  high_pool_matrix <- matrix(nrow = num_pop, ncol= num_var_sites)
  low_pool_matrix <- matrix(nrow = num_pop, ncol= num_var_sites)
  for (i in 1:num_pop){
    sample <- get_samples(sim$pt[[i]])
    high_pool_matrix[i, ] <- pool_freq(sim$gt_list[[i]], sample$hi, variant_sites)
    low_pool_matrix[i, ] <- pool_freq(sim$gt_list[[i]], sample$lo, variant_sites)
  }
  return(list(hi = high_pool_matrix,
              lo = low_pool_matrix))
}
