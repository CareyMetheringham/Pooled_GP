#' Simulate the test population
#'
#' @param es vector containing estimated effect sizes - eg from training pop
#' @param af vector containing allelic frequency - eg from training pop
#' @param h2 heritability
#' @param test_ind number of individuals in test population
#'
#' @return
#' @export
#'
#' @examples
#' train <- sim_training_pops(10, 100, 100, 0.5)
#' sim_test_pop(train$es, train$af)
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
#'
#' @param num_ind number of individuals
#' @param af vector of allelic frequency
#'
#' @return matrix of genotypes: 0, 1, 2
#' @export
#'
#' @examples
#' sim_test_gt(100, generate_allelic_freqency(100))
sim_test_gt <- function(num_ind, af){
  num_sites <- length(af)
  test_pop_gt <- t(generate_population(num_ind, num_sites, af))
  sim_snp_id <- paste("snp", 1:num_sites, sep = "_")
  sim_ind_id <- paste("ind", 1:num_ind, sep = "_")
  rownames(test_pop_gt) <- sim_snp_id
  colnames(test_pop_gt) <- sim_ind_id
  return(test_pop_gt)
}
