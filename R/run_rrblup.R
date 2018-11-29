#' Find Difference in Allele Frequency from Provinence Mean
#'
#' @param gt_matrix matrix of allele frequences
#' @param prov_list list of provinances
#'
#' @return matrix of deviation from prov means
#' @export
#'
#' @examples
#' get_af_diff(matrix(rnorm(50), 10), c(rep("A", 2), rep("B", 3)))
#' need to check this as appears something odd going on
get_af_diff <- function(gt_matrix, prov_list){
  prov <- rank(prov_list)
  unique_prov <- unique(prov)
  number_of_prov <- length(unique_prov)
  diff_matrix <- gt_matrix
  for (i in 1:number_of_prov){
    current_prov <- unique_prov[i]
    pick <- gt_matrix[, prov == current_prov]
    prov_mean_freq <- rowSums(pick) / ncol(pick)
    freq_diff <- pick - prov_mean_freq
    # store results in correct columns of diff_matrix
    for (j in colnames(pick)){
      diff_matrix[, j] <- freq_diff[, j]
    }
  }
  return(diff_matrix)
}

#' Mixed Solve For Major and Minor Alleles
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' mixed_solve_both(produce_sim_data(10, 100, 1000))
mixed_solve_both <- function(data){
  freq_diff_mia <- data$mia
  freq_diff_maa <- data$maa
  mia_fit <- mixed.solve(data$y, t(freq_diff_mia), SE = TRUE)
  maa_fit <- mixed.solve(data$y, t(freq_diff_maa), SE = TRUE)
  return(list(mia = mia_fit,
              maa = maa_fit,
              snps = data$snp_id))
}
