#functions to estimate breeding values in a test population

#need estimated effect sizes for loci and gt matrix for test population

get_ebv <- function(ees_table, gt_matrix){
  num_ind <- ncol(gt_matrix)
  alt_matrix <- 2 - gt_matrix
  ebv <- list()
  for (i in 1:num_ind){
    effect_mia <- ees_table$EES.MIA * gt_matrix[, i]
    effect_maa <- ees_table$EES.MAA * alt_matrix[, i]
    ebv[[i]] <- sum(effect_maa, effect_mia)
  }
  return(unlist(ebv))
}

#' est_bv <- testPop$genotype %*% (estEffect * 100) #why multiply by 100?
#' #' Calculate Estimated Breeding Value
#' #' @param gt_subset
#' #' @param my_ees
#' #' @param ind_sample_names
#' #' @param snp_header
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' calculate_ebv <- function(gt, mia_ees, maa_ees, sample_names, snp_header){
#'   gt_rev <- 2 - gt
#'   ebv_mia <- apply(as.matrix(gt), 2, allele_by_ees)
#'   ebv_maa <- apply(as.matrix(gt_rev), 2, allele_by_ees2)
#'   combined_ebv <- ebv_mia + ebv_maa
#'   ebv <- colSums(combined_ebv)
#'   my_ebv <- data.frame(sample_names, ebv)
#'   colnames(my_ebv) <- c("Individual", snp_header)
#'   return(my_ebv)
#' }
