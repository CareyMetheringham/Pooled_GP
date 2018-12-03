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

#' Mixed Solve Without Prov Effects Major and Minor Alleles
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' mixed_solve_both(produce_sim_data(10, 100, 1000))
mixed_solve_both_af_diff <- function(data){
  freq_diff_mia <- get_af_diff(data$mia, data$prov)
  freq_diff_maa <- get_af_diff(data$maa, data$prov)
  mia_fit <- mixed.solve(data$y, t(freq_diff_mia), SE = TRUE)
  maa_fit <- mixed.solve(data$y, t(freq_diff_maa), SE = TRUE)
  return(list(mia = mia_fit,
              maa = maa_fit,
              snps = data$snp_id))
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

#' Write EES Results to File
#'
#' @param wd working directory
#' @param fit model fit output to use
#'
#' @return a table printed to file ees.table in named working dir
#' # SNP	MAA.EES	MAA.EES.SE	MIA.EES	MIA.EES.SE
#' @export
#'
#' @examples
#' create_ees_table(mixed_solve_both(produce_sim_data(10, 100, 100)))
create_ees_table <- function(fit){
  ees_table <-
    data.frame(fit$snps, fit$mia$u, fit$mia$u.SE, fit$maa$u, fit$maa$u.SE)
  colnames(ees_table) <-
    c("SNP", "EES.MIA", "EES.MIA.SE", "EES.MAA", "EES.MAA.SE")
  return(ees_table)
}
