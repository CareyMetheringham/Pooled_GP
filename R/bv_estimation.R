#functions to estimate breeding values in a test population

#' Calculate EBV from GT and ES
#'
#' @param ees_table data table containing results of rrblup
#' @param gt_matrix matrix of 0, 1, 2 and NA
#'
#' @return ebv vector
#' @export
#'
#' @examples
get_ebv <- function(ees_table, gt_matrix){
  gt_matrix <- gt_matrix[order(rownames(gt_matrix)),]
  ees_table <- ees_table[order(ees_table$SNP), ]
  alt_matrix <- 2 - gt_matrix
  effect_mia <-  t(gt_matrix) %*% ees_table$EES.MIA
  effect_maa <-  t(alt_matrix) %*% ees_table$EES.MAA
  ebv <- effect_mia + effect_maa
  return(as.vector(ebv))
}
