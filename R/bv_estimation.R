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
