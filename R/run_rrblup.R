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

#test scripts need to check lengths

#remove prov effects

#run rrblup with or without repeats

#take input from simulation or pools_rc !!

#produce sanity check boxplots

#print output to files??




# # 6. Get rrBLUP params and check lengths
# pheno <- pool_info$Health
# prov <- rank(pool_info$Prov)
#
# if ( length(pheno) != length(prov)){
#   stop("Pool Number Mismatch - check input files")
# }
#
# # 7. Find divergence from provinence means
# maa_freq_diff <- remove_prov_effects(maa_dec, prov)
# mia_freq_diff <- remove_prov_effects(mia_dec, prov)
#
# # 8. Run rrBLUP on major and minor alleles
# glmm_maa <- mixed.solve(pheno, t(maa_freq_diff), X = prov, SE = TRUE)
# glmm_mia <- mixed.solve(pheno, t(mia_freq_diff), X = prov, SE = TRUE)
# ees_table <-
#   data.table(
#     snp_names,
#     glmm_maa$u,
#     glmm_maa$u.SE,
#     glmm_mia$u,
#     glmm_mia$u.SE,
#     major_allele,
#     minor_allele
#   )
# rownames(ees_table) <- snp_names
# colnames(ees_table) <-
#   c("SNP",
#     "MAA.EES",
#     "MAA.EES.SE",
#     "MIA.EES",
#     "MIA.EES.SE",
#     "MAJOR",
#     "MINOR")
#
# # 9. Print output to file - table and error results etc
# print("Major Allele")
# print(paste("Vu:", glmm_maa$Vu, sep = " "))
# print(paste("Ve:", glmm_maa$Ve, sep = " "))
# print(paste("beta:", glmm_maa$beta, sep = " "))
# print(paste("beta.SE:", glmm_maa$beta.SE, sep = " "))
# print("Minor Allele")
# print(paste("Vu:", glmm_mia$Vu, sep = " "))
# print(paste("Ve:", glmm_mia$Ve, sep = " "))
# print(paste("beta:", glmm_mia$beta, sep = " "))
# print(paste("beta.SE:", glmm_mia$beta.SE, sep = " "))
#
# # 10. Produce boxplot of EBV in training populatio - Sanity Check
# # 1 should be lower than 2
# #mia_ees <- glmm_mia$u
# #maa_ees <- glmm_maa$u
# #ebv_train <-
# #  calculate_ebv(
# #    gt = mia_dec,
# #    mia_ees = glmm_mia$u,
# #    maa_ees = glmm_maa$u,
# #    sample_names = pool_info$Syspop,
# #    snp_header = "EBV"
# #  )
# #training <- data.frame(pool_info, ebv_train)
# #create_boxplot(wd, "training_boxplot.pdf", "Estimated BV in Training Population", training, 5)
#
# #head(training)
# #TG1 <- training[training$Health == 1, ]
# #TG2 <- training[training$Health == 2, ]
# #pdf(paste(wd, "training_boxplot.pdf", sep = "/"))
# #boxplot(TG1$EBV, TG2$EBV,
# #        main = "Estimated BV in Training Population",
# #        ylab = "Estimated Breeding Value")
# #dev.off()
#
# write.table(
#   ees_table,
#   file = paste(wd, "ees.table", sep = "/"),
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )
