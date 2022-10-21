#Genomic prediction from data object

#Libraries
library("argparse")
library("rrBLUP")

######## INPUT ARGS ###########
#Define input arguments and parse
parser <-
  ArgumentParser(description = 'This simulates population for GP')
parser$add_argument('--sim_training', '-tr', help = 'Simulated training data object')
parser$add_argument('--sim_test', '-te', help = 'Simulated test data object')
xargs <- parser$parse_args()

######## FUNCTIONS ###########

#' Run rrBLUP as loop
#' @param data object containing data from training population
#' @param X Boolean - should attempt be made to fit X - default FALSE
#' @param Diff Should the difference from Group2 mean be used in place of raw values - default = TRUE
#' @param both Should both minor and msjor alleles be used - default = TRUE
#' @param rep How many times should model fitting be repeated
#' @return list containing vectors with estimated effect sizes and snp_id
rrblup_loop <-
  function(data,
           X = FALSE,
           Diff = TRUE,
           both = TRUE,
           rep = 10) {
    #need to add other outputs to the data structure?
    mia_u_df <- data.frame()
    maa_u_df <- data.frame()
    mia_u_SE_df <- data.frame()
    maa_u_SE_df <- data.frame()
    for (i in 1:rep) {
      freq_diff_mia <- get_af_diff(data$mia, data$prov)
      freq_diff_maa <- get_af_diff(data$maa, data$prov)
      mia_fit <- mixed.solve(data$y, t(freq_diff_mia), SE = TRUE)
      maa_fit <- mixed.solve(data$y, t(freq_diff_maa), SE = TRUE)
    }
    
    mia_u_df <- rbind(mia_u_df, mia_fit$u)
    maa_u_df <- rbind(maa_u_df, maa_fit$u)
    mia_u_SE_df <- rbind(mia_u_SE_df, mia_fit$u.SE)
    maa_u_SE_df <- rbind(maa_u_SE_df, maa_fit$u.SE)
  }
#get the averages
mean_mia_u <- colMeans(mia_u_df)
mean_maa_u <- colMeans(maa_u_df)
mean_mia_u_SE <- colMeans(mia_u_SE_df)
mean_maa_u_SE <- colMeans(maa_u_SE_df)
names(mean_mia_u) <- data$snp_id
names(mean_maa_u) <- data$snp_id
#return mean
return(list(
  mia = list(u = mean_mia_u, u.SE = mean_mia_u_SE),
  maa = list(u = mean_maa_u, u.SE = mean_maa_u_SE),
  snps = data$snp_id
))
}

#' Find Difference in Allele Frequency from Provinence Mean
#' @param gt_matrix matrix of allele frequences
#' @param prov_list list of provinances
#' @return matrix of deviation from prov means
get_af_diff <- function(gt_matrix, prov_list) {
  prov <- rank(prov_list)
  unique_prov <- unique(prov)
  number_of_prov <- length(unique_prov)
  diff_matrix <- gt_matrix
  for (i in 1:number_of_prov) {
    current_prov <- unique_prov[i]
    pick <- gt_matrix[, prov == current_prov]
    prov_mean_freq <- rowSums(pick) / ncol(pick)
    freq_diff <- pick - prov_mean_freq
    # store results in correct columns of diff_matrix
    for (j in colnames(pick)) {
      diff_matrix[, j] <- freq_diff[, j]
    }
  }
  return(diff_matrix)
}

######## RUN ###########

training_data <- load(xargs$)
test_data <- load(xargs$)

fit_rrblup <- rrblup_loop(training_data)
ees_table <- create_ees_table(fit_rrblup)


matched_ees <- match_snps_in_ind(ees_table, test_data$gt)
matched_gt <- get_gt_subset(matched_ees$SNP, test_data$gt)
print("Estimate breeding values")
ebv <- get_ebv(matched_ees, matched_gt)
print("Correlation of EBV and true BV in test population")
print(cor(ebv, as.vector(test_data$bv)))
plot(as.vector(test_data$bv),
     ebv,
     xlab = "Input Breeding Value",
     ylab = "Estimated Breeding Value")
print("Correlation of EBV and observed phenotype in test population")
print(cor(ebv, as.vector(test_data$ph)))
plot(as.vector(test_data$ph),
     ebv,
     xlab = "Phenotypic Value",
     ylab = "Estimated Breeding Value")
}