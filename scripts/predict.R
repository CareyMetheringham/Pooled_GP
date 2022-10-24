#Genomic prediction from data object

#Libraries
library("argparse")
library("rrBLUP")

######## INPUT ARGS ###########
#Define input arguments and parse
parser <-
  ArgumentParser(description = 'This simulates population for GP')
parser$add_argument('--data', '-d', help = 'RData object')
parser$add_argument('--reps', '-r', help = 'rrBLUP repeats')
parser$add_argument('--out', '-o', help = 'Output')
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
  function(data, rep = 10) {
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

#' Write EES Results to Fileless 
#' @param wd working directory
#' @param fit model fit output to use
create_ees_table <- function(fit){
  ees_table <-
    data.frame(fit$snps, fit$mia$u, fit$mia$u.SE, fit$maa$u, fit$maa$u.SE)
  colnames(ees_table) <-
    c("SNP", "EES.MIA", "EES.MIA.SE", "EES.MAA", "EES.MAA.SE")
  rownames(ees_table) <- fit$snps
  return(ees_table)
}

######## RUN ###########
#load in the pool data object
load(xargs$data)

print(head(poolData))

fit_rrblup <- rrblup_loop(poolData, rep = xargs$reps)
ees_table <- create_ees_table(fit_rrblup)

#Write output as csv table
write.csv(ees_table, file = xargs$out)