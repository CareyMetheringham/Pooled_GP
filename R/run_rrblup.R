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
mixed_solve_both_af_diff_X <- function(data){
  freq_diff_mia <- get_af_diff(data$mia, data$prov)
  freq_diff_maa <- get_af_diff(data$maa, data$prov)
  mia_fit <- mixed.solve(data$y, X = rank(data$prov), t(freq_diff_mia), SE = TRUE)
  maa_fit <- mixed.solve(data$y, X = rank(data$prov), t(freq_diff_maa), SE = TRUE)
  return(list(mia = mia_fit,
              maa = maa_fit,
              snps = data$snp_id))
}


#' Mixed Solve Without Prov Effects Major and Minor Alleles
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' mixed_solve_both_af_diff(produce_sim_data(10, 100, 1000))
mixed_solve_both_af_diff <- function(data){
  freq_diff_mia <- get_af_diff(data$mia, data$prov)
  freq_diff_maa <- get_af_diff(data$maa, data$prov)
  mia_fit <- mixed.solve(data$y, t(freq_diff_mia), SE = TRUE)
  maa_fit <- mixed.solve(data$y, t(freq_diff_maa), SE = TRUE)
  return(list(mia = mia_fit,
              maa = maa_fit,
              snps = data$snp_id))
}

#' Mixed Solve Without Prov Effects Major and Minor Alleles
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' mixed_solve_both_af_diff(produce_sim_data(10, 100, 1000))
ml_both_af_diff <- function(data){
  freq_diff_mia <- get_af_diff(data$mia, data$prov)
  freq_diff_maa <- get_af_diff(data$maa, data$prov)
  mia_fit <- mixed.solve(data$y, t(freq_diff_mia), SE = TRUE, method = "ML")
  maa_fit <- mixed.solve(data$y, t(freq_diff_maa), SE = TRUE, method = "ML")
  return(list(mia = mia_fit,
              maa = maa_fit,
              snps = data$snp_id))
}

#' Mixed Solve Without Prov Effects Major and Minor Alleles
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' mixed_solve_both_af_diff(produce_sim_data(10, 100, 1000))
ml_both_af_diff_X <- function(data){
  freq_diff_mia <- get_af_diff(data$mia, data$prov)
  freq_diff_maa <- get_af_diff(data$maa, data$prov)
  mia_fit <- mixed.solve(data$y, t(freq_diff_mia), X = rank(data$prov), SE = TRUE, method = "ML")
  maa_fit <- mixed.solve(data$y, t(freq_diff_maa), X = rank(data$prov), SE = TRUE, method = "ML")
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

#' Mixed Solve For Major and Minor Alleles
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' mixed_solve_both(produce_sim_data(10, 100, 1000))
mixed_solve_both_with_X <- function(data){
  freq_diff_mia <- data$mia
  freq_diff_maa <- data$maa
  mia_fit <- mixed.solve(data$y, t(freq_diff_mia), X = rank(data$prov), SE = TRUE)
  maa_fit <- mixed.solve(data$y, t(freq_diff_maa), X = rank(data$prov), SE = TRUE)
  return(list(mia = mia_fit,
              maa = maa_fit,
              snps = data$snp_id))
}

rrblup_loop <- function(data, X = TRUE, Diff = TRUE, both = TRUE, rep = 10){
  #need to add other outputs to the data structure?
  mia_u_df <- data.frame()
  maa_u_df <- data.frame()
  mia_u_SE_df <- data.frame()
  maa_u_SE_df <- data.frame()
  for (i in 1:rep){
    freq_diff_mia <- get_af_diff(data$mia, data$prov)
    freq_diff_maa <- get_af_diff(data$maa, data$prov)
    mia_fit <- mixed.solve(data$y, X = rank(data$prov), t(freq_diff_mia), SE = TRUE)
    maa_fit <- mixed.solve(data$y, X = rank(data$prov), t(freq_diff_maa), SE = TRUE)
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
  #return mean
  return(list(mia = list(u = mean_mia_u, u.SE = mean_mia_u_SE),
              maa = list(u = mean_maa_u, u.SE = mean_maa_u_SE),
              snps = data$snp_id
              ))
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

#' Find Difference from Provenance Mean Frequency
#' @param maa_dec
#' @param prov
#'
#' @return
#' @export
#'
#' @examples
remove_prov_effects <- function(dec, prov){
  uni_prov <- unique(prov)
  prov_diff <- dec
  for (i in 1:length(uni_prov)){
    current_prov <- uni_prov[i]
    correct_prov <- prov == current_prov
    pick <- dec[, correct_prov == TRUE]
    av <- rowSums(pick) / ncol(pick)
    diff <- (pick - av)
    for (j in colnames(diff)){
      prov_diff[, j] <- diff[, j]
    }
  }
  return(prov_diff)
}
