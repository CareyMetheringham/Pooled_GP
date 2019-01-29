#demo pipeline: simulation
gppool_demo <- function(  n_pop = 10,
                          n_ind = 1000,
                          n_site = 1000,
                          h2 = 0.3,
                          MAF = 0.01,
                          threshold = 0.2){
  training_data <-
    produce_sim_data(n_pop, n_ind, n_site, h2, MAF, threshold)
  print(paste(
    "Simulate ",
    n_pop,
    " populations each containing ",
    n_ind,
    " individuals",
    sep = ""
  ))
  print(paste("Use ", n_site, " sites with MAF of ", MAF, sep = ""))
  print(paste(
    "Produces ",
    length(training_data$snp_id),
    " varient sites",
    sep = ""
  ))
  fit_rrblup <- rrblup_loop(training_data)
  ees_table <- create_ees_table(fit_rrblup)
  print("Simulate test population")
  test_data <- sim_test_pop(training_data$es, training_data$af, h2, 200)
  matched_ees <- match_snps_in_ind(ees_table, test_data$gt)
  matched_gt <- get_gt_subset(matched_ees$SNP, test_data$gt)
  print("Estimate breeding values")
  ebv <- get_ebv(matched_ees, matched_gt)
  print("Correlation of EBV and true BV in test population")
  print(cor(ebv, as.vector(test_data$bv)))
  plot(as.vector(test_data$bv), ebv, xlab = "Input Breeding Value", ylab = "Estimated Breeding Value")
  print("Correlation of EBV and observed phenotype in test population")
  print(cor(ebv, as.vector(test_data$ph)))
  plot(as.vector(test_data$ph), ebv, xlab = "Phenotypic Value", ylab = "Estimated Breeding Value")
}

#' Demo using example dataset
#'
#' @param training_snps
#' @param test_snps
#'
#' @return print output to terminal and create boxplot
#' @export
#'
#' @examples
gppool_data_demo <- function(training_snps = 10, test_snps = 6){
  gwas_hits <- "./extdata/test.gwas"
  pools_rc_file <- ("./extdata/test.pool_rc")
  pop_info <- fread("./extdata/test.pool_info")
  ind_info_file <- "./extdata/test.ind_info"
  pool_data <- read_in_pools_rc(pools_rc_file, pop_info, gwas_hits, training_snps)
  fit_rrblup <- rrblup_loop(pool_data)
  ees_table <- create_ees_table(fit_rrblup)
  ind_gt <- read_gt_table("./extdata/test.gt")
  ind_fix <- read_fix_table("./extdata/test.fix")
  matched <- match_and_subset(ees_table, ind_gt, ind_fix, pool_data, test_snps)
  ebv <- get_ebv(matched$ees, matched$gt)
  ind_info <- read_ind_info(ind_info_file)
  accuracy <- calculate_accuracy(create_ebv_table(ind_info, ebv))
  print(paste("Accuracy:", accuracy))
  correlation <- calculate_correlation(create_ebv_table(ind_info, ebv))
  print(paste("Correlation:", correlation))
  ggplot(create_ebv_table(ind_info, ebv), aes(x = group, group = Group, y =EBV)) + geom_boxplot()
}
