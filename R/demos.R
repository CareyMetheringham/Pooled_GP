#demo pipeline
#print out results at intervals

gppool_demo <- function(){
  n_pop <- 10
  n_ind <- 1000
  n_site <- 1000
  h2 <- 0.3
  MAF <- 0.01
  threshold <- 0.2
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
  fit_rrblup <- mixed_solve_both(training_data)
  ees_table <- create_ees_table(fit_rrblup)
  test_data <- sim_test_pop(training_data$)
  return(ees_table)
}

gppool_data_demo <- function(){
  gwas_hits <- "./extdata/example_100_hits.gwas"
  pools_rc_files <- find_pools_rc("./extdata/Pools_RC")
  pop_info_file <- "./extdata/example_pop_data.csv"
  pool_data <- read_in_pools_rc(pools_rc_files, pop_info_file, gwas_hits, 50)
  fit_rrblup <- mixed_solve_both(pool_data)
  ees_table <- create_ees_table(fit_rrblup)
  ind_gt <- read_gt_table("./extdata")
  ind_fix <- read_fix_table("./extdata")
  matched <- match_and_subset(ees_table, ind_gt, ind_fix, pool_data, 20)
  ebv <- get_ebv(matched$ees, matched$gt)
  plot(ebv)
}
