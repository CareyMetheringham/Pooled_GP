#demo pipeline
#print out results at intervals

gppool_demo <- function(){
  num_populations <- 10
  ind_in_each <- 1000
  sites_to_simulate <- 1000
  heritability <- 0.3
  training_populations <-
    pool_sim(num_populations, ind_in_each, sites_to_simulate, heritability)
  MAF <- 0.01
  threshold <- 0.2
  pools <- get_pools(training_populations, MAF, threshold)
}

gppool_data_demo <- function(){
  top_gwas_hits <- get_hits_from_file("./extdata/example_100_hits.gwas", 50)
  pools_rc_files <- find_pools_rc("./extdata/Pools_RC")
  pop_info_file <- "./extdata/example_pop_data.csv"
  pool_data <- read_in_pools_rc()
  snps_to_use <- find_top_snps(pools_rc_files, top_gwas_hits, pop_info_file)
  major_allele_freq <- get_allele_freq(snps_to_use, pop_info_file, "major")
  minor_allele_freq <- get_allele_freq(snps_to_use, pop_info_file, "minor")
  return(major_allele_freq)
}
