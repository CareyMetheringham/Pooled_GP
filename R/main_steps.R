#demo pipeline
#print out results at intervals

gppool_demo <- function(){
  num_populations <- 10
  ind_in_each <- 1000
  sites_to_simulate <- 1000
  heritability <- 0.3
  MAF <- 0.01
  threshold <- 0.2
  training_data <-
    produce_sim_data(num_populations, ind_in_each, sites_to_simulate, heritability, MAF, threshold)
  #issue - mia contains NaN!!
}

gppool_data_demo <- function(){
  gwas_hits <- "./extdata/example_100_hits.gwas"
  pools_rc_files <- find_pools_rc("./extdata/Pools_RC")
  pop_info_file <- "./extdata/example_pop_data.csv"
  pool_data <- read_in_pools_rc(pools_rc_files, pop_info_file, gwas_hits, 50)
  return(pool_data)
}
