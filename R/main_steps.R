#demo pipeline
#print out results at intervals

gppool_demo <- function(){
  num_populations <- 10
  ind_in_each <- 1000
  sites_to_simulate <- 1000
  heritability <- 0.3
  training_populations <-
    pool_sim(num_populations, ind_in_each, sites_to_simulate, heritability)
}


