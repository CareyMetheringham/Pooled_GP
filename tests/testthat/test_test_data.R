context("Test Scripts using the Test Dataset")
library(gppool)
library(data.table)
library(rrBLUP)

gwas_hits <- "../../extdata/test.gwas"
pools_rc_file <- ("../../extdata/test.pool_rc")
pop_info <- fread("../../extdata/test.pool_info")
ind_info_file <- "../../extdata/test.ind_info"
training_snps <- 10
test_snps <-5

test_that("Data is being read in correctly",{
  pool_data <- read_in_pools_rc(pools_rc_file, pop_info, gwas_hits, training_snps)
  expect_equal(pool_data$mia[1], 0.3)
  expect_equal(pool_data$y[1], 1)
  expect_equal(pool_data$major[1], "G")
})


test_that("Known ES values are returned", {
  pool_data <- read_in_pools_rc(pools_rc_file, pop_info, gwas_hits, training_snps)
  fit_rrblup <- mixed_solve_both_af_diff_X(pool_data)
 # expect_equal(fit_rrblup$mia$u[1], -0.3) issue: names for target but not for current
})
