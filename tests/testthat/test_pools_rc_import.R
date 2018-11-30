context("Test Ability to import Pool-rc Data")
library(gppool)
library(data.table)

#Example files to use in tests
pool_dir <- "../../extdata/Pools_RC"
pool_info <- "../../extdata/example_pop_data.csv"
gwas <- "../../extdata/example_100_hits.gwas"

test_that("Pools rc data object as expected",{
  expect_named(read_in_pools_rc(find_pools_rc(pool_dir), pool_info, gwas, 10),
               c("y", "prov", "maa", "mia", "snp_id"))
})

test_that("gt is a matrix",{
  expect_that(read_in_pools_rc(find_pools_rc(pool_dir), pool_info, gwas, 10)$mia, is_a("matrix"))
  expect_that(read_in_pools_rc(find_pools_rc(pool_dir), pool_info, gwas, 10)$maa, is_a("matrix"))
})

test_that("y is an integer",{
  expect_that(read_in_pools_rc(find_pools_rc(pool_dir), pool_info, gwas, 10)$y, is_a("integer"))
})

test_that("Length of y equals length of prov",{
  data <- read_in_pools_rc(find_pools_rc(pool_dir), pool_info, gwas, 10)
  expect_length(data$y, length(data$prov))
})

test_that("Nrows in gt matrix = number of snps",{
 data <- read_in_pools_rc(find_pools_rc(pool_dir), pool_info, gwas, 10)
  expect_length(data$snp_id, nrow(data$mia))
  expect_length(data$snp_id, nrow(data$maa))
})

test_that("Nrows in gt matrix = length y",{
  data <- read_in_pools_rc(find_pools_rc(pool_dir), pool_info, gwas, 10)
  expect_length(data$y, ncol(data$maa))
  expect_length(data$y, ncol(data$mia))
})
