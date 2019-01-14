context("Test Ability to Import Pool-rc Data")
library(gppool)
library(data.table)

#Example files to use in tests
pool_rc <- "../../extdata/test.pool_rc"
pool_info <- "../../extdata/test.pool_info"
gwas <- "../../extdata/test.gwas"

test_that("Pools rc data object as expected", {
  expect_named(read_in_pools_rc(pool_rc, fread(pool_info), gwas, 10),
               c("y", "prov", "maa", "mia", "snp_id", "major", "minor"))
})

test_that("gt is a matrix", {
  expect_that(read_in_pools_rc(pool_rc,
                               fread(pool_info), gwas, 10)$mia,
              is_a("matrix"))
  expect_that(read_in_pools_rc(pool_rc,
                               fread(pool_info), gwas, 10)$maa,
              is_a("matrix"))
})

test_that("y is an integer", {
  expect_that(read_in_pools_rc(pool_rc,
              fread(pool_info), gwas, 10)$y,
              is_a("integer"))
})

test_that("Length of y equals length of prov", {
  data <- read_in_pools_rc(pool_rc, fread(pool_info), gwas, 10)
  expect_length(data$y, length(data$prov))
})

test_that("Nrows in gt matrix = number of snps", {
 data <- read_in_pools_rc(pool_rc, fread(pool_info), gwas, 10)
  expect_length(data$snp_id, nrow(data$mia))
  expect_length(data$snp_id, nrow(data$maa))
})

test_that("Nrows in gt matrix = length y", {
  data <- read_in_pools_rc(pool_rc, fread(pool_info), gwas, 10)
  expect_length(data$y, ncol(data$maa))
  expect_length(data$y, ncol(data$mia))
})

test_that("maa freq + mia freq = 1", {
  data <- read_in_pools_rc(pool_rc, fread(pool_info), gwas, 10)
  expect_equal((data$maa[1,1] + data$mia[1,1]), 1)
})
# This needs generalisation
