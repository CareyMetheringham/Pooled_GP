context("Test Results of Population Simulation")
library(gppool)

test_that("Length of y equals length of prov",{
  expect_equal(length(produce_sim_data(1, 10, 10)$y), length(produce_sim_data(1, 10, 10)$prov))
})

test_that("Nrows in gt matrix = number of snps",{
  sim_data <- produce_sim_data(1, 10, 10)
  expect_equal(nrow(sim_data$mia), length(sim_data$snp_id))
  expect_equal(nrow(sim_data$maa), length(sim_data$snp_id))
})

test_that("Nrows in gt matrix = length y",{
  expect_equal(length(produce_sim_data(1, 10, 10)$y), ncol(produce_sim_data(1, 10, 10)$maa))
  expect_equal(length(produce_sim_data(1, 10, 10)$y), ncol(produce_sim_data(1, 10, 10)$mia))
})

#test range of values in maa and mia
#test that only two values are found in y
