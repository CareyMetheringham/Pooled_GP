context("Test MCMC Functions")
library(gppool)
library(data.table)

test_that("Varience Calcs Correct",{
  expect_equal(calc_var_e(0.1, 1), 0)
  expect_equal(calc_var_e(0, 0.5), 0)
})
