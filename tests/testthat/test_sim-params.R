context("Pool-seq Simulation")
library(gppool)

test_that("Length allelicF is equal to nLoci",{
  expect_equal(length(generate_allelic_freqency(1)), 1)
  expect_equal(length(generate_allelic_freqency(10)), 10)
  expect_equal(length(generate_allelic_freqency(1000)), 1000)
   expect_equal(length(generate_allelic_freqency(1000000)), 1000000)
})

test_that("Length es is equal to nLoci",{
  expect_equal(length(generate_effect_size(1)), 1)
  expect_equal(length(generate_effect_size(10)), 10)
  expect_equal(length(generate_effect_size(1000)), 1000)
  expect_equal(length(generate_effect_size(1000000)), 1000000)
})
