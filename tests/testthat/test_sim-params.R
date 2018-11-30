context("Test Parameters used in Population Simulation")
library(gppool)

test_that("Length allelicF is equal to snps",{
  expect_length(generate_allelic_freqency(1), 1)
  expect_length(generate_allelic_freqency(1000), 1000)
})

test_that("Length es is equal to snps",{
  expect_length(generate_effect_size(1), 1)
  expect_length(generate_effect_size(1000), 1000)
})

test_that("Length loci is equal to snps",{
  expect_length(create_loci(1, generate_allelic_freqency(1)), 1)
  expect_length(create_loci(10000, generate_allelic_freqency(10000)), 10000)
})

test_that("Simulated effect size has expected mean of 0",{
  expect_equal(mean(generate_effect_size(100)), 0, tolerance = 0.1)
  expect_equal(mean(generate_effect_size(1000000)), 0, tolerance = 0.01)
})

test_that("Simulated effect size has expected sd of 0.1",{
  expect_equal(sd(generate_effect_size(100)), 0.1, tolerance = 0.1)
  expect_equal(sd(generate_effect_size(1000000)), 0.1, tolerance = 0.01)
})

test_that("Simulated effect sizes are numeric",{
  expect_that(generate_effect_size(10), is_a("numeric"))
})

test_that("Simulated allelic frequencies are numeric",{
  expect_that(generate_allelic_freqency(10), is_a("numeric"))
})
