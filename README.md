# Gppool

Gppool is an R based pipeline for performing genomic prediction on pool-seq data, using rrBLUP.

The pipeline was a designed and tested on data on ash dieback susceptibility in Fraxinus excelsior.

Included are functions fortesting the performance of the gppool pipeline on a simulated dataset.

This is a beta version of the pipeline

## Getting Started

..

### Prerequisites

The gppool pipeline requires the following R packages to be installed:


rrBLUP

vcfR

data.table

plyr

ggplot2

## Demos and Examples

An example pipeline for using the gppool scripts is included in gppool_example_pipeline

The function gppool_demo() runs a pipeline on a simulated dataset

The funtion gppool_data_demo() runs a test pipeline on a small example dataset 

## Running Tests

Test functions can be run in R using the devtools::test() function

Included tests:

Test Parameters used in Population Simulation

Test Results of Population Simulation

Test Ability to Import Pool-rc Data

Test the Correction of Major and Ref Mismatches

Test Scripts using the Test Dataset

Test MCMC Functions

## License 

