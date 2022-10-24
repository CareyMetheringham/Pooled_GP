An snakemake and R based pipeline for performing genomic prediction on pool-seq data, using rrBLUP, as a follow up to a genome wide association study (GWAS).

The pipeline was a designed and tested on data on ash dieback susceptibility in Fraxinus excelsior.

The pipeline takes the pools_rc file produced by Popoolation2 as input for the pooled training population, and a vcf file of SNPs as input for the test individuals. 

Included are functions for testing the performance of the gppool pipeline on a simulated dataset.

This is a beta version of the pipeline - currently undergoing revisions 

## Getting Started

Paths to the input files are listed in config.yaml. Currently paths are set to example data, update paths to run on own data.

#Run Full Analysis

#Estimate Breeding Values from existing set of Effect Sizes

Once the model has been trained for a dataset, the table of effect sizes can be used to retest the effect size estimates on addditional datasets

### Prerequisites

Conda

Snakemake 

The pipeline installs the following R packages:

rrBLUP

vcfR

data.table

plyr

ggplot2


## License 

Apache License, Version 2.0
