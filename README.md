An snakemake and R based pipeline for performing genomic prediction on pool-seq data, using rrBLUP, as a follow up to a genome wide association study (GWAS).

The pipeline was a designed and tested on data on ash dieback susceptibility in Fraxinus excelsior.

The pipeline takes the pools_rc file produced by Popoolation2 as input for the pooled training population, and a vcf file of SNPs as input for the test individuals. 

This is a beta version of the pipeline - currently undergoing revisions 

## Getting Started

Paths to the input files are listed in config.yaml. Currently paths are set to example data, update paths to run on own data.

## Input Files

#pool_rc

#pool info file

#list of snps to use

#individual info file

#individual genotype matrix

#individual fix file

### Prerequisites

Conda

Snakemake 

The pipeline installs the following R packages:

rrBLUP

vcfR

data.table

## Running Snakemake pipeline

Set up conda enviroment with snakemake (see above)

conda activate snakemake
snakemake --use-conda

## Simulations

Scripts for simulation under construction

## License 

Apache License, Version 2.0
