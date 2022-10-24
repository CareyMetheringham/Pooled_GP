An snakemake and R based pipeline for performing genomic prediction on pool-seq data, using rrBLUP, as a follow up to a genome wide association study (GWAS).

The pipeline was a designed and tested on data on ash dieback susceptibility in Fraxinus excelsior.

The pipeline takes the pools_rc file produced by Popoolation2 as input for the pooled training population, and a vcf file of SNPs as input for the test individuals. 

This is a beta version of the pipeline - currently undergoing revisions 

## Getting Started

Paths to the input files are listed in config.yaml. Currently paths are set to example data, update paths to run on own data.

## Input Files

#pool_rc
Ouput file from Popoolation2

#pool info file
File containing info on the pool conditions

#list of snps to use
A list of snps on which to run GP - i.e. a list of gwas hits

#individual info file
File containing info on test individuals

#individual genotype matrix
A matrix containing the genotypes of test individuals

## Output Files

effect size tables written out to data/output/effect_sizes.txt
gebv tables written to data/output/gebv.txt

### Prerequisites

Snakemake: https://snakemake.readthedocs.io/en/stable/

The pipeline installs the following R packages:

argparse

rrBLUP

vcfR

data.table

## Running Snakemake pipeline

Set up conda enviroment with snakemake (see above)

conda activate snakemake
snakemake --use-conda

## Simulations

Scripts for simulation are being updated

## License 

Apache License, Version 2.0
