#Simulate object for testing GP prediction

#Libraries
library("argparse")
#Source additional functions
source("/cluster/ggs_lab/cmetheringham001/GP/R/simulate_test_pop.R")
source("/cluster/ggs_lab/cmetheringham001/GP/R/simulate_training_pop.R")

######## INPUT ARGS ###########
#Define input arguments and parse
parser <-
  ArgumentParser(description = 'This simulates population for GP')
parser$add_argument('--n_pop', '-p', help = 'Number of populations')
parser$add_argument('--n_ind', '-i', help = 'Number of individuals')
parser$add_argument('--n_site', '-s', help = 'Number of sites')
parser$add_argument('--h2', '-h2', help = 'Heritability')
parser$add_argument('--maf', '-m', help = 'MAF')
parser$add_argument('--threshold', '-t', help = 'Threshold')
parser$add_argument('--out1', '-o1', help = 'Simulated training data')
parser$add_argument('--out2', '-o2', help = 'simulated test data')
xargs <- parser$parse_args()

######## FUNCTIONS ###########

#' Generate the distribution of effect sizes and allelic frequency
#' To be used for training pools and test population
#' @param nLoci total number of sites in simulation 
#' @return a list containing effectSizes and allelicF - both vectors of length nLoci
#' @examples
#' get_sim_params(nLoci=100)
get_sim_params <- function(nLoci) {
  #the effect size of haplotypes at these loci follow a normal distribution
  es <- rnorm(nLoci, 0, .1)
  #the frequency of alleles varies between loci and follows a beta distribution
  allelicF <- rbeta(nLoci, .1, .1) #gives ~40% variable sites
  return(list(effectSizes = es,
              allelicF = allelicF))
}


#' Create Loci from Binomial Distribution 
#' create_loci(100,rbeta(100, .1, .1))
create_loci <- function(nLoci, allelicF){
  loci <- rbinom(nLoci, 2, allelicF)
  return(loci)
}

#' Create the Genotype and Phenotype of a Single Population
#' create_pop(es=rnorm(100),nInd=10,nLoci=100,allelicF=rbeta(nLoci, .1, .1),h2=0.5)
create_pop <- function(es, nInd, nLoci, allelicF, h2) {
  print(nInd)
  print(nLoci)
  lociP <- matrix(nrow = nInd, ncol = nLoci)
  for (j in 1:nInd) {
    #each loci can take one of three values (0-2)
    lociP[j, ] <- create_loci(nLoci, allelicF)
  }
  #calculate breeding values
  bv <- (lociP %*% es)
  #normalise breeding values
  bv <- scale(bv)
  #heritability factor used to calculate varience
  var <- calculate_varience( bv = bv, h2 = h2 ) 
  envVar <- rnorm( length(bv), 0, var)
  #Store values in lists
  return(list(
    geno = lociP,
    bv = bv,
    pheno = (bv + envVar),
    var = var
  ))
}

#' Create a gentype matrix from the geno_list
#' create_geno_matrix(geno_list=list(rep(1,10)),nPop=1,nInd=1,nLoci=10)
create_geno_matrix <- function(geno_list, nPop, nInd, nLoci){
  genotypes <- matrix(nrow = nPop * nInd, ncol = nLoci)
  k <- 1
  for (i in 1:nPop) {
    geno <- geno_list[[i]]
    genotypes[ k:(nInd + k - 1), ] <- geno
    k <- k + nInd
  }
  return(genotypes)
}

#' Indentify Variant and Invaraiant Sites 
#' find_varient_sites(genotypes=(matrix(rbinom(10 * 5, 1, 0.5), ncol = 5, nrow = 10)), 2, 5, 5)
find_varient_sites <- function(genotypes, nPop, nInd, nLoci, MAF) {
  sums <- colSums(genotypes)
  #minumum nInd of times allele can occur
  min_count <- (nPop * nInd * MAF)
  max_count <- (nPop * nInd * (2 - (MAF)))
  fixed <- (sums <= min_count | sums >= max_count)
  
  return(list(
    fixed = (sums <= min_count | sums >= max_count),
    vLoci = nLoci - sum(fixed)
  ))
}

#' Select the Individuals of Extreme Phenotype from each Population
get_extreme_pheno <- function(pheno_list,i,cutoff){
  topInd <-
    pheno_list[[i]] > quantile(pheno_list[[i]], (1 - cutoff))
  lowInd <- pheno_list[[i]] < quantile(pheno_list[[i]], cutoff)
  return(list(topInd=topInd, lowInd=lowInd))
}

#' Get Frequency of Alleles in the Pools
pool_freq <- function(geno_list, i, sample, fixed){
  freq <- colMeans((geno_list[[i]]) [sample, ])
  var_freq <- freq[fixed == FALSE]
  return(var_freq)
}

#' Create a Matrix containing High and Low Individuals
create_hilo_matrix <-
  function(nPop,
           variant,
           pheno_list,
           geno_list,
           cutoff) {
    
    #find the difference in allelicF between the extremes
    hi <- (matrix(nrow = nPop, ncol = variant$vLoci))
    lo <- (matrix(nrow = nPop, ncol = variant$vLoci))
    
    for (i in 1:nPop) {
      #get the high and low individuals 
      samp <- get_extreme_pheno(pheno_list, i, cutoff)
      
      #find the pooled frequency of the variant sites 
      hi[i, ] <- pool_freq(geno_list, i, samp$topInd, variant$fixed)
      lo[i, ] <- pool_freq(geno_list, i, samp$lowInd, variant$fixed)
    }
    return(list(hi=hi,
                lo=lo))
  }

#' Create a Matrix containing High and Low Individuals
create_hilo_ind_matrix <-  function(nPop,
                                    nInd,
                                    variant,
                                    pheno_list,
                                    geno_list,
                                    cutoff) {
  
  #create matrix of high and low individuals
  hiInd <- matrix(nrow = nPop * nInd * cutoff, ncol = variant$vLoci)
  loInd <- matrix(nrow = nPop * nInd * cutoff, ncol = variant$vLoci)
  #set counter 
  k <- 1
  for (i in 1:nPop) {
    #get the high and low individuals 
    samp <- get_extreme_pheno(pheno_list, i, cutoff)
    
    #store high and low individuals in the matrices
    hiInd[k:(nInd*cutoff + k - 1), ] <- (geno_list[[i]])[samp$topInd, variant$fixed==FALSE]
    loInd[k:(nInd*cutoff + k - 1), ] <- (geno_list[[i]])[samp$lowInd, variant$fixed==FALSE]
    k <- k + nInd*cutoff
  }
  return(list(hiInd=hiInd,
              loInd=loInd))
}

#' Generate Simulated Data for pools from multiple populations
#' using parameters generated in get_sim_params
#' @param nLoci number of sites in each genome under examination
#' @param nPop number of populations
#' @param nInd number of individuals in each population
#' @param h2 heritability factor of the trait
#' @param cutoff pecentage of the population in each extreme - i.e. cutoff=0.5 divides the population at the midpoint
#' @param MAF  minimum allele frequency allowed - eg 0.01 = 1%
simulate_training_pools <-
  function(simParams,
           nLoci = 1000,
           nPop = 10,
           nInd = 1000,
           h2 = 0.5,
           cutoff = 0.5,
           MAF = 0.05) {
    
    #set up lists to fill with values
    geno_list <- list()
    bv_list <- list()
    pheno_list <- list()
    varience_list <- list()
    
    for (i in 1:nPop) {
      population <- create_pop(es=simParams$effectSizes, nInd=nInd, nLoci=nLoci, allelicF=simParams$allelicF, h2=h2)
      geno_list[[i]] <- population$geno
      bv_list[[i]] <-population$bv
      pheno_list[[i]] <- population$pheno
      varience_list[[i]] <- population$var
    }
    
    #create matrix for genotypes
    genotypes <- create_geno_matrix(geno_list, nPop, nInd, nLoci)
    
    #look only at varient sites
    variant <- find_varient_sites(genotypes, nPop, nInd, nLoci,MAF)
    
    #find the difference in allele frequency between the extremes
    HiLo <- create_hilo_matrix(nPop, variant,pheno_list, geno_list, cutoff)
    HiLoInd <- create_hilo_ind_matrix(nPop, nInd, variant, pheno_list, geno_list,cutoff)
    
    mean_high <- colMeans(HiLo$hi)
    mean_low <- colMeans(HiLo$lo)
    
    obsFreq <- colSums(genotypes) / (nInd * nPop) / 2
    
    return(
      list(
        nPop = nPop,
        effectSize = simParams$effectSize[variant$fixed == FALSE],
        mHigh = mean_high,
        mLow = mean_low,
        High = HiLo$hi,
        Low = HiLo$lo,
        obsFreq = obsFreq[variant$fixed == FALSE],
        numVarSites = variant$vLoci,
        fixedSites = variant$fixed,
        Varience = varience_list,
        Genotype = geno_list,
        Phenotype = pheno_list,
        BV =bv_list,
        HiInd = HiLoInd$hiInd,
        LoInd = HiLoInd$loInd
      )
    )
  }

#' Simulate Test Population 
simulate_test_pop <- function(simParams, nLoci, nInd, h2, fixed, cutoff) {
  lociP <- matrix(nrow = nInd, ncol = nLoci)
  for (j in 1:nInd) {
    #each loci can take one of three values (0-2)
    lociP[j, ] <- create_loci(nLoci, simParams$allelicF)
  }
  #calculate breeding values
  popBV <- (lociP %*% simParams$effectSizes)
  #heritability factor used to calculate enviormental varience
  env <-
    rnorm(length(popBV), 0, calculate_varience(popBV, h2))
  
  phenotype <- popBV + env
  
  topInd <- phenotype > quantile(phenotype, (1 - cutoff))
  lowInd <- phenotype < quantile(phenotype, cutoff)
  
  #calculate the difference in alleles between the two pools
  High <- lociP[topInd, ]
  Low <- lociP[lowInd, ]
  
  return(
    list(
      genotype = lociP[, fixed == FALSE ],
      phenotype = phenotype,
      breedingValue = popBV,
      highInd = High,
      lowInd = Low
    )
  )
}

######### CREATE SIMULATAIONS ##########

#Calculate effect size and allelic frequency 
simParams <- get_sim_params(xargs$n_site)

#Create a training data set
simData <- simulate_training_pools(simParams,
           nLoci = as.numeric(xargs$n_site),
           nPop = as.numeric(xargs$n_pop),
           nInd = as.numeric(xargs$n_ind),
           h2 = as.numeric(xargs$h2),
           cutoff = as.numeric(xargs$threshold),
           MAF = as.numeric(xargs$maf)) 
# Save training data object
save(simData, file = xargs$out1)
   
#simulate a test population of individuals
testPop <-
      simulate_test_pop(simParams, nLoci = xargs$n_site, nInd = 150, h2 = xargs$h2, simData$fixedSites, cutoff = xargs$threshold)
# Save test data object
save(testPop, file = xargs$out2)