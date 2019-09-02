#Simulation and testing of GP from pool-seq data
#Carey Metheringham 

#requires library rrBLUP to perform modelling 
library(rrBLUP)

check_input_values <- function(paramList){
  if (paramList$h2 > 1){
    print("Warning! h2 > 1")
  }
}

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

#' Calculate Enviromental Variance of BV from h2
#' @param bv 
#' @param h2 
#' @return 
#' @export
#' @examples
#' calculate_varience(bv=rnorm(100),h2=0.5)
calculate_varience <- function(bv, h2){
  Var <- sd(bv) * sqrt( (1 - h2) / h2)
  return(Var)
}

#' Create Loci from Binomial Distribution 
#' @param nLoci 
#' @param allelicF 
#' @return
#' @export
#' @examples
#' create_loci(100,rbeta(100, .1, .1))
create_loci <- function(nLoci, allelicF){
  loci <- rbinom(nLoci, 2, allelicF)
  return(loci)
}

#' Create the Genotype and Phenotype of a Single Population
#' @param es 
#' @param nInd 
#' @param nLoci 
#' @param allelicF 
#' @param h2 
#'
#' @return
#' @export
#'
#' @examples
#' create_pop(es=rnorm(100),nInd=10,nLoci=100,allelicF=rbeta(nLoci, .1, .1),h2=0.5)
create_pop <- function(es, nInd, nLoci, allelicF, h2) {
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
#' Could potentially be simplified - needs test function 
#' @param geno_list 
#' @param nPop 
#' @param nInd 
#' @param nLoci 
#'
#' @return
#' @export
#'
#' @examples
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
#'
#' @param genotypes # a matrix 
#' @param nPop 
#' @param nInd 
#' @param nLoci 
#'
#' @return #list containing fixed sites (boolean vector) and number of varient sites (int)
#' @export
#'
#' @examples
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
#' @param pheno_list 
#' @param i 
#' @param cutoff 
#' @return
#' @export
#' @examples
get_extreme_pheno <- function(pheno_list,i,cutoff){
  topInd <-
    pheno_list[[i]] > quantile(pheno_list[[i]], (1 - cutoff))
  lowInd <- pheno_list[[i]] < quantile(pheno_list[[i]], cutoff)
  return(list(topInd=topInd, lowInd=lowInd))
}

#' Get Frequency of Alleles in the Pools
#' @param geno_list 
#' @param i 
#' @param sample 
#' @param fixed 
#' @return
#' @export
#' @examples
pool_freq <- function(geno_list, i, sample, fixed){
  freq <- colMeans((geno_list[[i]]) [sample, ])
  var_freq <- freq[fixed == FALSE]
  return(var_freq)
}

#' Create a Matrix containing High and Low Individuals
#' @param nPop 
#' @param nInd 
#' @param variant 
#' @param pheno_list 
#' @param geno_list 
#' @return
#' @export
#' @examples
#' 
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
#'
#' @return
#' @export
#'
#' @examples
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
      population <- create_pop(es=simParams$effectSizes, nInd, nLoci, allelicF=simParams$allelicF, h2)
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

#' Calculate Expected Mean 
#' @param p 
#' @param effect 
#' @return
#' @export
#' @examples
#' calculate_expected_mean(0.05, 0.1)
calculate_expected_mean <- function(p, effect){
  q <- 1 - p
  exptMean	<- (2 * p ^ 2 + 2 * p * q) * effect
  return(exptMean)
}

#' Calculate Expected Variance 
#' @param p 
#' @param effect 
#' @param exptMean 
#' @return
#' @export
#' @examples
calculate_expected_var <- function(p, effect, exptMean){
  q <- 1 - p
  exptVar	<- (p ^ 2  * (2 * effect - exptMean) ^ 2 +
                2 * p * q * (effect - exptMean) ^ 2 +
                q ^ 2 * exptMean ^ 2)
  return(exptVar)
}

#' Calaculate parameter matrix from a set of effect sizes
#' @param simData Simulated data object
#' @param effect Vector of effect sizes - from normal distribution
#' @return
#' @export
#' @examples
get_parameters <- function(simData, effect) {
  
  #mean observed frequency across the two pools
  p <- (simData$mHigh + simData$mLow) / 4 
  
  #calculate expected mean and variance 
  exptMean <- calculate_expected_mean(p, effect)
  exptVar <- calculate_expected_var(p, effect, exptMean)
  
  #calculate the combined effects over loci
  combMean <- sum(exptMean)
  combVar	<- sum(exptVar)
  
  # calculate the values excluding the focal locus
  exclMean	<- combMean	- exptMean
  exclVar	<- combVar	- exptVar
  
  return(
    list(
      combMean = combMean,
      combVar = combVar,
      exclMean = exclMean,
      exclVar = exclVar,
      obsFreq = p
    )
  )
}

#' Use rrBLUP to estimate effect sizes 
#' From pooled data 
#' @param simData 
#' @return
#' @export
#' @examples
run_rrBLUP <- function(simData, nPop) {
  y <- c(rep(1,nPop),rep(0,nPop))
  HiLo <- rbind(simData$High, simData$Low)
  #use rrBLUP to find most likley effect sizes to use as start point 
  glmm_start <- mixed.solve(y,HiLo)
  startEffects <-glmm_start$u
  return(startEffects)
}

#' Get the adjusted r squared for the checkpoint
#' @param simData
#' @param estEffect
#' @return adjusted r squared value
get_r_squared <- function(simData, estEffect) {
  mod3 <- lm(simData$effectSize ~ estEffect)
  mod3_summary <- summary(mod3)
  return(mod3_summary$adj.r.squared)
}

#' Calculat expected enviromental varience 
#' @param Vg 
#' @param h2 
#' @return
#' @export
#' @examples
getVe <- function(Vg,h2){
  Ve <- (Vg*(1-h2))/h2 #calculate Ve based on assumed h2
  return(Ve)
}

#' Calculate the probability of the allele occuring in the upper pool
#' @param f = observed frequency of allele
#' @param m  = estimated average effect of all other alleles
#' @param e = effect of allele
#' @param sd = standard deviation 
#' @return
#' @export
#' @examples
#' calc_pHi(f=0.5,m=0,e=0,sd=1) 
#' calc_pHi(f=1,m=0,e=0,sd=1) <- should be 1????????
#' calc_pHi(f=0.5,m=0,e=0.5,sd=1)
calc_pHi <- function(f,m,e,sd){
  q=1-f #q = observed frequency of alternate allele
  pHi <- f * pnorm(m + 2 * e, 0.5, sd) + q * pnorm(m + e, 0.5, sd)
  return(pHi)
}

#' Calculate probability of allele occuring in lower pool ...
#' @param f = observed frequency of allele
#' @param m  = estimated average effect of all other alleles
#' @param e = effect of allele
#' @param sd = standard deviation 
#' @return
#' @export
#' @examples
#' calc_pLo(f=0.5, m=0, e=0, sd=1) 
#' calc_pLo(f=1, m=0, e=0, sd=1)
#' calc_pLo(f=0.5, m=0, e=0.5, sd=1)
calc_pLo <- function(f, m, e, sd){
  q <- 1 - f #q = observed frequency of alternate allele
  pLo <-  f * pnorm(m + 2 * e, paramList$cutoff, sd) + q * pnorm(m + e, paramList$cutoff, sd)
  return(pLo)
}

find_posterior <- function(y, n, pHi){
  #liklihood function
  liklihood <- dbinom(round(y), round(n), pHi) #does not take non int values
  #diffuse prior 
  prior <- rbeta(length(pHi), 1, 1)
  post <- prior * liklihood
  return(post)
}

#' Run an MCMC to itteratively improve effect size estimates 
#' @param chainLength
#' @param simData
#' @param effect
#' @param nInd
#' @param nPop
#' @return
#' @export
#' @examples
run_mcmc <-
  function(paramList,
           simData,
           effect) {
    checkpoints <- seq(from = 1,
                       to = chainLength,
                       length.out = checkpoint)
    rValues <- list()
    c <- 1
    for (i in 1:paramList$chainLength) {
      params <-
        get_parameters(simData, effect) #generate parameters for the starting vector of random effect sizes
      #select one locus at random
      l <- sample(1:simData$numVarSites, 1) 
      #vector of possible effect sizes
      e <- -100:100 / 200
      
      f <- params$obsFreq[[l]] # observed frequency of the allele
      m <- params$exclMean[[l]] #estimated average effect of all other alleles
      Vg <- params$exclVar[[l]]  #genetic varience of all other alleles
      
      #calculate Ve based on Vg and assumed h2
      Ve <- getVe(Vg, paramList$h2) 
      
      sdCombined <- sqrt(Vg + Ve) # need to add enviromental varience
      
      #calculate probability of the allele being in the upper/lower pool
      pHi <- calc_pHi(f, m, e, sd=sdCombined)
      pLo <- calc_pLo(f, m, e, sd=sdCombined)  #<----------------------------------#NOT USED!!!
      
      #observed frequency of the allele in the upper pool
      obsHi <- simData$mHigh[[l]] / 2
      
      #n = number of times the allele is observed
      n <- f * paramList$nPop * paramList$nInd * 4 
      y <- obsHi * paramList$nPop * paramList$nInd * 2
      
      #find the posterior 
      post <- find_posterior(y, n, pHi)
      #store maximum value
      effect[l] <- e[which.max(post)]
      
      if (i %in% checkpoints) {
        r <- get_r_squared(simData, effect)
        rValues[[c]] <- r
        info <- paste("Checkpoint", c, "- Iteration", i, "- r2 =", r, sep = " ")
        print(info)
        c <- c + 1
      }
    }
    return(list(estEffect = effect / paramList$nLoci,
                rValues = rValues))
  }

#' Create a vector of phenotype classification
#' @param nPop 
#' @param nInd 
#' @param cutoff 
#' @return
#' @export
#' @examples
#' y <- create_y_vector(10,100,0.2)
create_y_vector <- function(nPop, nInd, cutoff){
  y <- c(rep(1, nPop * nInd * cutoff), rep(0, nPop * nInd * cutoff))
  return(y)
}

#' rrBLUP for individual data 
#' @param simData 
#' @return
#' @export
#' @examples
get_rrBLUP_ind <- function(simData, paramList) {
  y <-
    create_y_vector(nPop = paramList$nPop,
                    nInd = paramList$nInd,
                    cutoff = paramList$cutoff)
  H <- simData$HiInd
  L <- simData$LoInd
  HiLo <- rbind(H, L)
  #use rrBLUP to find most likley effect sizes to use as start point 
  glmm_start <- mixed.solve(y, HiLo)
  startEffects <- glmm_start$u
  return(startEffects)
}

#' Simulate Test Population 
#' @param simParams
#' @param nLoci
#' @param nInd
#' @param h2
#' @param fixed
#' @return a single test population of nInd individuals
#' @export
#' @examples
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

#' Test Predictive Power in Test Population
#' @param testPop 
#' @param estEffect 
#' @param cutoff 
#' @return
#' @export
#' @examples
test_prediction <- function(testPop, estEffect, cutoff){
  est_bv <- testPop$genotype %*% (estEffect * 100) #why multiply by 100? 
  known_bv <- testPop$breedingValue
  
  pdf( "bv_vs_bv.pdf" )
  plot(est_bv ~ known_bv,
       main = "Predictive Power in Test Population",
       xlab = "Simulated Breeding Value",
       ylab = "Estimated Breeding Value")
  dev.off()
  
  mod <- lm(known_bv ~ est_bv)
  mod_summary <- summary(mod)
  r2_bv <- mod_summary$adj.r.squared
  
  predTopIndBV <- est_bv > quantile(est_bv, (1 - cutoff))
  predLowIndBV <- est_bv < quantile(est_bv, cutoff)
  
  trueBVTop <- known_bv[predTopIndBV]
  trueBVLow <- known_bv[predLowIndBV]
  bv <- c(trueBVTop, trueBVLow) / 4
  state <- c(rep("High", length(trueBVTop)), rep("Low", length(trueBVLow)))
  testData <- data.frame(state, bv)
  
  pdf("test_boxplot.pdf")
  plot <- boxplot(testData$bv ~ testData$state,
                  main = "Testing",
                  ylab = "Simulated Breeding Value",
                  xlab = "Classification Based on Estimated Breeding Value")
  dev.off()
  
  return(list(testData = testData,
              R2 = r2_bv,
              testBoxPlot = plot))
}

#' Test for Significant Difference
#'Students paired T-test
#' @param a 
#' @param b 
#' @return
#' @export
#' @examples
is_sig_diff <- function(a, b){
  if (length(a) == length(b)){
    test_result <- t.test(a, b, paired = TRUE)
    print( paste("p-value:", test_result$p.value, sep = " "))
    if (test_result$p.value > 0.05){
      print("No significant difference")
    }
    if (test_result$p.value <= 0.05){
      print("Significant difference!!")
    }
  }
  else{
    print("Error: inputs into test are not of equal length")
  }
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
pooling_simulation <- function(h, t, reps) {
  result_tab <-
    data.frame(
      loci = integer(),
      h2 = numeric(),
      threshold = numeric(),
      var_sites = integer(),
      ind_cor = numeric(),
      pool_cor = numeric(),
      ttest_bv = numeric(),
      ttest_es = numeric()
    )
  line <- 1
  
  for (rep in 1:reps){
    #result_tab[line,] <- rep(0,7)
    
    #need to phase out these references - proper contained functionalisation
    chainLength <- 10 ^ 3
    checkpoint <- 10
    
    paramList <- list(
      nLoci = 100,
      h2 = h,
      cutoff = t,
      nPop = 10,
      nInd = 100,
      MAF = 0.01,
      chainLength = 10 ^ 3,
      checkpoint = 10
    )
    
    check_input_values(paramList)
    
    print("SIMULATION")
    simParams <- get_sim_params(nLoci = paramList$nLoci)
    
    #generate the simulated data
    simData <-
      simulate_training_pools(
        simParams,
        nLoci = paramList$nLoci,
        nPop = paramList$nPop,
        nInd = paramList$nInd,
        h2 = paramList$h2,
        cutoff = paramList$cutoff,
        MAF = paramList$MAF
      )
    print(paste("Simulated", paramList$nPop, "populations with", paramList$nInd, "individuals each", sep = " "))
    print(paste("Total loci = ", paramList$nLoci, sep = " "))
    print(paste("Number of varient loci:", simData$numVarSites, sep = " "))
    print(paste("Captured heritability: h2 = ", paramList$h2, sep = " "))
    
    #simulate a test population of individuals
    testPop <-
      simulate_test_pop(simParams, paramList$nLoci, nInd = 150, paramList$h2, simData$fixedSites, paramList$cutoff)
    
    ######################################
    print("INDIVIDUAL PREDICTION")
    start <- Sys.time()
    totalInd <- paramList$nPop * paramList$nInd
    print(paste(totalInd * paramList$cutoff * 2, "Individuals used out of", totalInd, sep = " "))
    rrBLUPind <- get_rrBLUP_ind(simData, paramList)
    indR <- get_r_squared(simData, rrBLUPind)
    
    testResultsInd <- test_prediction(testPop, rrBLUPind, paramList$cutoff)
    print( paste("Predictive r2 in test population (rrBLUP):", testResultsInd$R2, sep = " "))
    
    ########################################
    print("POOL PREDICTION")
    print(paste("Total of", paramList$nPop * 2, "pools containing", paramList$cutoff * paramList$nInd, "individuals each", sep = " " ))
    
    effectStart <- run_rrBLUP(simData, paramList$nPop)
    testResults1 <- test_prediction(testPop, effectStart, paramList$cutoff)
    print( paste("Predictive r2 in test population (rrBLUP):", testResults1$R2, sep = " "))
    
    ttest_bv <- t.test(testResultsInd$testData$bv, testResults1$testData$bv, paired = TRUE)
    ttest_es <- t.test(rrBLUPind, effectStart, paired = TRUE)
    
    result_tab[line,]<- data.frame(paramList$nLoci, paramList$h2, paramList$cutoff, simData$numVarSites, testResultsInd$R2, testResults1$R2, ttest_bv$p.value, ttest_es$p.value)
    
    line<-line+1
  }
  
}

##################RUN################

result_tab<- data.frame(loci=integer(), h2=numeric(), threshold=numeric(), var_sites=integer(),ind_cor=numeric(),pool_cor=numeric(),ttest_bv=numeric(),ttest_es=numeric())
line <- 1
h_list<- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
t_list<- c(0.1, 0.2, 0.3, 0.4, 0.5)

for (h in h_list){
  for (t in t_list){
    for (rep in 1:10){
      #result_tab[line,] <- rep(0,7)
      
      #need to phase out these references - proper contained functionalisation
      chainLength <- 10 ^ 3
      checkpoint <- 10
      
      paramList <- list(
        nLoci = 100,
        h2 = h,
        cutoff = t,
        nPop = 10,
        nInd = 100,
        MAF = 0.01,
        chainLength = 10 ^ 3,
        checkpoint = 10
      )
      
      check_input_values(paramList)
      
      print("SIMULATION")
      simParams <- get_sim_params(nLoci = paramList$nLoci)
      
      #generate the simulated data
      simData <-
        simulate_training_pools(
          simParams,
          nLoci = paramList$nLoci,
          nPop = paramList$nPop,
          nInd = paramList$nInd,
          h2 = paramList$h2,
          cutoff = paramList$cutoff,
          MAF = paramList$MAF
        )
      print(paste("Simulated", paramList$nPop, "populations with", paramList$nInd, "individuals each", sep = " "))
      print(paste("Total loci = ", paramList$nLoci, sep = " "))
      print(paste("Number of varient loci:", simData$numVarSites, sep = " "))
      print(paste("Captured heritability: h2 = ", paramList$h2, sep = " "))
      
      #simulate a test population of individuals
      testPop <-
        simulate_test_pop(simParams, paramList$nLoci, nInd = 150, paramList$h2, simData$fixedSites, paramList$cutoff)
      
      ######################################
      print("INDIVIDUAL PREDICTION")
      start <- Sys.time()
      totalInd <- paramList$nPop * paramList$nInd
      print(paste(totalInd * paramList$cutoff * 2, "Individuals used out of", totalInd, sep = " "))
      rrBLUPind <- get_rrBLUP_ind(simData, paramList)
      indR <- get_r_squared(simData, rrBLUPind)
      
      testResultsInd <- test_prediction(testPop, rrBLUPind, paramList$cutoff)
      print( paste("Predictive r2 in test population (rrBLUP):", testResultsInd$R2, sep = " "))
      
      ########################################
      print("POOL PREDICTION")
      print(paste("Total of", paramList$nPop * 2, "pools containing", paramList$cutoff * paramList$nInd, "individuals each", sep = " " ))
      
      effectStart <- run_rrBLUP(simData, paramList$nPop)
      testResults1 <- test_prediction(testPop, effectStart, paramList$cutoff)
      print( paste("Predictive r2 in test population (rrBLUP):", testResults1$R2, sep = " "))
      
      ttest_bv <- t.test(testResultsInd$testData$bv, testResults1$testData$bv, paired = TRUE)
      ttest_es <- t.test(rrBLUPind, effectStart, paired = TRUE)
      
      result_tab[line,]<- data.frame(paramList$nLoci, paramList$h2, paramList$cutoff, simData$numVarSites, testResultsInd$R2, testResults1$R2, ttest_bv$p.value, ttest_es$p.value)
      
      line<-line+1
    }}}

write.table(result_tab, file="sim_loop_results.table", row.names = FALSE, quote = FALSE)