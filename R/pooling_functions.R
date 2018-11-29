#' Indentify Variant and Invaraiant Sites
#'
#' @param gt matrix of genotypes for all populations
#' @param MAF minor allele frequency cutoff
#'
#' @return #boolean vector where TRUE is a varient site
#' @export
#'
#' @examples
#' find_varient_sites(sim_training_pops(10, 100, 100, 0.5)$gt_matrix, 0.01)
find_varient_sites <- function(gt, MAF) {
  sums <- colSums(gt)
  total_ind <- nrow(gt)
  #minumum nInd of times allele can occur
  min_count <- (total_ind * MAF)
  max_count <- (total_ind * (2 - (MAF)))
  variant <- (sums >= min_count & sums <= max_count)
  return(variant)
}

get_pools <- function(sim, MAF = 0.01){
  num_var_sites <- sum(find_varient_sites(sim$gt_matrix, MAF))
  num_pop <- length(sim$bv)
  high_pool_matrix <- matrix(nrow = num_pop, ncol= num_var_sites)
  for (i in 1:num_pop){
    sample <- get_samples()
  }
  #low_pool_matrix <-
  return(high_pool_matrix)
}

get_samples <- function(phenotypes, threshold = 0.5){

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

#end product - data structure for rrBLUP - ydiff, prov, gt_freq_matrix
#count varient sites = sum(v, na.rm = TRUE)

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


#geno_matrix <- matrix(nrow = pop * ind, ncol = sites)



# find varient sites

# create a hilo matrix type = ind, type=pool

# get observed frequencies

#get means of high and low?





# function(simParams,
#          nLoci = 1000,
#          nPop = 10,
#          nInd = 1000,
#          h2 = 0.5,
#          cutoff = 0.5,
#          MAF = 0.05) {
#
#   #set up lists to fill with values
#   geno_list <- list()
#   bv_list <- list()
#   pheno_list <- list()
#   varience_list <- list()
#
#   for (i in 1:nPop) {
#     population <- create_pop(es=simParams$effectSizes, nInd, nLoci, allelicF=simParams$allelicF, h2)
#calculate breeding values
# bv <- (lociP %*% es)
# #normalise breeding values
# bv <- scale(bv)
# #heritability factor used to calculate varience
# var <- calculate_varience( bv = bv, h2 = h2 )
# envVar <- rnorm( length(bv), 0, var)
# #Store values in lists
# return(list(
#   geno = lociP,
#   bv = bv,
#   pheno = (bv + envVar),
#   var = var
# ))
#     geno_list[[i]] <- population$geno
#     bv_list[[i]] <-population$bv
#     pheno_list[[i]] <- population$pheno
#     varience_list[[i]] <- population$var
#   }
#
#   #create matrix for genotypes
#   genotypes <- create_geno_matrix(geno_list, nPop, nInd, nLoci)
#
#   #look only at varient num_sites
#   variant <- find_varient_num_sites(genotypes, nPop, nInd, nLoci,MAF)
#
#   #find the difference in allele frequency between the extremes
#   HiLo <- create_hilo_matrix(nPop, variant,pheno_list, geno_list, cutoff)
#   HiLoInd <- create_hilo_ind_matrix(nPop, nInd, variant, pheno_list, geno_list,cutoff)
#
#   mean_high <- colMeans(HiLo$hi)
#   mean_low <- colMeans(HiLo$lo)
#
#   obsFreq <- colSums(genotypes) / (nInd * nPop) / 2
