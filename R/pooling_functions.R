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
