#create data object using a subset - not quite working yet

make_rerun_object <- function(pools_rc_file, info_file, ees_table, subset_size = 10){
  sorted_ees <- order_by_ees(ees_table)
  top_ees <- head(sorted_ees, subset_size)
  rc <- fread(pools_rc_file)
  snp_names <- get_snp_id(rc)
  snps_to_use <-  subset(rc, snp_names %in% top_ees$SNP)
  maa_freq_d <- fraction_to_decimal(get_allele_freq(snps_to_use, info_file, "major"))
  mia_freq_d <- fraction_to_decimal(get_allele_freq(snps_to_use, info_file, "minor"))
  info <- fread(info_file)
  colnames(info) <- c("Sample", "Group", "Group2")
  y <- info$Group
  prov <- info$Group2
  major <- get_allele(snps_to_use$allele_states, "major")
  minor <- get_allele(snps_to_use$allele_states, "minor")
return(list(
  y = y,
  prov = prov,
  maa = maa_freq_d,
  mia = mia_freq_d,
  snp_id = get_snp_id(snps_to_use),
  major = major,
  minor = minor
))
}



# # 1. Get Starting pools_rc file
# snps_to_use <- read.table(
#   file = paste(wd, "snps_used.table",
#                sep = "/"),
#   header = TRUE
# )
# # 2.b. Read in info on pools
# pool_info <- fread(info)
# print(head(pool_info))
# pool_col_names <- get_pool_colnames(pool_info)
#
# # 2. Get the estimated effect sizes from the first run
# prior_effect_sizes <- read.table(
#   file = paste(wd, "ees.table", sep = "/"),
#   sep = "\t",
#   header = TRUE
# )
# head(prior_effect_sizes)
#
# # 3. Pick a subset of the top X sites e.g. 100
# sorted_ees_prior <-
#   prior_effect_sizes[order((prior_effect_sizes$MIA.EES) ^ 2, decreasing = TRUE),]
# head_ees <- head(sorted_ees_prior, 100)
# head(head_ees)
#
# # # 4. create a pool_rc table for the subset - match
# snp_id <- paste(snps_to_use$contig, snps_to_use$pos, sep = "_")
# subset_snps <- subset(snps_to_use, snp_id %in% head_ees$SNP)
# head(subset_snps)
#
# # ######################
# #
# # # 3. Get major and minor alleles
# major_allele <- get_allele(subset_snps$allele_states, "major")
# minor_allele <- get_allele(subset_snps$allele_states, "minor")
# #
# # # 4. Get columns containing minor and major allele counts
# maa_freq <- subset_snps[, pool_col_names$maa]
# mia_freq <- subset_snps[, pool_col_names$mia]
# snp_names <- paste(subset_snps$contig, subset_snps$pos, sep = "_")
# #
# # # 5. Convert fractions to decimals
# maa_dec <- fraction_to_decimal(as.data.frame(maa_freq))
# mia_dec <- fraction_to_decimal(as.data.frame(mia_freq))
# #
# # # 6. Get rrBLUP params and check lengths
# pheno <- pool_info$Health
# prov <- rank(pool_info$Prov)
# #
# # if ( length(pheno) != length(prov)){
# #   stop("Pool Number Mismatch - check input files")
# # }
# #
# # # 7. Find divergence from provinence means
# maa_freq_diff <- remove_prov_effects(maa_dec, prov)
# mia_freq_diff <- remove_prov_effects(mia_dec, prov)
# #
