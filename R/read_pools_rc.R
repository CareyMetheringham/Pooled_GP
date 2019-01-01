#' Create Lists of Colnames for pools_rc
#' @param pool_info_file file containing a table - need to make more flexible
#'
#' @return list of colnames for pools_rc file
#' @export
#'
#' @examples
#' get_pool_colnames(fread("./extdata/example_pop_data.csv"))
get_pool_colnames <- function(pool_info){
  colnames(pool_info) <- c("Sample", "Group", "Group2")
  mia_names <- paste("mia", pool_info$Sample, sep = "_")
  maa_names <- paste("maa", pool_info$Sample, sep = "_")
  pool_col_names <-
    c(
      "contig",
      "pos",
      "rc",
      "allele_count",
      "allele_states",
      "deletion_sum",
      "snp_type",
      "major_alleles.maa.",
      "minor_alleles.mia.",
      maa_names,
      mia_names)
  return(list(
    all = pool_col_names,
    maa = maa_names,
    mia = mia_names
  ))
}

#' Find all .pools_rc files in a directory
#'
#' @param pool_loc path to a directory containing wanted files
#'
#' @return list of pools_rc files
#' @export
#'
#' @examples
#' find_pools_rc("./extdata/Pools_RC")
find_pools_rc <- function(pool_loc){
  pool_files <- list.files(path = pool_loc,
                           pattern = "pools_rc.*", full.names = TRUE)
  return(pool_files)
}

#' Pull Top GWAS Hits from pools.rc files
#' @param pool_files list of pools_rc files - can be 1
#' @param top_hits table of snps to use
#' @param pool_info_file info on pools inc groupings and phenotype
#'
#' @return table with the top snps as listed in hits file
#' @export
#'
#' @examples
#' find_top_snps(
#'fread("./extdata/Pools_RC/pools_rc.1"),
#'get_hits_from_file("./extdata/example_100_hits.gwas", 10),
#'fread("./extdata/example_pop_data.csv"))
find_top_snps <- function(pools_rc, hits, pool_info) {
  #get the column names
  pool_col_names <- get_pool_colnames(pool_info)$all
    colnames(pools_rc) <- pool_col_names
    pools_rc$pos <- as.character(pools_rc$pos)
    top_snps <- merge(hits, pools_rc, by = c("contig", "pos"))
  return(top_snps)
}

#' Get SNPs of Interest from GWAS Results
#'
#' @param p_val_file file containing data: contig pos snp_id p-value
#' @param num_hits number of snps needed
#'
#' @return subset of input with lowest p-value
#' @export
#'
#' @examples
#' get_hits_from_file("./extdata/example_100_hits.gwas", 10)
get_hits_from_file <- function(p_val_file, num_hits){
  hits <- fread(file = p_val_file)
  colnames(hits) <- c("contig", "pos", "snp", "P")
  hits$pos <- as.character(hits$pos)
  hits <- hits[order(hits$P), ]
  top_hits <- head(hits, num_hits)
  return(top_hits)
}

#' Get the Major/Minor Allele
#' @param alleles column containing major minor e.g A/T
#' @param which which allele to get
#'
#' @return either major or minor allele
#' @export
#'
#' @examples
#' get_allele("A/T", "major")
get_allele <- function(alleles, which){
  if (which == "major") {
    pick <- 1
  }
  if (which == "minor") {
    pick <- 2
  }
  allele <-
    sapply(strsplit(as.character(alleles), "/"), function(alleles) {
      alleles <- alleles[pick]
    })
  return(allele)
}

#' Get the Major/Minor Allele Frequency
#' @param pool_rc table of pools_rc data
#' @param pool_info_file file with info on training population
#' @param which which allele to get
#'
#' @return table with freqency of either major or minor alleles in each pool
#' @export
#'
#' @examples
#'get_allele_freq(fread("./extdata/example_pool_rc"), fread("./extdata/example_pop_data.csv"), "minor"
get_allele_freq <- function(pool_rc, pool_info, which){
  pool_col_names <- get_pool_colnames(pool_info)
  if (which == "major") {
    freq_data <- pool_rc[, pool_col_names$maa, with = FALSE]
  }
  if (which == "minor") {
    freq_data <- pool_rc[, pool_col_names$mia, with = FALSE]
  }
  return(as.data.frame(freq_data))
}

#' Transform allele frequencies from fraction to decimal
#' @param frac_data# a data frame containing only the fractions
#'
#' @return frequency of allele
#' @export
#'
#' @examples
#' fraction_to_decimal(get_allele_freq(fread("./extdata/example_pool_rc"), "./extdata/example_pop_data.csv", "minor"))
fraction_to_decimal <- function(frac_data) {
  dec_data <- data.frame()
  for (i in 1:length(frac_data)) {
    temp <-
      sapply(strsplit(as.character(frac_data[, i]), "/"), function(frac_data) {
        frac_data <- as.numeric(frac_data)
        frac_data[1] / frac_data[2]
      })
    dec_data <- rbind(dec_data, temp)
  }
  dec_data <- t(dec_data)
  colnames(dec_data) <- colnames(frac_data)
  return(dec_data)
}

#' Get Joined SNP ID from pools_rc format
#'
#' @param pool_rc file
#'
#' @return list of snps
#' @export
#'
#' @examples
#' get_snp_id(fread("./extdata/example_pool_rc"))
get_snp_id <- function(pool_rc){
  snp_names <- paste(pool_rc$contig, pool_rc$pos, sep = "_")
  return(snp_names)
}

#' Get Data Structure from pools_rc
#'
#' @param dir
#' @param info
#' @param gwas
#' @param hit_num
#'
#' @return
#' @export
#'
#' @examples
#'read_in_pools_rc(find_pools_rc("./extdata/Pools_RC"), fread("./extdata/example_pop_data.csv"), "./extdata/example_100_hits.gwas", 10)
read_in_pools_rc <- function(pools_rc_file, info, gwas, hit_num){
  colnames(info) <- c("Sample", "Group", "Group2")
  top_gwas_hits <- get_hits_from_file(gwas, hit_num)
  snps_to_use <- find_top_snps(fread(pools_rc_file), top_gwas_hits, info)
  snp_names <- get_snp_id(snps_to_use)
  maa_freq <- get_allele_freq(snps_to_use, info, "major")
  mia_freq <- get_allele_freq(snps_to_use, info, "minor")
  maa_freq_d <- fraction_to_decimal(maa_freq)
  mia_freq_d <- fraction_to_decimal(mia_freq)
  y <- info$Group
  prov <- info$Group2
  major <- get_allele(snps_to_use$allele_states, "major")
  minor <- get_allele(snps_to_use$allele_states, "minor")
  return(list(
    y = y,
    prov = prov,
    maa = maa_freq_d,
    mia = mia_freq_d,
    snp_id = snp_names,
    major = major,
    minor = minor
  ))
}
