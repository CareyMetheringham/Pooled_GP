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
    c("contig",
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
  hits <- hits[order(as.numeric(paste(hits$P))), ]
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
#' @param frac_data a data frame containing only the fractions
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
#' @param pools_rc_file input file pools_rc
#' @param info table containing data on the grouping of pools <- needs to be read in
#' @param gwas file containing snps and p-values from gwas
#' @param hit_num number of snps to be used to calculate ees
#'
#' @return list data read in from pools_rc
#' @export
#'
#' @examples
#'read_in_pools_rc("./extdata/test.pool_rc", fread("./extdata/test.pool_info"), "./extdata/test.gwas", 10)
read_in_pools_rc <- function(pools_rc_file, info, gwas, hit_num){
  colnames(info) <- c("Sample", "Group", "Group2")
  top_gwas_hits <- get_hits_from_file(gwas, hit_num)
  snps_to_use <- find_top_snps(fread(pools_rc_file, stringsAsFactors = FALSE), top_gwas_hits, info)
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
    snp_id = snp_names, # use as rownames
    major = major,
    minor = minor
  ))
}


#' Exclude one named Group2 grouping
#' #can only exclude one at a time
#need to include tests for equal length
#'
#' @param info data table with 3 columns
#' @param pool_data data output by read_in_pools_rc
#' @param group_name a group in column 3 (Group2) of the info table
#'
#' @return pool_data excluding the named group
#' @export
#'
#' @examples
exclude_group2 <- function(info, pool_data, group_name){
  new_y = pool_data$y[info[, 3] != group_name]
  new_prov = pool_data$prov[info[, 3] != group_name]
  new_maa = pool_data$maa[, info[, 3] != group_name]
  new_mia = pool_data$mia[, info[, 3] != group_name]
  return(list(
    y = new_y,
    prov = new_prov,
    maa = new_maa,
    mia = new_mia,
    snp_id = pool_data$snp_id,
    major = pool_data$major,
    minor = pool_data$minor
  ))
}

#' Get random subset pools_rc
#'
#' @param pools_rc_file input file pools_rc
#' @param info table containing data on the grouping of pools <- needs to be read in
#' @param gwas file containing snps and p-values from gwas
#' @param hit_num number of snps to be used to calculate ees
#'
#' @return list of a random subset of data read in from pools_rc
#' @export
#'
#' @examples
#'read_in_pools_rc("./extdata/test.pool_rc", fread("./extdata/test.pool_info"), 10)
read_pools_rc_random <- function(pools_rc_file, info, hit_num){
  colnames(info) <- c("Sample", "Group", "Group2")
  random_hits <- get_random(fread(pools_rc_file, stringsAsFactors = FALSE), hit_num)
  snps_to_use <- find_top_snps(fread(pools_rc_file, stringsAsFactors = FALSE), random_hits[, 1:2], info) #ERROR HERE!!
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
    snp_id = snp_names, # use as rownames
    major = major,
    minor = minor
  ))
}

#' Get SNPs of Interest from GWAS Results
#'
#' @param pools_rc data frame
#' @param num_hits number of snps needed
#'
#' @return random subset of rows from pools_rc
#' @export
#'
#' @examples
#' get_random(fread("./extdata/test.pool_rc"), 10)
get_random <- function(pools_rc, num_hits){
  my_sample <- pools_rc[sample(nrow(pools_rc), num_hits),]
  return(my_sample)
}
