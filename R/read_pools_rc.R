#' Create Lists of Colnames for pools_rc
#' @param pool_info_file file containing a table with ...
#'
#' @return list of colnames for pools_rc file
#' @export
#'
#' @examples
#' get_pool_colnames("./extdata/example_pop_data.csv")
get_pool_colnames <- function(pool_info_file){
  pool_info <- read.table(pool_info_file, sep = "\t")
  colnames(pool_info) <- c("Sample", "Group1", "Group2")
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
#' @examples find_pools_rc("./extdata/Pools_RC")
find_pools_rc <- function(pool_loc){
  pool_files <- list.files(path = pool_loc,
                           pattern = "pools_rc.*", full.names = TRUE)
  return(pool_files)
}

#' Pull Top GWAS Hits from pools.rc files
#' @param pool_files
#' @param hits_file
#' @param pool_col_names
#'
#' @return table with the top snps as listed in hits file
#' @export
#'
#' @examples
find_top_snps <- function(pool_files, hits_file, pool_info_file) {
  #get the column names
  pool_col_names <- get_pool_colnames(pool_info_file)$all
  hits <- read.table(file = hits_file, sep = "\t")
  colnames(hits) <- c("contig", "pos", "snp", "P")
  topsnp_list <- list()
  for ( i in 1:length(pool_files)) {
    pools_rc <- read.table(file = pool_files[i], sep = "\t")
    colnames(pools_rc) <- pool_col_names
    my_top_snps <- merge(hits, pools_rc, by = c("contig", "pos"))
    topsnp_list[[i]] <- my_top_snps
  }
  #join the top tables
  top_snp <- rbindlist(topsnp_list)
  return(top_snp)
}

#' Transform allele frequencies from fraction to decimal
#' @param frac_data# a data frame containing only the fractions
#'
#' @return frequency of allele
#' @export
#'
#' @examples
fraction_to_decimal <- function(frac_data) {
  dec_data <- data.frame()
  for (i in 1:length(frac_data)) {
    temp <-
      sapply(strsplit(as.character(frac_data[, i]), "/"), function(frac_data) {
        frac_data <- as.numeric(frac_data)
        frac_data[1] / frac_data[2]
      })
    dec_data <- rbind(temp, dec_data)
  }
  dec_data <- t(dec_data)
  colnames(dec_data) <- colnames(frac_data)
  return(dec_data)
}
