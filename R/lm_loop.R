#' Get P values for each SNP
#'
#' @param pools_rc a table containing read counts in each pool for snps of interest
#' @param pop_info_file the location of the file containing pool info
#'
#' @return a table containing p_values for each snp
#' @export
#'
#' @examples
#' get_snp_pvals(fread("./extdata/example_pool_rc"), fread("./extdata/example_pop_data.csv"))
get_snp_pvals <- function(pools_rc, pop_info){
  df <- cbind(paste(pools_rc$contig, pools_rc$pos, sep = "_"), pools_rc)
  colnames(df)[1] <- "SNP"
  p_values_from_glm <- data.frame(SNP = character(), P = numeric())
  for (snp in 1:nrow(df)){
    table <- as.data.frame(create_lm_table(df[snp, ], pop_info))
    mod <- glm(cbind(major, minor) ~ prov + health, data = as.data.frame(table), family = "binomial")
    p_val <- (coefficients(summary(mod)))[14,4]
    p_values_from_glm <- rbind(p_values_from_glm, cbind(df$SNP[snp], p_val))
  }
  return(p_values_from_glm)
}

#' Create a Data Table for GLM
#'
#' @param pool_rc_line a line from the pool-rc table + SNP  column!!
#' @param info table contain pop info
#'
#' @return a table for use in glm
#' @export
#'
#' @examples
#' create_lm_table(cbind("SNP_N", fread("./extdata/example_pool_rc")), fread("./extdata/example_pop_data.csv"))
create_lm_table <- function(pool_rc_line, info){
  major <- unlist(get_allele_count(get_allele_freq(pool_rc_line, info, "major")))
  minor <- unlist(get_allele_count(get_allele_freq(pool_rc_line, info, "minor")))
  prov <- info$Prov
  health <- info$Health
  lm_table <- cbind(major, minor, prov, health)
  return(lm_table)
}

#' Get Allele Count
#'
#' @param allele_cols selected cols from pool_rc table
#'
#' @return the allele count from those cols
#' @export
#'
#' @examples
#' get_allele_count(c("1/5", "3/7", "2/13"))
get_allele_count <- function(allele_cols){
  allele_split <- lapply(allele_cols, split_allele_count)
  return(allele_split)
}

#' Split Allele count Cell
#'
#' @param allele_count alele count in format "1/2
#'
#' @return allele count
#' @export
#'
#' @examples
#' split_allele_count("10/20")
split_allele_count <- function(allele_count){
  split <- strsplit(as.character(allele_count), "/")[[1]]
  return(split[1])
}
