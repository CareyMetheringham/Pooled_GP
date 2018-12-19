#' Get P values for each SNP
#'
#' @param pools_rc
#' @param pop_info_file
#'
#' @return
#' @export
#'
#' @examples
#' pools_rc <- find_top_snps(find_pools_rc("./extdata/Pools_RC"), get_hits_from_file("./extdata/example_100_hits.gwas", 10), "./extdata/example_pop_data.csv")
#' get_snp_pvals(pools_rc, "./extdata/example_pop_data.csv")
get_snp_pvals <- function(pools_rc, pop_info_file){
  info <- fread(pop_info_file)
  df <- cbind(paste(pools_rc$contig, pools_rc$pos, sep = "_"), pools_rc)
  colnames(df)[1] <- "SNP"
  p_values_from_glm <- data.frame(SNP = character(), P = numeric())
  for (snp in 1:nrow(df)){
    table <- as.data.frame(create_lm_table(df[snp, ], info))
    mod <- glm(cbind(major, minor) ~ prov + health, data = as.data.frame(table), family = "binomial")
    p_val <- (coefficients(summary(mod)))[14,4]
    p_values_from_glm <- rbind(p_values_from_glm, cbind(df$SNP[snp], p_val))
  }
  return(p_values_from_glm)
}

#' Create a Data Table for GLM
#'
#' @param pool_rc_line
#' @param info
#'
#' @return
#' @export
#'
#' @examples
create_lm_table <- function(pool_rc_line, info){
  major <- unlist(get_allele_count(pool_rc_line[, 13:43]))
  minor <- unlist(get_allele_count(pool_rc_line[, 44:74]))
  prov <- info$Prov
  health <- info$Health
  lm_table <- cbind(major, minor, prov, health)
  return(lm_table)
}

#' Get Allele Count
#'
#' @param allele_cols
#'
#' @return
#' @export
#'
#' @examples
get_allele_count <- function(allele_cols){
  allele_split <- lapply(allele_cols, split_allele_count)
  return(allele_split)
}

#' Split Allele count Cell
#'
#' @param allele_count
#'
#' @return
#' @export
#'
#' @examples
split_allele_count <- function(allele_count){
  split <- strsplit(as.character(allele_count), "/")[[1]]
  return(split[1])
}
