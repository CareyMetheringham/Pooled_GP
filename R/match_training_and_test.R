#' Find SNPs Occuring in Training and Test
#'
#' @param ees_table
#' @param gt_ind
#'
#' @return list of snps which occur in training and test
#' @export
#'
#' @examples
#' match_snps_in_ind(fread("./extdata/test.ees_table"), read_gt_table("./extdata/test.gt"))
match_snps_in_ind <- function(ees_table, gt_ind){
  use_snps <- subset(ees_table, ees_table$SNP %in% rownames(gt_ind))
  return(use_snps)
}

#' Sort EES Table by MIA^2
#'
#' @param ees_table
#'
#' @return snps ordered by estimated effect size
#' @export
#'
#' @examples
#' order_by_ees(fread("./extdata/test.ees_table")) - test that top > bottom
order_by_ees <- function(ees_table){
  ordered_matches <-
    ees_table[order( (ees_table$EES.MIA) ^ 2, decreasing = TRUE), ]
  return(ordered_matches)
}

#' Find and Fix Major Allele Mismatches
#'
#' @param ees_table
#' @param gt_table
#' @param fix_table
#' @param pop_data
#'
#' @return gt table corrected for any mismatched minor and major alleles
#' @export
#'
#' @examples
#' ind_fix <- read_fix_table("./extdata/test.fix")
#' ind_gt <- read_gt_table("./extdata/test.gt")
#' test_data <- read_in_pools_rc("./extdata/test.pool_rc", fread("./extdata/test.pool_info"), "./extdata/test.gwas", 10)
#' fix_allele_mismatch(fread("./extdata/test.ees_table"), ind_gt, ind_fix, test_data)
fix_allele_mismatch <- function(ees_table, gt, fix, major_snp){
  ees_and_maa <- merge(ees_table, major_snp, by = "SNP") #need to ensure order is correct
  match_snps <- merge(ees_and_maa, fix, by = "SNP")
  corrected_gt <- correct_gt(gt, match_snps)
  return(corrected_gt)
}

#' Correct Non Matching Alleles - consider using gsub instead
#' @param gt
#' @param match_snps
#'
#' @return swapped gt
#' @export
#'
#' @examples
correct_gt <- function(gt, match_snps){
  swapped_gt <- data.frame()
  for (i in 1:nrow(gt)){
    my_row <- rownames(gt)[i]
    if (my_row %in% match_snps$SNP){
      fix_line <- match_snps[match_snps$SNP == my_row, ]
      gt_line <- gt[i, ]
      if (fix_line$MAJOR[1] != fix_line$REF[1]){
        gt_line[gt[i,] == 2] <- 0
        gt_line[gt[i,] == 0] <- 2
      }
      swapped_gt <- rbind(swapped_gt, gt_line)
    }
  else{next}
  }
  swapped_gt[is.na(swapped_gt)] <- 0
  rownames(swapped_gt) <- rownames(gt)
  colnames(swapped_gt) <- colnames(gt)
  return(swapped_gt)
}

#' Get a Subset of SNPs with Largest EES
#' @param ees_table
#' @param subset_size
#'
#' @return subset of snps selected by ees
#' @export
#'
#' @examples
#' get_ees_subset(fread("./extdata/test.ees_table"), 5)
get_ees_subset <- function(ees_table, subset_size){
  ordered_by_ees <- order_by_ees(ees_table)
  ees_subset <- head(ordered_by_ees, subset_size)
  return(ees_subset)
}

#' Subset GT Matrix by SNP list
#'
#' @param snp_list
#' @param gt
#'
#' @return subset of gt table
#' @export
#'
#' @examples
#' ind_gt <- read_gt_table("./extdata/test.gt")
#' get_gt_subset(get_ees_subset(fread("./extdata/test.ees_table"), 5)$SNP, ind_gt)
get_gt_subset <- function(snp_list, gt){
  gt_subset <- subset(gt, row.names(gt) %in% snp_list)
  return(gt_subset)
}

#' Match Test & Training, Correct & Subset
#' @param ees_table
#' @param gt
#' @param fix
#' @param pool_data
#' @param subset_size
#'
#' @return list containing subsets of the ees and gt tables
#' @export
#'
#' @examples
#' ind_gt <- read_gt_table("./extdata/test.gt")
#' ind_fix <- read_fix_table("./extdata/test.fix")
#' test_data <- read_in_pools_rc("./extdata/test.pool_rc", fread("./extdata/test.pool_info"), "./extdata/test.gwas", 10)
#' match_and_subset(fread("./extdata/test.ees_table"), ind_gt, ind_fix, test_data, 5)
match_and_subset <- function(ees_table, gt, fix, pool_data, subset_size){
  major_and_snp <- data.table(pool_data$snp_id, pool_data$major)
  colnames(major_and_snp) <- c("SNP", "MAJOR")
  match_snps <- match_snps_in_ind(ees_table, gt)
  subset_ees <- get_ees_subset(match_snps, subset_size)
  gt_subset <- get_gt_subset(subset_ees$SNP, gt)
  corrected_mismatch <-
    fix_allele_mismatch(subset_ees, gt_subset, fix, major_and_snp)

  return(list(gt = corrected_mismatch,
              ees = subset_ees))
}
