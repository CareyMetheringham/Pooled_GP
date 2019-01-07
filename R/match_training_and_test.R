#' Find SNPs Occuring in Training and Test
#'
#' @param ees_table
#' @param gt_ind
#'
#' @return
#' @export
#'
#' @examples
#' match_snps_in_ind(fread("./extdata/test.ees_table"), read_gt_table("./extdata/test.gt"))
match_snps_in_ind <- function(ees_table, gt_ind){
  use_snps <- subset(ees_table, ees_table$SNP %in% gt_ind$SNP)
  return(use_snps)
}

#' Sort EES Table by MIA^2
#'
#' @param ees_table
#'
#' @return
#' @export
#'
#' @examples
order_by_ees <- function(ees_table){
  ordered_matches <-
    ees_table[order( (ees_table$EES.MIA) ^ 2, decreasing = TRUE), ]
  return(ordered_matches)
}

#' Find and Fix Major Allele Mismatches - NOT WORKING!!!
#'
#' @param ees_table
#' @param gt_table
#' @param fix_table
#' @param pop_data
#'
#' @return
#' @export
#'
#' @examples
fix_allele_mismatch <- function(ees_table, gt_table, fix_table, pop_data){
  ees_and_maa <- data.frame(ees_table, pop_data$major)
  colnames(ees_and_maa)[6] <- "MAJOR"
  match_by_snp <- merge(ees_and_maa, fix_table, by = "SNP")
  corrected_gt <- correct_gt(gt_table, match_by_snp)
  return(corrected_gt)
}

#' Correct Non Matching Alleles
#' @param gt
#' @param match_snps
#'
#' @return
#' @export
#'
#' @examples
correct_gt <- function(gt, match_snps){
  swapped_gt <- gt
  for (i in 1:nrow(gt)){
    my_row <- rownames(gt)[i]
    if (my_row %in% row.names(match_snps)){
      if (match_snps[my_row, "MAJOR"] != match_snps[my_row, "REF"]){
        # Mismatch - swap values
        swapped_gt[i, ] <- sub("0", "NO", gt[i, ])
        swapped_gt[i, ] <- sub("2", "0", gt[i, ])
        swapped_gt[i, ] <- sub("NO", "2", gt[i, ])
      }
    }
  }
  swapped_gt[is.na(swapped_gt)] <- 0
  return(swapped_gt)
}

#' Get a Subset of SNPs with Largest EES
#'
#' @param ees_table
#' @param subset_size
#'
#' @return
#' @export
#'
#' @examples
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
#' @return
#' @export
#'
#' @examples
get_gt_subset <- function(snp_list, gt){
  gt_subset <- subset(gt, row.names(gt) %in% snp_list)
  return(gt_subset)
}

#' Match Test & Training, Correct & Subset
#'
#' @param ees_table
#' @param gt
#' @param fix
#' @param pool_data
#' @param subset_size
#'
#' @return
#' @export
#'
#' @examples
match_and_subset <- function(ees_table, gt, fix, pool_data, subset_size){
  corrected_mismatch <-
    fix_allele_mismatch(ees_table, gt, fix, pool_data)
  match_snps <- match_snps_in_ind(ees_table, corrected_mismatch)
  subset_snps <- get_ees_subset(match_snps, subset_size)
  gt_subset <- get_gt_subset(subset_snps$SNP, corrected_mismatch)
  return(list(gt = gt_subset,
              ees = subset_snps))
}
